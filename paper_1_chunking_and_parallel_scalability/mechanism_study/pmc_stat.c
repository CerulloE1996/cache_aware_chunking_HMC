/*
 * pmc_stat: counts CPU hardware events (user space only) for one command and every thread/process it starts,
 * using the kernel's perf_event_open interface directly (no perf binary or root access needed).
 *
 * Usage:  pmc_stat <output.csv> <label> -- <command> [arguments...]
 *
 * Appends one CSV line when the command ends: label, wall_seconds, exit_status, max_rss_kb, then one column per event.
 * Snapshots: the command can send SIGUSR1 to pmc_stat (its PID is in the environment variable PMC_STAT_PID) to append a line
 * "<label>:snapshot_<k>" with the running totals at that moment, so a window (e.g. one sampling call) is the difference of two
 * snapshots and excludes start-up work such as loading R and packages.
 * Events (AMD Family 19h, Zen 3 and Zen 4; PMCx044 = "any data cache fills by data source"):
 *   cycles, instructions                         generic hardware events
 *   fills_local_L2   = PMCx044 umask 0x01          L1 data cache fills served by the core's own L2
 *   fills_local_L3   = PMCx044 umask 0x02          fills served by the same CCX (its shared L3, or another core's L2)
 *   fills_other_CCX  = PMCx044 umask 0x04 | 0x10   fills served by a cache in another CCX (same or other node)
 *   fills_DRAM       = PMCx044 umask 0x08 | 0x40   fills served by DRAM or IO (near or far)
 * Six events = the six AMD core counters, so nothing is multiplexed (values are still scaled if the kernel had to).
 */
#define _GNU_SOURCE
#include <linux/perf_event.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <signal.h>

struct event_definition { const char *name; uint32_t type; uint64_t config; };

static const struct event_definition events[] = {
    { "cycles",          PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES   },
    { "instructions",    PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS },
    { "fills_local_L2",  PERF_TYPE_RAW,      0x0144 },
    { "fills_local_L3",  PERF_TYPE_RAW,      0x0244 },
    { "fills_other_CCX", PERF_TYPE_RAW,      0x1444 },
    { "fills_DRAM",      PERF_TYPE_RAW,      0x4844 },
};
#define N_EVENTS ((int)(sizeof(events) / sizeof(events[0])))

static volatile sig_atomic_t snapshot_requested = 0;
static void on_sigusr1(int signal_number) { (void)signal_number; snapshot_requested = 1; }

static void write_line(const char *output_file, const char *label, double wall_seconds, int exit_status, long max_rss_kb, const int *fds) {
    FILE *out = fopen(output_file, "a");
    if (!out) { perror("fopen"); return; }
    fprintf(out, "%s,%.6f,%d,%ld", label, wall_seconds, exit_status, max_rss_kb);
    for (int i = 0; i < N_EVENTS; ++i) {
        uint64_t values[3] = { 0, 0, 0 };
        double scaled = -1.0;
        if (fds[i] >= 0 && read(fds[i], values, sizeof(values)) == (ssize_t)sizeof(values)) {
            scaled = (values[2] > 0) ? (double)values[0] * (double)values[1] / (double)values[2] : 0.0;
        }
        fprintf(out, ",%.0f", scaled);
    }
    fprintf(out, "\n");
    fclose(out);
}

static long perf_event_open_syscall(struct perf_event_attr *attr, pid_t pid, int cpu, int group_fd, unsigned long flags) {
    return syscall(__NR_perf_event_open, attr, pid, cpu, group_fd, flags);
}

int main(int argc, char **argv) {
    if (argc < 5 || strcmp(argv[3], "--") != 0) {
        fprintf(stderr, "usage: %s <output.csv> <label> -- <command> [arguments...]\n", argv[0]);
        return 2;
    }
    const char *output_file = argv[1];
    const char *label = argv[2];
    int go_pipe[2];
    if (pipe(go_pipe) != 0) { perror("pipe"); return 2; }
    struct sigaction action;
    memset(&action, 0, sizeof(action));
    action.sa_handler = on_sigusr1;          /* no SA_RESTART: wait4 returns EINTR so the snapshot is taken promptly */
    sigaction(SIGUSR1, &action, NULL);
    char parent_pid[32];
    snprintf(parent_pid, sizeof(parent_pid), "%d", (int)getpid());
    pid_t child = fork();
    if (child == 0) {
        setenv("PMC_STAT_PID", parent_pid, 1);
        signal(SIGUSR1, SIG_DFL);
        close(go_pipe[1]);
        char go;
        if (read(go_pipe[0], &go, 1) != 1) _exit(126);
        execvp(argv[4], &argv[4]);
        perror("execvp");
        _exit(127);
    }
    close(go_pipe[0]);
    int fds[N_EVENTS];
    for (int i = 0; i < N_EVENTS; ++i) {
        struct perf_event_attr attr;
        memset(&attr, 0, sizeof(attr));
        attr.size = sizeof(attr);
        attr.type = events[i].type;
        attr.config = events[i].config;
        attr.disabled = 1;
        attr.enable_on_exec = 1;
        attr.inherit = 1;
        attr.exclude_kernel = 1;
        attr.exclude_hv = 1;
        attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
        fds[i] = (int)perf_event_open_syscall(&attr, child, -1, -1, 0);
        if (fds[i] < 0) fprintf(stderr, "pmc_stat: cannot open %s: %s\n", events[i].name, strerror(errno));
    }
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    if (write(go_pipe[1], "g", 1) != 1) { perror("write"); return 2; }
    close(go_pipe[1]);
    int status = 0, n_snapshots = 0;
    struct rusage usage;
    memset(&usage, 0, sizeof(usage));
    for (;;) {
        pid_t done = wait4(child, &status, 0, &usage);
        if (snapshot_requested) {
            snapshot_requested = 0;
            clock_gettime(CLOCK_MONOTONIC, &end);
            double t = (double)(end.tv_sec - start.tv_sec) + 1e-9 * (double)(end.tv_nsec - start.tv_nsec);
            char snapshot_label[512];
            snprintf(snapshot_label, sizeof(snapshot_label), "%s:snapshot_%d", label, ++n_snapshots);
            write_line(output_file, snapshot_label, t, -2, 0, fds);
        }
        if (done == child) break;
        if (done < 0 && errno != EINTR) { perror("wait4"); break; }
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    double wall_seconds = (double)(end.tv_sec - start.tv_sec) + 1e-9 * (double)(end.tv_nsec - start.tv_nsec);
    write_line(output_file, label, wall_seconds, WIFEXITED(status) ? WEXITSTATUS(status) : -1, usage.ru_maxrss, fds);
    return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}
