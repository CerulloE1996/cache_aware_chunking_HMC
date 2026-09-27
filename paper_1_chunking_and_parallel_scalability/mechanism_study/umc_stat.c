/*
 * umc_stat: counts DRAM traffic at the memory controllers (AMD Zen 4 "amd_umc_<k>" uncore PMUs, system-wide) while a command
 * runs. Event 0x0a = CAS commands; rdwrmask 1 = reads, 2 = writes; each CAS moves one 64-byte line.
 * Needs kernel.perf_event_paranoid <= 0 (system-wide counters).
 *
 * Usage:  umc_stat <output.csv> <label> -- <command> [arguments...]
 *
 * Like pmc_stat, the command can send SIGUSR1 (PID in the environment variable UMC_STAT_PID) to append a snapshot line
 * "<label>:snapshot_<k>" with the running totals; a final line is written when the command ends.
 * Columns: label, seconds_since_start, exit_status, DRAM_read_bytes, DRAM_write_bytes.
 */
#define _GNU_SOURCE
#include <linux/perf_event.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>
#include <signal.h>

#define MAX_UMC 32
static int read_fds[MAX_UMC], write_fds[MAX_UMC], n_umc = 0;
static volatile sig_atomic_t snapshot_requested = 0;
static void on_sigusr1(int s) { (void)s; snapshot_requested = 1; }

static long perf_event_open_syscall(struct perf_event_attr *attr, pid_t pid, int cpu, int group_fd, unsigned long flags) {
    return syscall(__NR_perf_event_open, attr, pid, cpu, group_fd, flags);
}

static int open_umc(int type, uint64_t config) {
    struct perf_event_attr attr;
    memset(&attr, 0, sizeof(attr));
    attr.size = sizeof(attr);
    attr.type = (uint32_t)type;
    attr.config = config;
    int fd = (int)perf_event_open_syscall(&attr, -1, 0, -1, 0);
    if (fd < 0) fprintf(stderr, "umc_stat: cannot open type %d config 0x%lx: %s\n", type, (unsigned long)config, strerror(errno));
    return fd;
}

static double total(const int *fds) {
    double sum = 0.0;
    for (int k = 0; k < n_umc; ++k) {
        uint64_t value = 0;
        if (fds[k] >= 0 && read(fds[k], &value, sizeof(value)) == (ssize_t)sizeof(value)) sum += (double)value;
    }
    return 64.0 * sum;
}

static void write_line(const char *file, const char *label, double t, int status) {
    FILE *out = fopen(file, "a");
    if (!out) return;
    fprintf(out, "%s,%.6f,%d,%.0f,%.0f\n", label, t, status, total(read_fds), total(write_fds));
    fclose(out);
}

int main(int argc, char **argv) {
    if (argc < 5 || strcmp(argv[3], "--") != 0) { fprintf(stderr, "usage: %s <output.csv> <label> -- <command> [arguments...]\n", argv[0]); return 2; }
    for (int k = 0; k < MAX_UMC; ++k) {
        char path[256];
        snprintf(path, sizeof(path), "/sys/bus/event_source/devices/amd_umc_%d/type", k);
        FILE *f = fopen(path, "r");
        if (!f) break;
        int type = -1;
        if (fscanf(f, "%d", &type) != 1) type = -1;
        fclose(f);
        read_fds[n_umc] = open_umc(type, 0x10a);
        write_fds[n_umc] = open_umc(type, 0x20a);
        ++n_umc;
    }
    struct sigaction action;
    memset(&action, 0, sizeof(action));
    action.sa_handler = on_sigusr1;
    sigaction(SIGUSR1, &action, NULL);
    char parent_pid[32];
    snprintf(parent_pid, sizeof(parent_pid), "%d", (int)getpid());
    struct timespec start, now;
    clock_gettime(CLOCK_MONOTONIC, &start);
    write_line(argv[1], "start_marker_ignore", 0.0, -3);
    pid_t child = fork();
    if (child == 0) {
        setenv("UMC_STAT_PID", parent_pid, 1);
        signal(SIGUSR1, SIG_DFL);
        execvp(argv[4], &argv[4]);
        perror("execvp");
        _exit(127);
    }
    int status = 0, n_snapshots = 0;
    for (;;) {
        pid_t done = waitpid(child, &status, 0);
        if (snapshot_requested) {
            snapshot_requested = 0;
            clock_gettime(CLOCK_MONOTONIC, &now);
            char label[512];
            snprintf(label, sizeof(label), "%s:snapshot_%d", argv[2], ++n_snapshots);
            write_line(argv[1], label, (double)(now.tv_sec - start.tv_sec) + 1e-9 * (double)(now.tv_nsec - start.tv_nsec), -2);
        }
        if (done == child) break;
        if (done < 0 && errno != EINTR) break;
    }
    clock_gettime(CLOCK_MONOTONIC, &now);
    write_line(argv[1], argv[2], (double)(now.tv_sec - start.tv_sec) + 1e-9 * (double)(now.tv_nsec - start.tv_nsec),
               WIFEXITED(status) ? WEXITSTATUS(status) : -1);
    return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}
