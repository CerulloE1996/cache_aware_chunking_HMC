/*
 * power_log: logs the CPU package energy (RAPL "power/energy-pkg", system-wide) once per second until killed.
 * Columns: unix_time, package_energy_joules (cumulative since start).
 * Needs kernel.perf_event_paranoid <= 0 (system-wide counters).
 *
 * Usage:  power_log <output.csv>
 */
#define _GNU_SOURCE
#include <linux/perf_event.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <errno.h>

static long perf_event_open_syscall(struct perf_event_attr *attr, pid_t pid, int cpu, int group_fd, unsigned long flags) {
    return syscall(__NR_perf_event_open, attr, pid, cpu, group_fd, flags);
}

static double read_number(const char *path) {
    FILE *f = fopen(path, "r");
    double value = 0.0;
    if (f) {
        char buffer[128] = {0};
        if (fgets(buffer, sizeof(buffer), f)) {
            char *equals = strchr(buffer, '=');
            value = strtod(equals ? equals + 1 : buffer, NULL);
        }
        fclose(f);
    }
    return value;
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <output.csv>\n", argv[0]); return 2; }
    const char *dir = "/sys/bus/event_source/devices/power";
    char path[512];
    snprintf(path, sizeof(path), "%s/type", dir);
    int type = (int)read_number(path);
    snprintf(path, sizeof(path), "%s/events/energy-pkg", dir);
    uint64_t config = (uint64_t)read_number(path);
    snprintf(path, sizeof(path), "%s/events/energy-pkg.scale", dir);
    double scale = read_number(path);
    struct perf_event_attr attr;
    memset(&attr, 0, sizeof(attr));
    attr.size = sizeof(attr);
    attr.type = (uint32_t)type;
    attr.config = config;
    int fd = (int)perf_event_open_syscall(&attr, -1, 0, -1, 0);
    if (fd < 0) { fprintf(stderr, "power_log: cannot open power/energy-pkg: %s\n", strerror(errno)); return 1; }
    FILE *out = fopen(argv[1], "a");
    if (!out) return 1;
    fprintf(out, "unix_time,package_energy_joules\n");
    for (;;) {
        uint64_t value = 0;
        struct timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        if (read(fd, &value, sizeof(value)) == (ssize_t)sizeof(value))
            fprintf(out, "%.3f,%.6f\n", (double)now.tv_sec + 1e-9 * (double)now.tv_nsec, scale * (double)value);
        fflush(out);
        sleep(1);
    }
}
