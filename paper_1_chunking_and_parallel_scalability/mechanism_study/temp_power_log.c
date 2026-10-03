/*
 * temp_power_log: one process that logs, once per second, the current phase label (read from a phase file),
 * Tctl and the twelve per-CCD temperatures (k10temp hwmon), the CPU package energy (RAPL "power/energy-pkg",
 * system-wide) and the cumulative pages swapped in/out (/proc/vmstat). No child processes are started, so the log
 * keeps running under heavy load. It stops when the phase file is removed.
 * Columns: unix_time, phase, Tctl, Tccd1..Tccd12 (degrees C), package_energy_joules, pswpin, pswpout.
 * Needs kernel.perf_event_paranoid <= 0 (system-wide counters).
 *
 * Usage:  temp_power_log <output.csv> <phase file> <hwmon directory>
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

static long perf_event_open_syscall(struct perf_event_attr *attr, pid_t pid, int cpu, int group_fd,
                                    unsigned long flags) {
    return syscall(__NR_perf_event_open, attr, pid, cpu, group_fd, flags);
}

static double read_number(const char *path) {
    FILE *f = fopen(path, "r");
    double value = -1.0;
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

static int read_phase(const char *path, char *phase, size_t size) {
    FILE *f = fopen(path, "r");
    if (!f) return 0;
    if (!fgets(phase, (int)size, f)) phase[0] = '\0';
    fclose(f);
    phase[strcspn(phase, "\r\n")] = '\0';
    return 1;
}

static void read_swap(long long *pswpin, long long *pswpout) {
    FILE *f = fopen("/proc/vmstat", "r");
    char key[64];
    long long value;
    *pswpin = -1; *pswpout = -1;
    if (!f) return;
    while (fscanf(f, "%63s %lld", key, &value) == 2) {
        if (strcmp(key, "pswpin") == 0) *pswpin = value;
        if (strcmp(key, "pswpout") == 0) *pswpout = value;
    }
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <output.csv> <phase file> <hwmon directory>\n", argv[0]);
        return 2;
    }
    const char *power_dir = "/sys/bus/event_source/devices/power";
    char path[512];
    snprintf(path, sizeof(path), "%s/type", power_dir);
    int type = (int)read_number(path);
    snprintf(path, sizeof(path), "%s/events/energy-pkg", power_dir);
    uint64_t config = (uint64_t)read_number(path);
    snprintf(path, sizeof(path), "%s/events/energy-pkg.scale", power_dir);
    double scale = read_number(path);
    struct perf_event_attr attr;
    memset(&attr, 0, sizeof(attr));
    attr.size = sizeof(attr);
    attr.type = (uint32_t)type;
    attr.config = config;
    int energy_fd = (int)perf_event_open_syscall(&attr, -1, 0, -1, 0);
    if (energy_fd < 0) fprintf(stderr, "temp_power_log: cannot open power/energy-pkg: %s\n", strerror(errno));
    FILE *out = fopen(argv[1], "a");
    if (!out) return 1;
    fprintf(out, "unix_time,phase,Tctl");
    for (int k = 1; k <= 12; ++k) fprintf(out, ",Tccd%d", k);
    fprintf(out, ",package_energy_joules,pswpin,pswpout\n");
    char phase[256];
    while (read_phase(argv[2], phase, sizeof(phase))) {
        struct timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        snprintf(path, sizeof(path), "%s/temp1_input", argv[3]);
        fprintf(out, "%.3f,%s,%.3f", (double)now.tv_sec + 1e-9 * (double)now.tv_nsec, phase,
                read_number(path) / 1000.0);
        for (int k = 3; k <= 14; ++k) {
            snprintf(path, sizeof(path), "%s/temp%d_input", argv[3], k);
            fprintf(out, ",%.3f", read_number(path) / 1000.0);
        }
        uint64_t energy = 0;
        double joules = -1.0;
        if (energy_fd >= 0 && read(energy_fd, &energy, sizeof(energy)) == (ssize_t)sizeof(energy))
            joules = scale * (double)energy;
        long long pswpin, pswpout;
        read_swap(&pswpin, &pswpout);
        fprintf(out, ",%.3f,%lld,%lld\n", joules, pswpin, pswpout);
        fflush(out);
        sleep(1);
    }
    fclose(out);
    return 0;
}
