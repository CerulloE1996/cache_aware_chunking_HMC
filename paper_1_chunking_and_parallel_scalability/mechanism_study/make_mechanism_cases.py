#!/usr/bin/env python3
"""Writes the mechanism experiment list for one machine: cases_HPC.csv or cases_Laptop.csv.
Columns: label, experiment, algorithm, N, chunks, chains, threads_per_chain, n_iter, cpus ("all" = no pinning)."""
import csv, sys

device = sys.argv[1]
n_iter = {("MD", 50000): 20, ("MD", 10000): 50, ("AD", 50000): 10, ("AD", 10000): 30}
bayesmvp_chunks = {50000: [1, 4, 25, 100, 500], 10000: [1, 4, 25, 100]}
stan_chunks = {50000: [10, 50, 250, 500], 10000: [4, 25, 100]}
rows = []

def add(experiment, algorithm, N, chunks, chains, tpc, cpus, placement):
    family = "MD" if algorithm.startswith("MD") else "AD"
    label = f"{experiment}_{algorithm}_N{N}_chunks{chunks}_chains{chains}_tpc{tpc}_{placement}"
    rows.append([label, experiment, algorithm, N, chunks, chains, tpc, n_iter[(family, N)], cpus])

if device == "HPC":
    chains_E1 = [1, 8, 48, 96, 180]
    one_ccd, spread8 = "0-7", "0,8,16,24,32,40,48,56"
    smt16, two_ccd16 = "0-7,96-103", "0-15"
else:
    chains_E1 = [1, 4, 8, 16]
    one_ccd, spread8 = None, None
    physical8 = "0,2,4,6,8,10,12,14"

## E1: chunk count x number of concurrent chains (one thread per chain), BayesMVP first
for N in (50000, 10000):
    for chunks in bayesmvp_chunks[N]:
        for chains in chains_E1:
            add("E1", "MD_BayesMVP", N, chunks, chains, 1, "all", "unpinned")

if device == "HPC":
    ## E2: 8 chains sharing ONE chiplet's 32 MB L3 vs spread one per chiplet (each with its own L3)
    for placement, cpus in (("oneCCD", one_ccd), ("spread8CCD", spread8)):
        for chunks in (1, 25, 500):
            add("E2", "MD_BayesMVP", 50000, chunks, 8, 1, cpus, placement)
        add("E2", "AD_Stan", 50000, 1, 8, 1, cpus, placement)
        for chunks in (10, 250):
            add("E2", "AD_Stan_tape_chunked", 50000, chunks, 8, 1, cpus, placement)

## E3: within-chain parallelism, 1 chain x 8 threads, WCP-only (chunks = 8) vs chunks + WCP; HPC: one chiplet vs spread
wcp_joint_chunks = {("MD", 50000): 200, ("MD", 10000): 100, ("AD", 50000): 250, ("AD", 10000): 100}
for N in (50000, 10000):
    for algorithm in ("MD_BayesMVP_WCP", "AD_Stan_WCP"):
        family = "MD" if algorithm.startswith("MD") else "AD"
        for chunks in (8, wcp_joint_chunks[(family, N)]):
            if device == "HPC":
                for placement, cpus in (("oneCCD", one_ccd), ("spread8CCD", spread8)):
                    add("E3", algorithm, N, chunks, 1, 8, cpus, placement)
            else:
                add("E3", algorithm, N, chunks, 1, 8, "all", "unpinned")

## E4: SMT. HPC: 16 chains on one chiplet's 8 cores (SMT) vs 16 chains on two chiplets (no SMT); 8 chains on one chiplet is in E2.
##     Laptop: 8 chains on the 8 physical cores (one thread per core) vs 16 chains on all 16 threads (E1).
for chunks, algorithm in ((1, "MD_BayesMVP"), (25, "MD_BayesMVP"), (500, "MD_BayesMVP"), (1, "AD_Stan"), (250, "AD_Stan_tape_chunked")):
    if device == "HPC":
        add("E4", algorithm, 50000, chunks, 16, 1, smt16, "SMT_oneCCD")
        add("E4", algorithm, 50000, chunks, 16, 1, two_ccd16, "noSMT_twoCCD")
    else:
        add("E4", algorithm, 50000, chunks, 8, 1, physical8, "physical_cores_only")

## E1 for Stan last: plain Stan (one chunk) and tape chunking
for N in (50000, 10000):
    for algorithm, chunk_list in (("AD_Stan", [1]), ("AD_Stan_tape_chunked", stan_chunks[N])):
        for chunks in chunk_list:
            for chains in chains_E1:
                add("E1", algorithm, N, chunks, chains, 1, "all", "unpinned")

with open(f"cases_{device}.csv", "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n"); w.writerow(["label", "experiment", "algorithm", "N", "chunks", "chains", "threads_per_chain", "n_iter", "cpus"]); w.writerows(rows)
print(device, len(rows), "cases")
