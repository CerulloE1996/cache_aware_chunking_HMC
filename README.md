<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
# Cache-aware chunking dramatically improves parallel scaling for HMC with autodiff and manual gradients: Application to multivariate probit model in Stan and NicoStan/BayesMVP
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

Enzo Cerullo¹, Olivia Carter², Hayley E Jones³, Tim Lucas¹, Nicola J. Cooper¹, Alex J. Sutton¹

¹ Biostatistics Research Group, Division of Public Health & Epidemiology, School of Medical Sciences, University of Leicester, Leicester, UK

² Queen's Veterinary School Hospital, Cambridge, University of Cambridge, UK

³ Population Health Sciences, Bristol Medical School, University of Bristol, UK

[Abstract](#abstract) ·
[Key results](#key-results) ·
[Repository contents](#repository-contents) ·
[Reproducing the results](#reproducing-the-results) ·
[Related software](#related-software) ·
[How to cite](#how-to-cite)


<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
## Abstract
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

We show that cache-aware chunking - a technique used in high-performance computing, but novel in the context of MCMC sampling - can yield large improvements in parallel scaling.
Applied to the LC-MVP model implemented in NicoStan+BayesMVP, on a 96-core server, chunking was up to ~17.5× faster, and on a consumer laptop it was up to ~11.9× faster, with the optimal chunk count depending on dataset size and the allocation of threads.

These gains are due to a shift in the computational bottleneck: without chunking, the `lp_grad()` evaluation is memory-bandwidth-bound, leaving no room for additional threads to contribute.
With chunking, more of the data fits within the L2/L3 caches, shifting the computation to a more compute-bound regime.
This has an important consequence for simultaneous multithreading (SMT): HMC workloads typically see no benefit - or even degradation - from SMT, whereas chunked NicoStan+BayesMVP achieved SMT efficiency gains of ~20% - 29% on a 96-core AMD server CPU (AMD EPYC 9654), and ~9% - 22% on an 8-core laptop (AMD Ryzen 5800H) - both without 3D V-Cache.

Our results suggested three computational profiles - bandwidth-bound (plain Stan, and NicoStan+BayesMVP without chunking), compute-bound (NicoStan+BayesMVP with chunking), and latency-bound (Mplus) - with predictable consequences for parallel scaling and SMT behaviour.

These parallel scaling results should be interpreted alongside raw efficiency: NicoStan+BayesMVP with chunking is approximately [TBD]× and [TBD]× faster than Mplus and Stan, respectively, in time to achieve a target ESS (at N = 10,000; final values to follow).

Looking ahead, CPUs with larger L3 caches - such as AMD's 3D V-Cache line - would likely amplify these benefits.
Furthermore, tape chunking improved the parallel scaling of the Stan model in this study; similar strategies applied to the Stan math C++ library may improve parallel scaling for general Stan models.

The full BayesMVP algorithm and benchmarks, and acceleration of arbitrary Stan models via NicoStan (with preliminary speed-ups of ~1.5-15×), will be described in two upcoming papers.


<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
## Key results
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

![Figure 1](results/figures/Figure_N_chunks_pilot_study_plot_1_n_threads_SMT_vs_no_SMT_cache_lines.png)

*Figure 1: Mean NicoStan+BayesMVP sampling time (seconds) against the number of chunks (log scale), for N = 500, 2,500, 10,000 and 50,000. Direct comparison of HPC (solid lines) and laptop (dashed lines), without SMT (HPC: 96 threads, laptop: 8 threads; blue lines) and with SMT (HPC: 180 threads, laptop: 16 threads; orange lines), with error bars showing the standard deviation across 4 runs. The vertical lines show the number of chunks at which one chunk (1,608 bytes per individual) first fits within the L3 or L2 cache per active thread, in the colour and line type of the corresponding data line; to the left of an L3 line, one chunk no longer fits within the L3 cache per active thread.*

![Figure 2](results/figures/Figure_ps2_plot_2_adj_scalability_markers.png)

*Figure 2: Normalised parallel scaling for each algorithm, computed as: (N_chains/second) × T_0, where T_0 is the time of each algorithm's own one-chain, one-thread reference at the same number of iterations (for WCP, the corresponding implementation without WCP); i.e., the speed-up over the one-chain, one-thread reference, as defined in the paper (with the same reference as the paper's parallel efficiency tables). Note: this normalisation removes absolute speed differences between algorithms, allowing direct comparison of how well each algorithm's performance scales with additional threads, but relative to its own baseline. Top half: HPC (1-180 threads; WCP allocations from 8 up to 128 or 176 threads). Bottom half: Laptop (1-16 threads; WCP allocations from 4 threads). Higher values indicate greater throughput. The dashed grey line shows perfect linear scaling (i.e., speed-up = number of threads). The dotted vertical lines mark 96 threads (HPC) and 8 threads (laptop), above which SMT is used. Note: these results reflect parallelisation behaviour only, and do not account for differences in effective sample size (ESS) per second between algorithms; a method that scales perfectly here may still require more total computation to achieve equivalent posterior precision.*


<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
## Repository contents
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

| Folder | Contents |
| --- | --- |
| `paper_1_chunking_and_parallel_scalability/` | The R scripts for the benchmarks, figures and tables (see [Reproducing the results](#reproducing-the-results)), the Stan models (`stan_models/`: the unpartitioned model, `LC_MVP_bin_PartialLog_v5.stan`; the chunked model used for the Stan + chunking configurations, `LC_MVP_bin_PartialLog_v5_chunked.stan`; and the model used for the Stan tape-chunking and WCP configurations, `LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan`), and the saved study summaries (`paper_1_computational_outputs/`, `burnin_outputs/`) |
| `paper_1_chunking_and_parallel_scalability/mechanism_study/` | The counter programs (C), run scripts, case lists, and the measured counts and times for the hardware-counter profiling study (Experiment 4) |
| `0_utilities/` | Shared R functions, including the simulation of the COVID-19-based LC-MVP datasets and the Stan compilation settings |
| `1_appendix_pilot_studies/` | The Mplus and Stan pilot-study results, and the NicoStan+BayesMVP pilot-study functions, used for the absolute efficiency comparison (Experiment 3) |
| `results/data/` | CSV files with every measured configuration (e.g., `measured_cases.csv`, `configurations.csv`, `scaling.csv`, `wcp_chunk_search.csv`) |
| `results/figures/`, `results/tables/` | The figures and tables in the paper and its supplementary material |


<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
## Reproducing the results
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

The benchmarks are run from R, using our R packages [NicoStan](https://github.com/CerulloE1996/NicoStan) and [BayesMVP](https://github.com/CerulloE1996/BayesMVP) (see their installation instructions), as well as [BridgeStan](https://roualdes.us/bridgestan/latest/) and [cmdstanr](https://mc-stan.org/cmdstanr/).
The Mplus comparisons also require [Mplus](https://www.statmodel.com/) (commercial software) and the MplusAutomation R package.
The other R packages used are listed in `0_utilities/shared_configs/load_R_packages.R`.

The folders in this repository follow the same structure as our analysis folder; hence, to run the scripts, set `algorithm_study_dir` (and the other paths at the top of each script) to the location of this repository on your computer.
Note that the Stan models were compiled with the AMD AOCC compiler on Linux; the compiler paths and flags are set in `R_fns_alg_paper_1_chunking_WCP_par_scaling.R`.

| Script | What it does |
| --- | --- |
| `alg_paper_1_chunking_WCP_par_scaling.R` | The main sampling benchmark: chunking, WCP and parallel scaling for NicoStan+BayesMVP, the Stan model (via NicoStan) and Mplus, for N = 500, 2,500, 10,000 and 50,000, on the HPC and the laptop |
| `alg_paper_1_figures_tables.R` | The figures and tables, from the saved study summaries (no sampling) |
| `ps_1_burnin_optimizing_N_chunks_and_WCP_threads.R` | The burn-in benchmark (number of chunks and WCP threads for NicoStan+BayesMVP's burn-in) |
| `alg_paper_1_burnin_report.R` | The burn-in figures and tables, from the saved burn-in results |
| `alg_paper_1_experiment_3_table.R` | The absolute efficiency table (time to a target ESS for NicoStan+BayesMVP, Stan and Mplus), from the pilot-study results |
| `mechanism_study/run_mechanism_experiments.sh`, `mechanism_study/analyse_mechanism_study.R`, `mechanism_study/make_paper_figure_exp4.R` | The hardware-counter profiling study (Experiment 4): the counter runs, their summaries, and the DRAM bandwidth figure |

Note that the complete saved run outputs are not included here because of their size (over 200 GB).
Hence, `alg_paper_1_experiment_3_table.R` also needs the saved outputs of our NicoStan+BayesMVP pilot study (or a re-run of that study) for its NicoStan+BayesMVP rows.


<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
## Related software
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

- [NicoStan](https://github.com/CerulloE1996/NicoStan) ([website](https://cerulloe1996.github.io/NicoStan/)): Efficient MCMC for Stan models, with advanced between-chain adaptation and diffusion-pathspace HMC.
- [BayesMVP](https://github.com/CerulloE1996/BayesMVP): Accelerated multivariate probit models using NicoStan.


<!-- ------------------------------------------------------------------------------------------------------------------------------- -->
## How to cite
<!-- ------------------------------------------------------------------------------------------------------------------------------- -->

```bibtex
@misc{Cerullo_2026_cache_aware_chunking,
  title  = {Cache-aware chunking dramatically improves parallel scaling for {HMC} with autodiff and manual gradients: Application to multivariate probit model in {Stan} and {NicoStan}/{BayesMVP}},
  author = {Cerullo, Enzo and Carter, Olivia and Jones, Hayley E. and Lucas, Tim and Cooper, Nicola J. and Sutton, Alex J.},
  year   = {2026}
}
```
