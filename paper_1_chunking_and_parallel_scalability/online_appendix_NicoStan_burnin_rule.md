<!-- ----------------------------------------------------------------------------------------------------------------------------- -->
# Online appendix: the partial sums of `reduce_sum_static()`, and NicoStan's automatic burn-in configuration for tape-chunked Stan models
<!-- ----------------------------------------------------------------------------------------------------------------------------- -->

Online appendix of the paper "Cache-based chunking dramatically improves parallel scaling for HMC, for autodiff and manual gradients: Application to multivariate probit in Stan and NicoStan/BayesMVP".
It gives (i) the number of partial sums which Stan's `reduce_sum_static()` actually runs for a requested $N_{\text{chunks}}$, and (ii) the rule which NicoStan uses to choose the burn-in $N_{\text{threads/chain}}$ and $N_{\text{chunks}}$ of Stan models whose likelihood is summed by `reduce_sum_static()`, together with our in-sample check of this rule against the burn-in timings of experiment 1 (E1; Part VI).
Section, table and equation numbers refer to the paper.

## 1. The partial sums which `reduce_sum_static()` runs

Stan's `reduce_sum_static()` passes the range of the $N$ individuals to `tbb::parallel_reduce()` with TBB's `simple_partitioner`, with a grainsize of $\lceil N/N_{\text{chunks}} \rceil$ individuals (i.e., the `chunk_size` which our Stan models pass to it).
The partitioner halves each range - the first half holding $\lfloor \text{size}/2 \rfloor$ individuals - for as long as it holds more individuals than the grainsize; hence, the number of partial sums which are run is the number of leaves of this recursive halving, which is not necessarily the requested $N_{\text{chunks}}$.

More specifically, for every $N_{\text{chunks}}$ which we tested, the number of partial sums was $N_{\text{chunks}}$ rounded up to a power of two (see table 4 of the paper for the candidate $N_{\text{chunks}}$ at each $N$; e.g., $N_{\text{chunks}} = 10$, $25$, $250$ and $1000$ ran as $16$, $32$, $256$ and $1024$ partial sums).
For instance, for WCP-only (i.e., with the requested $N_{\text{chunks}} = N_{\text{threads/chain}}$), $N_{\text{threads/chain}} = 6$, $11$ or $12$, $22$ or $24$, and $44$ or $45$ ran $8$, $16$, $32$ and $64$ partial sums, respectively, at every $N$; hence, `AD_Stan_WCP` ran one partial sum per thread only when $N_{\text{threads/chain}}$ was a power of two.
Throughout the paper, $N_{\text{chunks}}$ gives these partial sums for `AD_Stan_tape_chunked`, `AD_Stan_WCP` and `AD_Stan_WCP_chunking`; the container-chunked `AD_Stan_chunked` loops over its blocks itself, and so ran exactly the requested $N_{\text{chunks}}$, as did NicoStan+BayesMVP.

Note that, where two requested $N_{\text{chunks}}$ ran as the same partial sums (e.g., $N_{\text{chunks}} = 44$ and $50$, which both ran as $64$ partial sums), they were the same configuration, timed twice.
In the figures and tables of E1 Part III, we used their mean, whereas the other figures show both timings; and, for Stan's own burn-in (via cmdstanr) at $N = 10{,}000$, the near-tie between WCP-only and chunking + WCP with $4 \times 44$ (table 17 of the paper) was between two such timings.

The number of partial sums is not a power of two in general; for instance, $N_{\text{chunks}} = 500$ at $N = 500$ would run $500$ partial sums (one per individual), and $N_{\text{chunks}} = 5000$ at $N = 10{,}000$ would run $5904$ (neither of which we tested).
However, a request of $N_{\text{chunks}} = 2^k$ runs exactly $2^k$ partial sums whenever each holds at least three individuals (i.e., $\lceil N/2^k \rceil \ge 3$).
We checked this by counting the leaves for every $N \le 12{,}000$ (and for $N = 25{,}000$, $50{,}000$ and $100{,}000$) with every $2^k \le N$; in all $7891$ of these $151{,}676$ cases with fewer partial sums than requested, each partial sum held at most two individuals.

The function `fn_number_of_partial_sums_run_by_reduce_sum_static()` (in [`R_fn_number_of_partial_sums_run_by_reduce_sum_static.R`](R_fn_number_of_partial_sums_run_by_reduce_sum_static.R)) counts the partial sums in this way; the figures and tables of the paper which give the Stan $N_{\text{chunks}}$ as partial sums were made with:

- [`alg_paper_1_E1_Part_III_Stan_figures_tables_with_partial_sums_run_by_reduce_sum_static.R`](alg_paper_1_E1_Part_III_Stan_figures_tables_with_partial_sums_run_by_reduce_sum_static.R) (E1 Part III);
- [`alg_paper_1_Stan_chunk_search_figures_with_partial_sums_run_by_reduce_sum_static.R`](alg_paper_1_Stan_chunk_search_figures_with_partial_sums_run_by_reduce_sum_static.R) (the E1 Part IV chunk searches);
- [`alg_paper_1_burnin_report_Stan_with_partial_sums_run_by_reduce_sum_static.R`](alg_paper_1_burnin_report_Stan_with_partial_sums_run_by_reduce_sum_static.R) and [`alg_paper_1_burnin_report_cmdstanr_with_partial_sums_run_by_reduce_sum_static.R`](alg_paper_1_burnin_report_cmdstanr_with_partial_sums_run_by_reduce_sum_static.R) (the E1 Part VI burn-in figures of the Stan model via NicoStan and via cmdstanr).

Each of these reads the saved runs only (no sampling); with the second argument `requested`, they draw the requested $N_{\text{chunks}}$ instead, which reproduces the earlier figures exactly.
The re-scored evaluation of the $N_{\text{chunks}}$ rule for tape-chunked Stan (table 26 of the paper) is in [`rule_evaluation_Stan_partial_sums/`](rule_evaluation_Stan_partial_sums/).

## 2. NicoStan's burn-in rule

NicoStan (`fn_compute_burnin_n_threads_WCP_and_num_chunks()`) chooses the burn-in within-chain threads ($T = N_{\text{threads/chain}}$) and $N_{\text{chunks}}$ of a Stan model which uses `reduce_sum_static()` from $M_{\text{bytes/chain}}$ (see equation 7a of the paper) and the CPUs available to the R process.
More specifically, with $H$ the number of hardware threads available, $N_{\text{L3}}$ the number of L3 caches which they belong to, $H_{\text{L3}}$ the largest number of them which share one L3 cache, $S_{\text{L3}}$ the size of this cache and $N_{\text{chains}}$ the number of burn-in chains:

- $T$ is the largest power of two which is at most $\min\left(\lfloor H / N_{\text{chains}} \rfloor,\; H_{\text{L3}}\right)$ (and at least $1$);
- $N_{\text{active/L3}} = \min\left(H_{\text{L3}},\; \lceil N_{\text{chains}} \cdot T / N_{\text{L3}} \rceil\right)$, and $S_{\text{target}} = S_{\text{L3}} / N_{\text{active/L3}}$ (i.e., the Stan expression of equation 7g of the paper, with the burn-in layout);
- $N_{\text{chunks}} = 2^{\left\lceil \log_2 \max\left(\lceil M_{\text{bytes/chain}} / S_{\text{target}} \rceil,\; T \cdot j\right) \right\rceil}$, at most $2^{\lfloor \log_2 N \rfloor}$, where $j$ is the number of partial sums per thread.

In other words, each burn-in chain gets the hardware threads available per chain, but no more than those of one L3 cache (on the local-HPC, $N_{\text{threads/chain}} = 16$ was the fastest at $N \ge 10{,}000$, with $24$ within run-to-run noise; see figure 19 and table 17 of the paper), and each partial sum fits within its thread's share of that L3 cache.
Since $T$ and $N_{\text{chunks}}$ are both powers of two, every thread then evaluates the same number of partial sums (see section 1 above).
For a Stan model, NicoStan times the rule's configurations for $j = 1$, $2$ and $4$, at $T$ and at $T/2$ threads per chain, within its own burn-in, and then uses the fastest.

## 3. In-sample check against the E1 burn-in timings

Table A1 compares the rule's choice with the fastest configuration with $N_{\text{burn\_chains}} = 4$ in the E1 burn-in benchmark of the Stan model via NicoStan (see section 3.2.6 and table 17 of the paper; with the $N_{\text{chunks}}$ given as the partial sums which were run, and $M_{\text{bytes/chain}} = 19{,}152$ bytes per individual).
Note that this is an in-sample check, since the rule was derived from these same timings.

**Table A1.** NicoStan's burn-in rule vs. the fastest timed configuration with $N_{\text{burn\_chains}} = 4$, given as $N_{\text{threads/chain}} \times N_{\text{chunks}}$; $j$ is the number of partial sums per thread of the rule's choice, and the run-to-run noise is the coefficient of variation over the repeats of the fastest configuration.

| Machine | $N$ | $M_{\text{bytes/chain}}$ (MiB) | Fastest timed | Rule | $j$ | Rule's configurations timed | Rule slower than fastest (%) | Run-to-run noise (%) |
|---|---:|---:|---|---|---:|---:|---:|---:|
| laptop | $500$ | 9.132 | $4 \times 16$ | $4 \times 16$ | 1 | 3 | 0 | 1.1 |
| | $2500$ | 45.662 | $4 \times 64$ | $4 \times 64$ | 1 | 6 | 0 | 5.0 |
| | $10{,}000$ | 182.648 | $4 \times 256$ | $4 \times 256$ | 1 | 6 | 0 | 13.3 |
| | $50{,}000$ | 913.239 | $4 \times 1024$ | $4 \times 1024$ | 1 | 6 | 0 | 18.4 |
| local-HPC | $500$ | 9.132 | $8 \times 16$ | $8 \times 16$ | 2 | 2 | 0 | 4.5 |
| | $2500$ | 45.662 | $8 \times 32$ | $8 \times 32$ | 4 | 3 | 0 | 5.3 |
| | $10{,}000$ | 182.648 | $16 \times 64$ | $16 \times 64$ | 1 | 6 | 0 | 13.2 |
| | $50{,}000$ | 913.239 | $16 \times 256$ | $16 \times 256$ | 1 | 6 | 0 | 18.1 |

The rule selected the fastest timed configuration in all eight cases.
However, it is important to note that, on the local-HPC at $N \le 2500$, $N_{\text{threads/chain}} = 8$ was the largest tested (see the caption of table 17 of the paper); hence, only the rule's configurations with $T/2 = 8$ threads per chain were timed there (i.e., $2$ and $3$ of them, at $N = 500$ and $2500$, respectively).
