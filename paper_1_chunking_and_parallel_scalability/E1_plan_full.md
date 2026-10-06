<!-- ----------------------------------------------------------------------------------------------------------------------------- -->
# Experiment 1 (E1): full plan
<!-- ----------------------------------------------------------------------------------------------------------------------------- -->

Full plan of Experiment 1 (E1) of the paper "Cache-based chunking dramatically improves parallel scaling for HMC, for autodiff and manual gradients: Application to multivariate probit in Stan and NicoStan/BayesMVP" (section 3.1 of the paper gives a shortened version).
Section, table and equation numbers refer to the paper.

## 3.1 E1: Plan

E1 selects the best configuration of each implementation, which the later experiments then use.

### 3.1.1 E1 Plan; Part I: NicoStan+BayesMVP, chunking only: chunking implementations and determining the optimal $N_{\text{chunks}}$

#### The BayesMVP chunking implementations (standard vs. chunking)

Part I assessed two BayesMVP configurations (see section 3.2.1 for results):

- (i) `MD_BayesMVP`: The BayesMVP baseline has $N_{\text{chunks}} = 1$ and $N_{\text{threads/chain}} = 1$.

- (ii) `MD_BayesMVP_chunking`: The same as `MD_BayesMVP` (i.e., $N_{\text{threads/chain}} = 1$), but with multiple sequential observation-level chunks (i.e., $N_{\text{chunks}} > 1$).

#### Candidate $N_{\text{chunks}}$ and iterations for the BayesMVP implementations

We estimated the cost of each candidate $N_{\text{chunks}}$ by running sampling-only HMC for a fixed $N_{\text{iter}}$ with $10$ leapfrog steps ($L = 10$) per iteration, and taking the mean adjusted wall-clock time over four independent runs (see table 3, and section 2.3 for more details).

**Table 3 of the paper.** Candidate $N_{\text{chunks}}$ values, and number of HMC iterations ($N_{\text{iter}}$) without WCP,
for NicoStan/BayesMVP, for each dataset size $N$
(see table 4 for the Stan implementations,
and table 18 for the $N_{\text{iter}}$ with WCP).
We also varied $N_{\text{chains}}$ (and hence threads), namely:
$N_{\text{chains}} \in \{1,\; 2,\; 4,\; 8,\; 16,\; 32,\; 64,\; 96,\; 180\}$ on the HPC,
and $N_{\text{chains}} \in \{1,\; 2,\; 4,\; 8,\; 16\}$ on the laptop.
$L$ is the number of leapfrog steps.
On the local-HPC, we also tested $N_{\text{chunks}} = 6$ and $8$ at $N{=}500$ and $N{=}2500$
for the one-chain, one-thread reference runs
(see section 2.3).
On the laptop, we also tested $N_{\text{chunks}} = 8$ at $N{=}500$ and $N{=}2500$ for these reference runs.

| $N$ | $N_{\text{chunks}}$ | $N_{\text{iter}}$/$L$ (NicoStan+BayesMVP) |
|---|---|---|
| $500$ | $\{1,\; 2,\; 4,\; 10\}$ | $400/10$ |
| $2500$ | $\{1,\; 2,\; 4,\; 10,\; 25\}$ | $80/10$ |
| $10{,}000$ | $\{1,\; 4,\; 10,\; 25,\; 50,\; 100\}$ | $80/10$ |
| $50{,}000$ | $\{1,\; 4,\; 10,\; 25,\; 50,\; 100,\; 250,\; 500,\; 1000\}$ | $50/10$ |

### 3.1.2 E1 Plan; Part II: NicoStan+BayesMVP with within-chain parallelism (WCP-only vs. chunking + WCP vs. chunking-only)

#### The BayesMVP WCP configurations (WCP-only vs. chunking + WCP)

Part II assessed two BayesMVP WCP configurations (see section 3.2.2 for results):

- (i) `MD_BayesMVP_WCP`: The "WCP-only" configuration, with $N_{\text{threads/chain}} > 1$ and one chunk per thread (i.e., $N_{\text{chunks}} = N_{\text{threads/chain}}$).

- (ii) `MD_BayesMVP_WCP_chunking`: The "chunking + WCP" configuration, in which each thread can evaluate several smaller chunks, each of which can fit within the CPU cache (see section 6.2 for evidence on this).

#### BayesMVP chunk and WCP selection

For both chunked BayesMVP configurations (i.e., `MD_BayesMVP_chunking` and `MD_BayesMVP_WCP_chunking`), we selected the fastest $N_{\text{chunks}}$ separately for each device, $N$, $N_{\text{chains}}$ and $N_{\text{threads/chain}}$.

In E2 (see sections 4.1 and 2.3), for each $N_{\text{threads}}$ we use the fastest combination of $N_{\text{chains}}$, $N_{\text{threads/chain}}$, and $N_{\text{chunks}}$ ($N_{\text{threads}} = N_{\text{chains}} \cdot N_{\text{threads/chain}}$) for each configuration in table 1; hence, WCP-only ($N_{\text{chunks}} = N_{\text{threads/chain}}$) and chunking + WCP ($N_{\text{chunks}} \ge N_{\text{threads/chain}}$) are selected separately.

For `MD_BayesMVP_WCP_chunking`, the candidate $N_{\text{chunks}}$ were those in table 3, restricted to $N_{\text{chunks}} \ge N_{\text{threads/chain}}$, together with $N_{\text{chunks}} = N_{\text{threads/chain}}$. Furthermore, the candidate $N_{\text{threads/chain}}$ for the BayesMVP WCP configurations (i.e., those of the allocations, $N_{\text{chains}} \times N_{\text{threads/chain}}$, in table 10) were also used for the Stan model and Mplus (see sections 3.1.4 and 3.1.5). Furthermore, on the HPC, we also used $N_{\text{chains}} \times N_{\text{threads/chain}} = 32 \times 2$, $48 \times 2$, $64 \times 2$, $90 \times 2$, $24 \times 4$, $32 \times 4$ and $45 \times 4$ at every $N$ (only the $\times 2$ allocations for Mplus); for these, the timed runs used fewer iterations (except for BayesMVP at $N \le 2500$, and for Mplus at $N{=}50{,}000$), and their times were scaled to the same $N_{\text{iter}}$ (i.e., $T_{\mathrm{adj}}$ in equation 1, multiplied by $N_{\text{iter}}/n_l$).

### 3.1.3 E1 Plan; Part III: Stan model (via NicoStan), chunking only: chunking implementations and determining the optimal $N_{\text{chunks}}$

#### The Stan chunking implementations (standard vs. container-chunking vs. tape-chunking)

Part III assessed three Stan variations (see section 3.2.3 for results):

- (i) `AD_Stan`: The plain (or baseline) Stan implementation - this uses our standard, unpartitioned, and vectorised LC-MVP Stan model; hence, each $\text{lp\_grad}{\left(\cdot\right)}$ evaluation is made over the full set of observations.

- (ii) `AD_Stan_chunked`: This is the "naive" autodiff-chunking (i.e., "container-chunked") implementation, which retains the same parameterisation and $\text{lp\_grad}{\left(\cdot\right)}$ as `AD_Stan`, but evaluates the likelihood in a sequential loop over observation blocks of size $\lceil N/N_{\text{chunks}}\rceil$. With $N_{\text{threads/chain}} = 1$, it tests whether reducing the size of the `stan::math::var` containers alone changes autodiff behaviour without dividing the autodiff tape (e.g., via `reduce_sum()`).

- (iii) `AD_Stan_tape_chunked`: The "tape-chunked" implementation uses Stan's `reduce_sum_static()` function (Stan Development Team, 2024; Weber, 2020) with $N_{\text{threads/chain}} = 1$. Here, partitioning acts on the autodiff computation itself, rather than only on the `stan::math::var` containers; furthermore, the static partitioning is deterministic, with at most $\lceil N/N_{\text{chunks}}\rceil$ observations assigned to each partial sum. Note that, since `reduce_sum_static()` repeatedly halves the data until each partial sum has at most $\lceil N/N_{\text{chunks}}\rceil$ observations, the actual number of partial sums is the smallest power of two which is $\ge N_{\text{chunks}}$ (e.g., $16$ for $N_{\text{chunks}} = 10$, and $64$ for $N_{\text{chunks}} = 50$); hence, for `AD_Stan_tape_chunked` - and for `AD_Stan_WCP` and `AD_Stan_WCP_chunking`, which use the same `.stan` file - $N_{\text{chunks}}$ is the requested (rather than the actual) number of chunks.

#### Candidate $N_{\text{chunks}}$ and iterations for the Stan implementations

For Stan, we used the NicoStan+BayesMVP cost procedure (see section 3.1.1) with $L = 1$ (see table 4).

**Table 4 of the paper.** Candidate $N_{\text{chunks}}$ values for the three chunked-Stan implementations
(i.e., `AD_Stan_chunked`, `AD_Stan_tape_chunked`,
and `AD_Stan_WCP_chunking` [see section 3.1.4]),
and the number of HMC iterations ($N_{\text{iter}}$) for the three Stan implementations without WCP
(i.e., `AD_Stan`, `AD_Stan_chunked` and `AD_Stan_tape_chunked`)
(see table 18 for the two Stan WCP implementations),
for each dataset size $N$.
Note that, for the plain Stan implementation (`AD_Stan`), $N_{\text{chunks}} = 1$ is the standard LC-MVP
Stan model.
$L$ is the number of HMC leapfrog steps.
On the local-HPC, for the one-chain, one-thread reference runs of `AD_Stan_tape_chunked`
(see section 2.3),
we also tested $N_{\text{chunks}} = 6$ and $8$ at $N{=}500$,
and $N_{\text{chunks}} = 2$, $6$ and $8$ at $N{=}2500$.
On the laptop, we also tested $N_{\text{chunks}} = 8$ at $N{=}500$,
and $N_{\text{chunks}} = 2$ and $8$ at $N{=}2500$.

| $N$ | $N_{\text{chunks}}$ (Stan [via NicoStan]) | $N_{\text{iter}}$/$L$ (Stan [via NicoStan]) |
|---|---|---|
| $500$ | $\{2,\; 4,\; 10\}$ | $200/1$ |
| $2500$ | $\{4,\; 10,\; 25,\; 50,\; 100,\; 250\}$ | $40/1$ |
| $10{,}000$ | $\{4,\; 10,\; 25,\; 50,\; 100,\; 250,\; 500,\; 1000\}$ | $20/1$ |
| $50{,}000$ | $\{10,\; 25,\; 50,\; 100,\; 250,\; 500,\; 1000,\; 2000,\; 5000\}$ | $10/1$ |

### 3.1.4 E1 Plan; Part IV: Stan model (via NicoStan) with within-chain parallelism (WCP-only vs. chunking + WCP vs. tape-chunking only)

#### The Stan WCP implementations (WCP-only vs. chunking + WCP)

Part IV assessed two Stan WCP variations (see section 3.2.4 for results):

- (i) `AD_Stan_WCP`: This is the Stan WCP implementation, which uses the same `.stan` model file as `AD_Stan_tape_chunked` (i.e., the serial counterpart of both Stan WCP models), but with $N_{\text{threads/chain}} > 1$ and one chunk per thread (i.e., "WCP-only"; $N_{\text{chunks}} = N_{\text{threads/chain}}$).

- (ii) `AD_Stan_WCP_chunking`: This is the "chunking + WCP" Stan implementation, which is the same as `AD_Stan_WCP` but permits more than one chunk per thread (i.e., $N_{\text{chunks}} \ge N_{\text{threads/chain}}$).

#### Stan chunk and WCP selection

For `AD_Stan_WCP_chunking`, candidate $N_{\text{chunks}}$ were the values in table 4 (see section 3.1.3), restricted to $N_{\text{chunks}} \ge N_{\text{threads/chain}}$, together with $N_{\text{chunks}} = N_{\text{threads/chain}}$; other selection followed BayesMVP (see section 3.1.2).

### 3.1.5 E1 Plan; Part V: Mplus: fitted model, multi-chain timing and iteration controls ("BITERATIONS" vs. "FBITERATIONS")

The standard Mplus (Muthén and Muthén, 1998-2017; Hallquist and Wiley, 2018) configuration (`Mplus_standard`) allocates one processor (that is, one thread - not necessarily a single physical CPU core) to each chain ($N_{\text{threads/chain}} = 1$), whilst the WCP configuration (`Mplus_WCP`) allocates multiple processors within each chain ($N_{\text{threads/chain}} \ge 2$; see section 3.1.2).

Note that, apart from experiment 3 (see section 5.1), we used the single-population LC-MVP model throughout this paper, since Mplus - at the time of writing - only supports WCP for single-population models. However, the two-population LC-MVP model is more appropriate for the datasets used in this paper (which were simulated from two populations; see section 2.1); hence, we used it in experiment 3 (see section 5.1).

Mplus has two iteration-control settings (Muthén and Muthén, 1998-2017):

- (i) `BITERATIONS`: the requested number of iterations is an upper limit, since Mplus stops the chains early once its convergence criterion (`BCONVERGENCE`) is met; hence, we set `BCONVERGENCE=0` (i.e., never met).

- (ii) `FBITERATIONS`: the requested number of iterations is fixed (i.e., there is no convergence check); however, Mplus 8.10 (the version which we used) only runs whole blocks of 100 iterations (i.e., $100 \cdot \lfloor N_{\text{iter}}/100 \rfloor$ iterations per chain, and at least 100, with one exception; for the numbers of iterations which Mplus actually ran, see section 3.2.5).

We tested both settings since, for the same $N_{\text{iter}}$, they do not necessarily lead to the same run time (e.g., if one setting involves additional computations); the results (see section 3.2.5) helped us determine which setting to use for the cross-algorithm/software comparison (see section 4.2).

Each Mplus configuration was run twice: an overhead run with one requested iteration and a timed run with $N_{\text{iter}}$ (see table 18; $n_l < N_{\text{iter}}$ for $N_{\text{chains}} \ge 32$, $N_{\text{threads/chain}} = 2$ and $N \le 10{,}000$). Then, using the notation of equation 1, the adjusted Mplus time is:

$$
T_{\mathrm{adj}} = \frac{T_l - T_s}{n_l^{\mathrm{obs}}} \cdot N_{\text{iter}}
$$
(equation 4 of the paper)

where $n_l^{\mathrm{obs}}$ is the number of iterations Mplus actually ran in the timed run; $N_{\text{iter}}$ is the requested number, which can differ from $n_l^{\mathrm{obs}}$; hence, scaling by $N_{\text{iter}}$ ensures a common iteration budget.

### 3.1.6 E1 Plan; Part VI: Optimal burn-in configurations for NicoStan+BayesMVP and the Stan model (WCP vs. chunking vs. both)

#### NicoStan+BayesMVP

During NicoStan's burn-in (which has options for ChEES-R-HMC (Sountsov and Hoffman, 2022), SNAPER-HMC (Sountsov and Hoffman, 2022), or ChEES-HMC (Hoffman et al, 2021)), the chains advance together by only a single iteration before the shared HMC adaptation quantities are updated. Hence, allocating more threads within each chain can reduce the burn-in time, even when that same configuration is less favourable for the post-burn-in sampling phase.

We measured seconds per burn-in iteration with the HMC settings (step size $\epsilon$, metric $M$, and mean trajectory length $\hat\tau$) held fixed, so every configuration performed the same computations on average.

Burn-in candidates were $N_{\text{chunks}}$ from table 3, restricted to $N_{\text{chunks}} \ge N_{\text{threads/chain}}$, plus $N_{\text{chunks}} = N_{\text{threads/chain}}$; candidate $N_{\text{threads/chain}}$ were $N_{\text{threads/chain}} = 1$ and those in section 3.1.2. We used $N_{\text{burn\_chains}} \in \{4, 8, 16\}$ on the local-HPC and $\{4, 8\}$ on the laptop, timing 200, 40, 10, and 10 iterations at $N = 500$, $2500$, $10{,}000$, and $50{,}000$, respectively, with $3$ repetitions each. Note that this burn-in study only measures computational cost, i.e., the seconds per burn-in iteration of each configuration, with the HMC settings held fixed (see below); this is the same kind of measurement for NicoStan and for Stan (see the Stan via NicoStan and Stan via cmdstanr paragraphs in section 3.1.6). More specifically, since NicoStan pools its adaptation across the burn-in chains, $N_{\text{burn\_chains}}$ is a genuine choice for NicoStan (i.e., for both NicoStan+BayesMVP and the Stan model via NicoStan), and hence we measured the cost per burn-in iteration for each $N_{\text{burn\_chains}}$; however, it is important to note that this study does not show which $N_{\text{burn\_chains}}$ adapts best - i.e., how $N_{\text{burn\_chains}}$ affects the number of burn-in iterations needed - which is beyond the scope of this paper. In contrast, since Stan's chains adapt independently, more burn-in chains cannot reduce the number of burn-in iterations which Stan needs (whatever their cost); hence, for Stan's own burn-in (via cmdstanr), we only measured the cost of the default $N_{\text{burn\_chains}} = 4$ (see section 3.1.6). Note that this is a property of Stan's algorithm, rather than a result of this study.

At each iteration, NicoStan draws the trajectory length uniformly, such that: $\tau \sim \text{Uniform}(0, 2\hat\tau)$, with mean $\hat\tau = L\epsilon$ (where we fixed $L = 20$ and $\epsilon{=}10^{-5}$; i.e., around 20 leapfrog steps on average). This draw was shared across the burn-in chains (i.e., $\tau_{k} = \tau$ for the $k$-th chain), so that trajectory length is not a source of computational variation between chains.

For each device, $N$, and $N_{\text{burn\_chains}}$, we selected the configuration with the lowest mean seconds per burn-in iteration across repetitions. Furthermore, for each device and $N$, table 17 gives the configuration with the largest speed-up over all $N_{\text{burn\_chains}}$ (each relative to $N_{\text{threads/chain}} = 1$ at the same $N_{\text{burn\_chains}}$); note that this compares the cost per burn-in iteration only, not how well each $N_{\text{burn\_chains}}$ adapts. We use two references at the same $N_{\text{burn\_chains}}$: (i) $N_{\text{chunks}} = 1$ and $N_{\text{threads/chain}} = 1$, for the combined chunking and WCP improvement; and (ii) the best measured $N_{\text{chunks}}$ with $N_{\text{threads/chain}} = 1$, for the additional WCP improvement after chunking (see section 3.2.6 for results).

#### Stan via NicoStan

For Stan, we repeated the NicoStan+BayesMVP burn-in design using tape-chunking, with WCP when $N_{\text{threads/chain}} > 1$ (see section 3.1.4).

The $N_{\text{burn\_chains}}$, HMC settings, candidate $N_{\text{threads/chain}}$ and selection of the best configuration were as for NicoStan+BayesMVP; however, the candidate $N_{\text{chunks}}$ were those in table 4, restricted to $N_{\text{chunks}} \ge N_{\text{threads/chain}}$ (see section 3.1.4), together with $N_{\text{chunks}} = N_{\text{threads/chain}} > 1$; furthermore, at $N = 10{,}000$ and $50{,}000$, we timed fewer iterations (5 and 3, vs. 10). Since these exclude $N_{\text{chunks}} = 1$, we only used reference (ii) (i.e., the best measured $N_{\text{chunks}}$ with $N_{\text{threads/chain}} = 1$).

#### Stan via cmdstanr

In contrast to NicoStan, Stan itself (which we ran via cmdstanr (Gabry et al, 2021)) has no between-chain adaptation, i.e., each chain adapts its step size and metric independently; hence, the burn-in chains of Stan's NUTS sampler do not advance together by a single iteration (as NicoStan's do), and the fastest burn-in configuration for the Stan model via NicoStan need not be the fastest for Stan's own burn-in. Furthermore, running more than $N_{\text{burn\_chains}} = 4$ burn-in chains in Stan gives little benefit, but still slows the burn-in down. We therefore also benchmarked Stan's own burn-in with $N_{\text{burn\_chains}} = 4$, to inform the Stan burn-in configuration of E3 (see section 5.1) and of future benchmarks; we compared `AD_Stan`, `AD_Stan_tape_chunked`, `AD_Stan_WCP` and `AD_Stan_WCP_chunking`, with $N_{\text{threads/chain}} \in \{1, 2, 4, 8, 16, 22, 44\}$ on the local-HPC and $\{1, 2, 4\}$ on the laptop, and the candidate $N_{\text{chunks}}$ of the Stan via NicoStan study above. More specifically, we used the same step size ($\epsilon{=}10^{-5}$), timed iterations and repetitions, and a maximum tree depth of $4$ (i.e., $15$ leapfrog steps per iteration, the closest to the mean of $20$ above); we then selected the configuration with the lowest seconds per burn-in iteration of the slowest chain, since the burn-in only finishes once its slowest chain does.

## References

- Gabry J, Češnovar R, Bales B, Morris M, Popov M, Lawrence M (2021). CmdStanR: A lightweight interface to Stan for R users. R package. https://mc-stan.org/cmdstanr/
- Hallquist MN, Wiley JF (2018). MplusAutomation: An R Package for Facilitating Large-Scale Latent Variable Analyses in Mplus. Structural Equation Modeling.
- Hoffman M, Radul A, Sountsov P (2021). An Adaptive-MCMC Scheme for Setting Trajectory Lengths in Hamiltonian Monte Carlo. Proceedings of The 24th International Conference on Artificial Intelligence and Statistics, PMLR. https://proceedings.mlr.press/v130/hoffman21a.html
- Muthén LK, Muthén BO (1998-2017). Mplus User's Guide (8th edition). https://www.statmodel.com/download/MplusUserGuideVer_8.pdf
- Sountsov P, Hoffman MD (2022). Focusing on Difficult Directions for Learning HMC Trajectory Lengths.
- Stan Development Team (2024). Stan User's Guide: Reduce-Sum. https://mc-stan.org/docs/stan-users-guide/reduce-sum.html
- Weber S (2020). New Within-Chain Parallelisation in Stan 2.23. https://statmodeling.stat.columbia.edu/2020/05/05/easy-within-chain-parallelisation-in-stan/



