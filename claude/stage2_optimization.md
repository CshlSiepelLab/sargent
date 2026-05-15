# Stage 2 Optimization: Efficient Trace Estimation for Segment Weights

## Problem

In Stage 2 of the SARGENT variational algorithm, the update for each segment
responsibility $\rho_{sm}$ depends on the REML-style score

$$\frac{\partial \mathcal{L}}{\partial w_s^2} \approx \tfrac{1}{2}(A_s - B_s),$$

where

$$A_s = \tilde{\mathbf{r}}^\top \mathbf{G}_s \tilde{\mathbf{r}}, \qquad
B_s = \operatorname{tr}(\mathbf{V}^{-1} \mathbf{G}_s).$$

The current implementation computes $B_s$ by solving
$\mathbf{V}\mathbf{z}_i = \mathbf{G}_s \boldsymbol\epsilon_i$ via CG for each
Hutchinson probe $i$ and each segment $s$, then computing
$\boldsymbol\epsilon_i^\top \mathbf{z}_i$.

**Cost:** $S \times N_{\text{probe}}$ CG solves per outer iteration. Each CG
solve internally requires $k_{\text{CG}}$ MVPs with $\mathbf{V}$ (cost $O(Sn)$
each). Total per-iteration cost is

$$O(S^2 \cdot n \cdot N_{\text{probe}} \cdot k_{\text{CG}}),$$

which is intractable for $S$ in the thousands and impossible for genome-wide
analyses with $S$ in the millions.

## Fix 1: Reuse Hutchinson probes across segments (priority)

Use the cyclic property of the trace to move $\mathbf{V}^{-1}$ off
$\mathbf{G}_s$:

$$B_s = \mathbb{E}_\epsilon\!\left[(\mathbf{V}^{-1}\boldsymbol\epsilon)^\top
\mathbf{G}_s \boldsymbol\epsilon\right] \approx
\frac{1}{N_{\text{probe}}} \sum_{i=1}^{N_{\text{probe}}}
\mathbf{z}_i^\top \mathbf{G}_s \boldsymbol\epsilon_i,
\qquad \mathbf{z}_i \equiv \mathbf{V}^{-1}\boldsymbol\epsilon_i.$$

**Key observation:** the vectors $\{\mathbf{z}_i\}$ do not depend on $s$.
Precompute them **once per outer iteration** via $N_{\text{probe}}$ CG solves of
the form $\mathbf{V}\mathbf{z}_i = \boldsymbol\epsilon_i$. The per-segment work
then reduces to $N_{\text{probe}}$ MVPs with $\mathbf{G}_s$, which is exactly
the operation made cheap by the tree-based MVP algorithm
(Lee et al. 2025; Zhu et al. 2026) at $O(n)$ per MVP.

### Revised per-iteration cost accounting

| Step                                              | Cost                                              |
|---------------------------------------------------|---------------------------------------------------|
| $N_{\text{probe}}$ CG solves for $\mathbf{V}^{-1}\boldsymbol\epsilon_i$ | $O(N_{\text{probe}} \cdot S n \cdot k_{\text{CG}})$ |
| 1 CG solve for $\tilde{\mathbf{r}} = \mathbf{V}^{-1}\mathbf{r}$ | $O(S n \cdot k_{\text{CG}})$ |
| All $A_s$: one $\mathbf{G}_s$-MVP per segment     | $O(S n)$                                          |
| All $B_s$: $N_{\text{probe}}$ $\mathbf{G}_s$-MVPs per segment | $O(S n \cdot N_{\text{probe}})$           |
| **Total**                                         | $O(S n \cdot N_{\text{probe}} \cdot k_{\text{CG}})$ |

Net improvement: factor of $S$ over the original scheme. This is the standard
trick used in BOLT-REML for multi-component variance estimation.

### Implementation note

This is a one-loop refactoring: pull the CG solve out of the segment loop and
replace the per-segment CG solve with a per-segment $\mathbf{G}_s$-MVP. The
loop structure becomes:

```
for outer iteration:
    solve V z_i = epsilon_i  for i = 1..N_probe         # O(N_probe) CG solves
    solve V r_tilde = r                                  # 1 CG solve
    for s = 1..S:
        A_s = r_tilde^T G_s r_tilde                      # 1 G_s-MVP
        B_s = mean_i [ z_i^T G_s epsilon_i ]             # N_probe G_s-MVPs
        update rho_sm using (A_s - B_s)
```

## Fix 2: Tree-sequence batching for genome-wide $S$

For $S$ in the millions, even $O(Sn \cdot N_{\text{probe}})$ in $\mathbf{G}_s$
MVPs is significant. The tskit-style tree sequence representation provides a
further speedup: adjacent local trees differ by a small number of SPR-like
edits, so the full sequence of MVPs $\{\mathbf{G}_s \mathbf{v}\}_{s=1}^S$ can
be computed incrementally in $O(Tn)$ where $T$ is the total number of distinct
edges across the tree sequence (sub-linear in $S$).

Lee et al. 2025 and Zhu et al. 2026 implement this. If the current SARGENT code
computes each $\mathbf{G}_s \mathbf{v}$ independently from a per-segment
double-centered distance matrix, the incremental form is the next optimization
to pursue after Fix 1.

## Fix 3: Constant-$\bar{B}$ approximation (analog of $\bar{u}$)

In the polygenic regime $\mathbf{V} \approx c \mathbf{I}$, so
$B_s \approx \operatorname{tr}(\mathbf{G}_s) / c$, and
$\operatorname{tr}(\mathbf{G}_s)$ is trivial (sum of diagonal entries). More
generally, $B_s / \operatorname{tr}(\mathbf{G}_s)$ should be approximately
constant across segments under polygenicity.

**Calibration procedure:**

1. Sample $\sim 30$ calibration segments.
2. Compute $B_s$ exactly via Hutchinson on each.
3. Fit $B_s \approx \bar{c} \cdot \operatorname{tr}(\mathbf{G}_s)$.
4. Use $\bar{c} \cdot \operatorname{tr}(\mathbf{G}_s)$ for the remaining
   segments.

This mirrors the $\bar{u}$ calibration already used for SNPs (main.tex,
lines 1242–1252) and admits the same justification. Validate empirically by
plotting exact vs. approximated $B_s$ on a real simulation before deploying.

## Fix 4: Stochastic and lazy updates

Most $w_s^2$ will be in the null mixture component, so fresh $B_s$ values are
not needed for every segment on every iteration.

- **Mini-batch updates.** Sample a random subset of segments (e.g., 10%) per
  outer iteration; update only those $\rho_{sm}$. Standard stochastic
  coordinate ascent; converges to the same fixed point.
- **Lazy refresh.** Track when each segment's $\rho_{sm}$ last changed
  materially; refresh $B_s$ less often for segments stuck in the null
  component.

## Variance considerations for $N_{\text{probe}}$

With Rademacher probes, the variance of the Hutchinson estimator scales with

$$\|\mathbf{V}^{-1/2} \mathbf{G}_s \mathbf{V}^{-1/2}\|_F^2.$$

For segments where $\mathbf{G}_s$ has a few large eigenvalues (rare,
high-impact local trees), $N_{\text{probe}}$ may need to be larger than for
SNP $u_j$ calibration. Typical BOLT-REML values are in the 10–30 range, but
the right value depends on the spectrum of $\mathbf{G}_s$ in real ARGs. Check
per-segment Monte Carlo SE on realistic simulations before committing to a
small $N_{\text{probe}}$.

## Recommended implementation order

1. **Fix 1 first.** One-line refactor of the trace loop. Validate the
   correctness of $B_s$ against the current per-segment scheme on a small ARG.
   Benchmark wall-clock improvement.
2. **Then Fix 3** (constant-$\bar{B}$) if $S$ is large enough that the
   $O(SN_{\text{probe}})$ MVP cost is still painful. Validate by comparing
   exact and approximated $B_s$ values.
3. **Fix 2** (tree-sequence batching) for biobank-scale claims.
4. **Fix 4** (stochastic/lazy updates) as further optimization if needed.

## Benchmark suggestion

When the refactor is complete, write a benchmark script that times Stage 2 on
a small simulated ARG ($n = 300$, $S = 100$) before and after each fix. This
gives concrete numbers for the manuscript and catches regressions if any
approximation introduces noticeable error.
