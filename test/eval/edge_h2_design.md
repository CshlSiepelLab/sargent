# Per-edge posterior contribution to phenotypic variance (Block 2.7.1)

Date: 2026-06-04 PM (design v1).
Author: Claude (laptop session) for Adam.
Backs the v6 reframe's primary SARGENT-distinctive application (see
`~/newwriting/arg-mapping/REFRAME_NOTES.md` § Evolutionary analyses;
`ROADMAP.md` § Block 2.7).

## Motivation

The naïve per-edge variance contribution `Δh²_b = ℓ_b · ŵ²_s` is the
*prior* expected contribution under the fitted model.  Since ŵ²_s is
constant across all edges in tree s, the per-edge rank within a tree is
just by branch length, and the cross-tree rank is by `ŵ²_s × ℓ_b`.
This reduces to "ranked by ŵ²_s with ties broken by branch length."
It doesn't use the data.

We want a **posterior** edge contribution that conditions on the actual
phenotype `y`.  Edges with large posterior contribution are those that
(a) have non-trivial prior support AND (b) carry data signal — i.e.,
edges whose descendants have *consistently signed* BLUP residuals.
These are the edges that separate phenotypically distinguishable groups.

The same calculation, with different aggregation, produces three
SARGENT-distinctive readouts:

1. Per-edge ranked list (top-K, annotated by time / population / position
   / length) — the "linkage-style" output Adam proposed.
2. Per-age-bin aggregation — a posterior sharpening of Block 2.1 v2's
   WHEN axis (replaces per-tree-mean-age binning with per-edge-age
   binning, weighted by the data).
3. Per-population-bin aggregation — pre-split vs post-split branch h²
   contribution under structured demographies; backs the cross-
   ancestry story directly.

## SARGENT's model and per-edge decomposition

Verified against `src/lib/trees.c::tree_cov_matvec` and
`src/lib/variational.c::v_matvec`.

The phenotype model is
> y  =  Xβ + u + ε,    ε ~ N(0, σ² I),    u = Σ_s √(σ² · ŵ²_s) · g_s,
> g_s ~ N(0, G_s),    V := Cov(y) = σ² (I + Σ_s ŵ²_s · G_s).

The per-tree kinship `G_s` is the diploid analog of pairwise shared
branch length:
> G_s[i, j]  =  sum of branch lengths shared between diploid
>               individuals i and j down to MRCA(i, j) in tree s.

Decomposing `tree_cov_matvec`'s post-/pre-order traversal in edge form
(see Appendix A for the algebraic check):
> G_s  =  (1/2) · M_s · diag(ℓ) · M_sᵀ,
> M_s[i, e]  =  c_i^{u(e)}  ∈ {0, 1, 2},

where `c_i^{u(e)}` is the number of individual i's haplotypes
descending from edge e's child node `u(e)`, and `ℓ` is the per-edge
branch-length vector for tree s.  Edges in tree s contribute additively
to the per-tree covariance with a factor of 1/2 from SARGENT's diploid
hap-splitting convention.

The implied **per-edge prior** on a Brownian-style coefficient is then
> α_e^s  ~  N(0, σ² · ŵ²_s · ℓ_e / 2),

and the per-tree random effect can be written
> √(σ² ŵ²_s) g_s  =  Σ_{e in tree s} α_e^s · m_e,
> m_e[i]  =  c_i^{u(e)}.

**Centering check.** Sargent's per-tree REML EM uses raw `G_s` (no
centering); `center_kinship` in `predict.py` is applied only inside
the BLUP arms (GBLUP / SARGENT-U) as the Block 1.3 fix for the rank-1
constant-offset bias.  Our edge calculation therefore works with raw
`G_s` (and raw `α_e`) — no centering correction needed.

## The posterior calculation

Define the BLUP residual
> z  :=  V⁻¹ (y − Xβ).

`z` is already computed inside sargent's CG solve (the .gall
preconditioner makes it cheap), and `predict.py` uses an equivalent
form for the SARGENT prediction arm.  All per-edge / per-tree /
per-bin posterior contributions below reduce to a single
descendant-sum of `z` per edge, plus bookkeeping for the aggregation.

### Per-edge posterior mean

By the standard linear-Gaussian conditional formula (corrected against
the in-code parameterisation `V = σ² I + Σ ŵ²_s · G_s`, so per-tree
`Var(u_s) = ŵ²_s · G_s` with **no extra σ² factor**; verified by the
self-consistency identity in `test/eval/edge_h2.py`, max residual 9e-15
on a Block 1.2 cell):

> E[α_e^s | y]  =  (ŵ²_s · ℓ_e / 2) · m_eᵀ · z
>                =  (ŵ²_s · ℓ_e / 2) · Σ_{i: c_i^{u(e)} > 0} c_i^{u(e)} · z_i
>                =  ŵ²_s · ℓ_e · S_z[u(e)],

where `S_z[u] = ½ · m_eᵀ · z` follows the matvec post-order accumulator
convention.

**Implementation note**: the Python helpers `streaming_combined_kinship`
(in `test/eval/predict.py`) and `hap_to_indiv_2d` together return *half*
the matvec-scale `G_s`, because the helper uses `0.5 · bl` at every
branch AND the `hap_to_indiv_2d` collapse contributes another `0.5`.
`edge_h2.py` accordingly doubles the Python-side kinship before
assembling `V` (and again when computing the self-consistency check's
expected per-tree contribution).  Document this here because future
extensions that reuse the Python helpers will need to do the same.

### Per-edge posterior contribution to predicted phenotype

The edge's contribution to the predicted phenotype vector is
> u_e^posterior  =  m_e · E[α_e^s | y]
> u_e^posterior[i]  =  c_i^{u(e)} · E[α_e^s | y].

Its squared norm (the "variance footprint" of edge e) is
> ‖u_e^posterior‖²  =  (E[α_e^s | y])² · Σ_i (c_i^{u(e)})².

(For c_i ∈ {0, 1, 2}: Σ_i c_i² equals the number of heterozygous
descendants plus 4× the number of homozygous descendants in subtree(u(e)).)

### Per-edge fraction of phenotypic variance (the ranking statistic)

> Δh²_e^posterior  :=  ‖u_e^posterior‖² / (n · Var(y)).

This is the headline statistic for the ranked-edge list.  Properties:

- **Data-driven**: depends on `z` through the descendant sum.
- **Differentiates edges with the same ŵ²_s · ℓ_e**: edges whose
  descendants happen to have aligned z values rank above edges with
  heterogeneous z.
- **Cheap**: descendant sums of z can be computed for ALL edges in a
  tree in a single O(n) pre-order pass (Appendix B); same shape as the
  `tree_cov_matvec` algorithm itself.
- **Self-consistent**: summed over all edges in a tree,
  `Σ_e u_e^posterior = u_s^posterior` (the per-tree posterior
  contribution).  Summed over all trees and edges, reproduces
  `Σ_s u_s^posterior = E[u | y]`.  Useful sanity check.

### Per-tree posterior contribution (Adam: "does this track ŵ²_s?")

The per-tree analog
> u_s^posterior  =  σ² · ŵ²_s · G_s · z,
> Δh²_s^posterior  =  ‖u_s^posterior‖² / (n · Var(y))

is cheap (matvec already implemented in sargent) but **largely tracks
ŵ²_s**: at REML convergence ŵ²_s itself is defined by a quantity of the
form `z^T G_s z / tr(G_s G⁻¹)`, so the rank order would essentially
reproduce ŵ²_s.  Distinction is only at the second-decimal level
(diffuse vs aligned z patterns).  **Not recommended as a separate
manuscript figure** — included here for completeness and as a sanity-
check sum target.

### Per-age-bin posterior contribution (sharper WHEN-axis, with SNP/ARG split)

Bin edges by their age (child-node age or midpoint).  For each bin b
compute the SNP, ARG, and total contributions **separately**:

> u_b^SNP[i]    =  Σ_{e ∈ bin b} c_i^{u(e)} · β_e
> u_b^ARG[i]   =  Σ_{e ∈ bin b} c_i^{u(e)} · E[α_e^s | y]
> u_b^total[i] =  Σ_{e ∈ bin b} c_i^{u(e)} · ( β_e + E[α_e^s | y] )
>
> Δh²_b^SNP   =  ‖u_b^SNP‖²  / (n · Var(y))
> Δh²_b^ARG   =  ‖u_b^ARG‖²  / (n · Var(y))
> Δh²_b^total =  ‖u_b^total‖² / (n · Var(y))

Each of these aggregates across many trees and many edges, so none
track ŵ²_s.  The per-age-bin readout is the posterior, per-edge
sharpening of Block 2.1 v2's per-tree-mean-age binning, now resolved
into the typed (SNP) and untyped (ARG-only) components.

**Manuscript role**: replaces the current `figs2/fig5_pBC.png` Panel A
(WHEN axis) with this sharper version.  The new panel can show
stacked bars (SNP bottom, ARG top) per age bin, with the height
giving total Δh² and the partition showing how much of the
contribution at each timescale is typed vs untyped.

**Hypothesis worth testing**: under standard demographic models,
*younger* branches should carry a larger fraction of their Δh² in the
typed-SNP component (recent mutations are well-tagged by common typed
SNPs), while *older* branches should be more ARG-dominated (deep
diversity, including the rare-variant fraction, is captured by the ARG
kinship and not the typed panel).  If this holds, it's a direct
quantitative validation of the "ARG captures untyped variance"
manuscript claim, at the resolution of evolutionary time.

### Per-edge fixed-effect (SNP) contribution

Under the infinite-sites assumption, every typed SNP `j` arose from a
single mutation at one position on one branch of one local tree.
`tskit` records this directly: each mutation has a `node` field that
identifies the child node of the branch on which it sits.

For SNP `j` arising on branch `e`, the genotype `x_{ij}` is exactly the
descendant count `c_i^{u(e)}`: individual `i` carries the variant iff
at least one of its haplotypes descends from `u(e)`, with copies
matching the descendant count.

So the BayesR fixed-effect contribution decomposes cleanly across
branches.  Define the **per-edge cumulative SNP effect**
> β_e  :=  Σ_{j : SNP j sits on edge e}  β̂_j
> (β̂_j from sargent's `.beta.tsv`; un-typed SNPs are excluded by construction)

Then the per-edge contribution to predicted phenotype via the fixed
effect is
> u_e^SNP[i]  =  c_i^{u(e)} · β_e,
> ‖u_e^SNP‖²  =  β_e² · Σ_i (c_i^{u(e)})²,
> Δh²_e^SNP   =  ‖u_e^SNP‖² / (n · Var(y)).

This uses **exactly the same descendant-count structure** as the ARG
contribution `Δh²_e^ARG`.  Per-edge, the two are on the same basis and
can be summed or contrasted.

**Identifiability check** (a key advantage over Block 2.2).  Multiple
typed SNPs sitting on the same branch are perfectly correlated by
construction (they're co-inherited by all descendants of the branch).
BayesR's mixture prior assigns the signal to ONE of them with the rest
≈ 0, but the **sum β_e is well-determined** — it's just the joint
posterior expectation of the cumulative effect at branch e.  Per-edge
SNP attribution is therefore clean, in contrast with the Block 2.2
ARG-per-tree collinearity issue.  Different edges *within the same tree*
have different descendant sets (different positions in the tree), so
BayesR has signal to discriminate.

### Per-edge total contribution: SNP + ARG

Since `Δh²_e^SNP` and `Δh²_e^ARG` use the same descendant-count
structure, the total is the squared norm of the combined contribution:

> u_e^total[i]  =  c_i^{u(e)} · ( β_e + E[α_e^s | y] )
> ‖u_e^total‖² =  ( β_e + E[α_e^s | y] )² · Σ_i (c_i^{u(e)})²
> Δh²_e^total  =  ‖u_e^total‖² / (n · Var(y)).

(NOT just `Δh²_e^SNP + Δh²_e^ARG`; the cross-term `2·β_e·E[α_e|y]` is
non-trivial.  Edges with aligned SNP and ARG contributions have
`Δh²_e^total > Δh²_e^SNP + Δh²_e^ARG`; edges with opposing contributions
have `Δh²_e^total < Δh²_e^SNP + Δh²_e^ARG`.  Both readouts are
manuscript-relevant.)

**Manuscript composition for per-edge ranked list**: for each top-K
edge report the triple `(SNP, ARG, total)` and the cross-term sign.
A stacked-bar visualization (SNP component on bottom, ARG above)
shows how typed and untyped contributions combine at the
branch level.  No SNP-based or summary-stats method can produce this.

### Per-population-bin posterior contribution (cross-ancestry decomposition with SNP/ARG split)

In structured demographies (Block 1.1 v3 onward), bin edges by the
population identity of their child node:
- **Pre-split** edges: shared across all populations.
- **Post-split, pop A specific**: in the ancestral A lineage post-split.
- **Post-split, pop B specific**: in the ancestral B lineage post-split.
- (For admixed test individuals: also "introgressed from archaic.")

Same formula structure as per-age-bin:
> Δh²_{pop b}^SNP   =  ‖Σ_{e ∈ pop bin b} c_·^{u(e)} · β_e‖²    / (n · Var(y))
> Δh²_{pop b}^ARG  =  ‖Σ_{e ∈ pop bin b} c_·^{u(e)} · E[α_e|y]‖² / (n · Var(y))
> Δh²_{pop b}^total =  ‖Σ_{e ∈ pop bin b} c_·^{u(e)} · (β_e + E[α_e|y])‖² / (n · Var(y))

The fraction of total h² in each bin is a **direct decomposition of
phenotypic variance by evolutionary origin**, separated by whether the
typed panel captured it (SNP component) or not (ARG component).  Backs
the cross-ancestry portability story quantitatively.

**Hypotheses to test under Block 1.1 v3**:
- High pop-A-specific Δh² fraction → low portability A → B predicted.
- High pre-split fraction → high portability across the F_ST ladder.
- For admixed-test individuals (AmericanAdmixture ASW/MXL/PEL):
  Δh²-contribution split between AFR-ancestry, EUR-ancestry, NAT-
  ancestry branches → matches each individual's known admixture
  proportions?  Direct AoU pitch.
- **NEW** with the SNP/ARG split: at high F_ST cross-pop, does the SNP
  contribution from post-split pop-A-specific branches port to pop B?
  Mechanistic prediction: NO (typed-SNP effects don't transfer
  cross-pop, per Block 1.1 v2's β-liability finding), while the ARG
  contribution from the same branches *does* port via shared kinship
  structure.  Direct mechanistic validation of the cross-ancestry
  result.

## Algorithm and cost

### Single-pass per-tree descendant sums

For each tree s in {top-K by ŵ²_s}, compute in one O(n) post-order pass:
> S_z[u]      :=  Σ_{leaves ℓ ∈ subtree(u)} 0.5 · z[indiv(ℓ)]
> S_z2[u]     :=  Σ_{leaves ℓ ∈ subtree(u)} (c_indiv(ℓ)^u)² · 0.25  *(or alternative*
                  Σ_i (c_i^u)² accumulator at node u as integers)*

(The factor 0.5 mirrors the sargent matvec; `S_z[u]` is exactly the
analog of `S[u]` in `tree_cov_matvec` with `z` substituted for `in`.)

Then for each edge `e` with child node `u(e)`:
> E[α_e | y]  =  (σ² · ŵ²_s · ℓ_e) · S_z[u(e)]    (factor of 1/2 absorbed)
> ‖u_e^posterior‖²  =  (E[α_e | y])² · S_cnt2[u(e)]   *(or equivalent)*

Total cost per tree: O(n + n_nodes_in_tree).  Across top-K trees:
O(K · n).  For typical K = top-100 trees by ŵ²_s and n = 1000:
~10⁵ operations.  Negligible vs sargent's main fit.

### Implementation plan

`test/eval/edge_h2.py` — single-file Python helper:

```python
def edge_h2(ts, ws2, beta, beta_snp_to_mutation, z, sigma2, var_y,
            top_k_trees=None, age_bins=None, pop_assignments=None,
            top_k_edges=None):
    """
    Compute per-edge posterior contributions to phenotypic variance,
    separated into SNP and ARG components.

    Args:
        ts: tskit tree sequence
        ws2: per-tree variance weights (sargent output, length S)
        beta: BayesR posterior means (sargent <outroot>.beta.tsv), length p
        beta_snp_to_mutation: array of length p mapping SNP j -> mutation id in ts
            (built from case.bim + ts.mutations cross-reference)
        z: BLUP residual V^{-1} (y - X beta), length n
        sigma2: sargent fitted sigma2
        var_y: phenotypic variance
        top_k_trees: int — restrict to trees with the K largest ws2 (None = all)
        age_bins: list of (lo, hi) tuples for age-bin aggregation
        pop_assignments: dict {node_id: pop_label} for per-population aggregation
        top_k_edges: int — keep only the K largest-Δh²_total edges in the returned list

    Returns:
        ranked_edges: pandas DataFrame with columns
            tree_idx, edge_child_node, length, child_age, n_descendants,
            beta_e, E_alpha,
            dh2_snp, dh2_arg, dh2_total,
            position_lo, position_hi, pop_label
        age_summary: dict {bin: {snp, arg, total}}
        pop_summary: dict {pop: {snp, arg, total}}
    """
```

**Key building block**: a per-edge `beta_e` accumulator.  For each
tree s, iterate `ts.mutations()` restricted to mutations whose position
falls in tree s's interval; for each mutation `m`, find the SNP index
in the BED panel (via `m.site` and `case.bim` cross-reference), look up
`beta[j]`, and accumulate `beta_e[m.node] += beta[j]`.  Untyped
mutations are skipped (they have no BayesR estimate).

### Validation

1. **Self-consistency** (both ARG and SNP sides):
   - **ARG side**: summed over all edges in a tree, per-edge `u_e^ARG`
     must equal the per-tree `u_s^ARG = σ² ŵ²_s G_s z`.  Test on a small
     simulation.  Summed over all trees, `Σ u_s^ARG = E[u | y]`.
   - **SNP side**: summed over all edges, `Σ_e u_e^SNP[i] = Σ_j x_{ij}·β̂_j`
     (the BayesR fixed-effect prediction).  Test against
     `predict_bayesR` in `predict.py`.
   - **Total side**: `Σ_e u_e^total[i]` = sargent's full prediction `X β̂ + E[u|y]`.
2. **Truth-side recovery on a known-architecture simulation**:
   - Block 1.2's sparse-3-block architecture has known true active blocks.
     The top-K edges by Δh²_total should concentrate in those blocks.
     Quantify with active-block AUC at the edge level.
   - Pure-SNP simulation (h²_arg=0): Δh²_e^SNP should dominate; Δh²_e^ARG
     should be near zero at convergence.
   - Pure-ARG simulation (h²_snp=0): Δh²_e^ARG should dominate;
     Δh²_e^SNP near zero (β̂≈0 from BayesR).
3. **Cross-check vs Block 2.1 v2 age-binning**:
   - Aggregating per-edge to age bins should produce numbers comparable to
     Block 2.1 v2's per-tree-mean-age binning at the same age cuts.
     Discrepancy = the sharpening (per-edge vs per-tree-mean-age) plus
     the SNP component (new).
4. **Diffuse-control**:
   - On a panmictic, uniform-ŵ² simulation with all-typed SNPs, no age
     bin or population bin should dominate; SNP and ARG components
     should be roughly co-distributed across bins.

## Manuscript composition

Three figures or panels follow from this calculation, each carrying
the SNP / ARG / total decomposition:

| readout | what it shows | manuscript role |
|---|---|---|
| Ranked-edge list (top-K) | top edges by posterior Δh²_total, annotated time/pop/length/position; per-edge stacked SNP+ARG bars | Fig 6 or supplementary; SARGENT-distinctive primary contribution |
| Per-age-bin Δh² (SNP / ARG / total) | fraction of h² by coalescent time, partitioned typed vs untyped | Fig 5 Panel A (replaces / sharpens existing WHEN panel from Block 2.1 v2) |
| Per-population-bin Δh² (SNP / ARG / total) | fraction of h² by evolutionary origin × typed/untyped | Fig 3 Panel C (or sup); backs cross-ancestry portability story with a mechanistic decomposition |

All three are produced by the same single calculation (descendant
sums of `z` and `β_e` per edge); only the aggregation step differs.
Total marginal cost beyond sargent fit is ~1 day of implementation +
validation.

## Open questions and possible extensions

1. **Posterior variance per edge (for credible intervals)**.
   `Var(α_e | y) = σ² ŵ²_s ℓ_e / 2 · (1 − σ² ŵ²_s ℓ_e / 2 · m_eᵀ V⁻¹ m_e)`.
   The `m_eᵀ V⁻¹ m_e` quadratic form is O(n²) per edge with the .gall
   eigendecomposition; for the top-K edges in the ranked list this is
   tractable.  Useful for honest reporting (CI bands on Fig 6).
2. **Leave-one-tree-out per-edge contribution.**  Replace `z` with
   `z_{−s} := V_{−s}⁻¹ (y − Xβ − E[u_{−s} | y])` for each tree s.
   Sherman-Morrison-style rank-1 update keeps the cost manageable.
   Yields a sharper per-tree attribution (each edge is credited only for
   what *its* tree contributes uniquely, not double-counted across trees).
   Worth implementing for the top-K edges; may not change the ranking
   much in practice but is the cleaner statistic.
3. **Belief propagation**.  Not needed for posterior means (the
   descendant-sum is already a single tree pass).  May be useful if we
   want posterior covariances between pairs of edges — but this isn't
   currently a manuscript need.
4. **Edge "position" annotation**.  For a tree s spanning genomic
   coordinates [a_s, b_s], every edge in tree s shares the same
   position interval.  For position-stratified analyses, edge entries
   can be assigned to position bins by their tree's interval.  This
   would let us cross-tab "position × time" — does a particular genomic
   region's h² contribution come from recent or ancient ancestry?
5. **`top_k_edges` filtering vs full output**.  At biobank scale, the
   full edge list can be hundreds of millions of entries.  For
   manuscript figures we only need the top-K by Δh²; for aggregation
   (age bins, populations) we need the full sum but not the per-edge
   detail.  Implementation streams to both outputs in one pass.

## Implementation timeline

- Day 1: write `edge_h2.py`; cross-check posterior-mean calculation
  against sargent's BLUP prediction (the sum-over-edges identity).
- Day 2: apply to a representative Block 1.2 v1 cell; produce the
  per-edge ranked list + per-age-bin aggregation; cross-check vs
  Block 2.1 v2 binning.
- Day 3 (gated on Block 1.1 v3 landing): apply to a structured-
  demography cell; produce per-population-bin aggregation; draft
  Fig 5 Panel A replacement + Fig 6 supplementary.

Posterior-variance / leave-one-tree-out extensions: separate ~2-day
follow-ups, only if needed for review.

---

## Appendix A — per-edge decomposition of tree_cov_matvec

Reading `src/lib/trees.c::tree_cov_matvec` (the production
implementation):

```
Post-order: S[u] = Σ_{leaves ℓ ∈ subtree(u)} 0.5 · in[indiv(ℓ)]
           = 0.5 · Σ_i c_i^u · in[i]            (c_i^u = #haps of i in subtree(u))
Pre-order:  acc[u] = acc[parent(u)] + ℓ(u-up-edge) · S[u]
Read-off:   out[i] = acc[hap_a(i)] + acc[hap_b(i)]
```

Each edge `e` (child node u, length ℓ_e) contributes
`acc[u]` += `ℓ_e · S[u]` and that contribution propagates downward.
Tracking it forward to `out`:

```
contribution_e[i] = ℓ_e · S[u] · c_i^u
                  = ℓ_e · 0.5 · (Σ_j c_j^u · in[j]) · c_i^u
                  = (ℓ_e / 2) · c_i^u · (m_eᵀ · in),       m_e[j] := c_j^u(e)
```

Summing across edges:
```
G_s · in  =  Σ_e (ℓ_e / 2) · m_e · (m_eᵀ · in)
G_s        =  (1/2) · M_s · diag(ℓ) · M_sᵀ                □
```

## Appendix B — single-pass per-tree descendant sums

To get `S_z[u]` (the descendant sum of `z` at each node u) for all
nodes in a tree in O(n_nodes):

```
post-order over nodes u:
    if u is a leaf (sample haplotype h):
        S_z[u] = 0.5 · z[indiv(h)]
    else:
        S_z[u] = Σ_{children c of u} S_z[c]
```

(Identical structure to the `S[u]` computation in `tree_cov_matvec`;
substitute `z` for `in`.)

Per-edge posterior means then read out directly:
```
for each edge e with child u, length ℓ_e:
    E[α_e | y] = σ² · ŵ²_s · ℓ_e · S_z[u]      (factor of 1/2 absorbed
                                                 into the prior scaling)
```

Per-edge squared-norm requires `S_cnt2[u] := Σ_i (c_i^u)²`, which is
also a post-order accumulator with leaf init `0.25` (per haplotype)
and internal accumulation
`S_cnt2[u] = Σ_{children c} S_cnt2[c] + (cross-terms within u)`.
The cross-terms collapse cleanly for diploid: see implementation.
