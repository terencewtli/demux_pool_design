# RESULTS — what we know, and what shows it

Current validated findings only. No history: how each result was arrived at (including claims
that were retracted along the way) is in `JOURNAL.md` / `JOURNAL.archive.md`; what has finished
running is in `PROGRESS.md`; pipeline bugs are in `docs/pipeline_issues.md`.

**Basis for every number below:** the chr22-matrix `ambisim/` grid, rescored at **n=128 pools**
(105 ordinary + 23 adversarial) on 2026-09-19 (`notebooks/ambisim/01e-01h`, `results/*.csv`).
This supersedes the n=88 run: 40 of the 128 were demuxed after that snapshot and act as a
replication set. The genome-wide `ambisim_final/` grid will supersede these numbers when it
lands (see `PROGRESS.md`). Last updated 2026-09-19.

**Standing rule:** every pool-level correlation is reported with the 23 adversarial pools
excluded, and with Spearman alongside Pearson. The adversarial pools sit at the low end of
`min_dist` and are high-leverage; two whole-sample claims have already been retracted because of
them (`JOURNAL.md` status board, `docs/CONVENTIONS.md`).

---

## 0. Read this first: the two margins are not the same column

`DIFF.LLK.BEST.NEXT` is what this project called "the LL-gap" from the start. **It is not the
donor-discrimination margin.** For a droplet demuxlet calls SNG, `NEXT.GUESS` is a *doublet*
model (6,705 of 6,706 in a spot-checked pool), so that column is numerically identical to
`DIFF.LLK.SNG.DBL` — it is the **singlet-vs-doublet** margin.

The **donor-vs-donor** margin is `SNG.BEST.LLK − SNG.NEXT.LLK`, carried as `sng_gap` from 01e on.
Medians across 128 pools: `ll_gap` 6.7 / `sng_gap` 53.9 (GEX), 19.4 / 135.6 (ATAC); ρ=0.81–0.88.

Sections 2 and 4 below were rewritten on 2026-09-19 because they had been reported on the wrong
column. Anything citing "LL-gap" from before that date is a statement about doublet calling.

## 1. Singlet accuracy is mostly not a genetics measurement

Singlet accuracy is **0.799–0.862 (GEX), 0.903–0.926 (ATAC)** across 128 pools. Among the 105
ordinary pools, `min_dist` explains ~0% of it (GEX r=−0.087, ATAC r=+0.108, both n.s.), and the
same null holds on the 40 held-out pools (r=−0.225 / +0.149, n.s.).

The reason is that the metric multiplies two decisions with wildly different error rates:

```
accuracy = P(called SNG | true singlet) x P(right donor | called SNG)
gex:  0.8410 = 0.8419 x 0.99898
atac: 0.9170 = 0.9170 x 0.99995
```

Across ordinary pools `sd(accuracy)=0.0107` and `sd(P called SNG)=0.0107`, but
`sd(P donor|SNG)=0.0004`. **All** the pool-to-pool spread is the singlet/doublet call. The true
singlets demuxlet loses are ambient-heavy (median ambient 0.35 vs 0.16 for the ones it keeps),
and ambient was fixed at 20% by simulation design and never swept.

So "singlet accuracy" is dominated by an axis donor genetics barely moves, and buries a
donor-identity error ~200x smaller. **Do not use it as the outcome for a design question.**

## 2. Donor distance does predict the donor margin — and `mean_dist` beats `min_dist`

On `sng_gap`, ordinary pools, n=105 (Pearson / Spearman):

| metric | GEX | ATAC |
|---|---|---|
| `min_dist` | +0.341 / +0.349 | +0.363 / +0.352 |
| **`mean_dist`** | **+0.706 / +0.717** | **+0.728 / +0.731** |
| `kl_min` | +0.543 / +0.528 | +0.601 / +0.572 |

`mean_dist` is roughly twice `min_dist`, and the relationship is identical in every coverage bin
(ρ=+0.69 to +0.74 from <1,000 to >8,000 informative SNPs) — a uniform multiplicative scaling, not
something concentrated in weak droplets. **Caveat (§6): `greedy_maxmin` and `greedy_maxmean`
share 5–6/8 donors, so "mean beats min" needs `ambisim_new/` to be clean.**

## 3. The margin converts to accuracy only below a threshold, and pools sit above it

P(right donor | called SNG) against `sng_gap`, pooled over all 128 pools:

| sng_gap | GEX | ATAC |
|---|---|---|
| 2–3 | 0.817 | 0.829 |
| 3–5 | 0.929 | 0.929 |
| 5–7.5 | 0.980 | 0.975 |
| 7.5–10 | 0.996 | 0.999 |
| 10–15 | 0.999 | 0.998 |
| **≥20** | **1.000** | **1.000** |

Saturation at `sng_gap ≈ 20` in both modalities. Pool means are 58.5–84.2 (GEX) — about 3x above
it. `sng_gap` is a **sufficient statistic**: within an `sng_gap` bin, P(right donor) agrees to
the 4th decimal across `random` / `greedy_maxmin` / `greedy_maxmean` / `ancestry_balanced` /
`adversarial_mindist`. Pools differ in *where they sit on one shared curve*, not in the curve.

What design moves is the fraction of droplets below threshold: `frac(sng_gap<20)` is 0.112
(ordinary) vs 0.165 (adversarial) in GEX, tracking `mean_dist` at ρ=−0.693 and `min_dist` at
ρ=−0.304.

In cells: **4.8 wrong donors per 10,000** in the 10 highest-`mean_dist` ordinary pools vs **12.9**
in the 10 lowest. ~8 cells per 10,000 is the whole prize for distance optimization.

## 4. The design effect does *not* fade with ambient RNA

`mean_dist` vs `sng_gap` by ambient bin, ordinary pools: ρ=+0.71 / +0.72 / +0.72 / +0.71 (GEX,
0–10% through 30–40%) and +0.73 / +0.74 / +0.72 / +0.73 (ATAC). Flat.

The previous version of this section said the benefit fades to n.s. by ≥40% ambient. That was
measured on `ll_gap` (the doublet margin, §0) and is true of *that* quantity. It is not true of
the donor margin. Sequencing depth remains fixed by design (~20,021 RNA / ~20,240 ATAC reads,
identical median across pools).

## 5. Pool composition shifts the singlet/doublet operating point (new, 2026-09-19)

The simulation makes 900 true doublets per pool alongside 8,100 true singlets. Those had never
been scored. They should have been:

- **P(true doublet correctly called DBL) vs `mean_dist`**, ordinary pools: **ρ=+0.600 (GEX,
  p<10⁻⁴)**, ρ=+0.437 (ATAC, p<10⁻⁴). Stronger than anything in §2–§3 for singlets.
- **Within-pool control** — rank each droplet's own two donors by distance inside its own pool
  (same depth, same ambient, same pool): doublets from the closest 20% of pairs are caught
  0.9471 of the time vs 0.9645 for the farthest 20% (GEX), monotone across all five quintiles.

And the two sides of the boundary move in **opposite** directions with distance:

| | GEX | ATAC |
|---|---|---|
| `mean_dist` ~ P(SNG \| true singlet) | ρ=−0.148 (n.s.) | ρ=+0.450 (p<10⁻⁴) |
| `mean_dist` ~ P(DBL \| true doublet) | ρ=+0.600 (p<10⁻⁴) | ρ=+0.437 (p<10⁻⁴) |

In GEX the two partly cancel, which is why net singlet accuracy (§1) looks flat. **Pool
composition changes where demuxlet draws the singlet/doublet line, not how often it is right.**

Mechanistically this is close to expected — a mixture of two similar donors should resemble one
donor. What is new is the sign relative to the singlet metric, which is what dissolves §1's null.

## 6. Real relatives: a pool-wide cost, not a localized one

Ancestry-matched (the family arms are EUR, so the all-universe ordinary mean is the wrong
control — it inflated this to 3.4x in an earlier pass), GEX:

| group | pools | min_dist | sng_gap | wrong donor /10k |
|---|---|---|---|---|
| ordinary EUR_only | 12 | 106.9 | 66.7 | 11.2 |
| **family EUR_only** | 3 | **74.9** | **45.9** | **32.8** |
| ordinary EUR_AFR | 12 | 107.0 | 78.0 | 7.0 |
| **family EUR_AFR** | 3 | **74.9** | **51.2** | **26.0** |

A **2.9x** relatedness penalty. Family pools sit *inside* the ordinary `mean_dist` range
(102–108) but far below its `min_dist` floor (74.9 vs 90.7–120.8): relatedness acts through the
closest pair, not the average.

**The deficit is localized to the related donors** (GEX). Relatedness is derived from the 1000G
pedigree per pool, not assumed: each of the 3 reps is a different CEPH pedigree (1420, 1463,
1362), each a couple plus all four of their parents, so 6/8 donors have a first-degree relative
in-pool and 2 are unrelated fill.

| donors in family pools | sng_gap | wrong donor /10k |
|---|---|---|
| with a first-degree relative in-pool (6/8) | 38.9 | **36.6** |
| unrelated fill, same pools (2/8) | 77.3 | **8.0** |
| ordinary EUR_only, for reference | 66.7 | 11.2 |

A **4.6x** within-pool gap, and **62.7% of wrong calls name a first-degree relative of the true
donor** (57.1% in ATAC). The unrelated fill donors are not harmed at all — 8.0/10k, if anything
better than ordinary EUR pools. This is a pairwise confusion between relatives, not a pool-wide
loss of diversity.

**The relatedness tested is parent-offspring only** (~0.5 IBD), 4 such pairs per pool. No sibs,
no second-degree, no gradient. That is the single biggest gap in this result.

**Full cost, per 10,000 loaded droplets (GEX):**

| | singlets called DBL/AMB | wrong donor | doublets missed | total |
|---|---|---|---|---|
| ordinary EUR_only | 1,421.0 | 8.5 | 47.5 | 1,477 |
| family EUR_only | 1,539.8 | 24.4 | 60.4 | 1,625 |

**~148 extra bad droplets per 10,000**, of which only ~16 are donor misassignment. A relatedness
check before pooling is cheap and worth doing; the cost of missing one is not catastrophic.

`adversarial_mindist` (17 pools, unrelated, worst-case-selected) shows a weaker version:
`sng_gap` 65.9 vs 72.3, 12.0 wrong donors/10k vs 8.8, 1,530 total bad vs 1,477.

## 7. Mixing ancestries gives distance room to matter

- Per single-ancestry universe (n=8–10 pools each): noisy, sign-inconsistent — underpowered.
- Multi-ancestry universes span a much wider `min_dist` range (16.5–24.8 vs 9.4–16.1), which is
  where the §2 correlations get their leverage.

Restricting to one ancestry (as real-data comparisons do) removes power rather than sharpening
signal.

## 8. Min vs mean distance cannot be separated with the current strategies

`greedy_maxmin` and `greedy_maxmean` pick a median 5–6/8 identical donors per universe/rep:

| source | min_dist vs. mean_dist r |
|---|---|
| real 10x pools (n=63) | 0.14 (n.s.) |
| 20,000 random 8-donor draws, `EUR_only` | 0.536 |
| greedy pools, chr22 matrix | 0.69–0.84 |
| `greedy_maxmin`, genome-wide matrix (n=45) | 0.992 |

So every maxmin-vs-maxmean comparison is not a test of worst-case vs average-case design. It does
**not** affect §1, §3 or §5 (which are continuous regressions or within-pool contrasts), but it
does gate §2's "mean beats min". Mechanism and the fix: `docs/minmean_entanglement.md`. The
pilot (`ambisim_new/`, 4 `EUR_only` pools) is demuxed but **not yet analyzed**.

---

## Overall read

Curating unrelated donors by genetic distance does not meaningfully move demultiplexing outcome.
The effect is real and well-powered on the donor margin (`mean_dist` ρ≈0.71, replicated, holds at
every depth and ambient level), but the margin is ~3x above the threshold where it converts to
accuracy, so the payoff is **~8 cells per 10,000**. Avoiding close relatives is worth ~148 droplets per
10,000 pool-wide — but the cost is concentrated: for the affected donors themselves it is 36.6 vs
8.0 wrong calls per 10,000, a 4.6x hit that their pool-mates do not share. That is a screening
decision, not an optimization one, and only parent-offspring (~0.5 IBD) has been tested.

The finding worth building on is §5: pool composition shifts demuxlet's singlet/doublet operating
point, with closer pools systematically under-calling doublets (ρ=+0.60, confirmed within-pool).
That reframes the null instead of merely reporting it, and it makes a falsifiable prediction for
the real data that sequencing depth cannot explain — see `PROGRESS.md` → Next steps.

## Where the data lives

| what | path |
|---|---|
| per-droplet scores, n=128, with coverage + call type | `csv/ambisim/droplet_scores_cov.csv` (422 MB, not mirrored) |
| per-pool accuracy, both margins, all design metrics | `results/pool_margins_n128.csv`, `results/pool_summary_n128.csv` |
| donor-margin threshold curve | `results/sng_gap_accuracy_curve.csv` |
| doublet detection vs within-pool pair distance | `results/doublet_pair_distance.csv` |
| ambient-stratified donor margin | `results/ambient_stratified_sng_gap.csv` |
| relatedness, ancestry-matched | `results/family_ancestry_matched.csv` |
| full cell accounting per 10,000 droplets | `results/cell_accounting_per10k.csv` |
| coverage-matched design effect | `results/coverage_matched_design_allmetrics.csv` |
| nominated donor lists (chr22 matrix) | `results/nominated_pools_n8.tsv` |
