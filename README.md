# Pool design and genotype-based demultiplexing

Does the genetic composition of a multiplexed single-cell pool (which donors you
choose to pool together) affect how well genotype-based demultiplexing (demuxlet)
can tell them apart? This repo holds the whole-genome ambisim simulation arm of
that question: 8-donor pools built from 1000 Genomes individuals via several
donor-selection strategies, run through a demuxlet-only pipeline (GEX and ATAC
pileup + call, no cell-typing), scored against ambisim's own ground truth.

Companion real-data work (the original motivation — ~40 real 10x multiome
pools) lives in the parent project, not duplicated here; `notebooks/vcf_metrics/`
and `notebooks/real_data/` are copied in only because their results are directly
load-bearing for interpreting the simulation.

**For pipeline bugs, methodology caveats, and full derivations behind the numbers
below, see `NOTES.md`. For the session-by-session history of how this project got
here, see `JOURNAL.md`.**

## Pool construction

All pools draw n=8 donors from a "universe" of unrelated 1000 Genomes
individuals — either single-ancestry (`EUR_only`, `AFR_only`, `EAS_only`,
`SAS_only`, `AMR_only`) or multi-ancestry (`EUR_AFR`, `EUR_EAS`, `all_3_major`,
`all_5`). See `scripts/ambisim/lib/pool_nomination.py`. Every non-random
strategy runs the same greedy loop — start with one random donor, repeatedly
add whichever remaining candidate optimizes a pairwise objective against the
donors already chosen — with the objective swapped out:

| strategy | objective |
|---|---|
| `random` | no optimization; baseline |
| `greedy_maxmin` | maximize the *minimum* pairwise genetic distance in the pool |
| `greedy_maxmean` | maximize the *mean* pairwise genetic distance |
| `greedy_maxkl` | like `greedy_maxmin` but using KL divergence instead of Euclidean distance (single-ancestry only) — **buggy, treat as unreliable, see `NOTES.md`** |
| `ancestry_balanced` | proportional sampling across sub-populations, no distance optimization (multi-ancestry only) |
| `adversarial_mindist` | mirror of `greedy_maxmin` — *minimizes* the pool's spread (worst-case floor) |
| `adversarial_family` / `adversarial_family_mixed` | a real 1000G pedigree (~50% IBD), padded with unrelated fill — a non-synthetic worst case |

Each cell (universe × strategy) has 3 replicates; the full design is 132 pools
(`txt/pool_experiments.txt`). **As of 2026-09-14, 88/132 pools have completed
demuxlet on both modalities** (full core grid rep1-3 + a growing adversarial
subset — 15 of 132 pools so far).

## Results (n=88 pools: 73 ordinary + 15 adversarial)

Every correlation here is reported for the 73 ordinary (unrelated-donor:
`random`/`greedy_max*`/`ancestry_balanced`) pools with the 15 adversarial pools
excluded, because the adversarial pools are high-leverage outliers that inflate
whole-sample correlations without a real trend behind them — see `NOTES.md`
for why that check matters and the full old-vs-corrected numbers.

**Singlet accuracy: 0.799–0.860 (GEX), 0.903–0.926 (ATAC) across all 88 pools —
narrow, and not explained by design.** Among the 73 ordinary pools, `min_dist`
explains ~0% of that variance (r≈0, n.s., both modalities). The one accuracy
difference that *is* real: pools built from actual close relative pairs
(`adversarial_family`/`_mixed`, 6 pools, real 1000G pedigrees) score measurably
lower than everything else (Mann-Whitney p<0.001 both modalities, ~1-2 accuracy
points lower). `adversarial_mindist` (unrelated donors, worst-case-selected)
shows a smaller version of the same gap (p=0.002-0.005, mean diff as small as
+0.004 ATAC) despite `min_dist` values overlapping ordinary pools — thin
evidence, not yet a settled result.

**Classifier margin (LL-gap) is the one design signal that holds up.**
`min_dist` vs. mean LL-gap, ordinary pools only: r=0.31 (GEX), r=0.31 (ATAC),
both p<0.01 — smaller than the naive whole-sample r=0.66-0.73 but real, and it
replicates the original n=32 finding (r=0.51-0.56 there) in direction and
significance. `mean_dist` still outpredicts `min_dist` at the bottleneck-donor
level (r≈0.66 vs. r≈0.30-0.40) among ordinary pools — same as the n=32
snapshot; the earlier report of a reversal was a leverage artifact (`NOTES.md`).

**Design does not rescue high-ambient ("bad") droplets — the opposite.**
`min_dist` vs. LL-gap, split by ambient RNA bin (ordinary pools only): strongest
at low ambient (0-10%: r=0.31 both modalities, p<0.01), weakening monotonically
to non-significant by ≥40% ambient. Design's margin benefit is a cushion in
already-decent conditions, not a rescue once a droplet is already in trouble.
Per-pool sequencing depth was checked as a possible confound alongside ambient
fraction and ruled out — it's fixed by simulation design (median identical
across all 88 pools).

**Restricting to a single ancestry (as real-data comparisons do) doesn't
sharpen the signal — mixing ancestries does.** Per-single-ancestry-universe
correlations (n=8-10 pools each) are noisy and directionally inconsistent.
Pooling all single-ancestry universes (n=44) gives a *weaker* margin
correlation (r=0.20, n.s.) than the full ordinary-pool set; pooling the 4
multi-ancestry universes instead (n=29) gives the *strongest* of any subset
(r=0.56-0.59, p<0.005) — because multi-ancestry pools have the widest
within-universe `min_dist` range, giving the metric more room to vary at all.
Restricting to one ancestry mostly just removes power (details in `NOTES.md`).

## Commentary — does pool design matter?

**Narrower than "yes" or "no."** Curating unrelated donors by genetic distance
— any strategy tested here, vs. `random` — has not been shown to move raw
demuxlet accuracy, at n=88 with real statistical power, in any ambient
condition, in any ancestry subset. It does move classifier margin, reliably
and reproducibly, but that margin advantage does not extend to the droplets
that would need it most (high ambient). The one thing that clearly does move
raw accuracy is avoiding actual close relatives — which argues for a basic
relatedness/pedigree check before pooling, not for a sophisticated genetic-
distance optimization algorithm over otherwise-unrelated donors.

That is a real, useful, if less exciting than hoped, result: it says the
field's implicit assumption that clever donor curation improves demultiplexing
accuracy is not supported here, while giving a precise, quantifiable account
of what a VCF *can* still tell you pre-experiment (fragility/margin) and what
it apparently can't (a guaranteed accuracy win) for anything short of literal
relatedness. Two things could still change this before calling it final:

- 44/132 pools aren't simulated yet (`NEXT_STEPS.md`), including most of the
  adversarial floor — more `adversarial_family`/`_mixed` replicates would say
  whether the relatedness effect (n=6 pools currently) is as robust as its
  p-values suggest.
- The downsampling/degradation robustness check (`md/analyses.md` Tier 3) —
  does the margin advantage ever convert to an accuracy advantage under
  stress, for *ordinary* (non-relative) pools — hasn't been run in simulation
  yet, and is probably the highest-leverage remaining experiment for settling
  whether this is a real accuracy story or a margin-only one.

## Next steps

See `NEXT_STEPS.md` for the concrete pool list and qsub notes.
