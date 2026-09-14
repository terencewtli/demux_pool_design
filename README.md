# Pool design and genotype-based demultiplexing

Does the genetic composition of a multiplexed single-cell pool (which donors you
choose to pool together) affect how well genotype-based demultiplexing (demuxlet)
can tell them apart? This repo holds the whole-genome ambisim simulation arm of
that question: 8-donor pools built from 1000 Genomes individuals via several
donor-selection strategies, run through a demuxlet-only pipeline (GEX and ATAC
pileup + call, no cell-typing), scored against ambisim's own ground truth.

Companion real-data work (the original motivation — ~40 real 10x multiome
pools) lives in the parent project, not duplicated here; `notebooks/vcf_metrics/`
and `notebooks/real_data/` below are copied in only because their results are
directly load-bearing for interpreting the simulation.

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
| `greedy_maxkl` | like `greedy_maxmin` but using KL divergence between genotype distributions instead of Euclidean distance (single-ancestry universes only) |
| `ancestry_balanced` | proportional sampling across sub-populations, no distance optimization (multi-ancestry universes only) |
| `adversarial_mindist` | mirror of `greedy_maxmin` — *minimizes* the pool's spread (worst-case floor) |
| `adversarial_family` / `adversarial_family_mixed` | a real 1000G pedigree (~50% IBD), padded with unrelated fill — a non-synthetic worst case |

Each cell (universe × strategy) has 3 replicates; the full design is 132 pools
(`txt/pool_experiments.txt`). **As of this snapshot (2026-09-14), 88/132 pools
have completed demuxlet on both modalities** — the full core grid (rep1-3) plus
a growing set of adversarial pools across several universes. This is the first
snapshot with adversarial pools in the data, which changes several conclusions
below relative to the n=32, rep1-only snapshot this section used to describe.

## `greedy_maxkl` is buggy/weird — treat all its results as unreliable for now

**Status: partially fixed, still broken, deprioritized rather than resolved.**
Do not use `greedy_maxkl` pools (including the regenerated ones described
below) as a valid "good design" data point until this is actually resolved.

`nominate_greedy_maxkl`'s inner KL computation added its numerical-stability
epsilon *before* taking the complement probability (`1 - p`), so any
homozygous-alt site (dosage=2 → p=1.0) produced `1 - p < 0` and
`log(negative) = NaN`, silently corrupting the greedy argmax every time it hit
a common homozygous site. This is why the 3 `greedy_maxkl` pools originally
simulated scored at or below `random` on their own target metric (`kl_min`) —
the nomination loop wasn't actually optimizing anything. Fixed in
`pool_nomination.py` (confirmed: no more NaN/divide-by-zero warnings post-fix).

**But re-running the (fixed) nomination across all 5 single-ancestry
universes shows the bug wasn't the whole story.** Re-nominated canonical
(best-of-5-seed) `kl_min` values — AFR_only 1.738, AMR_only 1.702, EAS_only
1.498, EUR_only 1.639, SAS_only 1.799 — are **still below that universe's
`random` baseline in every single case** (random: 1.851, 1.725, 1.675, 1.773,
1.880 respectively), just by a smaller margin than the pre-fix version. So
`greedy_maxkl`, even after the NaN fix, still doesn't demonstrably achieve
what its name claims (a pool that maximizes minimum pairwise KL divergence).
The search-time objective (`kl_to_pool`) and the scoring function
(`pool_kl_divergence`) were checked and now use consistent epsilon-handling,
so this isn't simply a scoring/search mismatch — root cause not identified.
Possibly a genuine property of greedy max-min search being more prone to
local-optimum traps for a KL-shaped objective than for Euclidean distance
(`greedy_maxmin`, using the same greedy structure, reliably beats `random` by
a wide margin) — or a subtler remaining bug not yet isolated.

**Decision: deprioritized, not resolved.** Given `kl_min`/`min_dist` are
r=0.999 collinear in this design space anyway, `greedy_maxmin` already serves
the "good design by genetic separation" comparator role, and the adversarial
pools are a far higher-value use of remaining compute for the paper's actual
question. Donor lists were regenerated with the fixed code
(`txt/donors/*greedy_maxkl*rep1.txt`, backups of the pre-fix
`nominated_pools_n8.tsv`/`simulation_pools_n8.tsv` kept alongside), but **the
3 already-simulated `greedy_maxkl` pools (AMR_only, EUR_only, SAS_only, rep1)
have deliberately not been deleted/re-simulated** — there was no point
spending ~20h/pool of compute regenerating a strategy that still isn't
working as intended. Treat every `greedy_maxkl` result in this repo
(pre- or post-fix) as unreliable until this is actually root-caused.

## Results so far (n=88 pools, full rep1-3 core grid + partial adversarial floor)

**Headline change from the n=32/rep1-only snapshot: adding adversarial pools
(the first pools in this design to contain an actually-close donor pair)
flipped raw accuracy from a null result to a real, significant design signal,
and flipped the min_dist-vs-mean_dist bottleneck comparison to match the
real-data result.** Both were explicitly predicted as likely once adversarial
pools landed (see the old point-3 reasoning kept below) — this snapshot is the
first evidence for or against that prediction, and it came out in favor.

**1. Raw singlet accuracy now shows a real, significant design signal.**
min_dist vs. accuracy: r=+0.349, p=0.001 (GEX); r=+0.564, p<0.001 (ATAC).
kl_min is an even stronger predictor: r=+0.399 (GEX), r=+0.668 (ATAC), both
p<0.001. (Previously: r=0.062 p=0.74 GEX, r=0.435 p=0.013 ATAC, n=32,
rep1-only, no adversarial pools.) See
`notebooks/ambisim/01b_pool_level_summary.ipynb`.

**2. Classifier margin (LL-gap) still shows a real, significant design
signal, now even stronger.** Pool-mean LL-gap vs. `min_dist`: r=0.657 (GEX),
r=0.729 (ATAC), both p<0.001; low-margin-call fraction (t100) vs. `min_dist`:
r=-0.507 (GEX), r=-0.679 (ATAC), both p<0.001.
(`notebooks/ambisim/01d_bottleneck_donor.ipynb`)

**3. `min_dist`/`kl_min` (worst-pair geometry) now beats `mean_dist`/`kl_mean`
(pool-average geometry) at the bottleneck-donor level — matching the real-data
result, reversing the n=32 snapshot.** Bottleneck LL-gap: min_dist r=+0.732
vs. mean_dist r=+0.593 (GEX); min_dist r=+0.765 vs. mean_dist r=+0.561 (ATAC).
This confirms the leading hypothesis from the previous snapshot: those 32
pools were all built from *unrelated* 1000G donors with no genuinely close
pair, so `min_dist` had no real edge case to detect yet; the adversarial pools
now landing are exactly what creates one. At the pool-*average* level the
picture is more mixed: GEX still slightly favors `mean_dist` (r=0.729 vs.
r=0.655), while ATAC now favors `min_dist` (r=0.728 vs. r=0.689) — so the
reversal is real and modality-dependent, not a clean sweep. Bottleneck
*true-donor accuracy* (as opposed to LL-gap) also picked up a real ATAC signal
that wasn't there before (r=+0.389, p<0.001; GEX r=+0.175, p=0.102, still
n.s.) — the first sign the margin effect is starting to convert into an
accuracy effect at the hardest-donor level, at least for ATAC. Caveat: only a
partial adversarial set has landed so far (44/132 pools not simulated at
all) — treat this reversal as a real, replicated-across-many-universes
signal, not yet the final word.

**4. Ambient RNA contamination still swamps design signal at high
contamination, but the crossover point moved with more data.** min_dist vs.
LL-gap p5 correlation: r=0.606/0.709 (GEX/ATAC) at low ambient (<10%),
dropping to r=0.046 (n.s., GEX) / r=0.352 (p=0.001, still significant, ATAC)
at high ambient (≥30%) — GEX design signal now collapses fully at high
ambient while ATAC retains a real, if weaker, effect even there.
(`notebooks/ambisim/01c_droplet_ambient_stratified_llgap.ipynb`)

**4b. The GEX/ATAC accuracy split from the n=32 snapshot has mostly resolved
now that adversarial pools are in the mix.** Previously GEX's pooled
accuracy-vs-`min_dist` correlation looked null everywhere, ambient bin or not
(only one of 5 bins marginally significant). Now: GEX accuracy vs. min_dist is
significant in 3 of 5 ambient bins (0-10%: r=0.567; 10-20%: r=0.539; 20-30%:
r=0.424; all p<0.001), losing significance only once ambient hits ≥30% where
accuracy is already collapsing toward its noise floor for every pool
regardless of design. ATAC is similar or stronger in the same low/mid bins
(10-20%: r=0.545; 20-30%: r=0.548; 30-40%: r=0.437, all p<0.001) but, unlike
GEX, keeps a bit of signal into higher ambient and loses it instead at the
very lowest ambient bin (0-10%: r=0.099, n.s. — likely a ceiling effect,
mean_acc=0.999 there, essentially no room for design to move the needle).
Net read: the earlier "GEX shows no design signal, period" conclusion was an
artifact of not yet having pools with a real close pair to differentiate on,
not a true modality difference — GEX and ATAC now look qualitatively similar
(design matters at low/mid ambient, washes out at the high-ambient floor),
with ATAC's edge being that it degrades more gracefully at the top end.

## Commentary — does pool design matter?

**Updated 2026-09-14: yes, on this evidence.** The n=32/rep1-only snapshot's
honest read was "not yet shown on accuracy, but real and robust on margin" —
that was explicitly conditioned on the adversarial floor not existing yet.
It now partially does (a growing subset of adversarial pools across several
universes, 88/132 pools total), and accuracy moved: min_dist/kl_min now
predict raw singlet accuracy directly (point 1 above), not just margin, and
the min_dist-over-mean_dist bottleneck reversal predicted in the old point 3
happened exactly as hypothesized. This is meaningfully different from a naive
"designed pools are better" story, though — the reversal is bottleneck-level
and ATAC-favoring more than pool-average and GEX-uniform, and the effect is
strongest for genuinely close pairs (adversarial) rather than a smooth
gradient across the "well-designed-ish" unrelated-donor space tested before.

What would still change this picture:
- 44/132 pools (mostly more adversarial reps and remaining rep2/rep3 of the
  core grid) haven't been simulated at all yet — see `NEXT_STEPS.md`. The
  correlations above are already many-universe and highly significant, but a
  bigger n, especially more adversarial_family/adversarial_family_mixed
  pedigree pools, would tighten the accuracy-side estimates (currently the
  weaker of the two effects, especially for GEX bottleneck accuracy, still
  n.s. at r=+0.175 p=0.10) and clarify whether the pool-average GEX result
  (still mean_dist-favoring) is a stable finding or will also flip with more
  adversarial data.
- The one comparison that could still test whether margin differences convert
  into accuracy differences under stress — downsampling/degradation
  robustness (`md/analyses.md` Tier 3) — hasn't been run in simulation yet.

The publication framing from the old snapshot ("a precise account of what a
VCF alone can and can't tell you pre-experiment") still holds, just with the
answer to "can it predict accuracy" now trending toward "yes, once the pool
actually contains a hard pair" rather than "no." That is a more interesting
and more publishable result than either a clean null or a naive "designed
pools are always better" — it says pool curation matters specifically insofar
as it avoids close pairs, which is a testable, actionable claim rather than a
vague diversity heuristic.

## Next steps — highest-value pools to simulate next

See `NEXT_STEPS.md` for the concrete pool list and qsub notes.
