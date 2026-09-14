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

**Correction (2026-09-14, same day as first written): the initial pass through
this section overstated things by reporting whole-sample correlations without
checking whether they were driven by a handful of leverage points. They were.**
Every correlation below is now reported both on the full n=88 sample and with
the 15 adversarial pools excluded (73 pools: `random`/`greedy_max*`/
`ancestry_balanced`, all unrelated 1000G donors) — that split is the load-
bearing check in this whole section, not a footnote.

**1. Raw singlet accuracy: the "significant design signal" is entirely an
adversarial-pool artifact, not a general trend.** Full sample: min_dist vs.
accuracy r=+0.349 p=0.001 (GEX), r=+0.564 p<0.001 (ATAC). Excluding the 15
adversarial pools: r=-0.011 p=0.93 (GEX), r=+0.101 p=0.40 (ATAC) — nothing
left in either modality. GEX's full-sample correlation isn't even robust to
rank transform (Spearman rho=0.205, p=0.056, n.s.), a sign the Pearson number
was leverage from a few points, not a monotonic trend. **The honest
1-pool-level plot is a floor, not a slope**: 6 pools built from real 1000G
pedigrees (`adversarial_family`/`adversarial_family_mixed`, ~50% IBD relative
pairs) sit at accuracy 0.81–0.83 (GEX) / 0.90–0.91 (ATAC), measurably below
everyone else (mean 0.843/0.918; Mann-Whitney p=0.0008 GEX, p=0.0001 ATAC).
The 9 `adversarial_mindist` pools (unrelated donors chosen to minimize genetic
spread) show a smaller but still real gap (p=0.0016 GEX, p=0.0045 ATAC, mean
diff only +0.013/+0.004) despite `min_dist` values that heavily overlap the
ordinary `random`/`greedy_max*` pools — so even that gap isn't well explained
by `min_dist` as a continuous metric; it looks more like a relatedness/IBD
effect the Euclidean distance only partly captures. **Effect size, even
including every adversarial pool**: the fitted slope over the full observed
min_dist range predicts only a 1.2–1.9 percentage point accuracy swing.
Among the 73 non-adversarial pools specifically — the comparison that answers
"does curating unrelated donors by genetic distance move raw accuracy" — the
answer is still no, same as the n=32 snapshot.
See `notebooks/ambisim/01b_pool_level_summary.ipynb`.

**2. Classifier margin (LL-gap) is the one result that survives the
adversarial-pool check.** Full sample: min_dist vs. ll_gap_mean r=0.657 (GEX),
r=0.729 (ATAC). Excluding adversarial pools: r=+0.314 p=0.007 (GEX), r=+0.311
p=0.007 (ATAC) — smaller than the full-sample number (as expected, since the
adversarial pools are real high-leverage points, not noise) but still real
and significant among ordinary unrelated-donor pools alone. This replicates
the n=32 finding (r=0.51–0.56 there) — the point estimate is a bit lower here
(0.31), plausibly the n=32 estimate running high by chance rather than the
effect actually weakening, but the direction and significance hold up with
more than double the pools. `low_margin_frac_t100` vs. min_dist also survives
(r=-0.271 p=0.02 GEX, r=-0.349 p=0.003 ATAC, non-adversarial only). This is
still the best-supported result in the whole repo.
(`notebooks/ambisim/01d_bottleneck_donor.ipynb`)

**3. The "min_dist now beats mean_dist" bottleneck reversal does NOT survive
excluding adversarial pools — it's the same artifact as point 1.** Full
sample: bottleneck LL-gap vs. min_dist r=+0.732 (GEX) / r=+0.765 (ATAC) vs.
mean_dist r=+0.593 / r=+0.561 — min_dist ahead, as originally reported.
Excluding the 15 adversarial pools: min_dist drops to r=+0.295 (GEX) /
r=+0.271 (ATAC) while mean_dist actually *strengthens* to r=+0.658 / r=+0.666
— mean_dist wins again, same direction and similar magnitude as the n=32
snapshot (which had no adversarial pools at all). **This result is reversed
from what was first written: `min_dist`'s apparent real-data-like advantage
was the adversarial pools acting as leverage points on the correlation, not a
property that emerged once "any" close pair existed in the design.** It may
still be true that `min_dist` matters specifically once relatedness is real
(the adversarial_family pools are the closest analog to the real-data
pedigrees), but that's a different, narrower claim than "min_dist beats
mean_dist in general," which is what the first pass claimed.

**4. Ambient-stratified accuracy correlations (previously reported as
resolving the GEX/ATAC split): also mostly an adversarial-pool artifact.**
Redone excluding adversarial pools: GEX accuracy-vs-min_dist is significant in
only 1 of 5 ambient bins (10-20%: r=0.264, p=0.024 — down from 3/5 bins in the
full-sample version), and ATAC likewise drops to 1 of 5 (10-20%: r=0.343,
p=0.003, down from 3/5). This is essentially the same pattern the n=32
snapshot already reported ("only one bin marginally significant") — the
"split resolved" conclusion in the first pass was wrong.
(`notebooks/ambisim/01c_droplet_ambient_stratified_llgap.ipynb`)

## Commentary — does pool design matter?

**Still "not yet shown on accuracy, but real and robust on margin" — the
n=32 conclusion holds up, it just needed to survive an adversarial-pool
leverage check first.** An earlier version of this section (written the same
day, before that check) said adding adversarial pools had flipped accuracy to
a real signal and reversed the min_dist-vs-mean_dist comparison. Both claims
turned out to be entirely driven by the 15 adversarial pools acting as
leverage points on the correlation — see the corrections in "Results so far"
above. Among the 73 non-adversarial pools (the actual "does curating unrelated
donors help" question), accuracy is still null and `mean_dist` still beats
`min_dist`, same as n=32.

**What is genuinely new in this snapshot:** pools built from real close
relative pairs (`adversarial_family`/`adversarial_family_mixed`, ~50% IBD)
show a real, replicated accuracy *and* margin deficit relative to every other
strategy (Mann-Whitney p<0.001 both modalities) — small in absolute terms
(1–2 accuracy points) but the first direct evidence in this whole project that
margin differences can convert into an accuracy difference at all, under a
condition extreme enough to produce one. `adversarial_mindist` (unrelated
donors, worst-case-selected) shows a smaller version of the same thing despite
`min_dist` values that overlap ordinary `random` pools, suggesting the
Euclidean `min_dist` metric doesn't fully capture whatever adversarial
selection or relatedness is actually detecting.

So the honest framing is narrower than "does pool design matter": **curating
unrelated donors by genetic distance (any of the strategies tested —
`greedy_max*`, `ancestry_balanced`, vs. `random`) has not been shown to move
raw accuracy, even now at n=88 with real statistical power. Avoiding actual
close relatives clearly does.** Whether that second finding generalizes into
something a `min_dist`-style score can predict *before* you know a pool
contains relatives, or whether it only shows up once you already know to flag
relatedness directly (e.g. from pedigree/IBD data, not simulated Euclidean
genotype distance), is unresolved — the `adversarial_mindist` result is
suggestive but weak enough (p=0.002–0.005, mean diff as small as +0.004 ATAC)
that it needs more replicates before leaning on it.

What would change this picture:
- 44/132 pools (mostly more adversarial reps and remaining rep2/rep3 of the
  core grid) haven't been simulated at all yet — see `NEXT_STEPS.md`. More
  `adversarial_family`/`adversarial_family_mixed` replicates specifically
  would tell us whether the relatedness effect is as robust as its p-values
  suggest at n=6 pools, or partly a small-sample fluke of exactly which
  pedigrees got drawn.
- The one comparison that could test whether margin differences convert into
  accuracy differences under stress even for *unrelated* donors — downsampling
  /degradation robustness (`md/analyses.md` Tier 3) — hasn't been run in
  simulation yet. This is arguably higher-value than more adversarial
  replicates for answering the paper's actual question.
- **Any new correlation reported in this repo going forward should be checked
  against the adversarial-pool subset before being written up as a finding**
  — this snapshot is the second time (after the n=32 min_dist/mean_dist
  reasoning) that a plausible-sounding pattern turned out to hinge entirely on
  a small, structurally different subset of pools rather than a trend across
  the design space.

## Next steps — highest-value pools to simulate next

See `NEXT_STEPS.md` for the concrete pool list and qsub notes.
