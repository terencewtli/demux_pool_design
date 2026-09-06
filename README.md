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
(`txt/pool_experiments.txt`). **As of this snapshot, 32/132 pools have
completed demuxlet on both modalities** — all rep1 of the core grid (no
adversarial pools have been simulated yet).

## Bug found and fixed during this analysis

`nominate_greedy_maxkl`'s inner KL computation added its numerical-stability
epsilon *before* taking the complement probability (`1 - p`), so any
homozygous-alt site (dosage=2 → p=1.0) produced `1 - p < 0` and
`log(negative) = NaN`, silently corrupting the greedy argmax every time it hit
a common homozygous site. This is why the 3 `greedy_maxkl` pools already
simulated scored at or below `random` on their own target metric (`kl_min`) —
the nomination loop wasn't actually optimizing anything. Fixed in
`pool_nomination.py` (verified: post-fix, `greedy_maxkl` now converges to the
same pool as `greedy_maxmin` on a smoke test, as expected given `kl_min` and
`min_dist` are r=0.999 collinear in this data — see below). **The 3 already-run
`greedy_maxkl` pools (AMR_only, EUR_only, SAS_only, rep1) reflect the buggy
selection and should be re-nominated and re-simulated**, not treated as valid
`greedy_maxkl` data points.

## Results so far (n=32 pools, rep1 only, no adversarial floor yet)

**1. Raw singlet accuracy shows no detectable design signal.** Range is
0.809–0.859 (GEX) / 0.905–0.924 (ATAC) across pools — real variation, but
uncorrelated with pool genetics: accuracy vs. `min_dist` is r=0.062 (p=0.74,
GEX) and r=0.435 (p=0.013, ATAC). See `notebooks/ambisim/01b_pool_level_summary.ipynb`.

**2. Classifier margin (LL-gap) shows a real, significant design signal that
accuracy misses.** Pool-mean LL-gap vs. `min_dist`: r=0.51 (GEX), r=0.56
(ATAC), both p<0.005; low-margin-call fraction vs. `min_dist`: r=-0.37 (GEX),
r=-0.61 (ATAC, p<0.001). This holds at every aggregation level tried,
including the per-pool worst-behaved donor ("bottleneck") specifically
(`notebooks/ambisim/01d_bottleneck_donor.ipynb`) — even the bottleneck
donor's *accuracy* shows no design correlation (r=-0.02 to +0.25, all n.s.),
while its LL-gap does (r=0.56-0.77).

**3. `mean_dist`/`kl_mean` (pool-average geometry) consistently outpredicts
`min_dist`/`kl_min` (worst-pair geometry) in this data** — at both the
pool-average and bottleneck-donor level. This is the opposite of what the real
10x pool data suggested (`notebooks/vcf_metrics/01a_expected_ll.ipynb`:
`min_kl` r=0.40 vs. `mean_kl` r=0.13 vs. observed LL-gap;
`notebooks/real_data/01h_bottleneck_donor.ipynb`: bottleneck LL-gap vs.
`min_dist` r=0.35 vs. `mean_dist` r=-0.19). Ruled out: this isn't an
aggregation-level artifact (the bottleneck-matched test still favors mean).
Leading hypothesis: these 32 pools are all built from *unrelated* 1000G
donors, so no pool yet contains a genuinely hard (close) pair — `min_dist`
has no real edge case to detect yet. `min_dist`'s real-data advantage likely
only shows up once a pool contains an actually-close pair, which is exactly
what the adversarial pools (not yet run) are built to create.

**4. Ambient RNA contamination swamps design signal at high contamination.**
`min_dist` vs. LL-gap correlation is r≈0.48-0.53 (p<0.01) at low ambient
(<10%) and collapses to r≈0.07-0.17 (n.s.) at high ambient (≥30%) — accuracy
itself degrades to the same bad floor for every pool regardless of design once
ambient RNA is high enough (e.g. GEX: 0.99 at <10% ambient → 0.26 at ≥40%,
for every design tier alike). Design differentiation is a benefit you get in a
clean experiment, not a rescue for a noisy one.
(`notebooks/ambisim/01c_droplet_ambient_stratified_llgap.ipynb`)

**4b. Splitting the pooled *accuracy* correlation by ambient bin doesn't rescue
it for GEX, but partially does for ATAC.** If the pooled null GEX accuracy
result (point 1) were purely an artifact of averaging a real low-ambient
effect against a swamped high-ambient one, per-bin accuracy-vs-`min_dist`
correlations should look like the LL-gap pattern above. They don't, for GEX:
only one of 5 bins reaches even marginal significance (10-20%: r=0.365,
p=0.04 — doesn't survive correction for the 10 bin×modality tests run), and
higher-ambient bins trend slightly negative. ATAC is different: real signal
recovers in some non-ceiling bins (10-20%: r=0.593, p<0.001; 30-40%: r=0.443,
p=0.011). So this sharpens rather than resolves the GEX/ATAC split — ATAC's
pooled accuracy correlation (point 1) looks like a real, if ambient-band-
dependent, effect, while GEX's near-null result holds up across the whole
ambient range, not just at the ceiling. Plausible reason: ATAC pileups carry
far more informative SNPs/reads per cell (~3x larger raw LL-gaps than GEX),
giving design more room to move a thresholded call, not just the continuous
margin.

## Commentary — does pool design matter?

Not "no" — the honest read at this snapshot is **"not yet shown on accuracy,
but real and robust on margin, and the experiments that would actually settle
it haven't been run yet."** Three things keep this from being a clean null:

- Every pool tested so far comes from a fairly narrow, "well-designed-ish"
  slice of the space (1000G unrelated donors, single rep). The adversarial
  floor — the condition under which a real accuracy effect is most likely to
  appear — doesn't exist in the data yet.
- Margin (LL-gap) is a continuous, more sensitive statistic than a thresholded
  accuracy call; it showing an effect while accuracy doesn't is consistent
  with a ceiling effect (all pools currently easy enough that nobody's
  accuracy is being tested), not necessarily consistent with "the effect isn't
  real."
- The one comparison that could settle whether margin differences ever cash
  out into accuracy differences — downsampling/degradation robustness
  (`md/analyses.md` Tier 3) — hasn't been run in simulation yet.

If, after adversarial pools and a stress-test (lower depth / higher ambient)
are added, accuracy *still* doesn't move while margin keeps differentiating —
that's a genuinely useful, publishable negative result: it would directly
correct the likely-assumed-but-never-rigorously-tested practice of
ancestry/diversity-aware pool curation, with a precise account of what a VCF
alone *can* still tell you (fragility/margin, quantifiable pre-experiment) and
what it apparently can't (a guaranteed accuracy win). That is a more useful and
more publishable paper than a naive "designed pools are better," precisely
because most people assume the latter and have never checked it against
ground truth. But it's premature to call it null before the adversarial +
stress-test experiments are in.

## Next steps — highest-value pools to simulate next

See `NEXT_STEPS.md` for the concrete pool list and qsub notes.
