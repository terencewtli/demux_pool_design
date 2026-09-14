# Methodology notes and known pipeline issues

Detail that supports `README.md` but would clutter it. Read this when you want the
full reasoning behind a headline claim, or when a metric looks wrong and you want to
check it's not a known issue first. Session-by-session narrative (who did what, when)
is in `JOURNAL.md`, not here — this file is organized by topic, not by date.

## Robustness check: always test correlations against the adversarial-pool subset

**Why this file's numbers sometimes differ from an initial README pass, and why any
future finding should be checked the same way before being written up.**

The design includes 15 "adversarial" pools (6 `adversarial_family`/`adversarial_family_mixed`
pedigree pools + 9 `adversarial_mindist` pools) alongside 73 pools built from ordinary
unrelated-donor strategies (`random`, `greedy_max*`, `ancestry_balanced`). The adversarial
pools sit at meaningfully lower `min_dist` than everything else, which makes them high-
leverage points in any pool-level correlation — a handful of points at the low end of the
x-axis can produce a "significant" Pearson r even when the other 73 points show no relationship
at all. This bit twice in the 2026-09-14 session (see `JOURNAL.md`): an initial pass reported
whole-sample correlations for accuracy-vs-design and the min_dist-vs-mean_dist bottleneck
comparison without this check, and both claims evaporated once the adversarial pools were
excluded.

**Concretely, for every correlation in this repo, compute it twice**: once on the full
sample, once with `strategy.str.startswith('adversarial')` filtered out. If the excluded-
adversarial version isn't significant (or reverses sign), the full-sample number describes
"what happens when you include extreme cases," not "a trend across the design space" — say
so explicitly rather than reporting the full-sample r alone. Also check Spearman rank
correlation alongside Pearson; a Pearson r that isn't backed by a similar-sized Spearman rho
is another sign of leverage rather than a monotonic trend (this caught the GEX accuracy
result: Pearson r=0.349 p=0.001 but Spearman rho=0.205 p=0.056, n.s.).

The full point-by-point old-vs-corrected numbers from the 2026-09-14 session are in
`JOURNAL.md`'s correction note; the current, correct numbers are in `README.md`.

## `greedy_maxkl` is buggy/weird — treat all its results as unreliable

**Status: partially fixed, still broken, deprioritized rather than resolved.** Do not use
`greedy_maxkl` pools as a valid "good design" data point until this is root-caused.

`nominate_greedy_maxkl`'s inner KL computation added its numerical-stability epsilon
*before* taking the complement probability (`1 - p`), so any homozygous-alt site (dosage=2
→ p=1.0) produced `1 - p < 0` and `log(negative) = NaN`, silently corrupting the greedy
argmax every time it hit a common homozygous site. This is why the 3 `greedy_maxkl` pools
originally simulated scored at or below `random` on their own target metric (`kl_min`) — the
nomination loop wasn't actually optimizing anything. Fixed in `pool_nomination.py`
(confirmed: no more NaN/divide-by-zero warnings post-fix).

**But re-running the (fixed) nomination across all 5 single-ancestry universes shows the bug
wasn't the whole story.** Re-nominated canonical (best-of-5-seed) `kl_min` values — AFR_only
1.738, AMR_only 1.702, EAS_only 1.498, EUR_only 1.639, SAS_only 1.799 — are still below that
universe's `random` baseline in every single case (random: 1.851, 1.725, 1.675, 1.773, 1.880
respectively), just by a smaller margin than the pre-fix version. So `greedy_maxkl`, even
after the NaN fix, still doesn't demonstrably achieve what its name claims. The search-time
objective (`kl_to_pool`) and the scoring function (`pool_kl_divergence`) were checked and now
use consistent epsilon-handling, so this isn't simply a scoring/search mismatch — root cause
not identified. Possibly a genuine property of greedy max-min search being more prone to
local-optimum traps for a KL-shaped objective than for Euclidean distance (`greedy_maxmin`,
using the same greedy structure, reliably beats `random` by a wide margin) — or a subtler
remaining bug not yet isolated. Re-checked again with the rep2/rep3 data added 2026-09-14:
still underperforms `random` in every one of the 5 universes, no change to this status.

**Decision: deprioritized, not resolved.** Given `kl_min`/`min_dist` are r=0.999 collinear in
this design space anyway, `greedy_maxmin` already serves the "good design by genetic
separation" comparator role, and the adversarial pools are a far higher-value use of
remaining compute for the paper's actual question. Donor lists were regenerated with the
fixed code (`txt/donors/*greedy_maxkl*rep1.txt`, backups of the pre-fix
`nominated_pools_n8.tsv`/`simulation_pools_n8.tsv` kept alongside), but the 3 already-
simulated `greedy_maxkl` pools (AMR_only, EUR_only, SAS_only, rep1) have deliberately not
been deleted/re-simulated — there was no point spending ~20h/pool of compute regenerating a
strategy that still isn't working as intended.

## `demuxlet.sh`'s default memory (8GB) silently OOM-kills on larger pileups

The shared `demux_benchmark/template_demux/demuxlet.sh` template requests only
`h_data=4G, -pe shared 2` (8GB total). Whichever modality happens to have the larger
pileup for a given pool (not consistently GEX or ATAC — depends on the pool's realized
variant count, ~3.8M–6.9M observed) gets OOM-killed (`Killed`, no error propagated) partway
through loading the pileup, and `A02b_demuxlet_call.sh` fire-and-forgets both `qsub` calls
without checking their outcome — so this fails completely silently unless someone diffs the
GEX-vs-ATAC `.best` counts. Found and worked around 2026-09-14 (see `JOURNAL.md`) by
resubmitting affected calls with `qsub -l h_data=24G` (a few needed 64G). The shared template
itself wasn't edited (it's used elsewhere, e.g. the IGVF pilot pipeline) — if this recurs at
scale on a future batch, raising the template's default is probably worth it rather than
continuing to patch individual resubmissions after the fact.

## Extra checks run 2026-09-15: ambient×margin ("can design rescue bad cells?") and ancestry-restricted comparisons

Two follow-up questions, both answered directly against non-adversarial pools only (see the
robustness-check note above for why):

**Does design help more for high-ambient ("bad") droplets, or less?** Less. LL-gap vs.
`min_dist`, by ambient bin (non-adversarial pools only): the correlation is strongest at low
ambient (0-10%: r=0.31 GEX / 0.31 ATAC, both p<0.01) and weakens monotonically, vanishing by
40%+ ambient (r=-0.06 GEX n.s. / r=0.19 ATAC n.s.). So design's margin benefit is a cushion
available in already-decent conditions, not a rescue for droplets that are already in trouble
— the opposite of what "design saves bad cells" would predict. Per-pool sequencing depth
(`rna_nr_c + rna_nr_a`) was checked as a possible additional confound alongside ambient
fraction: it's essentially fixed by simulation design (median identical to the integer across
all 88 pools, ~20,021 RNA / ~20,240 ATAC reads), so no additional stratification by depth was
needed on top of the existing ambient-fraction bins.

**Does restricting to single-ancestry universes (where all real-data comparisons like
`AFR`/`EUR` panels are) sharpen the signal, on the theory that mixing ancestries dilutes it
across too much of the genetic continuum?** No — if anything, the opposite. Per-universe
correlations at n=8-10 pools each are noisy and directionally inconsistent (e.g. `AMR_only`
GEX ll_gap r=-0.75 p=0.02 vs. `SAS_only` r=+0.86 p=0.006 — opposite signs, likely both just
underpowered noise). Pooling all single-ancestry universes together (n=44) gives ll_gap
r=+0.20 (n.s.), weaker than the full non-adversarial sample's r=0.31. Pooling the 4
multi-ancestry universes instead (n=29: `EUR_AFR`, `EUR_EAS`, `all_3_major`, `all_5`) gives
the *strongest* margin correlation of any non-adversarial subset (ll_gap r=+0.56 GEX p=0.002,
r=+0.59 ATAC p=0.001) — because multi-ancestry pools have the widest within-universe
`min_dist` range (span 16.5-24.8 vs. 9.4-16.1 for single-ancestry universes), giving the
distance metric more room to vary at all. So the original "ancestry continuum dilutes the
effect" hypothesis was backwards: mixing ancestries is what creates enough distance variance
to detect a design effect in the first place; restricting to a single ancestry mostly just
removes power. Caveat: n=29 is still thin evidence, not a settled result.

## `min_dist` and `mean_dist` are far more entangled in our greedy-built pools than in real pools — a real methodology gap, not a null-result artifact

**This matters for interpreting every `min_dist` vs. `mean_dist` comparison in `README.md`,
including point 3 there.** `greedy_maxmin` (maximizes minimum pairwise distance) and
`greedy_maxmean` (maximizes mean pairwise distance) are meant to be two different optimization
targets, testing whether worst-case or average-case genetic separation matters more. In
practice they produce highly similar pools:

- **Donor overlap**: for the same (universe, rep), a median of 5-6 of 8 donors are the
  *identical* individuals chosen by both strategies (range 1-7/8 across all 15 single-ancestry
  and 12 multi-ancestry universe×rep combinations checked). Two causes: (1) both
  `nominate_greedy_*` functions in `pool_nomination.py` use `np.random.default_rng(seed)` with
  the *same seed* per rep, so they're guaranteed to start from the identical first donor before
  the greedy paths can diverge; (2) more fundamentally, a candidate that's far by the min-
  criterion from the current pool is usually also far by the mean-criterion when drawing from
  a single large, homogeneous, unrelated-donor panel, so the two greedy searches keep landing
  on the same outlying candidates regardless of the seed issue.
- **Metric correlation**: across the 88 simulated pools, `min_dist` vs. `mean_dist` r=0.69-0.72
  (non-adversarial), rising to r=0.84 in multi-ancestry universes specifically. `greedy_maxmin`
  boosts `mean_dist` by +10.4 over `random` — almost exactly what `greedy_maxmean` itself
  achieves (+10.9) — and `greedy_maxmean` boosts `min_dist` by +6.6, more than half of what
  `greedy_maxmin` achieves (+12.3). Neither strategy isolates its target axis.

**Real data shows these two metrics are close to independent, not correlated.** In the 63
real 10x pools (`notebooks/vcf_metrics/01a_expected_ll.ipynb` /
`notebooks/real_data/01h_bottleneck_donor.ipynb` source data,
`demux_benchmark/pool_design/csv/50_line/pool_genetic_dists.csv`): **min_dist vs. mean_dist
r=0.137 (p=0.29, n.s.)**. The mechanism is clean: `mean_dist` vs. each pool's EUR/AFR ancestry
ratio (`n_eur`, the only compositional axis that varies in this real cohort) is r=-0.91 —
essentially deterministic, since mean distance is a population-average quantity that tracks
aggregate ancestry composition almost perfectly. `min_dist` vs. that same ancestry ratio is
r=-0.16 (n.s.) — min_dist is set by whichever *specific pair* of individuals happens to be
closest, a local/idiosyncratic property of exactly who got pooled, independent of the
aggregate mix. Nothing forces these two kinds of quantity (a population average vs. an
extremal/order statistic) to move together — unless the construction procedure explicitly
makes them move together, which greedy point-repulsion does.

**Correction (2026-09-14, `ipynb/ambisim/00_pool_geometry_qc.ipynb`): an earlier draft of this
note claimed `random` pools decouple min_dist/mean_dist to r=-0.09, "matching real data almost
exactly."** That was based on only 7 actually-simulated `random`-strategy pools -- nowhere near
enough power to trust a correlation estimate (the same small-n trap flagged elsewhere in this
file). A properly powered check -- 20,000 genuine random 8-donor draws from the `EUR_only`
panel, distance-matrix lookups only, no simulation needed -- gives **r=0.536 (p<1e-100-ish, a
real effect)**, not zero. Random sampling does not fully decouple these metrics: if a random
draw happens to be broadly spread out, the closest pair also tends to be a bit farther apart,
just from the geometry of pairwise-distance order statistics, independent of any optimization.

The corrected three-tier picture: **real data r=0.14** (most decoupled) < **random draws
r=0.536** (a real, moderate, unavoidable baseline) < **greedy_maxmin/greedy_maxmean r=0.69-0.84**
(more entangled than even random). So greedy search still measurably over-entangles the two
metrics relative to pure random sampling -- the qualitative conclusion below is unchanged -- but
it is not true that "any unoptimized sampling" reproduces real data's independence; real data is
*more* decoupled than even i.i.d. random draws from a single large reference panel, plausibly
because real pool assignment isn't literally uniform sampling from one big homogeneous panel --
it's a smaller, fixed donor-line resource under study-specific ancestry-ratio constraints, which
apparently has different pairwise-distance geometry than a clean random draw. `ancestry_balanced`
(proportional ancestry sampling, no distance optimization) sits at r=0.64 (n=7, still
underpowered, not directly comparable to the well-powered numbers above). The correlation still
scales with how much explicit "push everyone apart" optimization pressure is applied — it is not
an intrinsic property of unrelated donors being "too similar to distinguish" (real data's min_kl
r=0.40 beating mean_kl r=0.13 against observed LL-gap already shows these metrics carry real,
different signal among unrelated donors when allowed to vary independently).

**What this means for the project, concretely:**
- The `greedy_maxmin` vs. `greedy_maxmean` comparisons throughout this repo were never a clean
  test of "does worst-case vs. average-case optimization matter" — they're a comparison of two
  algorithms that produce 50-60% the same donor list. This does not affect the broader
  `random`-vs-everything-else null result, which doesn't depend on separating min from mean.
- **A reviewer saying "min/mean dist are picking up the same donor structure" would be
  correct**, and the fix is not more replicates of the same greedy strategies — it's a
  different construction procedure. Plain random draws only partially help (r=0.536, not zero
  -- see the correction above), so the actual fix implemented is rejection sampling *on top of*
  random draws: generate a large pool of random candidates (20,000 per rep, cheap since it's
  distance-matrix lookups only), then explicitly keep the ones landing off the natural
  min/mean relationship -- e.g. among candidates with below-median `mean_dist`, pick the one
  with the highest `min_dist` ("high_min_low_mean"), and the mirror image
  ("low_min_high_mean"). See `scripts/ambisim_new/lib/generate_orthogonal_pools.py` and the
  4 `*_new` pools it produced (`EUR_only__highmin_lowmean_new__rep{1,2}`,
  `EUR_only__lowmin_highmean_new__rep{1,2}`) -- `00_pool_geometry_qc.ipynb` confirms these sit
  far off the random-draw relationship (e.g. `lowmin_highmean` pools: min_dist at the 0th
  percentile of random draws while mean_dist sits at the 75th-95th percentile), i.e. real
  disagreement between the two metrics that neither greedy search nor plain random sampling
  would produce on its own.
- This is a fixable methodology gap specific to how `min_dist`/`mean_dist` design pools are
  built, not evidence the project's underlying question is unanswerable or that the accuracy
  null result (`README.md`) is unsound — that result doesn't depend on min/mean separation.
