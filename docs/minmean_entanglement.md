# min_dist / mean_dist entanglement — methodology

Why `greedy_maxmin` and `greedy_maxmean` pools cannot answer "does worst-case or average-case
genetic separation matter more?", and how the orthogonal-sampling construction fixes it.
Headline numbers are in `RESULTS.md` §6; this file keeps the full derivation, including the
corrections made along the way. Moved verbatim from `NOTES.md` on 2026-09-18.

**This matters for interpreting every `min_dist` vs. `mean_dist` comparison in `RESULTS.md`,
including §2 and §6 there.** `greedy_maxmin` (maximizes minimum pairwise distance) and
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
  null result (`RESULTS.md` §1) is unsound — that result doesn't depend on min/mean separation.

**Methodology of the fix itself, and a real limitation found and patched in it (2026-09-14):**
`C(525,8) ≈ 1.36e17` possible `EUR_only` pools — the 20,000-draw search samples ~1.5e-13 of that
space, nowhere near exhaustive. It finds "reasonably extreme relative to a modest random
sample," not the global extreme, and the two quadrants converge very differently: re-running the
search at 1K/5K/20K/50K/200K draws shows `highmin_lowmean`'s achieved `min_dist` converges
smoothly and is already near-final at 20K (103.19 → 103.42 over 10x more draws — there's a
natural ceiling on how spread out 8 people can be). `lowmin_highmean` is heavier-tailed: a
diagnostic run (different seed than the actual pools) found `min_dist` jumping to 45.4 somewhere
around 20K draws and holding there through 200K — well below the real `adversarial_family`
pedigree pairs' 74-77 range, meaning naive extremal search on the "how close can two `unrelated`
people be" question risks surfacing a cryptic-relatedness pair the panel's `unrelated` flag
missed, not a genuinely coincidental close pair. **The actual generated pools are fine** — their
seeds (200001/200002, not the diagnostic seed) landed at min_dist=80.6/87.2, and I verified both
closest pairs (`HG00116`/`HG00120`, `HG00240`/`HG00238`) are flagged `unrelated=True` in `meta`
— but the script had no safeguard against this, so one was added: it now asserts the closest
pair in each selected pool is flagged unrelated and warns if their distance drops below 70
(approaching the real pedigree range), so a future rerun with more candidates or a different
universe can't silently turn a "control" pool into an accidental adversarial one.

Also worth being explicit about: the median-split design (below/above-median `mean_dist`, then
`max()`/`min()` `min_dist` within that half) is a reasonable default, not a proven-optimal one.
A stricter percentile cutoff (e.g. bottom quartile of `mean_dist`) would sharpen separation on
the constrained axis but shrinks the candidate pool available for the extremal search on the
target axis — extreme order statistics need a large sample to find a good tail value. This
tradeoff wasn't optimized, just chosen as the simplest 50/50 split that keeps the full 20,000
candidates available for the search.

**Why does random sampling correlate min_dist/mean_dist at all (r=0.54), given 525 candidates
to draw from?** Not a sample-size effect — checked directly. `min_dist` vs. the mean of the
*other 27 pairs* (excluding the min pair itself): r=0.494, barely lower than r=0.536 including
it. `min_dist` vs. the *second*-smallest pairwise distance in the same draw: r=0.72. This means
whether a random 8-person draw is "tight" or "spread out" is a property of the whole set, not
independent per-pair noise: population genetic-distance space isn't uniformly dense, so a draw
that happens to land 2 people in a locally dense pocket tends to have several other reduced
pairwise distances too (not just that one pair), which drags both the minimum and the mean down
together. A huge candidate panel doesn't prevent this — it's geometry, not scarcity.

**Are the ~54 existing `greedy_maxmin`/`greedy_maxmean` pools (across all 9 universes) now
useless?** No, but their role is narrower than originally intended. Still valid for: the
`random`-vs-any-optimized-strategy comparison (doesn't depend on separating min from mean) and
the continuous distance-vs-margin regression (their real, achieved `min_dist`/`mean_dist` values
still contribute genuine variance there — this is how the r≈0.31 non-adversarial LL-gap
correlation was estimated). Not valid for: isolating whether *worst-case* optimization matters
more than *average-case* optimization specifically — that comparison needs the two conditions
to actually differ, which they mostly don't. Don't discard this data; just don't cite it for
that one narrower claim.

**Recommended path forward (not yet decided/executed beyond the 4-pool pilot above):** don't
redo all ~54 pools under the new scheme preemptively — each is a full ambisim→cellranger-arc→
demuxlet run (~13x the pilot's compute). Sequence: (1) let the 4-pool `EUR_only` pilot land
first; (2) if `highmin_lowmean` and `lowmin_highmean` pools show a real, different LL-gap/
accuracy signal from each other, that's the first evidence in this project that `min_dist`
specifically (not distance in general) matters, and scaling up (more universes, plus the
`highmin_highmean`/`lowmin_lowmean` quadrants for a full 2x2 grid) is well justified; (3) if
the pilot shows nothing, that's still a cheap, useful result — it would mean the min/mean
distinction was never the story either, and the compute is better spent on the Tier-3
downsampling check, which is more likely to be decisive for the paper's actual accuracy
question.

**Postscript (2026-09-15): re-confirmed on the genome-wide matrix — same problem, slightly
worse, not a chr22 artifact.** Once the distance matrix was rebuilt genome-wide (LD-pruned,
459,999 SNPs — see `JOURNAL.md`), the obvious question was whether more SNPs would relieve the
entanglement. It doesn't: `greedy_maxmin` min_dist-vs-mean_dist correlation is **r=0.992** (n=45
genome-wide pools) vs. chr22's r=0.69-0.72; donor overlap between `greedy_maxmin`/
`greedy_maxmean` is still a median 5/8 identical donors. This is consistent with (not a
surprise given) the mechanistic explanation above — it's order-statistic geometry of pairwise
distances within a fixed candidate panel, not a SNP-count/statistical-power problem, so there
was no reason to expect genome-wide data to fix it, and it didn't.

The orthogonal-sampling fix itself (`generate_orthogonal_pools.py`) was re-run against the
genome-wide matrix and generalized from `EUR_only`-only to all 9 universes
(`generate_orthogonal_pools_genomewide.py`) — nomination is cheap (~65s total for 18
universe×rep combinations), so this doesn't change the simulation-cost gating above, only
which donor lists would be ready if/when that gate opens. One real bug caught in doing this:
the cryptic-relatedness safeguard's threshold (`70`) is meaningless on the new matrix's 1-IBS
scale (~0.10-0.29 vs. the old Euclidean-ish ~0-150) — copying it verbatim made **every** pool
trip the warning, which is itself the signal something's miscalibrated. Recalibrated against
this matrix's own known-truth pairs (1,205 real parent-child pairs at D=0.1056-0.1543, whole-
matrix 1st percentile 0.1984 — clean separation) to threshold=0.16. After recalibration, 7/36
pools warn, all in `lowmin_highmean` specifically — the same "heavier-tailed" quadrant already
flagged above, not a new problem. `EUR_only` (the universe with an actual pilot in flight) is
clean.

**The actual pilot question (does `highmin_lowmean` differ from `lowmin_highmean` in demux
accuracy) is still unanswered as of this postscript** — the `EUR_only` simulation is mid-
pipeline (pileup done, calling not yet run). Everything in this postscript is about whether the
*construction* method still works genome-wide (it does), not about whether the underlying
scientific question has been resolved (it hasn't).
