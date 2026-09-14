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
