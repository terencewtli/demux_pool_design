# Journal

Reverse-chronological. One entry per session that changes project state. This file — not a
proliferating set of dated `md/YYYYMMDD_*.md` files — is the single place a new session should
read first. See `README.md` for current project state/design and `NEXT_STEPS.md` for the pool
priority list; this file is *why* it got that way and what's still open, in the order it happened.

Entry template:

```
## YYYY-MM-DD — one-line summary

**State at start:** …
**Decisions made:** …
**Produced:** (files/commits)
**Open / next:** …
**If resuming, read:** (the one or two files that matter most)
```

---

## 2026-09-15 — genome-wide matrix landed (LD-pruned, superseding the 200k-reservoir attempt); min/mean entanglement re-confirmed (worse, not better) genome-wide; orthogonal-sampling fix extended to all 9 universes with a real threshold-calibration bug caught and fixed along the way

**State at start:** the 2026-09-14 evening entry left the genome-wide matrix rebuild
"in progress" using a 200,000-SNP reservoir sample (`build_genomewide_dist_matrix.py`, no LD
pruning) and `regenerate_pools_genomewide.py` written but not yet run. Between that entry and
this session, an **undocumented gap**: someone (this session found the job already running,
never wrote it) replaced the reservoir-sample approach with a proper LD-pruned PLINK pipeline
(`build_genomewide_dist_matrix_plink.py`: VCF→bed, `--indep-pairwise 50 5 0.2`, 1-IBS distance
on the pruned 459,999-SNP set) — the right fix (raw/unpruned SNPs over-weight whatever regions
happen to have long-range LD), but never journaled. This entry closes that gap retroactively
and covers everything from this session on top of it.

**Decisions made / found:**

1. **The PLINK matrix job OOM-killed twice, from two different, unrelated bugs** — worth
   recording precisely since both are non-obvious:
   - First OOM (silent, no traceback, log stops right after the `D shape` print, maxvmem 338G
     against a 128G request): assumed to be `sklearn.PCA(svd_solver='randomized')` on the
     3202×460k dosage array. Fixed by swapping in `flashpca` (compiled, genotype-specific PCA)
     on a separately-built pruned `.bed`/`.bim`/`.fam` — this fix alone dropped the PCA step
     from a 3+ hour crash to 3m37s.
   - **But the resubmitted job climbed past the same 338G crash point anyway** (caught by
     watching `qstat -j`, not assumed fixed just because PCA was) — the *actual* culprit was
     `pandas.read_csv` on the wide (3202-row × 460k-column) dosage `.raw` file: pandas' wide-
     dataframe path has per-column overhead that doesn't show up on a normal (many-rows,
     few-columns) table, regardless of how tight the requested dtype is. Fixed by parsing the
     file as plain text into a preallocated `float32` array via `np.fromstring` per line
     (`--output-missing-genotype 0` on the PLINK recode avoids `'NA'` string handling). Tested
     on the real file first (~3 min extrapolated, correct 0/1/2 values) before resubmitting.
   - Net result: `csv/designs/1kg_genomewide_dist_matrix.npz` (D, X, sample_ids, pcs) built
     successfully, 6m54s end-to-end once both bugs were fixed.

2. **Re-ran the min_dist/mean_dist entanglement check (NOTES.md's existing section) on the new
   genome-wide matrix — same problem, not better.** Using the fast-strategy pools generated
   against the new matrix (`nominated_pools_n8_genomewide.tsv`, 9 universes × 5 reps):
   `greedy_maxmin` min_dist-vs-mean_dist **r=0.992** (n=45) — *higher* than chr22's documented
   r=0.69-0.72, not lower; `greedy_maxmean` r=0.520; `random` r=0.210 (chr22's well-powered
   20k-draw random-pool estimate was r=0.536 — genome-wide random draws look somewhat more
   decoupled, though this specific n=45 estimate is far less powered than that 20k-draw check
   and shouldn't be over-read). Donor overlap between `greedy_maxmin`/`greedy_maxmean` pools:
   still a median 5/8 identical donors, same as chr22. **More SNPs did not fix the entanglement**
   — consistent with NOTES.md's own mechanistic explanation (order-statistic geometry of
   pairwise distances in a fixed candidate panel, not a SNP-count/power problem), now verified
   rather than assumed to carry over.

3. **Found and superseded a duplicate-effort collision.** Before finding
   `regenerate_pools_genomewide.py` (item 1's undocumented predecessor session had already
   written it, targeting the exact matrix filename this session also built), this session
   independently wrote a narrower `create_pools_genomewide.py` (fast strategies only, 31 of the
   full 132 universe×strategy combinations) and started a separate `greedy_maxkl`-only qsub job
   — both writing to the same `nominated_pools_n8_genomewide.tsv` path the complete script also
   targets. Caught by reading this file's own JOURNAL/NOTES before finalizing anything further
   (should have been the first move, not a later one). Resolution: `qdel`'ed the redundant
   `greedy_maxkl` job (14752676), left the narrower script in place with a superseded-by note
   rather than deleting it, and submitted `regenerate_pools_genomewide.py` instead (job 14752950,
   `scripts/ambisim/qsub/B04_regenerate_pools_genomewide.sh`, new — regenerates the full,
   original 132-pool grid's donor selections + metrics against the genome-wide matrix, copying
   `adversarial_family`/`adversarial_family_mixed` donor lists unchanged since their pedigree
   structure doesn't depend on which chromosome's distances are used).

4. **Extended `generate_orthogonal_pools.py`'s rejection-sampling fix to the genome-wide matrix,
   generalized from `EUR_only`-only to all 9 universes** (nomination is a cheap distance-matrix
   lookup — 20,000 draws/universe/rep, ~65s total for all 18 universe×rep combinations — so
   this is not the same cost decision as scaling up simulation, which NOTES.md's existing
   "recommended path forward" already gates behind the `EUR_only` pilot's outcome; that gate is
   unchanged, this is nomination only). New file:
   `scripts/ambisim_new/lib/generate_orthogonal_pools_genomewide.py`.
   - **Caught a real bug while doing this**: the cryptic-relatedness safeguard's "warn if
     closest pair is this close" threshold (`70`) was copied verbatim from the chr22 script,
     whose `D` was a raw Euclidean-ish genotype distance (real range ~0-150). The genome-wide
     matrix uses PLINK 1-IBS distance (real range ~0.10-0.29) — the old threshold fired on
     **18/18 pools** on first run, which is itself the tell (a 100% hit rate means the check is
     miscalibrated, not that every pool is secretly related). Recalibrated empirically rather
     than guessing a new number: 1,205 real 1000G parent-child pairs sit at
     D=0.1056–0.1543 (mean 0.1359, std 0.009) against a whole-matrix 1st percentile of 0.1984 —
     clean separation, no overlap — so the threshold is now `0.16`. Re-run: 7/36 pools flagged
     (all in the `lowmin_highmean` quadrant specifically, matching NOTES.md's already-documented
     "heavier-tailed" behavior for that quadrant at chr22 scale — not a new phenomenon).
     `EUR_only` (the universe with an actual simulation in flight, see below) is clean at both
     reps.

5. **Found the `EUR_only` orthogonal-sampling pilot (the 4 pools from the 2026-09-14 entry) is
   actively mid-simulation right now** (`ambisim_new/EUR_only__{highmin_lowmean,lowmin_highmean}
   _new__rep{1,2}/`, files timestamped today — ATAC pileup complete, demuxlet calling not yet
   run, no `.best` output exists yet). So the actual open question this whole entanglement
   thread exists to answer — **does `highmin_lowmean` show a different LL-gap/accuracy signal
   from `lowmin_highmean`** — is still unanswered; not run by this session, found in-flight from
   elsewhere. This is the single most important open item, not the genome-wide re-derivation.

**Produced:**
- `csv/designs/1kg_genomewide_dist_matrix.npz`, `1kg_genomewide_sample_meta.tsv` (LD-pruned,
  459,999 SNPs, 3202 samples)
- `scripts/ambisim/lib/build_genomewide_dist_matrix_plink.py` (supersedes
  `build_genomewide_dist_matrix.py`, kept not deleted)
- `scripts/ambisim/qsub/B02_build_genomewide_dist_matrix.sh` (updated to call the PLINK script)
- `scripts/ambisim/lib/create_pools_genomewide.py` + `run_greedy_maxkl_genomewide.py` +
  `scripts/ambisim/qsub/B03_greedy_maxkl_genomewide.sh` (superseded by item 3 above, kept with a
  superseded-by note, not deleted)
- `scripts/ambisim/qsub/B04_regenerate_pools_genomewide.sh` (new qsub wrapper for the
  pre-existing `regenerate_pools_genomewide.py`) — **job 14752950, running as of this entry**
- `scripts/ambisim_new/lib/generate_orthogonal_pools_genomewide.py` — produced
  `csv/designs/nominated_pools_new_genomewide.tsv`, `txt/donors_new_genomewide/*.txt` (36 pools:
  9 universes × 2 quadrants × 2 reps), `txt/pool_experiments_new_genomewide.txt`. **Nominated
  only, not simulated.**

**Open / next:**
1. **Check job 14752950's output** (`csv/designs/nominated_pools_n8_genomewide.tsv`,
   `csv/designs/genomewide_vs_chr22_donor_overlap.tsv`) once it finishes — this is the
   authoritative full-132-pool-grid regeneration against the genome-wide matrix, superseding
   anything from item 3's narrower script.
2. **Get the `EUR_only` orthogonal pilot's actual demux/LL-gap result** — item 5 above. This
   answers the real question ("is the min/mean distinction meaningful for demux accuracy, not
   just for pool construction"), which nothing in this entry does on its own.
3. **7 flagged pools from item 4** worth a quick human look before treating as controls
   (same category of check as the 2026-09-14 `HG00116`/`HG00120` etc. verification) — not
   necessarily wrong, `unrelated` flag still asserts true for all of them, just worth eyes-on
   given they're specifically in the quadrant most prone to this.
4. Per NOTES.md's existing gating: don't submit ambisim/cellranger-arc/demuxlet for the 36
   new genome-wide orthogonal pools (or the regenerated 132-pool grid) until item 2 resolves —
   nomination being cheap doesn't change that simulation is the expensive part this project has
   repeatedly had to ration.
5. `scripts/ambisim/lib/build_genomewide_dist_matrix.py` (reservoir-sample, unpruned) is now
   superseded by `build_genomewide_dist_matrix_plink.py` — not deleted, but don't use it for
   anything going forward.

**If resuming, read:** this entry, then `NOTES.md`'s "min_dist and mean_dist are far more
entangled" section (now has a genome-wide postscript), then check job 14752950's log
(`ambisim/logs/B04_regenerate_pools_genomewide.*`) and whether the `EUR_only` pilot's demuxlet
calling has finished.

**Addendum (later same day): job 14752950 crashed, not completed** — `pool_nomination.py`'s
`nominate()` dispatcher had no case for `adversarial_mindist` (the function
`nominate_greedy_mindist` existed but was never wired in), so the job died 111/132 pools in
with `ValueError: unknown strategy: adversarial_mindist`. Since the script only writes output
once at the very end of the loop, **nothing was saved** — `nominated_pools_n8_genomewide.tsv`
on disk was still the narrower `create_pools_genomewide.py` file from earlier, not a
regeneration, and `genomewide_vs_chr22_donor_overlap.tsv` didn't exist at all. Checked all 8
distinct strategies actually present in `ambisim/txt/pool_experiments.txt` against the
dispatcher before resubmitting (`adversarial_family`/`adversarial_family_mixed` are correctly
excluded via `FAMILY_STRATEGIES`, everything else already had a case) — added the missing
`adversarial_mindist` branch, resubmitted as **job 14753496**. Also submitted (held on the
still-running GEX/ATAC pileup jobs via `-hold_jid`) demuxlet calling for the `ambisim_n16`
(n=16 mini-experiment) and `ambisim_new` (EUR_only orthogonal pilot) pileups, and wrote the
full `scripts/ambisim_final/` pipeline (A00-A02b, mirroring `ambisim_new`'s structure, pointed
at the regenerated `txt/donors_genomewide/` donor lists) — **not submitted**, explicitly gated
in each script's header behind job 14753496 landing and the `EUR_only` pilot's actual
demux-accuracy result, per this entry's existing "Open / next" list.

---

## 2026-09-14 (evening, same day) — discovered the distance matrix is chr22-only; building genome-wide replacement; full 132-pool redo decided

**State at start:** user asked how to interpret the raw `min_dist`/`mean_dist` numbers (e.g. 80.6
vs. 106.1) from the orthogonal-sampling work above. Tracing this back to its source
(`ipynb/vcf_metrics/local/02a_distance_matrix.ipynb`) surfaced something not previously
documented anywhere in this repo: **every `min_dist`/`mean_dist`/`kl_min`/`kl_mean` value used
to select or describe a pool in this entire project — the full 132-pool grid, the n=16
experiment, the orthogonal-sampling pools — comes from a distance matrix built from chr22 only
(the smallest human autosome), further subsampled to 20,000 SNPs (MAF>=0.05).** Not genome-wide.
The source notebook's own comment: "for pool nomination, ranking is all that matters" — this was
a deliberate speed shortcut for interactive prototyping, never revisited once cluster compute was
available (`N_SITES = None` for "all sites" was noted as an option but not used).

Also surfaced: the `unrelated` flag used throughout (`universe_idx_for`, the orthogonal-sampling
safeguard) is `FatherID=='0' AND MotherID=='0'` — pedigree-founder status, not a formal kinship/
IBD threshold. It rules out documented 1000G trio/quad relationships, not cryptic relatedness in
general. Relevant given the orthogonal-sampling safeguard added earlier this session relies on it.

**Assessment (not "everything is wrong"):** population ancestry structure is strong and
consistent across chromosomes, so the ancestry-level findings (universe-level effects, the
margin/design correlation, the pedigree/family results) likely replicate under a genome-wide
matrix largely unchanged. What's actually at risk is donor-level precision: exactly which 8
people get selected as "most extreme" within a universe has more sampling noise from one
chromosome's 20k SNPs than from the whole genome, and extremal/rejection-sampling searches (the
orthogonal-sampling pools especially) are more exposed to cryptic relatedness the single-
chromosome, founder-only `unrelated` check wouldn't catch. The simulation pipeline itself
(ambisim/cellranger-arc/demuxlet) is completely unaffected — it's donor-list-agnostic.

**Decision (user's call, made explicitly):** rebuild the distance matrix genome-wide, regenerate
every pool's donor selection against it, then **redo the full 132-pool grid from scratch** (user
said 135; the actual grid size is 132, see the 2026-09-14 n=88 entry above) once current jobs
(n=16, orthogonal-sampling) finish. Flagged before proceeding: `/u/project/cluo` is at **95.7%**
capacity lab-wide (`myquota`, checked this session) and a full 132-pool resimulation is a
multi-week, heavy-disk undertaking (recall the 27/132-task 14h-wall-time-cap casualty wave from
the 2026-09-14 morning entry) — worth staging rather than launching all 132 simultaneously.

**Produced this session (data not yet regenerated — matrix build in progress as of this entry):**
- `scripts/ambisim/lib/build_genomewide_dist_matrix.py` + `scripts/ambisim/qsub/B02_build_genomewide_dist_matrix.sh`
  (submitted, job 14743297): reservoir-samples 200,000 sites (~10x the original N_SITES, ~2.9%
  of the 6,869,484 available) from the SAME already-filtered, genome-wide, common+biallelic VCF
  (`1000G.merged.common_biallelic.vcf.gz`, 24 contigs) that A00_prep_pool_vcfs.sh actually
  subsets for simulation — so the distance matrix and the real simulated genotypes now come from
  the same underlying variant set, which chr22-only did not. Two-pass extraction (cheap
  CHROM/POS listing, then targeted `bcftools -R` pull of just the sampled sites) rather than a
  slow single pass parsing full genotypes for all 6.87M sites.
- `pool_nomination.py`'s `load_design_data()` generalized to take `matrix_file`/`meta_file`
  arguments, defaulting to the ORIGINAL chr22 files — every existing caller is unaffected;
  pass the genome-wide filenames explicitly to opt in.
- `scripts/ambisim/lib/regenerate_pools_genomewide.py`: regenerates all 132 pools' donor
  selections against the new matrix, EXCEPT `adversarial_family`/`adversarial_family_mixed` (6
  pools) — their relatedness is real pedigree data independent of which chromosome's distances
  are used, and their "fill" donors are chosen by uniform random choice, not by `D`/`X` at all,
  so re-deriving them risks picking a different family for no scientific reason; their donor
  lists are copied unchanged, but their `min_dist`/`mean_dist`/`kl_*` metrics are still
  recomputed against the new matrix (reported for every pool regardless of selection method).
  For the other 126 pools, re-runs the same `nominate()` used originally with the same
  `seed = SEED + rep*100` convention, even for strategies whose own logic doesn't touch `D`/`X`
  (`random`, `ancestry_balanced`) — `universe_idx` is a list of integer positions into
  `sample_ids`, and the new genome-wide VCF's sample order isn't guaranteed to match the old
  chr22 VCF's, so re-deriving from scratch is safer than assuming identical output. Writes to
  new files (`nominated_pools_n8_genomewide.tsv`, `txt/donors_genomewide/`) for review — does
  NOT touch the existing `txt/donors/` or `nominated_pools_n8.tsv`.
- Not yet run: `regenerate_pools_genomewide.py` (waiting on the matrix build), any resimulation.

**Open / next:**
1. Once `B02_build_genomewide_dist_matrix` (job 14743297) finishes, run
   `regenerate_pools_genomewide.py` and review `genomewide_vs_chr22_donor_overlap.tsv` — it
   reports, per pool, how many of the 8 donors are identical to the old chr22-based selection.
   This tells you concretely how much actually changes before committing to resimulation.
2. Let the n=16 and orthogonal-sampling mini-experiments finish (both mid-pipeline, chr22-based
   — decide separately whether they're worth redoing genome-wide or are a big-enough effect/
   proof-of-concept either way).
3. Decide staging for the full 132-pool resimulation given the 95.7%-full shared quota — not
   recommended to launch all 132 simultaneously given the earlier 14h-wall-time crash wave.
4. `greedy_maxkl`'s "best-of-5-seed" canonical-selection ceremony (see the README/NOTES
   `greedy_maxkl` section) was NOT replicated in `regenerate_pools_genomewide.py` — it uses the
   plain `seed = SEED + rep*100` convention uniformly. Since `greedy_maxkl` is already flagged
   unreliable regardless of chr22-vs-genome-wide, this was judged not worth preserving, but flag
   if that reasoning should be revisited.

**If resuming, read:** this entry, then `qstat -u terencew | grep B02` for the matrix build
status, then `genomewide_vs_chr22_donor_overlap.tsv` once `regenerate_pools_genomewide.py` has
run, before deciding on resimulation scope.

---

## 2026-09-14 (later same day) — min_dist/mean_dist entanglement found, fixed, n=16 and orthogonal-sampling mini-experiments launched

**State at start:** n=88 results committed and pushed (see entry below). User asked a follow-up
about `min_dist` vs. `mean_dist` magnitude and whether the two objectives could be picking the
same donors "by coincidence."

**Investigation 1 — they do, substantially:**
- Donor overlap between `greedy_maxmin`/`greedy_maxmean` for the same universe/rep: median 5-6
  of 8 donors identical (checked all 15 single-ancestry + 12 multi-ancestry combinations).
  Two causes: both `nominate_greedy_*` use the same RNG seed per rep (guaranteed identical
  first donor), and more fundamentally both are the same greedy point-repulsion algorithm, so
  an outlier far by one criterion is usually far by the other too from the same candidate pool.
- `min_dist`/`mean_dist` correlate at r=0.69-0.84 across the 88 simulated pools (non-adversarial).
- Real 10x data (63 pools, `demux_benchmark/pool_design/csv/50_line/pool_genetic_dists.csv`):
  r=0.14 (n.s.) — far more independent. Mechanism confirmed: real `mean_dist` tracks each pool's
  EUR/AFR ancestry ratio almost deterministically (r=-0.91), `min_dist` doesn't (r=-0.16) — a
  population-average statistic and an extremal/order statistic have no reason to move together
  unless the construction procedure forces them to.

**Investigation 2 — an error in investigation 1, caught while building the QC notebook:** first
pass claimed `random` pools decouple the two metrics (r=-0.09) "matching real data almost
exactly," based on only 7 actually-simulated `random` pools — not enough power to trust. A
properly-powered check (20,000 genuine random draws, `00_pool_geometry_qc.ipynb`) gives
**r=0.536**, a real, moderate, unavoidable baseline correlation from pairwise-distance order
statistics, not zero. Corrected picture: real data r=0.14 < random draws r=0.54 < greedy
r=0.69-0.84. Same small-n-correlation trap flagged for the adversarial-pool leverage check
earlier — worth internalizing as a standing rule for this repo, not just a one-off fix.

**Decision and fix implemented:** built `scripts/ambisim_new/lib/generate_orthogonal_pools.py`
— rejection sampling on top of random draws (not a third greedy variant): draw 20,000 random
8-donor pools per rep (distance-matrix lookup only, no simulation, ~5s), then keep the ones
landing off the natural min/mean relationship (e.g. highest `min_dist` among below-median
`mean_dist` candidates, and the mirror image). Produced 4 pools in `EUR_only`
(`highmin_lowmean_new` / `lowmin_highmean_new`, 2 reps each, suffix `_new` per user request to
mark the new sampling scheme) in a separate `ambisim_new/` tree (mirrors the `ambisim_n16/`
pattern from the n=16 experiment below — keeps both mini-experiments from being silently
picked up by the n=8 grid's glob-based notebook discovery). `00_pool_geometry_qc.ipynb`
confirms real disagreement: `lowmin_highmean` pools sit at the 0th percentile of the random-draw
min_dist distribution while their mean_dist sits at the 75th-95th percentile.

**Also launched, same session: n=16 mini-experiment** (separate from the above — tests whether
donor *count* matters, not min/mean separation; see README/user request re:
https://elifesciences.org/reviewed-preprints/106769). 4 pools (`EUR_only`, `random`/
`greedy_maxmin`, 2 reps each, 16 donors) in `ambisim_n16/`. Pipeline fully parameterized
already (`nominate()` takes `n=`, `setup_drop.sh`/`run_ambisim.sh` take an arbitrary-length
donors file) — no new simulation code needed, just new arguments. `h_rt` set to 24h from the
start (not the n=8 grid's original, since-reverted 14h) and demuxlet call memory pre-set to
32GB (not the 8GB default that silently OOM-killed ~27/132 n=8 calls) to avoid rediscovering
both known failure modes for a 4-pool experiment.

**Adversarial pool walkthrough (recap for the record, numbers already in README point 1, ll_gap
numbers newly pulled this session):** family/family_mixed pools (n=6, real pedigrees) show BOTH
accuracy AND margin dropping together (GEX ll_gap_mean 7.03 vs. 9.92 rest, Mann-Whitney
p<0.0001; ATAC 18.96 vs. 27.44, p<0.0001) — unlike the ordinary min_dist-vs-margin comparisons
elsewhere in this repo, where margin moves but accuracy doesn't. `adversarial_mindist` (n=9,
unrelated donors, worst-case-selected) shows the same pattern more weakly (GEX ll_gap 8.93 vs.
9.92, p=0.018; ATAC 25.33 vs. 27.44, p=0.019). Read: ordinary distance-based curation moves
margin without moving accuracy (a small effect that doesn't cross the threshold to change a
discretized accuracy call); real relatedness moves both together, suggesting there is a
threshold and literal relatedness is large enough to cross it while ordinary design choices
aren't. This is the concrete sense in which "pool selection is not trivial" — it's not that
distance-optimization design doesn't matter at all, it's that avoiding actual close relatives is
the one selection decision demonstrated to move accuracy, and it's a different kind of decision
(a relatedness check) than what any of the distance-based strategies here implement.

**Open / next:**
1. Both mini-experiments (n=16, orthogonal-sampling) are mid-pipeline (ambisim/cr_arc/demuxlet
   not yet complete as of this entry) — check `qstat -u terencew` for `n16_*` and `new_*` jobs,
   then run demuxlet call and score against `ambisim_n16/`/`ambisim_new/` truth once ready.
   Neither has a results notebook yet — will need small dedicated analysis notebooks (not 01a-01d,
   which are wired to the main `ambisim/` tree and `pool_experiments.txt`).
2. Only 2 quadrants (4 pools) sampled for the orthogonal-sampling experiment so far
   (`highmin_lowmean`, `lowmin_highmean`) — the two "disagreement" cells. Could extend to
   `highmin_highmean`/`lowmin_lowmean` (roughly what greedy and adversarial_mindist already
   cover, respectively) for a fuller 2x2 grid if the disagreement cells show a real effect worth
   following up.
3. Real data is *more* decoupled (r=0.14) than even random simulated draws (r=0.54) — open
   question why, not yet investigated. Plausible: real pool assignment isn't i.i.d. sampling
   from one large homogeneous panel, it's a smaller fixed donor-line resource under
   study-specific ancestry-ratio constraints, which may have different pairwise-distance
   geometry. Worth a closer look if the orthogonal-sampling pools show a real effect and the
   project wants to push the simulation's geometry closer to real data's.

**If resuming, read:** `NOTES.md`'s "min_dist and mean_dist are far more entangled" section for
the full mechanism (including the correction), `00_pool_geometry_qc.ipynb` for the QC itself,
then `qstat` for both mini-experiments' status.

**Follow-up same session (user questions after the above landed):**
- **Why does random sampling correlate min_dist/mean_dist at all (r=0.54)?** Not underpowering
  — `min_dist` vs. mean of the *other 27 pairs* (excluding the min pair itself) is still r=0.494,
  and `min_dist` vs. the *second*-smallest distance in the same draw is r=0.72. A random draw's
  "tightness" is a property of the whole 8-person set, not per-pair noise — added to `NOTES.md`.
- **Per-donor localization of the adversarial-pool effect, split by pool's own distance
  matrix**: `adversarial_family` degradation is sharply localized to donors with a relative in
  the pool (GEX LL-gap 4.79 at-risk vs. 7.66 safe donor, same pool; ATAC 13.83 vs. 24.71).
  `adversarial_mindist`'s closest pair is barely distinguishable from the rest of its own pool
  (GEX 7.24 vs. 7.48; ATAC 23.14 vs. 23.16) — consistent with `nominate_greedy_mindist`
  uniformly compressing the whole pool's geometry (minimizes *max* distance to the pool at each
  greedy step) rather than embedding one bad pair. Added to `README.md` point 1.
- **Recommended path forward, not yet executed**: don't redo the ~54 existing
  `greedy_maxmin`/`greedy_maxmean` pools under the new orthogonal scheme preemptively (13x the
  4-pool pilot's compute) — let the pilot land first and decide from there. Full reasoning in
  `NOTES.md`.

---

## 2026-09-14 — root-caused the OOM-killed demuxlet calls, reran 01a-01d at n=88, first accuracy signal

**State at start:** 61/132 pools had both-modality demuxlet `.best` outputs (up from 32 on
2026-09-12, from the 96-pool array submitted that session working through the queue). A naive
`ls */demux/demuxlet/*/*.best | wc -l` gave 149, which overestimates readiness — it counts
lopsided pools where only one modality finished.

**Investigation:** 27 pools had exactly one modality's `.best` missing despite both pileups
(`.var.gz`) existing — the same silent-failure shape flagged but not root-caused for
`EUR_EAS__greedy_maxmean__rep1` on 2026-09-12. Root cause found this session: the shared
`demux_benchmark/template_demux/demuxlet.sh` template requests only `h_data=4G, -pe shared 2`
(8GB total). Whichever modality happens to have the larger pileup for a given pool (not
consistently GEX or ATAC — depends on the pool's realized variant count, ~3.8M-6.9M observed)
gets OOM-killed (`Killed`, no error propagated) partway through loading the pileup. `A02b`
fire-and-forgets both `qsub` calls and never checks their outcome, so this fails completely
silently unless someone diffs the GEX vs. ATAC `.best` counts.

**Decisions made:**
- Did not edit the shared `demuxlet.sh` template (used elsewhere, e.g. the IGVF pilot
  pipeline) — resubmitted the 27 affected (pool, modality) calls directly with
  `qsub -l h_data=24G` overrides instead. 3 of those (all sharing an unusually large ~6.87M-variant
  VCF) still OOM'd at 24G and needed `-l h_data=64G`. Worth raising `demuxlet.sh`'s default if this
  recurs at scale — 8GB is clearly too tight for these ambisim pileups in general, not just an
  edge case.
- Separately (pool_design infra, not this repo): found and reran the `AFR_only__greedy_maxkl__rep1`
  simulation blocked by a wrong-donor-samples VCF, killed a hung `cr_arc` job stuck 17.5h past its
  own wall-time deadline, and bumped `A01b_run_cr_arc.sh`'s `h_rt` 14h→24h after confirming (via
  the pipestance's own `_log` timestamp) that the 14h cap had killed 26 of 27 incomplete
  pipestances mid-run. Not part of this repo's history but relevant to why more pools are ready now.
- Added a readiness-reporting cell to `01a` (right after pool discovery) that diffs discovered
  `.best` files against the full `pool_experiments.txt` and prints ready/partial/not-started counts
  by name — makes the exact failure mode above visible on every rerun instead of requiring a manual
  GEX-vs-ATAC count diff to notice.

**Produced:**
- Reran `notebooks/ambisim/01a-01d` end to end at n=88/132 pools (up from 32) via
  `jupyter nbconvert --execute` using the `demux` kernel env (nbconvert itself isn't installed in
  `demux`; ran from `allcools`'s nbconvert pointed at the `demux` kernel).
- `README.md` "Results so far" and "Commentary" sections rewritten for n=88. Headline: adding
  adversarial pools flipped raw accuracy from null (n=32) to significant (min_dist r=+0.35 GEX /
  +0.56 ATAC, both p<0.005) and reversed the min_dist-vs-mean_dist bottleneck comparison to match
  the real-data direction (min_dist now wins), exactly as the n=32 snapshot's own hypothesis
  predicted. Full numbers in README.
- `greedy_maxkl` spot-checked again with the new rep2/rep3 data: still underperforms `random` on
  its own target metric (`kl_min`) in every one of the 5 single-ancestry universes. No change to
  its "unreliable, deprioritized" status.
- Updated `results/*.csv`, `results/nominated_pools_n8.tsv`, and `notebooks/ambisim/*.ipynb`
  (executed, with outputs) committed alongside.

**Open / next:**
1. 44/132 pools still not simulated at all (mostly further adversarial reps + remaining core-grid
   rep2/rep3) — see `NEXT_STEPS.md`. Rerun `01a-01d` again once more land; the readiness-reporting
   cell in `01a` will show exactly what's still missing without needing to re-derive it.
2. `demuxlet.sh`'s 8GB default is undersized for these pileups in general — if this OOM pattern
   recurs on the next batch, raise the template's default rather than continuing to patch
   individual `qsub -l h_data=` overrides after the fact.
3. Pool-average GEX still favors `mean_dist` over `min_dist` (unlike ATAC, which now favors
   `min_dist` at both bottleneck and pool-average level) — worth watching whether this flips too
   once more adversarial pools land, or is a stable GEX/ATAC difference.
4. Tier-3 downsampling/robustness check (`md/analyses.md`) — still the most direct way to test
   whether the margin advantage ever converts into an accuracy advantage under stress — still not
   run in simulation.

**If resuming, read:** `README.md` "Results so far" for current numbers, this entry for how n=88
was reached and what's still silently-failing-prone in the pipeline (item 2 above), and the
correction note immediately below for a mistake made and fixed the same day.

**Correction (same day, same session):** the first pass through `README.md`'s "Results so far"
reported whole-sample (n=88) correlations for accuracy-vs-design and the min_dist-vs-mean_dist
bottleneck comparison without checking whether they were driven by a handful of leverage points.
The user asked "which pools specifically contributed" and "is this a small effect size", which
prompted the check that should have been done before writing it up: excluding the 15 adversarial
pools from each correlation. Result — the "accuracy is now significant" and "min_dist now beats
mean_dist" claims both evaporate entirely without the adversarial pools (r drops to ~0 / reverses
back to favoring mean_dist), while the LL-gap margin correlation (the one result the n=32 snapshot
was actually confident about) survives the same check, just at a smaller point estimate (r≈0.31
non-adversarial vs. r≈0.66 full sample). The real new finding is narrower and more specific than
first written: pools built from actual close relative pairs (6 `adversarial_family`/
`_mixed` pools) show a small but statistically real accuracy and margin deficit
(Mann-Whitney p<0.001 both modalities) — not a general "genetic-distance-based design predicts
accuracy" result. README's "Results so far" and "Commentary" sections were rewritten in place to
reflect this; see the corrected numbers there rather than trusting a re-derivation from memory.
**Lesson for next session**: any pool-level correlation in this repo should be checked against the
adversarial-pool subset (`strategy.str.startswith('adversarial')`) before being reported as a
finding — this is the second time (after the original min_dist/mean_dist reasoning in the n=32
snapshot) that a structurally different small subgroup, not a trend across the design space, turned
out to be doing the work.

---

## 2026-09-12 — full remaining-pool array submitted (all 96 unsimulated pools + 3 stuck + 1 partial)

**State at start:** 32/132 pools had completed demuxlet on both modalities (README's "n=32,
rep1 only, no adversarial floor yet" snapshot). Re-checked actual on-disk state this session
rather than trusting that snapshot, since it's from 2026-09-06/07 and could be stale:
- All 132 pools already have a per-pool ground-truth VCF prepped (`ambisim/vcf/`) — including the
  8 adversarial pools `NEXT_STEPS.md` flagged as still needing `A00_prep_pool_vcfs.sh`. That step
  had already been run since `NEXT_STEPS.md` was written.
- 36 pools (lines 1-36) have ambisim simulation output; the other 96 (lines 37-132: all rep2/rep3
  replicates of the core grid, plus every adversarial pool beyond the 8 already VCF-prepped) have
  none yet.
- Of those 36, only 32 have finished demuxlet on both GEX and ATAC. 4 exceptions found:
  - Lines 1 (`AFR_only__greedy_maxmean__rep1`), 16 (`EAS_only__greedy_maxkl__rep1`), 36
    (`AFR_only__greedy_maxkl__rep1`) — matches `NEXT_STEPS.md`'s "3 stuck at cr_arc" list exactly;
    `cr_arc` was never run for these.
  - Line 25 (`EUR_EAS__greedy_maxmean__rep1`) — not previously flagged anywhere. `cr_arc` and both
    pileups completed, GEX demuxlet `.best` exists, but ATAC demuxlet call was never run (no
    `.best` in the ATAC dir despite the pileup `.var.gz` being present). Root cause not
    investigated — just backfilled the missing call.

**Decisions made:**
- User asked to launch the full remaining ~100 pools (not the ~11-pool NEXT_STEPS.md-prioritized
  subset) — an explicit choice to override that document's deferral of rep2/rep3 replicates and
  the general "greedy_maxkl is unreliable, don't spend more compute on it" guidance in `README.md`.
  Flagged both the disk-usage implication and that override before submitting; user confirmed
  proceeding with the full submission both times.
- `/u/project/cluo` is at 94.3% full (per `myquota`), ~34TB free lab-wide. Estimated the 96 new
  pools at ~8TB total (`submit_adversarial_rep1.sh`'s own ~83GB/pool estimate × 96) — within
  headroom but a real bite out of shared lab quota. User chose to proceed with the full
  submission rather than batching it down.
- SGE's `-t` flag only accepts one contiguous range, not a comma list — `-t 1,16,36` fails
  (`qsub: ERROR! -t option only allows one range specification`). Submitted lines 1/16/36 as three
  separate single-task array-job chains rather than one combined job.

**Produced:** No code changes. Submitted jobs (all via `qsub` from
`/u/project/cluo/terencew/claude/project_ideas/pool_design`):

| lines | stage(s) submitted | job IDs |
|---|---|---|
| 37-132 (96 new pools) | `A01a_run_simulation` → `A01b_run_cr_arc` → `A02a_demuxlet_pileup` (chained via `-hold_jid`) | 14732184 → 14732186 → 14732187 |
| 1 | `A01b_run_cr_arc` → `A02a_demuxlet_pileup` | 14732191 → 14732193 |
| 16 | `A01b_run_cr_arc` → `A02a_demuxlet_pileup` | 14732194 → 14732195 |
| 36 | `A01b_run_cr_arc` → `A02a_demuxlet_pileup` | 14732196 → 14732197 |
| 25 | `A02b_demuxlet_call` (ATAC only; GEX already had `.best` and was skipped by the script's own check) | 14732188 |

This `JOURNAL.md`, added and pushed to `github/demux_pool_design`.

**Open / next:**
1. **`A02a_demuxlet_pileup` only launches the actual pileup jobs (`A01a_popscle_gex`/
   `A01a_popscle_atac` per pool) and exits within seconds — it does not wait for them**, so
   `A02b_demuxlet_call` cannot be safely chained with `-hold_jid` on it (same caveat
   `submit_adversarial_rep1.sh` already documents). Once `qstat -u terencew` shows the
   `popscle_gex`/`popscle_atac` jobs have cleared for a given line range, submit the final stage:
   `qsub -t 37-132 scripts/ambisim/qsub/A02b_demuxlet_call.sh` and similarly `-t 1`, `-t 16`,
   `-t 36` (one at a time, same `-t` range-only limitation).
2. Check for failed tasks across this large a submission (96 ambisim + 96 cr_arc + 192 pileup
   jobs) before assuming the array finished cleanly — re-run failed task IDs individually rather
   than resubmitting the whole array, per `A01c_rerun_failed_pool.sh`'s existing pattern.
3. `greedy_maxkl` is still flagged unreliable in `README.md` (fixed NaN bug, still underperforms
   `random` on its own target metric, root cause not identified) — this submission includes its
   rep2/rep3 replicates (lines 43/44, 51/52, 75/76, 83/84, 107/108) anyway per the user's explicit
   "run everything" choice. Treat those results with the same skepticism the README already
   applies to the existing `greedy_maxkl` rep1 data.
4. Once the full array lands, re-run `notebooks/ambisim/01a-01d` (they auto-discover pools via
   glob) — now with adversarial pools and full replicate depth, this is the first point where
   `min_dist` vs. `mean_dist` and the accuracy-vs-margin question (README "Commentary" section)
   can actually be re-tested rather than hypothesized about.
5. Re-investigate why line 25's ATAC demuxlet call silently never ran despite its pileup
   completing — if it happened once outside any documented failure mode, it could recur across
   the 96-pool submission and wouldn't be caught by the "pileup exists" skip-check alone.

**If resuming, read:** this entry (job IDs above), then `qstat -u terencew` and
`ambisim/logs/A01a_run_simulation.14732184.*` / `A01b_run_cr_arc.*` / `A02a_demuxlet_pileup.*`
for job status, before `README.md` for full scientific context.
