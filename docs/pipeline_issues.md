# Pipeline issues — bugs, root causes, fixes

One entry per known failure, organized by topic rather than date. Check here first when a
number looks wrong or a job "succeeds" with nothing to show for it. The terse version of each
lives in `PROGRESS.md` → Recurring gotchas; the narrative of when each was found is in
`JOURNAL.archive.md`.

Status: **fixed** (code changed, verified) · **worked around** (avoided per run, root cause
remains) · **open** (not resolved) · **by design** (a constraint to work within).

---

## Silent job failures (the dangerous kind)

### demuxlet call OOM-killed with no error — worked around

- **Symptom:** a pool has both pileups (`.var.gz`) but only one modality's `.best`. First seen as
  `EUR_EAS__greedy_maxmean__rep1` (2026-09-12); 27/132 calls affected in the n=8 grid.
- **Cause:** the shared `demux_benchmark/template_demux/demuxlet.sh` requests `h_data=4G, -pe
  shared 2` (8 GB). Whichever modality has the larger pileup (~3.8M–6.9M variants; not
  consistently GEX or ATAC) is `Killed` while loading it. `A02b` fire-and-forgets both `qsub`s and
  never checks the outcome.
- **Fix:** resubmitted with `qsub -l h_data=24G` (3 pools sharing a ~6.87M-variant VCF needed
  64G). `ambisim_n16/` and `ambisim_final/` `A02b` pass `-l h_data=32G` from the start. The shared
  template is **not** edited (used by other pipelines).
- **Detection:** grep call logs for `Killed`; count pools with *both* `.best` files, never
  `ls */*.best | wc -l` (overcounts lopsided pools). `01a` has a readiness cell that prints
  ready/partial/not-started by name.

### Array task can't submit child jobs: "no submit host" — open

- **Symptom (2026-09-17, `ambisim_final` A02a):** 4 pools have cr_arc output but no pileup
  directory contents. The A02a array log shows
  `Unable to run job: denied: host "n7168" is no submit host.` and the task exits cleanly.
- **Cause:** `A02a`/`A02b` call `qsub` from inside an array task; some compute nodes are not SGE
  submit hosts. The script still prints `End:` and exits 0.
- **Affected:** `EUR_only__greedy_maxmin__rep1`, `SAS_only__random__rep1`,
  `EUR_only__random__rep2`, `SAS_only__greedy_maxmean__rep2`.
- **Fix:** resubmit those task IDs (skip-if-exists makes it safe). The same failure can hit
  `A02b`: after it runs, `grep -l "no submit host" ambisim_final/logs/A02b_*`. A durable fix would
  be to submit pileup/call jobs from the login node (a loop over the pool list) instead of from
  array tasks.

### popscle pileup killed at start — open

- `EUR_AFR__adversarial_mindist__rep3` GEX pileup (job 14780410, 2026-09-18 04:15, node n1888):
  exit 137 after 37 s, maxvmem 167 MB, log has only the start banner. Not OOM; looks like a node
  kill. ATAC for the same pool finished. Resubmit.

### Skip-if-exists lives in the shell script

`A02a` skips on `<sample>.pileup.var.gz`, `A02b` on `<sample>.best`. If an output path changes in
the tool/template, change the shell check too, or every task silently skips.

### A02a exits in seconds; it can't be `-hold_jid`'d — by design

`A02a` only *submits* `A01a_popscle_gex`/`A01a_popscle_atac` per pool and exits. Chaining `A02b`
with `-hold_jid` on `A02a` releases it before any pileup exists (A02b then skips every pool,
exit 0). Wait for `qstat` to clear the popscle jobs, or hold on those job *names*.

---

## Wall time and resources

### cellranger-arc 14 h wall-time cap — fixed

26 of 27 incomplete pipestances (2026-09-14) had been killed mid-run by `h_rt=14:00:00`
(confirmed from the pipestance `_log` timestamps). One cr_arc job hung 17.5 h past its own
deadline and was killed by hand. `A01b` now uses `h_rt=24:00:00`.

### SGE `-t` accepts one range only — by design

`qsub -t 1,16,36` fails (`-t option only allows one range specification`). Submit one job per
non-contiguous task ID.

### Genome-wide distance matrix OOM, twice — fixed

`build_genomewide_dist_matrix_plink.py` died silently at maxvmem 338G (128G requested), from two
unrelated causes:
1. `sklearn.PCA(svd_solver='randomized')` on the 3202 × 460k dosage array → replaced with
   `flashpca` on a pruned `.bed` (3m37s).
2. `pandas.read_csv` on the wide (3202 × 460k-column) PLINK `.raw` → per-column overhead
   regardless of dtype. Replaced by parsing lines into a preallocated `float32` array
   (`np.fromstring`), with `--output-missing-genotype 0` to avoid `NA` strings.

After both: 6m54s end to end. Lesson: don't assume the first fix was the right one — the
resubmitted job was watched (`qstat -j`) past the old crash point.

### Core dumps in `scripts/ambisim_final/qsub/` — open

Three ~9.4 GB `core.*` files (28 GB total, 2026-09-16 03:02–03:20) in the script directory — some
job launched with that directory as cwd crashed. Source not identified. Not needed for anything;
candidates for deletion by the user.

---

## Pool nomination / distance matrix

### Distance matrix was chr22-only — fixed (by resimulation)

Every `min_dist`/`mean_dist`/`kl_*` value used to *select* pools in `ambisim/`, `ambisim_n16/` and
`ambisim_new/` came from chr22 subsampled to 20,000 SNPs (MAF ≥ 0.05) —
`ipynb/vcf_metrics/local/02a_distance_matrix.ipynb`, a prototyping shortcut never revisited.
Replaced by a genome-wide LD-pruned matrix (`--indep-pairwise 50 5 0.2`, 459,999 SNPs, 1-IBS;
`csv/designs/1kg_genomewide_dist_matrix.npz`). Donor lists regenerated
(`regenerate_pools_genomewide.py`) and the grid resimulated as `ambisim_final/`. The simulation
itself was never affected — it is donor-list-agnostic. Ancestry-level conclusions are expected to
hold; donor-level extremal selections are what changes.

### `nominate()` dispatcher missing `adversarial_mindist` — fixed

`nominate_greedy_mindist` existed but was never wired into `pool_nomination.py`'s `nominate()`;
job 14752950 died at pool 111/132 with `ValueError: unknown strategy: adversarial_mindist`.
The script writes output only after the loop, so **nothing** was saved. Fixed, rerun as
14753496. Before resubmitting any nomination run, check every strategy in the task list against
the dispatcher.

### Cryptic-relatedness threshold on the wrong scale — fixed

The orthogonal-sampling safeguard's "closest pair too close" threshold (`70`) was copied from the
chr22 script, whose distance was Euclidean-ish (~0–150). On the 1-IBS matrix (~0.10–0.29) it
flagged 18/18 pools — a 100% hit rate means miscalibration. Recalibrated to **0.16** from the
matrix's own truth: 1,205 real parent–child pairs at D = 0.1056–0.1543 vs. whole-matrix 1st
percentile 0.1984. Now flags 7/36, all `lowmin_highmean` (the known heavy-tailed quadrant).

### `unrelated` means pedigree founder, not kinship-checked — by design

`unrelated` = `FatherID=='0' & MotherID=='0'`. It excludes documented 1000G trios/quads, not
cryptic relatedness. Extremal searches (adversarial, `lowmin_highmean`) are the most exposed; the
safeguard above is the backstop.

### Greedy strategies share a seed per rep — by design (documented)

Both `nominate_greedy_*` use `np.random.default_rng(seed)` with the same seed per rep, so
`greedy_maxmin` and `greedy_maxmean` start from the identical first donor. One of two causes of
their 5–6/8 donor overlap; see `docs/minmean_entanglement.md`.

### `greedy_maxkl` — partially fixed, still broken, deprioritized

- **Bug (fixed):** the KL computation added its epsilon *before* taking `1 - p`, so any
  homozygous-alt site (p=1.0) gave `log(negative) = NaN`, corrupting the greedy argmax.
- **Still broken:** after the fix, canonical best-of-5-seed `kl_min` is still below the `random`
  baseline in all 5 single-ancestry universes (AFR 1.738 vs. 1.851, AMR 1.702 vs. 1.725, EAS
  1.498 vs. 1.675, EUR 1.639 vs. 1.773, SAS 1.799 vs. 1.880); unchanged with rep2/rep3.
  Search objective and scoring function use consistent epsilon handling, so it's not a
  search/score mismatch. Root cause unknown.
- **Decision:** don't use `greedy_maxkl` pools as a "good design" data point. `kl_min` and
  `min_dist` are r=0.999 collinear, so `greedy_maxmin` covers the role. The 3 pre-fix simulated
  rep1 pools were not re-simulated; rep2/rep3 were simulated anyway on the user's "run
  everything" call (2026-09-12). `regenerate_pools_genomewide.py` uses plain
  `seed = SEED + rep*100`, not the best-of-5 ceremony.

### Wrong-donor VCF for `AFR_only__greedy_maxkl__rep1` — fixed

The pool's simulation was blocked by a VCF built from the wrong donor samples; rebuilt and rerun
(2026-09-14). This pool is still one of the 4 chr22-grid pools without demuxlet output.

### Duplicate-effort collision on `nominated_pools_n8_genomewide.tsv` — resolved

Two scripts (`create_pools_genomewide.py`, narrower; `regenerate_pools_genomewide.py`, complete)
targeted the same output path. The narrower one is kept with a superseded-by note; its
`greedy_maxkl` job (14752676) was `qdel`'ed. Read `JOURNAL.md` before writing new pipeline code.
