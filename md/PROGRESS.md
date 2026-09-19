# PROGRESS — operational state and next steps

Operational companion to `RESULTS.md` (findings) and `JOURNAL.md` (status board + log). This file
is the nitty-gritty: what has actually been produced, x/y counts, what is missing or stale, and
**the one place next steps live**. Update it whenever jobs land; keep it terse. Full bug write-ups
are in `docs/pipeline_issues.md`.

Last updated: 2026-09-19 (counts checked on disk, not from logs)

## Output inventory

Counts are pools with the output present for **both** modalities where applicable.

| tree | stage | count | status |
|---|---|---|---|
| `ambisim/` (chr22 matrix) | demuxlet `.best` | 128/132 | missing: `AFR_only__greedy_maxkl__rep1`, `EUR_EAS__ancestry_balanced__rep3`, `all_3_major__ancestry_balanced__rep3`, `EUR_AFR__adversarial_mindist__rep1` |
| `ambisim/` | scored in `01e-01i` | 128/132 | current; `RESULTS.md` is written from this run (2026-09-19) |
| `ambisim/` | scored in `01a-01d` | 88/132 | superseded by `01e-01i`; `droplet_scores.csv` kept, not overwritten |
| `ambisim_n16/` (n=16) | demuxlet `.best` | 4/4 | complete; **not analyzed** (no results notebook) |
| `ambisim_new/` (orthogonal pilot) | demuxlet `.best` | 3/4 | `EUR_only__highmin_lowmean_new__rep2` has no `.best`; **not analyzed** |
| `ambisim_final/` (genome-wide) | pool VCFs (`A00`) | 126/126 | complete; +6 family pools symlinked from `ambisim/` |
| `ambisim_final/` | cr_arc (`A01b`) | 118/126 | 8 pools with no cr_arc output (list below) |
| `ambisim_final/` | popscle pileup (`A02a`) | 113/126 | 113 both; 1 ATAC-only; 4 never submitted; 8 blocked on cr_arc |
| `ambisim_final/` | demuxlet call (`A02b`) | running | job **14809580** submitted 2026-09-18 22:26, 126-task array, `-tc 20`; 6/126 complete as of 01:00 on 09-19, no OOM |
| design | genome-wide distance matrix | — | `csv/designs/1kg_genomewide_dist_matrix.npz` (459,999 SNPs, 3,202 samples) |
| design | regenerated 132-pool donor lists | 132 | `txt/donors_genomewide/`, `nominated_pools_n8_genomewide.tsv` (job 14753496) |
| design | genome-wide orthogonal pools | 36 | `txt/donors_new_genomewide/` — nominated only, not simulated; 7 flagged by the relatedness safeguard |

`ambisim_final/` pools not yet through pileup:

- **No cr_arc output (rerun `A01b`, then `A02a`):** `AFR_only__greedy_maxmean__rep1`,
  `EUR_EAS__random__rep1`, `EAS_only__random__rep2`, `all_3_major__greedy_maxmin__rep2`,
  `all_5__random__rep2`, `AFR_only__adversarial_mindist__rep3`,
  `SAS_only__adversarial_mindist__rep1`, `EUR_AFR__adversarial_mindist__rep2`.
- **cr_arc done, pileup never submitted ("no submit host"):** `EUR_only__greedy_maxmin__rep1`,
  `SAS_only__random__rep1`, `EUR_only__random__rep2`, `SAS_only__greedy_maxmean__rep2`.
- **GEX pileup killed at start (job 14780410):** `EUR_AFR__adversarial_mindist__rep3`.

## Known-stale / do-not-use

- `RESULTS.md` numbers are from the chr22-matrix grid at n=128 — replace, don't average, once
  `ambisim_final/` is scored.
- **Any "LL-gap" number written before 2026-09-19** is `DIFF.LLK.BEST.NEXT`, which for an
  SNG-called droplet is the singlet-vs-doublet margin, not the donor margin (`RESULTS.md` §0).
  This includes `results/pool_summary.csv`, `results/bottleneck_donor_llgap.csv` and
  `results/droplet_ambient_stratified.csv`, and the `01a-01d` notebooks that produced them.
- `ambisim/` donor selections (all of them) were made on the chr22-only matrix — fine as a
  design-agnostic replicate set, not as "the" extremal pools.
- `greedy_maxkl` pools in every tree — unreliable nomination.
- `results/nominated_pools_n8.tsv` — chr22 matrix; genome-wide equivalent is
  `csv/designs/nominated_pools_n8_genomewide.tsv`.
- `scripts/ambisim/lib/build_genomewide_dist_matrix.py` (reservoir sample, unpruned) — superseded
  by `build_genomewide_dist_matrix_plink.py`.
- `scripts/ambisim/lib/create_pools_genomewide.py`, `run_greedy_maxkl_genomewide.py`,
  `B03_greedy_maxkl_genomewide.sh` — superseded by `regenerate_pools_genomewide.py` / `B04`.
- `scripts/ambisim_final/qsub/core.{42350,44550,46674}` — 28 GB of core dumps, no value.
- `ambisim/` BAMs for demuxlet-complete pools were removed by `pool_design/tmp.sh`
  (`ambisim/txt/demuxlet_complete_pools_20260915.txt`) — re-running pileup there requires re-running cr_arc.

## Recurring gotchas (cost time at least once — details in `docs/pipeline_issues.md`)

- demuxlet call needs ≥32G; `A02b` now asks for 40G as `h_data=4G` x `-pe shared 10`. The
  shared template's 8G default is killed silently.
- `A02b` runs demuxlet **inline** as of 2026-09-18; it no longer fans out child `qsub`s, so the
  "no submit host" trap below applies to `A02a` only. ~8 min/modality against a 1.3 GB pool VCF.
- `A02a` submits child jobs from array tasks → some nodes reject with "no submit host"; task
  still exits 0. Grep logs for it after every `A02a` run. (`A02b` no longer does this.)
- `A02a` exits in seconds; `-hold_jid` on it releases `A02b` too early (A02b then skips everything).
- Count readiness per pool with both `.best` files, never `ls */*.best | wc -l`.
- popscle child-job logs land in `scripts/<tree>/qsub/logs/`, not `<tree>/logs/`.
- SGE `-t` takes one range; submit non-contiguous task IDs separately.
- cr_arc needs `h_rt=24:00:00` (14 h killed 26/27 pipestances).
- Distance thresholds are matrix-specific: chr22 Euclidean ~0–150 vs. genome-wide 1-IBS ~0.10–0.29.
- `nominate()` scripts write output only at the end — check the strategy dispatcher first.
- Every pool-level correlation: adversarial-excluded + Spearman (`docs/CONVENTIONS.md`).

## Next steps (checked 2026-09-19)

**Operational — `ambisim_final/` (in order):**

1. Let job 14809580 finish, then count pools with **both** `.best` and grep the logs for `Killed`
   (`grep -ail killed ambisim_final/logs/A02b_*` — note demuxlet's own `--geno-error-offset`
   banner makes a naive `error` grep useless).
2. Resubmit `A02a` for the 5 pileup gaps (4 never-submitted + the GEX kill); skip-if-exists keeps
   it safe. Then `A02b` again for those.
3. Rerun `A01b` (cellranger-arc) for the 8 pools with no cr_arc output — **not submitted, not
   running as of 2026-09-19 01:30**; those 8 stay blocked until it is. Check their A01a/A01b logs
   for the cause first, then `A02a` → `A02b`. **Check disk first:** `/u/project/cluo` is at
   **99%** (10 T free of 546 T) and a cr_arc pipestance is ~83 GB/pool before BAM cleanup, so 8
   pools is ~660 GB. Stage them, and run `pool_design/tmp.sh`-style BAM cleanup on
   demuxlet-complete pools before starting.
4. Point `notebooks/ambisim/01e-01i` at `ambisim_final/` and rescore; rewrite `RESULTS.md` from
   that run (genome-wide numbers replace chr22 ones). `01e` has a `tree` variable; the pool
   list path also needs pointing at `ambisim/txt/pool_experiments_final126.txt`.

**Analyses (highest value first, reprioritized 2026-09-19):**

1. **Real-data test of the doublet result** (`RESULTS.md` §5) — the highest-value item, and the
   only one that puts a simulation finding in front of real data. In the ~40 reprogramming pools:
   compute each pool's `mean_dist`, and its *called* doublet rate against the rate predicted from
   loading concentration. Prediction: lower-`mean_dist` pools under-call doublets relative to
   loading. Not explainable by sequencing depth or wet-lab variation, which is exactly the
   "signal beyond coverage" the project needs.
2. **Orthogonal pilot** (`ambisim_new/`, 3/4 demuxed): does `highmin_lowmean` differ from
   `lowmin_highmean`? Gates `RESULTS.md` §2's "mean beats min" and §8. Needs a small dedicated
   notebook. Finish the missing rep2 call first.
3. Score `ambisim_final/` once it lands (operational 4) and confirm §1–§6 survive genome-wide
   donor selections. §1/§3/§5 are droplet-level and should not move; §2/§6 are pool-level and may.
4. **n=16** (`ambisim_n16/`, 4/4 demuxed): score and compare to n=8 `EUR_only` random/maxmin.
   More donors should push more droplets below the §3 threshold — a direct test of it.
5. More `adversarial_family`/`_mixed` replicates if the relatedness effect needs shoring up
   (n=3 per arm; `RESULTS.md` §6 is directionally clear but thin).
6. Eyes-on check of the 7 relatedness-flagged genome-wide orthogonal pools before any are
   simulated.
7. **Dropped:** a relatedness-gradient simulation (half-sibs / avuncular / cousins). `01i` puts
   all 994 donors on one nearest-pool-mate axis over which the error rate runs 4.1 -> 36.5 per
   10,000 monotonically, with the parent-offspring pairs at the bottom of the *same* curve. New
   degrees would interpolate between measured points. `RESULTS.md` §7 is the figure that answers
   a reviewer asking for the gradient.
8. **Dropped:** the read-downsampling stress test (was #2). The natural low-coverage tail already
   answers it — among ordinary pools `mean_dist` predicts accuracy in no coverage stratum, and
   projecting a 33x margin shrink through the shared §3 curve leaves mean donor accuracy above
   0.993. It would sharpen §6 (relatives), not rescue §1. Not worth the compute.
9. Unresolved, low priority: `greedy_maxkl` root cause; why ATAC's singlet call tracks `mean_dist`
   (ρ=+0.450) while GEX's does not (ρ=−0.148) — `RESULTS.md` §5.
