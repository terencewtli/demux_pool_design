# PROGRESS — operational state and next steps

Operational companion to `RESULTS.md` (findings) and `JOURNAL.md` (status board + log). This file
is the nitty-gritty: what has actually been produced, x/y counts, what is missing or stale, and
**the one place next steps live**. Update it whenever jobs land; keep it terse. Full bug write-ups
are in `docs/pipeline_issues.md`.

Last updated: 2026-09-18 (counts checked on disk, not from logs)

## Output inventory

Counts are pools with the output present for **both** modalities where applicable.

| tree | stage | count | status |
|---|---|---|---|
| `ambisim/` (chr22 matrix) | demuxlet `.best` | 128/132 | missing: `AFR_only__greedy_maxkl__rep1`, `EUR_EAS__ancestry_balanced__rep3`, `all_3_major__ancestry_balanced__rep3`, `EUR_AFR__adversarial_mindist__rep1` |
| `ambisim/` | scored in `01a-01d` | 88/132 | **stale** — notebooks last run at n=88; `RESULTS.md` numbers are from that run |
| `ambisim_n16/` (n=16) | demuxlet `.best` | 4/4 | complete; **not analyzed** (no results notebook) |
| `ambisim_new/` (orthogonal pilot) | demuxlet `.best` | 3/4 | `EUR_only__highmin_lowmean_new__rep2` has no `.best`; **not analyzed** |
| `ambisim_final/` (genome-wide) | pool VCFs (`A00`) | 126/126 | complete; +6 family pools symlinked from `ambisim/` |
| `ambisim_final/` | cr_arc (`A01b`) | 118/126 | 8 pools with no cr_arc output (list below) |
| `ambisim_final/` | popscle pileup (`A02a`) | 113/126 | 113 both; 1 ATAC-only; 4 never submitted; 8 blocked on cr_arc |
| `ambisim_final/` | demuxlet call (`A02b`) | 0/126 | **not submitted** — ready for the 113 |
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

- `RESULTS.md` numbers are from the chr22-matrix grid at n=88 — replace, don't average, once
  `ambisim_final/` is scored.
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

- demuxlet call needs ≥32G (`-l h_data=32G`); the shared template's 8G is killed silently.
- `A02a`/`A02b` submit from array tasks → some nodes reject with "no submit host"; task still
  exits 0. Grep logs for it after every run.
- `A02a` exits in seconds; `-hold_jid` on it releases `A02b` too early (A02b then skips everything).
- Count readiness per pool with both `.best` files, never `ls */*.best | wc -l`.
- popscle child-job logs land in `scripts/<tree>/qsub/logs/`, not `<tree>/logs/`.
- SGE `-t` takes one range; submit non-contiguous task IDs separately.
- cr_arc needs `h_rt=24:00:00` (14 h killed 26/27 pipestances).
- Distance thresholds are matrix-specific: chr22 Euclidean ~0–150 vs. genome-wide 1-IBS ~0.10–0.29.
- `nominate()` scripts write output only at the end — check the strategy dispatcher first.
- Every pool-level correlation: adversarial-excluded + Spearman (`docs/CONVENTIONS.md`).

## Next steps (checked 2026-09-18)

**Operational — `ambisim_final/` (in order):**

1. Submit demuxlet calls for the 113 ready pools:
   `cd scripts/ambisim_final/qsub && qsub A02b_demuxlet_call.sh` (skips pools lacking pileups).
   Afterwards grep `ambisim_final/logs/A02b_*` for `no submit host` and call logs for `Killed`.
2. Resubmit `A02a` for the 5 pileup gaps (4 never-submitted + the GEX kill); skip-if-exists keeps
   it safe. Then `A02b` again for those.
3. Rerun `A01b` for the 8 pools with no cr_arc output — check their A01a/A01b logs for the cause
   first. Then `A02a` → `A02b`.
4. Point `notebooks/ambisim/01a-01d` at `ambisim_final/` and score; rewrite `RESULTS.md` from that
   run (genome-wide numbers replace chr22 ones).

**Analyses (highest value first):**

1. **Orthogonal pilot** (`ambisim_new/`, 3/4 demuxed): does `highmin_lowmean` differ from
   `lowmin_highmean` in LL-gap/accuracy? First real test of min vs. mean distance. Needs a small
   dedicated notebook (01a-01d are wired to `ambisim/`). Finish the missing rep2 call first.
   Gate: if a real difference → scale the 36 genome-wide orthogonal pools and add the
   `highmin_highmean`/`lowmin_lowmean` quadrants; if not → spend the compute on (2).
2. **Downsampling / stress test** (`md/analyses.md` Tier 3): subsample reads/cells, re-run
   demuxlet, see whether narrow-margin pools lose accuracy first. The most direct test of whether
   margin ever converts to accuracy for ordinary (non-relative) pools.
3. **n=16** (`ambisim_n16/`, 4/4 demuxed): score and compare to n=8 `EUR_only` random/maxmin.
4. Eyes-on check of the 7 relatedness-flagged genome-wide orthogonal pools before any are
   simulated.
5. More `adversarial_family`/`_mixed` replicates if the relatedness effect (n=6) needs shoring up.
6. Why real pools decouple min/mean (r=0.14) more than random draws (r=0.54) — only if (1) shows
   the distinction matters.
7. Unresolved, low priority: `greedy_maxkl` root cause; whether GEX pool-average keeps favoring
   `mean_dist` while ATAC favors `min_dist`.
