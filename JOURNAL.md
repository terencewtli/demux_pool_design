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
