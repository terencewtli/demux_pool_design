# Conventions — how this project is operated

Read with `PROGRESS.md` (what has run) and `RESULTS.md` (what we know). This file is the
operating manual: the things a new session has to know before touching anything.

## Paths

- Working directory: `/u/project/cluo/terencew/claude/project_ideas/pool_design`
  (simulation trees `ambisim*/`, task lists `ambisim/txt/`, designs `csv/designs/`).
- This repo (`github/demux_pool_design/`) is a hand-maintained mirror: scripts, notebooks,
  results tables and all `.md` docs. There is no sync script — edit the docs here.
- Shared templates (not in this repo, used by other pipelines — don't edit):
  `/u/project/cluo/terencew/demux_benchmark/template_demux/{popscle_gex,popscle_atac,demuxlet}.sh`.
- 1000G source VCF: `demux_benchmark/pool_design/vcf/1000G/1000G.merged.common_biallelic.vcf.gz`.
- `/u/project/cluo` is shared and ~95% full; a full grid is ~83 GB/pool before BAM cleanup.
  Stage large submissions.

## Environments

- Notebook scoring (`notebooks/ambisim/01a-01d`): `demux` kernel. `nbconvert` isn't installed in
  `demux` — run `allcools`'s `jupyter nbconvert --execute` pointed at the `demux` kernel.
- Nomination / distance matrix: `allcools` (Python). `claude_env` is Node.js, not Python.

## Script layout and naming

`scripts/<tree>/{lib,qsub}/`, one tree per simulation directory (`ambisim`, `ambisim_final`,
`ambisim_n16`, `ambisim_new`). Per-pool pipeline:

| step | script | output |
|---|---|---|
| `A00` | `A00_prep_pool_vcfs.sh` | `<tree>/vcf/<pool>.vcf.gz` |
| `A01a` | `A01a_run_simulation.sh` | ambisim FASTQs + `drop_data_rand.txt` truth |
| `A01b` | `A01b_run_cr_arc.sh` | `<pool>/cr_arc/<sample>/outs/` |
| `A02a` | `A02a_demuxlet_pileup.sh` (submits popscle GEX + ATAC) | `<pool>/demux/demuxlet/{gex,atac}/<sample>.pileup.*` |
| `A02b` | `A02b_demuxlet_call.sh` (submits demuxlet GEX + ATAC) | `<pool>/demux/demuxlet/{gex,atac}/<sample>.best` |
| `B0*` | distance matrix / pool regeneration (`scripts/ambisim/qsub/`) | `csv/designs/` |

Sample name is always `20220928-IGVF-D0`. Pool IDs are `<universe>__<strategy>__rep<N>`.
Task lists: `ambisim/txt/pool_experiments.txt` (132, chr22 grid),
`pool_experiments_final126.txt` (`ambisim_final/`; the 6 family pools are symlinked from
`ambisim/`, not resimulated).

## Array-job pattern

`#$ -t 1-N`, `ID=$SGE_TASK_ID` (with a commented `# ID=1` for hand testing), skip-if-exists on
the step's own output, explicit `exit 0` for pools missing inputs, `time` before the main
command, logs to `<tree>/logs/<step>.$JOB_ID.$TASK_ID`. (The existing `ambisim_final/` scripts
lack the `# ID=1` line and `time`; add them to new scripts.) Child jobs submitted by `A02a`/`A02b`
log to the **cwd at submission** — for `ambisim_final/` that is
`scripts/ambisim_final/qsub/logs/A01a_popscle_{gex,atac}.$JOB_ID`, not `ambisim_final/logs/`.
Known traps in this pattern — `A02a`
can't be `-hold_jid`'d, child `qsub` can fail with "no submit host", demuxlet needs ≥32G — are in
`docs/pipeline_issues.md`.

## Analysis rules (each has cost us a retraction)

1. **Adversarial-leverage check.** Compute every pool-level correlation twice: full sample and
   with `strategy.str.startswith('adversarial')` excluded. If the excluded version isn't
   significant or flips sign, the full-sample number describes extreme cases, not a trend — say
   so. Two whole-sample claims were retracted on 2026-09-14 for failing this.
2. **Spearman alongside Pearson.** A Pearson r without a similar Spearman rho is leverage (GEX
   accuracy: Pearson r=0.349 p=0.001 vs. Spearman rho=0.205 p=0.056).
3. **No correlations from a handful of pools.** The `random`-pool "r=−0.09, matches real data"
   claim came from 7 pools; 20,000 draws gave r=0.536. Use a distance-matrix-only null (cheap)
   whenever the question is about pool geometry rather than demux outcome.
4. **Count readiness per pool, both modalities.** Never `ls */*.best | wc -l`.
5. **Don't cite greedy_maxmin vs. greedy_maxmean as a min-vs-mean test** (`RESULTS.md` §6).
6. **Treat `greedy_maxkl` as unreliable** (`docs/pipeline_issues.md`).

## Where things go

- New findings → `RESULTS.md` (current state only; rewrite in place, don't append history).
- Job counts, stale data, next steps → `PROGRESS.md`.
- A new bug or failure mode → `docs/pipeline_issues.md` (full entry) + one line in
  `PROGRESS.md` → Recurring gotchas.
- What happened this session, including anything retracted → a short dated entry in
  `JOURNAL.md` → Log, and move the claim on the status board.
