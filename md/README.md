# Pool design genotype-based demultiplexing

Does the genetic composition of a multiplexed single-cell pool (which donors you
choose to pool together) affect how well genotype-based demultiplexing (demuxlet)
can tell them apart? 

## Where to read what

| file | contents |
|---|---|
| `RESULTS.md` | **What we know.** Current validated findings only, organized by question. No history. |
| `PROGRESS.md` | **What has run and what to do next.** Per-tree x/y counts, stale / do-not-use data, recurring gotchas, next-session checklist. |
| `docs/pipeline_issues.md` | **Bugs and fixes.** Each known pipeline failure: symptom, root cause, fix, status. |
| `docs/minmean_entanglement.md` | Methodology deep-dive: why `greedy_maxmin`/`greedy_maxmean` can't separate min from mean distance, and the orthogonal-sampling fix. |
| `docs/CONVENTIONS.md` | Operating manual: environments, script naming, array-job pattern, the adversarial-leverage check every correlation must pass. |
| `JOURNAL.md` | Status board (verified / retracted / open) and a short dated log. |
| `JOURNAL.archive.md` | Full session-by-session record through 2026-09-15, verbatim. |

## Pool construction

All pools draw n=8 donors from a "universe" of unrelated 1000 Genomes
individuals — either single-ancestry (`EUR_only`, `AFR_only`, `EAS_only`,
`SAS_only`, `AMR_only`) or multi-ancestry (`EUR_AFR`, `EUR_EAS`, `all_3_major`,
`all_5`). We derived pools using the following set of metrics

| strategy | objective |
|---|---|
| `random` | no optimization; baseline |
| `greedy_maxmin` | maximize the *minimum* pairwise genetic distance in the pool |
| `greedy_maxmean` | maximize the *mean* pairwise genetic distance |
| `greedy_maxkl` | like `greedy_maxmin` but using KL divergence instead of Euclidean distance (single-ancestry only) — **buggy, treat as unreliable, see `docs/pipeline_issues.md`** |
| `ancestry_balanced` | proportional sampling across sub-populations, no distance optimization (multi-ancestry only) |
| `adversarial_mindist` | mirror of `greedy_maxmin` — *minimizes* the pool's spread (worst-case floor) |
| `adversarial_family` / `adversarial_family_mixed` | a real 1000G pedigree (~50% IBD), padded with unrelated fill — a non-synthetic worst case |

Each cell (universe × strategy) has 3 replicates; the full design is 132 pools
(`txt/pool_experiments.txt`).

## Simulation trees

| tree | distance matrix used to pick donors | status |
|---|---|---|
| `ambisim/` | chr22 only, 20k SNPs | 128/132 demuxed; source of every number in `RESULTS.md` (analyzed at n=88) |
| `ambisim_final/` | genome-wide, LD-pruned (459,999 SNPs, 1-IBS) | **in flight** — supersedes `ambisim/` once it lands |
| `ambisim_n16/` | chr22 | n=16 side experiment (does donor count matter?); demuxed, not yet analyzed |
| `ambisim_new/` | chr22 | orthogonal-sampling pilot (min and mean distance deliberately disagree); demuxed, not yet analyzed |

Pipeline per pool: `A00` prep pool VCF → `A01a` ambisim → `A01b` cellranger-arc → `A02a` popscle
pileup (GEX, ATAC) → `A02b` demuxlet call → `notebooks/ambisim/01a-01d` scoring.
