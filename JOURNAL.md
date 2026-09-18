# Journal

The record of *why* the current results can be trusted: which claims were checked, which were
retracted, and what happened when. Current findings are in `RESULTS.md`; operational state and
next steps are in `PROGRESS.md`; bug write-ups are in `docs/pipeline_issues.md`.

The full session-by-session record through 2026-09-15 (every wrong turn, job ID and decision) is
in `JOURNAL.archive.md`, verbatim. This file keeps only:

1. **Status board** — what has been checked and can be trusted, what's retracted, what's open.
2. **Log** — short dated entries, newest first. Long narrative goes in the archive, not here.

Labels: **[verified]** = rechecked directly against data, with the check described;
**[reported]** = produced by an earlier session, not independently rechecked since;
**[retracted]** = shown wrong, kept so nobody reuses it.

---

## 1. Status board (as of 2026-09-18)

### Findings

| claim | label | basis |
|---|---|---|
| Distance-based curation of unrelated donors does not move singlet accuracy | [verified] | n=73 ordinary pools, r≈0 both modalities, holds per ambient bin and ancestry subset (2026-09-14/15) |
| `min_dist` predicts LL-gap margin, r≈0.31 both modalities | [verified] | adversarial-excluded; replicates n=32 direction |
| Family pools lose accuracy (~1–2 pts) and margin (~30%), localized to related donors | [verified] | 6 pools, per-donor split on each pool's own matrix |
| `adversarial_mindist` has a weaker accuracy/margin deficit | [reported] | 9 pools, small effect; thin |
| Design's margin benefit fades with ambient RNA | [verified] | per-bin correlations; depth ruled out as confound |
| Multi-ancestry universes show the strongest margin signal (r≈0.56–0.59) | [reported] | n=29; thin |
| greedy_maxmin/maxmean pools share 5–6/8 donors; min/mean r=0.69–0.84 (chr22), 0.992 (genome-wide maxmin) | [verified] | donor overlap and r rechecked on the genome-wide matrix 2026-09-15 |
| Random 8-donor draws have min/mean r=0.536 | [verified] | 20,000 draws, `EUR_only` |
| Real 10x pools have min/mean r=0.14 | [reported] | 63 pools, `pool_genetic_dists.csv` |
| **Distance-based design predicts accuracy (r=+0.35 GEX / +0.56 ATAC)** | **[retracted]** | whole-sample n=88; vanishes with adversarial pools excluded (2026-09-14) |
| **`min_dist` beats `mean_dist` at the bottleneck donor** | **[retracted]** | same leverage artifact; ordinary pools favor `mean_dist` (and that comparison is itself confounded, `RESULTS.md` §6) |
| **Random pools decouple min/mean (r=−0.09), "matching real data"** | **[retracted]** | 7 pools; well-powered null is r=0.536 |
| **`greedy_maxkl` achieves higher `kl_min` than random** | **[retracted]** | below random in all 5 universes, before and after the NaN fix |

### Pipeline / data

| item | label | note |
|---|---|---|
| Genome-wide LD-pruned distance matrix | [verified] | built 2026-09-15; parent–child pairs cleanly separated from the rest |
| Regenerated 132-pool donor lists (genome-wide) | [reported] | job 14753496 after the dispatcher fix; used by `ambisim_final/` |
| `ambisim_final/` stage counts | [verified] | on-disk check 2026-09-18, see `PROGRESS.md` |
| demuxlet 32G memory is sufficient | [reported] | held for `ambisim_n16/`; unconfirmed at 126-pool scale |

### Open questions

1. Does min vs. mean distance matter for demux, once they're decoupled? (orthogonal pilot)
2. Does margin ever convert to accuracy under stress for ordinary pools? (downsampling)
3. Does donor count (n=16) change any of the above?
4. Do the n=88 chr22 findings hold on the genome-wide `ambisim_final/` grid?

Priorities and gating for all four: `PROGRESS.md` → Next steps.

---

## 2. Log (newest first)

Entry template:

```
## YYYY-MM-DD — one-line summary
**Did:** … **Found:** … **Changed on the status board:** … **Next:** (pointer into PROGRESS.md)
```

### 2026-09-18 — docs restructured; `ambisim_final/` pileups checked

**Did:** Split the docs into results / progress+next steps / documentation, following the
`asm_lr_hprc2` layout. `README.md` is now the project framing and a "where to read what" map.
Findings moved to `RESULTS.md`. `NEXT_STEPS.md` (stale since 2026-09-06) was replaced by
`PROGRESS.md` → Next steps. `NOTES.md` was split into `docs/pipeline_issues.md`,
`docs/minmean_entanglement.md` (verbatim) and `docs/CONVENTIONS.md`. The old `JOURNAL.md` is now
`JOURNAL.archive.md`.
**Found:** `ambisim_final/` was launched after the last archive entry (the A00 header says
"approved to run 2026-09-15"; the approval and submission weren't journaled). On disk: 118/126 have
cr_arc output, 113/126 have both pileups, 0 demuxlet calls. 4 pools never got pileups because
A02a's child `qsub` hit "no submit host". 1 GEX pileup was killed 37 s in. Also: the chr22 grid is
now at 128/132 demuxed (still analyzed only at 88), `ambisim_n16/` is 4/4 and `ambisim_new/` 3/4
demuxed, neither analyzed. Three 9.4 GB core dumps from 2026-09-16 sit in
`scripts/ambisim_final/qsub/`.
**Next:** `PROGRESS.md` → Next steps, operational 1–3.

### 2026-09-15 — genome-wide matrix; entanglement re-confirmed; `ambisim_final/` written

LD-pruned PLINK matrix built after two OOM fixes. min/mean entanglement is worse genome-wide
(r=0.992). The orthogonal-sampling fix was generalized to 9 universes (relatedness threshold
recalibrated 70 → 0.16). The `nominate()` dispatcher bug was fixed and the pools regenerated
(14753496). `ambisim_final/` was written, gated.
Archive: 2026-09-15 entry + addendum.

### 2026-09-14 (evening) — distance matrix found to be chr22-only; full grid redo decided

Every pool selection so far used chr22/20k SNPs. The user decided to rebuild genome-wide and
resimulate the grid. `unrelated` turns out to mean pedigree founder only. Archive: 2026-09-14
(evening).

### 2026-09-14 (later) — min/mean entanglement found; n=16 and orthogonal pilots launched

Retracted the random-pool r=−0.09 claim (r=0.536 at 20k draws). Built rejection-sampling
orthogonal pools. Localized the family-pool deficit to related donors. Archive: 2026-09-14
(later).

### 2026-09-14 — demuxlet OOM root-caused; n=88 analysis; two claims retracted same day

8G demuxlet template silently killed 27 calls. Reran 01a-01d at n=88. The whole-sample
accuracy and bottleneck claims were retracted after the adversarial-exclusion check; the
adversarial-leverage rule was adopted. Archive: 2026-09-14.

### 2026-09-12 — remaining 96 pools submitted

Full remaining grid launched on the user's "run everything" call (including `greedy_maxkl`
reps). Archive: 2026-09-12.
