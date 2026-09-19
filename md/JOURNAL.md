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

## 1. Status board (as of 2026-09-19)

### Findings

| claim | label | basis |
|---|---|---|
| Distance-based curation of unrelated donors does not move singlet accuracy | [verified] | n=105 ordinary pools at n=128; replicates on the 40 held-out pools (2026-09-19) |
| "Singlet accuracy" is ~entirely the SNG/DBL call, not donor identity | [verified] | accuracy = P(SNG)xP(donor\|SNG) = 0.8419x0.99898 GEX; sd 0.0107 vs 0.0004 across pools |
| `DIFF.LLK.BEST.NEXT` is the singlet-vs-doublet margin, not the donor margin | [verified] | NEXT.GUESS is a doublet model in 6705/6706 SNG droplets; == DIFF.LLK.SNG.DBL |
| Donor accuracy saturates at `sng_gap` ≥ 20, and `sng_gap` is a sufficient statistic | [verified] | shared curve across strategies to the 4th decimal; pool means 58-84 |
| Closer pools under-call true doublets | [verified] | P(DBL\|true doublet) ~ mean_dist rho=+0.600 GEX / +0.437 ATAC, n=105; monotone within-pool by donor-pair distance |
| Relatedness costs ~148 bad droplets per 10,000, ancestry-matched | [reported] | 3 family pools vs 12 ordinary EUR_only; 2.9x on wrong-donor rate |
| `min_dist` predicts the *doublet* margin, r≈0.32 both modalities | [verified] | adversarial-excluded, n=105; this is DIFF.LLK.BEST.NEXT, see above |
| `mean_dist` predicts the *donor* margin about 2x better than `min_dist` | [reported] | r=+0.706/+0.728 vs +0.341/+0.363, n=105; gated by the min/mean entanglement |
| Family pools lose accuracy and margin | [verified] | 6 pools; ancestry-matched 2.9x on wrong-donor rate |
| Family-pool deficit is localized to the related donors | [verified] | rechecked 2026-09-19 against the full pedigree: 36.6 vs 8.0 wrong donors/10k within the same pools; 62.7% of wrong calls name a first-degree relative |
| Relatedness tested is parent-offspring only (~0.5 IBD), n=3 pedigrees | [verified] | CEPH 1420/1463/1362, each a couple + all four parents; no sibs, no second-degree, no gradient |
| Relatedness is the bottom of one continuous nearest-pool-mate curve, not a separate mode | [verified] | 994 donor-level points, monotone over a 9x error range, no discontinuity at the related donors (`01i`) |
| `adversarial_mindist` has a weaker accuracy/margin deficit | [reported] | 9 pools, small effect; thin |
| **Design's margin benefit fades with ambient RNA** | **[retracted]** | true of the doublet margin only; on `sng_gap` it is flat (rho=+0.71/+0.72/+0.72/+0.71 GEX across ambient bins) (2026-09-19) |
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
| demuxlet 32G memory is sufficient | [reported] | `A02b` now requests 40G; no OOM in the first 13 tasks of job 14809580 |
| `A02b` runs demuxlet inline (one array, no child qsubs) | [verified] | rewritten 2026-09-18, dry-run against the old command line; job 14809580 completing ~8 min/modality |

### Open questions

1. Does the doublet-detection effect show up in the ~40 real reprogramming pools? (called
   doublet rate vs rate predicted from loading, regressed on `mean_dist`) — now the top priority.
2. Does min vs. mean distance matter for demux, once they're decoupled? (orthogonal pilot)
3. Does donor count (n=16) change any of the above?
4. Do the n=128 chr22 findings hold on the genome-wide `ambisim_final/` grid?
5. Why does ATAC's singlet call track `mean_dist` (rho=+0.450) while GEX's does not (-0.148)?

**Closed:** "does margin ever convert to accuracy under stress for ordinary pools?" — no, within
this regime. The natural low-coverage tail shows no `mean_dist`~accuracy signal in any stratum,
and a 33x margin shrink through the shared curve leaves donor accuracy above 0.993. The
downsampling run was dropped (`PROGRESS.md` → Next steps 7).

Priorities and gating: `PROGRESS.md` → Next steps.

---

## 2. Log (newest first)

Entry template:

```
## YYYY-MM-DD — one-line summary
**Did:** … **Found:** … **Changed on the status board:** … **Next:** (pointer into PROGRESS.md)
```

### 2026-09-19 — the two margins separated; doublets scored for the first time; n=88 -> n=128

**Did:** Rescored all 128 demuxed `ambisim/` pools carrying per-droplet coverage
(`NUM.SNPS`/`NUM.READS`), `DROPLET.TYPE` and `DIFF.LLK.SNG.DBL`, which `01a` had dropped
(`01e`). Re-ran the design correlations on the corrected axes with all three design metrics
(`01f`, `01g`). Scored the 900 true doublets per pool that no notebook had ever touched (`01h`).
Rewrote `A02b_demuxlet_call.sh` to run demuxlet inline and submitted it for `ambisim_final/`
(job 14809580).

**Found:** (1) `DIFF.LLK.BEST.NEXT` — "the LL-gap" everywhere in this repo — is the
singlet-vs-doublet margin, not the donor margin; the donor margin is `SNG.BEST.LLK -
SNG.NEXT.LLK`. Sections 2 and 4 of `RESULTS.md` had been written on the wrong column.
(2) Singlet accuracy is ~entirely `P(called SNG)`; donor identity contributes sd=0.0004 of
0.0107. (3) On the donor margin, `mean_dist` beats `min_dist` about 2:1, and the effect is flat
across coverage and ambient. (4) Donor accuracy saturates at `sng_gap` >= 20 and the curve is
shared across strategies, so pools differ only in where they sit on it — ~8 wrong donors per
10,000 separate the best and worst ordinary pools. (5) Closer pools systematically under-call
true doublets (rho=+0.600 GEX), monotone within-pool by donor-pair distance, and in GEX this runs
*opposite* to the singlet-call effect — which is what makes net accuracy look flat. (6) The
relatedness penalty is 2.9x ancestry-matched, not 3.4x — the earlier figure compared EUR family
pools against an all-universe ordinary mean. The deficit is localized to the related donors
(36.6 vs 8.0 wrong donors/10k inside the same pools; 62.7% of wrong calls name an actual
first-degree relative), confirming the 2026-09-14 claim. An intermediate version of this analysis
hardcoded one pedigree and so mislabelled reps 2 and 3 as unrelated, which briefly inverted the
result; relatedness is now derived from the pedigree file per pool.

(7) Scoring per donor rather than per pool puts all 994 donors on one axis — distance to nearest
pool-mate — over which the error rate moves monotonically from 4.1 to 36.5 per 10,000 with the
related donors at the bottom of the *same* curve. A relatedness-gradient simulation would
interpolate between measured points, so it is not worth running (`01i`).

**Changed on the status board:** one retraction ("margin benefit fades with ambient"), one claim
reworded from "LL-gap" to "doublet margin", seven new findings, one open question closed
(downsampling), one added as the new top priority (the real-data doublet test).

**Next:** `PROGRESS.md` → Next steps, analyses 1 (real-data doublet test) and operational 1
(let 14809580 finish). Writing can start on `RESULTS.md` §1/§3/§5/§7 without further simulation.

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
