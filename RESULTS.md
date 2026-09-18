# RESULTS — what we know, and what shows it

Current validated findings only. No history: how each result was arrived at (including claims
that were retracted along the way) is in `JOURNAL.md` / `JOURNAL.archive.md`; what has finished
running is in `PROGRESS.md`; pipeline bugs are in `docs/pipeline_issues.md`.

**Basis for every number below:** the chr22-matrix `ambisim/` grid, analyzed at **n=88 pools**
(73 ordinary + 15 adversarial) on 2026-09-14/15 (`notebooks/ambisim/01a-01d`, `results/*.csv`).
Since then the same grid has reached 128/132 demuxed pools but has **not been re-analyzed**. The
genome-wide `ambisim_final/` grid will supersede these numbers when it lands (see `PROGRESS.md`).
Last updated 2026-09-18.

**Standing rule:** every pool-level correlation is reported with the 15 adversarial pools
excluded, and with Spearman alongside Pearson. The adversarial pools sit at the low end of
`min_dist` and are high-leverage; two whole-sample claims have already been retracted because of
them (`JOURNAL.md` status board, `docs/CONVENTIONS.md`).

---

## The questions

1. Does genetic-distance-based donor curation improve demuxlet singlet accuracy?
2. Does it improve classifier margin (LL-gap between best and second-best donor)?
3. Does pooling actual relatives hurt, and how?
4. Does design help most where it's needed (high-ambient droplets)?
5. Does ancestry composition of the universe change the answer?
6. Can min vs. mean pairwise distance be separated as design targets?

---

## 1. Curating unrelated donors by distance does not move raw accuracy

Singlet accuracy is **0.799–0.860 (GEX), 0.903–0.926 (ATAC)** across all 88 pools — narrow. Among
the 73 ordinary pools, `min_dist` explains ~0% of accuracy variance (r≈0, n.s., both
modalities), in every ambient bin and every ancestry subset tested.

## 2. Distance does move margin — modestly

`min_dist` vs. mean LL-gap, ordinary pools: **r=0.31 (GEX), r=0.31 (ATAC), both p<0.01.** This
replicates the n=32 snapshot (r=0.51–0.56) in direction and significance at a smaller point
estimate; the whole-sample r=0.66–0.73 is inflated by adversarial leverage.

At the bottleneck-donor level, `mean_dist` outpredicts `min_dist` (r≈0.66 vs. r≈0.30–0.40) among
ordinary pools. **Do not read this as "mean beats min"** — see §6.

## 3. Real relatives hurt, locally; the singlet rate mostly survives

- **`adversarial_family`/`_mixed`** (6 pools, real 1000G pedigrees) score lower than everything
  else on accuracy (Mann-Whitney p<0.001 both modalities, ~1–2 points) *and* margin (GEX LL-gap
  7.03 vs. 9.92, ATAC 18.96 vs. 27.44, both p<0.0001).
- **The deficit is localized to the related donors.** Within the same pool, donors with a
  relative present: GEX LL-gap 4.79 vs. 7.66 for the unrelated fill donors; ATAC 13.83 vs. 24.71.
- **`adversarial_mindist`** (9 pools, unrelated donors, worst-case-selected) shows a weaker
  version (accuracy p=0.002–0.005 but mean diff as small as +0.004 ATAC; LL-gap GEX 8.93 vs. 9.92,
  p=0.018) with **no localized bad pair** (closest pair 7.24 vs. rest 7.48 GEX) — the greedy
  construction compresses the whole pool uniformly. Thin evidence; not settled.
- **Accuracy and margin move on different scales.** Family pools vs. ordinary: accuracy −2.4%
  (GEX) / −1.2% (ATAC); LL-gap −29.2% / −30.9%. 3/73 ordinary pools (GEX) and 1/73 (ATAC) already
  score at or below the family-pool mean accuracy from ordinary variation alone.

**Read:** multiplexing actual relatives mostly still works for raw singlet rate. A 30% margin drop
moves cells from "confidently correct" to "less-confidently correct"; the cost shows up under
downstream confidence filtering, not in the discrete call.

## 4. Design is a cushion in good droplets, not a rescue for bad ones

`min_dist` vs. LL-gap by ambient-RNA bin (ordinary pools): strongest at 0–10% ambient (r=0.31
both modalities, p<0.01), weakening monotonically to n.s. by ≥40% (r=−0.06 GEX, 0.19 ATAC).
Sequencing depth is not a confound: fixed by simulation design (~20,021 RNA / ~20,240 ATAC reads,
identical median across all 88 pools).

## 5. Mixing ancestries is what gives distance room to matter

- Per single-ancestry universe (n=8–10 pools each): noisy, sign-inconsistent (e.g. `AMR_only`
  GEX r=−0.75 vs. `SAS_only` r=+0.86) — underpowered.
- All single-ancestry universes pooled (n=44): LL-gap r=+0.20, n.s.
- Multi-ancestry universes pooled (n=29): **r=+0.56 GEX (p=0.002), +0.59 ATAC (p=0.001)** — the
  strongest of any subset, because multi-ancestry `min_dist` spans 16.5–24.8 vs. 9.4–16.1.

Restricting to one ancestry (as real-data comparisons do) removes power rather than sharpening
signal. n=29 is still thin.

## 6. Min vs. mean distance cannot be separated with the current strategies

`greedy_maxmin` and `greedy_maxmean` pick a median 5–6/8 identical donors per universe/rep, and
their `min_dist`/`mean_dist` correlate far more than in real pools:

| source | min_dist vs. mean_dist r |
|---|---|
| real 10x pools (n=63) | 0.14 (n.s.) |
| 20,000 random 8-donor draws, `EUR_only` | 0.536 |
| greedy pools, chr22 matrix | 0.69–0.84 |
| `greedy_maxmin`, genome-wide matrix (n=45) | 0.992 |

So every greedy-maxmin-vs-maxmean comparison (including §2's bottleneck result) is not a test of
worst-case vs. average-case design. It does **not** affect §1 (random vs. any optimized strategy)
or §2's continuous regression. Mechanism and the orthogonal-sampling fix:
`docs/minmean_entanglement.md`. The fix's pilot (`ambisim_new/`, 4 `EUR_only` pools) is demuxed but
**not yet analyzed** — the min-vs-mean question is open.

---

## Overall read

Curating unrelated donors by genetic distance — any strategy tested, vs. `random` — has not been
shown to move raw demuxlet accuracy at n=88, in any ambient condition or ancestry subset. It moves
classifier margin reliably but modestly, and not in the droplets that need it most. The one
selection decision shown to move accuracy is **avoiding close relatives** — an argument for a
relatedness check before pooling, not for distance optimization over unrelated donors.

What could still change this is listed in `PROGRESS.md` → Next steps (genome-wide grid, the
orthogonal pilot, n=16, and the downsampling stress test).

## Where the data lives

| what | path |
|---|---|
| per-pool summary (accuracy, LL-gap, distances) | `results/pool_summary.csv` |
| ambient-stratified LL-gap | `results/droplet_ambient_stratified.csv` |
| bottleneck donor | `results/bottleneck_donor_{accuracy,llgap}.csv` |
| pool geometry QC (random-draw baseline, orthogonal pools) | `results/pool_geometry_qc.csv`, `figures/00_*.pdf` |
| nominated donor lists (chr22 matrix) | `results/nominated_pools_n8.tsv`, `results/nominated_pools_new.tsv` |
