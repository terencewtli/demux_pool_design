# Priority pools to simulate next

Line numbers refer to `txt/pool_experiments.txt` (1-indexed; submit with
`qsub -t <N>` against the relevant stage script under `scripts/ambisim/qsub/`
— section 1 pools only need `A01b_run_cr_arc.sh` onward, section 2 pools need
the full pipeline from `A01a_run_simulation.sh` since they haven't been
simulated at all yet; `-t <N1>,<N2>,...` submits several at once).

## 1. Unstick the 3 rep1 core-grid cells that never got past cr_arc

These have a finished ambisim simulation (`drop_data_rand.txt` exists) but
`cr_arc` was never run — not a demuxlet problem, just needs
`A01b_run_cr_arc.sh` submitted:

```
1:  AFR_only__greedy_maxmean__rep1
16: EAS_only__greedy_maxkl__rep1
36: AFR_only__greedy_maxkl__rep1
```

(36 and, if re-nominated with the bug fix, also the already-run AMR_only/
EUR_only/SAS_only `greedy_maxkl` rep1 pools should be re-nominated first — see
README "Bug found and fixed" — before spending compute re-simulating them.)

## 2. Adversarial rep1 pools (highest scientific value)

These are the actual missing piece: every universe below already has full
`random`/`greedy_max*`/`ancestry_balanced` rep1 data, so adding one adversarial
pool per universe gives a direct, controlled worst-case comparison for the
first time, and is the only way to test whether `min_dist` ever wins over
`mean_dist` in this framework (see README point 3) — these pools are the first
ones expected to contain a genuinely close donor pair.

```
109: EUR_only__adversarial_family__rep1
112: EUR_only__adversarial_mindist__rep1
115: AFR_only__adversarial_mindist__rep1
118: EAS_only__adversarial_mindist__rep1
121: SAS_only__adversarial_mindist__rep1
124: AMR_only__adversarial_mindist__rep1
127: EUR_AFR__adversarial_mindist__rep1
130: EUR_AFR__adversarial_family_mixed__rep1
```

8 pools, rep1 only (replicates of the core grid can hold off per current
priority). `adversarial_family`/`adversarial_family_mixed` are especially
high-leverage: `md/analyses.md`'s own "new from claude code" note flags a
single validation figure — a monotone difficulty gradient from cross-ancestry
→ same-ancestry unrelated → first-degree relatives — as a cheap, high-value
addition to the paper pitch, and these 2 pools are exactly what that figure
needs.

## 3. Analysis (no new compute) to run in parallel

- Re-run `notebooks/ambisim/01a-01d` once (1) and (2) land — they auto-discover
  new pools via glob, no hardcoded pool list to edit.
- The Tier-3 downsampling/robustness check from `md/analyses.md` (subsample
  reads/cells, re-run demuxlet, see whether narrow-margin pools degrade first)
  is the single most direct way to test whether the margin advantage found so
  far ever converts into an accuracy advantage — worth prioritizing over
  additional pool replicates.
