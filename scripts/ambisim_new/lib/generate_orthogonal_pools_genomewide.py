#!/usr/bin/env python3
# Genome-wide re-run of generate_orthogonal_pools.py's rejection-sampling fix for the
# min_dist/mean_dist entanglement (see github/demux_pool_design/NOTES.md for the full writeup
# of why greedy_maxmin/greedy_maxmean don't isolate their target metric, and why plain random
# draws only partially decouple them).
#
# Re-confirmed on the genome-wide matrix before writing this (2026-09-15): greedy_maxmin
# min_dist-vs-mean_dist r=0.992 (n=45, even higher than chr22's r=0.69-0.72), greedy_maxmean
# r=0.520, random r=0.210 -- donor overlap between greedy_maxmin/greedy_maxmean still a median
# 5/8 identical donors. Same entanglement, same fix warranted.
#
# Same method as generate_orthogonal_pools.py (unchanged deliberately -- the approach already
# works, only the data source needed updating): draw N_CANDIDATES random pools per universe,
# median-split on mean_dist, take the min_dist-maximizing pool below the median
# ("highmin_lowmean") and the min_dist-minimizing pool above it ("lowmin_highmean"). Same
# cryptic-relatedness safeguard as the chr22 script (assert closest pair is flagged unrelated,
# warn if suspiciously close) -- BUT the "suspiciously close" threshold needed recalibrating,
# not copying verbatim: the chr22 script's D used raw Euclidean-ish genotype distance
# (threshold 70, real range ~0-150); this matrix uses PLINK 1-IBS distance (range ~0.10-0.29).
# A first run with the old threshold (70) fired on literally every pool -- caught immediately
# since a 100% warning rate is itself a sign the check is miscalibrated, not that every pool
# is secretly related. Recalibrated empirically against this matrix's own known-truth pairs:
# 1205 real 1000G parent-child pairs sit at D=0.1056-0.1543 (mean 0.1359, std 0.009), while the
# whole matrix's 1st percentile is 0.1984 -- clean separation, no overlap. Threshold set to
# 0.16 (just above the observed parent-child max), and only 1201/5,124,801 total pairs in the
# whole matrix fall below it -- almost exactly the count of known parent-child pairs, so this
# threshold has good specificity, not just a guess.
CRYPTIC_RELATEDNESS_THRESHOLD = 0.16
#
# DIFFERENCE from the original: generalized to all 9 universes, not just EUR_only. This is
# cheap to do for the NOMINATION step (distance-matrix lookups only, ~seconds/universe) --
# but per NOTES.md's own recommended path, actually SIMULATING (ambisim/cellranger-arc/
# demuxlet) more than the original EUR_only pilot should wait until that pilot shows a real
# LL-gap/accuracy difference between the two quadrants. This script only nominates + writes
# donor lists; it does not submit any simulation jobs.
#
# Writes:
#   csv/designs/nominated_pools_new_genomewide.tsv
#   txt/donors_new_genomewide/<pool>.txt
#   txt/pool_experiments_new_genomewide.txt

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / 'ambisim' / 'lib'))
from pool_nomination import load_design_data, universe_idx_for, pool_dist_metrics, score_pool

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design')
DESIGN_DIR = PROJDIR / 'csv' / 'designs'
DONOR_DIR = PROJDIR / 'txt' / 'donors_new_genomewide'
POOL_LIST = PROJDIR / 'txt' / 'pool_experiments_new_genomewide.txt'
NOMINATED_TSV = DESIGN_DIR / 'nominated_pools_new_genomewide.tsv'

UNIVERSE_NAMES = ['EUR_only', 'AFR_only', 'EAS_only', 'SAS_only', 'AMR_only',
                  'EUR_AFR', 'EUR_EAS', 'all_3_major', 'all_5']
N = 8
N_CANDIDATES = 20000
REPS = [1, 2]
SEED_BASE = 300000  # distinct from chr22's 200000 base and the n=8/n=16 grids


def sample_random_candidates(universe_idx, D, n, n_candidates, seed):
    rng = np.random.default_rng(seed)
    universe_arr = np.array(universe_idx)
    rows = []
    for _ in range(n_candidates):
        pool = rng.choice(universe_arr, n, replace=False)
        m = pool_dist_metrics(list(pool), D)
        rows.append((pool, m['min_dist'], m['mean_dist']))
    return rows


def main():
    D, X, sample_ids, meta = load_design_data(
        matrix_file='1kg_genomewide_dist_matrix.npz', meta_file='1kg_genomewide_sample_meta.tsv')

    rows = []
    pool_names = []
    DONOR_DIR.mkdir(parents=True, exist_ok=True)

    for universe in UNIVERSE_NAMES:
        universe_idx = universe_idx_for(universe, meta)
        print(f'\n--- {universe} (n={len(universe_idx)}) ---')

        for rep in REPS:
            seed = SEED_BASE + hash(universe) % 10000 + rep
            candidates = sample_random_candidates(universe_idx, D, N, N_CANDIDATES, seed)
            mean_dists = np.array([c[2] for c in candidates])
            median_mean = np.median(mean_dists)
            print(f'rep{rep}: median mean_dist={median_mean:.4f}, '
                  f'min_dist range {min(c[1] for c in candidates):.4f}-{max(c[1] for c in candidates):.4f}, '
                  f'mean_dist range {mean_dists.min():.4f}-{mean_dists.max():.4f}')

            below = [c for c in candidates if c[2] <= median_mean]
            above = [c for c in candidates if c[2] >= median_mean]
            best_hilo = max(below, key=lambda c: c[1])
            best_lohi = min(above, key=lambda c: c[1])

            for quadrant, (pool_arr, min_d, mean_d) in [
                ('highmin_lowmean', best_hilo), ('lowmin_highmean', best_lohi)
            ]:
                pool = list(pool_arr)

                # same cryptic-relatedness safeguard as the chr22-based script
                sub = D[np.ix_(pool, pool)]
                np.fill_diagonal(sub, np.inf)
                ci, cj = np.unravel_index(sub.argmin(), sub.shape)
                pair = (sample_ids[pool[ci]], sample_ids[pool[cj]])
                pair_dist = sub[ci, cj]
                for pid in pair:
                    is_unrelated = meta.loc[meta['SampleID'] == pid, 'unrelated'].iloc[0]
                    assert is_unrelated, \
                        f'{pid} in {universe}__{quadrant}_new_genomewide rep{rep} closest pair not flagged unrelated!'
                if pair_dist < CRYPTIC_RELATEDNESS_THRESHOLD:
                    print(f'  WARNING: {universe}__{quadrant} rep{rep} closest pair {pair} '
                          f'at {pair_dist:.4f} is below the calibrated parent-child range '
                          f'(0.1056-0.1543) -- likely cryptic relatedness, double check before '
                          f'using as a control.')

                pool_name = f'{universe}__{quadrant}_new_genomewide__rep{rep}'
                m = score_pool(pool, D, X)
                rows.append({'universe': universe, 'quadrant': quadrant, 'rep': rep,
                             'pool_name': pool_name, 'n_donors': N, 'seed': seed, **m})
                pool_names.append(pool_name)

                donor_ids = [sample_ids[i] for i in pool]
                (DONOR_DIR / f'{pool_name}.txt').write_text('\n'.join(donor_ids) + '\n')
                print(f'{pool_name}: min_dist={m["min_dist"]:.4f} mean_dist={m["mean_dist"]:.4f} '
                      f'(median_mean_dist_of_20k_draws={median_mean:.4f})')

    pd.DataFrame(rows).to_csv(NOMINATED_TSV, sep='\t', index=False)
    print(f'\nwrote {NOMINATED_TSV}')

    POOL_LIST.write_text('\n'.join(pool_names) + '\n')
    print(f'wrote {POOL_LIST} ({len(pool_names)} pools)')


if __name__ == '__main__':
    main()
