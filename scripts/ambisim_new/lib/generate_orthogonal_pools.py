#!/usr/bin/env python3
# "_new" pool sampling scheme: fixes the min_dist/mean_dist entanglement
# documented in NOTES.md (github/demux_pool_design) -- greedy_maxmin and
# greedy_maxmean are the same greedy point-repulsion algorithm applied to
# the same candidate panel, so they produce 50-60% the same donor list and
# their min_dist/mean_dist correlate at r=0.69-0.84 in our pools vs. r=0.14
# (independent) in real 10x data.
#
# Approach: rejection sampling on RANDOM draws, not greedy optimization.
# Random pools already decouple min_dist/mean_dist naturally (r=-0.09 in
# our own multi-ancestry universes, matching real data) -- so draw many
# random candidate pools, then explicitly keep the ones that land in the
# two "disagreement" quadrants a greedy search would never produce:
#   high_min_low_mean: high worst-case separation, but only typical/low
#                       average separation (greedy_maxmin's target with
#                       greedy_maxmean's opposite)
#   low_min_high_mean:  high average separation, but only typical/low
#                       worst-case separation (greedy_maxmean's target
#                       with greedy_maxmin's opposite)
# If demuxlet outcomes track min_dist specifically (not just "distance in
# general"), pools in these two quadrants should behave differently from
# each other even though a naive mean_dist-only or min_dist-only summary
# would call them similar -- exactly the test the greedy pools couldn't run
# because they never produced these combinations.
#
# n_candidates=20000 random draws per rep is a distance-matrix lookup only
# (no genotype simulation), so this is nearly free -- the actual simulation
# cost is the same 4 pools' worth as any other mini-experiment here.
#
# Writes:
#   csv/designs/nominated_pools_new.tsv   (design metrics + which quadrant)
#   txt/donors_new/<pool>.txt             (one donor ID per line)
#   txt/pool_experiments_new.txt          (the pool names, one per line)
#
# EUR_only chosen to match the n=8 baseline universe with the most existing
# data, so these can be compared directly against EUR_only random/
# greedy_maxmin/greedy_maxmean pools already simulated.

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / 'ambisim' / 'lib'))
from pool_nomination import load_design_data, universe_idx_for, pool_dist_metrics, score_pool

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design')
DESIGN_DIR = PROJDIR / 'csv' / 'designs'
DONOR_DIR = PROJDIR / 'txt' / 'donors_new'
POOL_LIST = PROJDIR / 'txt' / 'pool_experiments_new.txt'
NOMINATED_TSV = DESIGN_DIR / 'nominated_pools_new.tsv'

UNIVERSE = 'EUR_only'
N = 8
N_CANDIDATES = 20000
REPS = [1, 2]
QUADRANTS = ['highmin_lowmean', 'lowmin_highmean']
SEED_BASE = 200000  # distinct from both the n=8 grid (42, 542, 642) and the n=16 experiment (100000-100999)


def sample_random_candidates(universe_idx, D, n, n_candidates, seed):
    """Draw n_candidates random n-donor pools and score each. Distance-matrix
    lookups only -- no genotype simulation -- so this is cheap even at 20k."""
    rng = np.random.default_rng(seed)
    universe_arr = np.array(universe_idx)
    rows = []
    for _ in range(n_candidates):
        pool = rng.choice(universe_arr, n, replace=False)
        m = pool_dist_metrics(list(pool), D)
        rows.append((pool, m['min_dist'], m['mean_dist']))
    return rows


def main():
    D, X, sample_ids, meta = load_design_data()
    universe_idx = universe_idx_for(UNIVERSE, meta)
    print(f'{UNIVERSE} universe size: {len(universe_idx)} candidates')

    rows = []
    pool_names = []
    DONOR_DIR.mkdir(parents=True, exist_ok=True)

    for rep in REPS:
        seed = SEED_BASE + rep
        candidates = sample_random_candidates(universe_idx, D, N, N_CANDIDATES, seed)
        mean_dists = np.array([c[2] for c in candidates])
        median_mean = np.median(mean_dists)
        print(f'rep{rep}: {N_CANDIDATES} random candidates, median mean_dist={median_mean:.2f}, '
              f'min_dist range {min(c[1] for c in candidates):.1f}-{max(c[1] for c in candidates):.1f}, '
              f'mean_dist range {mean_dists.min():.1f}-{mean_dists.max():.1f}')

        below = [c for c in candidates if c[2] <= median_mean]
        above = [c for c in candidates if c[2] >= median_mean]

        # high_min_low_mean: among below-median mean_dist, maximize min_dist
        best_hilo = max(below, key=lambda c: c[1])
        # low_min_high_mean: among above-median mean_dist, minimize min_dist
        best_lohi = min(above, key=lambda c: c[1])

        for quadrant, (pool_arr, min_d, mean_d) in [
            ('highmin_lowmean', best_hilo), ('lowmin_highmean', best_lohi)
        ]:
            pool = list(pool_arr)

            # Safeguard added after a diagnostic run (not this one) found that pushing
            # N_CANDIDATES higher can surface a pair as close as min_dist~45 -- well below
            # the real adversarial_family pedigree pairs (~74-77) -- suggesting cryptic
            # relatedness the panel's 'unrelated' flag didn't catch, not genuine coincidence
            # among unrelated people. Verify the closest pair in THIS pool is (a) flagged
            # unrelated in meta and (b) not implausibly close, so a "control" pool doesn't
            # silently become an accidental adversarial one.
            sub = D[np.ix_(pool, pool)]
            np.fill_diagonal(sub, np.inf)
            ci, cj = np.unravel_index(sub.argmin(), sub.shape)
            pair = (sample_ids[pool[ci]], sample_ids[pool[cj]])
            pair_dist = sub[ci, cj]
            for pid in pair:
                is_unrelated = meta.loc[meta['SampleID'] == pid, 'unrelated'].iloc[0]
                assert is_unrelated, f'{pid} in {quadrant}_new rep{rep} closest pair is not flagged unrelated!'
            if pair_dist < 70:
                print(f'  WARNING: {quadrant}_new rep{rep} closest pair {pair} at {pair_dist:.1f} '
                      f'is close to the real adversarial_family pedigree range (~74-77) -- '
                      f'double check this isn\'t cryptic relatedness before using it as a control.')
            pool_name = f'{UNIVERSE}__{quadrant}_new__rep{rep}'
            m = score_pool(pool, D, X)
            rows.append({'universe': UNIVERSE, 'quadrant': quadrant, 'rep': rep,
                         'pool_name': pool_name, 'n_donors': N, 'seed': seed, **m})
            pool_names.append(pool_name)

            donor_ids = [sample_ids[i] for i in pool]
            (DONOR_DIR / f'{pool_name}.txt').write_text('\n'.join(donor_ids) + '\n')
            print(f'{pool_name}: min_dist={m["min_dist"]:.1f} mean_dist={m["mean_dist"]:.1f} '
                  f'(median_mean_dist_of_20k_draws={median_mean:.1f})')

    pd.DataFrame(rows).to_csv(NOMINATED_TSV, sep='\t', index=False)
    print(f'wrote {NOMINATED_TSV}')

    POOL_LIST.write_text('\n'.join(pool_names) + '\n')
    print(f'wrote {POOL_LIST} ({len(pool_names)} pools)')


if __name__ == '__main__':
    main()
