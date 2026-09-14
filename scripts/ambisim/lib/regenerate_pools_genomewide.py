#!/usr/bin/env python3
# Regenerate the full 132-pool grid's donor selections and design metrics
# using the genome-wide distance matrix (build_genomewide_dist_matrix.py),
# replacing the chr22-only, 20k-SNP matrix every pool in this project was
# originally selected against (github/demux_pool_design NOTES.md).
#
# NOT regenerated: adversarial_family / adversarial_family_mixed (6 pools).
# Their relatedness structure is real 1000G pedigree data -- independent of
# which chromosome's distance matrix is used -- and their "fill" donors are
# chosen by nominate_adversarial_family() via uniform random choice, not by
# D/X at all. Re-deriving them from scratch would risk picking a different
# family or fill donors for no scientific reason; existing donor lists are
# copied unchanged. Their score_pool() METRICS are still recomputed against
# the new matrix, since min_dist/mean_dist/kl_* are reported for every pool
# regardless of how it was selected.
#
# Everything else (random, greedy_maxmin, greedy_maxmean, greedy_maxkl,
# ancestry_balanced, adversarial_mindist -- 126 pools) is regenerated via a
# fresh nominate() call against the new matrix. Even strategies whose
# nomination logic doesn't directly use D/X (random, ancestry_balanced) are
# still re-run rather than assumed identical: universe_idx is a list of
# integer positions into sample_ids, and the new genome-wide VCF's sample
# order is not guaranteed to match the old chr22 VCF's, so re-deriving from
# scratch is safer than assuming reuse is equivalent.
#
# Writes (does NOT touch the existing ambisim/txt/donors/ or
# nominated_pools_n8.tsv -- this is for review before committing to
# resimulation):
#   csv/designs/nominated_pools_n8_genomewide.tsv
#   txt/donors_genomewide/<pool>.txt
#   csv/designs/genomewide_vs_chr22_donor_overlap.tsv  (diagnostic: how much
#     actually changed per pool)

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pool_nomination import load_design_data, universe_idx_for, nominate, score_pool, SEED

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design')
DESIGN_DIR = PROJDIR / 'csv' / 'designs'
OLD_DONOR_DIR = PROJDIR / 'txt' / 'donors'
NEW_DONOR_DIR = PROJDIR / 'txt' / 'donors_genomewide'
POOL_LIST = PROJDIR / 'ambisim' / 'txt' / 'pool_experiments.txt'
NOMINATED_TSV = DESIGN_DIR / 'nominated_pools_n8_genomewide.tsv'
OVERLAP_TSV = DESIGN_DIR / 'genomewide_vs_chr22_donor_overlap.tsv'

FAMILY_STRATEGIES = {'adversarial_family', 'adversarial_family_mixed'}
N = 8


def main():
    D, X, sample_ids, meta = load_design_data(
        matrix_file='1kg_genomewide_dist_matrix.npz', meta_file='1kg_genomewide_sample_meta.tsv')
    sid_to_idx = {s: i for i, s in enumerate(sample_ids)}
    print(f'loaded genome-wide matrix: {len(sample_ids)} samples, D shape {D.shape}')

    pool_names = [l.strip() for l in open(POOL_LIST)]
    print(f'{len(pool_names)} pools to process')

    NEW_DONOR_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    overlap_rows = []

    for pool_name in pool_names:
        universe, strategy, rep_str = pool_name.split('__')
        rep = int(rep_str.replace('rep', ''))
        old_donors = open(OLD_DONOR_DIR / f'{pool_name}.txt').read().split()

        if strategy in FAMILY_STRATEGIES:
            # Copied unchanged -- see module docstring for why.
            new_donors = old_donors
        else:
            universe_idx = universe_idx_for(universe, meta)
            seed = SEED + rep * 100
            pool_idx = nominate(strategy, universe_idx, D, X, meta, n=N, seed=seed)
            new_donors = [sample_ids[i] for i in pool_idx]

        # Recompute metrics against the new matrix for EVERY pool, including
        # the family ones (their donor identities are unchanged, but
        # min_dist/mean_dist/kl_* must be measured against the new matrix).
        new_idx = [sid_to_idx[d] for d in new_donors]
        m = score_pool(new_idx, D, X)

        rows.append({'universe': universe, 'strategy': strategy, 'rep': rep,
                     'pool': pool_name, **m})
        (NEW_DONOR_DIR / f'{pool_name}.txt').write_text('\n'.join(new_donors) + '\n')

        n_shared = len(set(old_donors) & set(new_donors))
        overlap_rows.append({'pool': pool_name, 'universe': universe, 'strategy': strategy,
                              'rep': rep, 'n_shared_of_8': n_shared,
                              'unchanged': n_shared == N})

    design_df = pd.DataFrame(rows)
    design_df.to_csv(NOMINATED_TSV, sep='\t', index=False)
    print(f'wrote {NOMINATED_TSV}')

    overlap_df = pd.DataFrame(overlap_rows)
    overlap_df.to_csv(OVERLAP_TSV, sep='\t', index=False)
    print(f'wrote {OVERLAP_TSV}')

    print()
    print('=== donor-list overlap summary (old chr22-based vs. new genome-wide-based) ===')
    print(overlap_df.groupby('strategy')['n_shared_of_8'].agg(['mean', 'min', 'max', 'count']))
    print()
    n_fully_unchanged = overlap_df['unchanged'].sum()
    print(f'{n_fully_unchanged}/{len(overlap_df)} pools have the IDENTICAL 8 donors under the new matrix')


if __name__ == '__main__':
    main()
