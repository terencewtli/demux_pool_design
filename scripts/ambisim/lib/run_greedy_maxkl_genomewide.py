#!/usr/bin/env python3
# The expensive half of create_pools_genomewide.py: greedy_maxkl over the genome-wide
# 459,999-SNP dosage matrix, for the 5 single-ancestry universes (same scope as the
# original 02b notebook restricted greedy_maxkl to -- KL_MAX_UNIVERSE_SIZE there kept it
# off the 4 multi-ancestry universes; pool_experiments.txt confirms greedy_maxkl was never
# run for EUR_AFR/EUR_EAS/all_3_major/all_5 either).
#
# Timed on this data before submitting (AMR_only, n=353): 116s for one 8-donor pool.
# Estimated total across 5 universes x 5 reps: ~70-90 min -- run via qsub
# (scripts/ambisim/qsub/B03_greedy_maxkl_genomewide.sh), not interactively.
#
# Appends to the SAME *_genomewide files create_pools_genomewide.py wrote (that script
# must be run first) -- never touches the chr22-based canonical files.

from pathlib import Path

import numpy as np
import pandas as pd

import pool_nomination as pn

DESIGN_DIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design/csv/designs')
AMBISIM_DIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim')

MATRIX_FILE = '1kg_genomewide_dist_matrix.npz'
META_FILE = '1kg_genomewide_sample_meta.tsv'

NOMINATED_TSV = DESIGN_DIR / 'nominated_pools_n8_genomewide.tsv'
SIMULATION_TSV = DESIGN_DIR / 'simulation_pools_n8_genomewide.tsv'
POOL_LIST_TXT = AMBISIM_DIR / 'txt' / 'pool_experiments_genomewide.txt'
DONORS_DIR = AMBISIM_DIR / 'txt' / 'donors'

UNIVERSE_NAMES = ['EUR_only', 'AFR_only', 'EAS_only', 'SAS_only', 'AMR_only']
STRATEGY = 'greedy_maxkl'
N_REPS = 5


def main():
    D, X, sample_ids, meta = pn.load_design_data(matrix_file=MATRIX_FILE, meta_file=META_FILE)
    print(f'D {D.shape}  X {X.shape}', flush=True)

    nominated = pd.read_csv(NOMINATED_TSV, sep='\t')
    already = nominated[nominated['strategy'] == STRATEGY]['universe'].unique()

    records = []
    for uname in UNIVERSE_NAMES:
        if uname in already:
            print(f'{uname}__{STRATEGY} already present, skipping', flush=True)
            continue
        uidx = pn.universe_idx_for(uname, meta)
        print(f'{uname} (n={len(uidx)}) starting {N_REPS} reps...', flush=True)
        for rep in range(N_REPS):
            rep_seed = pn.SEED + rep * 100
            pool = pn.nominate(STRATEGY, uidx, D, X, meta, seed=rep_seed)
            m = pn.score_pool(pool, D, X)
            records.append({'universe': uname, 'strategy': STRATEGY, 'rep': rep,
                             'pool': pool, **m})
            print(f'  rep {rep}: kl_min={m["kl_min"]:.4f} min_dist={m["min_dist"]:.4f}',
                  flush=True)

    if not records:
        print('nothing new to add', flush=True)
        return

    new_df = pd.DataFrame(records)
    updated = pd.concat([nominated, new_df], ignore_index=True)
    updated.to_csv(NOMINATED_TSV, sep='\t', index=False)
    print(f'appended {len(new_df)} rows to {NOMINATED_TSV}', flush=True)

    canonical = (new_df.sort_values('kl_min', ascending=False)
                 .groupby('universe').first().reset_index())

    sim = pd.read_csv(SIMULATION_TSV, sep='\t')
    new_sim_rows = []
    for _, row in canonical.iterrows():
        pidx = row['pool']
        for rank, i in enumerate(pidx):
            new_sim_rows.append({'universe': row['universe'], 'strategy': STRATEGY,
                                  'donor_rank': rank, 'SampleID': sample_ids[i],
                                  'Superpopulation': meta.loc[i, 'Superpopulation'],
                                  'min_dist': row['min_dist'], 'kl_min': row['kl_min'],
                                  'entropy_mean': row['entropy_mean']})
    updated_sim = pd.concat([sim, pd.DataFrame(new_sim_rows)], ignore_index=True)
    updated_sim.to_csv(SIMULATION_TSV, sep='\t', index=False)
    print(f'appended {len(new_sim_rows)} rows to {SIMULATION_TSV}', flush=True)

    pool_names = POOL_LIST_TXT.read_text().splitlines() if POOL_LIST_TXT.exists() else []
    for _, row in canonical.iterrows():
        pool_name = f"{row['universe']}__{STRATEGY}"
        donor_file = DONORS_DIR / f'{pool_name}_genomewide.txt'
        donor_file.write_text('\n'.join(sample_ids[i] for i in row['pool']) + '\n')
        pool_names.append(pool_name)
    POOL_LIST_TXT.write_text('\n'.join(pool_names) + '\n')
    print(f'{POOL_LIST_TXT} now has {len(pool_names)} pools', flush=True)


if __name__ == '__main__':
    main()
