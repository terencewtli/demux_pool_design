#!/usr/bin/env python3
# Regenerate pool nominations against the genome-wide, LD-pruned distance matrix
# (1kg_genomewide_dist_matrix.npz, 459,999 SNPs) instead of the old chr22-only,
# ~20k-SNP one -- mirrors ipynb/vcf_metrics/local/02b_pool_nomination_n8.ipynb
# cells 3-11, using pool_nomination.py's shared (already bug-fixed, see its
# nominate_greedy_maxkl EPS-ordering comment) functions rather than
# reimplementing them.
#
# Writes to NEW file names (*_genomewide.tsv/.txt), NOT overwriting
# nominated_pools_n8.tsv / simulation_pools_n8.tsv / pool_experiments.txt --
# those are what the already-submitted/planned ambisim simulation jobs
# reference against the OLD chr22 matrix, and switching the donor selection
# out from under them silently would invalidate anything already run against
# those pool definitions. Promote genomewide -> canonical explicitly, once
# the new pools are reviewed, rather than here.
#
# SCOPE: only the "main" strategies from 02b (random, greedy_maxmin,
# greedy_maxmean, ancestry_balanced, greedy_maxkl) across the 9 standard
# universes. Does NOT (yet) include the adversarial_mindist (02d) or
# adversarial_family_mixed (02e) strategies also present in the existing
# nominated_pools_n8.tsv -- those need their own ped/family-selection inputs
# re-examined against the new matrix before regenerating; out of scope here.
#
# FAST_STRATEGIES (use D only) run in this script directly, in ~seconds.
# greedy_maxkl (uses X, ~23x more SNPs than chr22 -> ~23x slower per step)
# is NOT run here -- see run_greedy_maxkl_genomewide.py + its qsub wrapper,
# submitted separately given the multi-hour cost measured on this data
# (116s for one 8-donor pool on the 353-donor AMR universe alone).

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

FAST_STRATEGIES = ['random', 'greedy_maxmin', 'greedy_maxmean', 'ancestry_balanced']
N_REPS = 5

UNIVERSE_NAMES = ['EUR_only', 'AFR_only', 'EAS_only', 'SAS_only', 'AMR_only',
                  'EUR_AFR', 'EUR_EAS', 'all_3_major', 'all_5']


def main():
    D, X, sample_ids, meta = pn.load_design_data(matrix_file=MATRIX_FILE, meta_file=META_FILE)
    print(f'D {D.shape}  X {X.shape}  n_samples {len(sample_ids)}')

    universes = {name: pn.universe_idx_for(name, meta) for name in UNIVERSE_NAMES}
    for name, idx in universes.items():
        print(f'{name:15s} n={len(idx)}')

    records = []
    for uname, uidx in universes.items():
        n_sp = meta.loc[uidx, 'Superpopulation'].nunique()
        for rep in range(N_REPS):
            rep_seed = pn.SEED + rep * 100
            for strategy in FAST_STRATEGIES:
                if strategy == 'ancestry_balanced' and n_sp == 1:
                    continue  # not meaningful for single-ancestry universes, same as 02b
                pool = pn.nominate(strategy, uidx, D, X, meta, seed=rep_seed)
                m = pn.score_pool(pool, D, X)
                records.append({'universe': uname, 'strategy': strategy, 'rep': rep,
                                 'pool': pool, **m})
        print(f'{uname}: {N_REPS} reps x fast strategies done')

    pool_df = pd.DataFrame(records)
    pool_df.to_csv(NOMINATED_TSV, sep='\t', index=False)
    print(f'\nwrote {len(pool_df)} rows to {NOMINATED_TSV}')

    ### canonical pool per (universe, strategy): highest kl_min, same rule as 02b cell 11
    canonical = (pool_df.sort_values('kl_min', ascending=False)
                 .groupby(['universe', 'strategy']).first().reset_index())

    rows = []
    for _, row in canonical.iterrows():
        pidx = row['pool']
        for rank, i in enumerate(pidx):
            rows.append({'universe': row['universe'], 'strategy': row['strategy'],
                         'donor_rank': rank, 'SampleID': sample_ids[i],
                         'Superpopulation': meta.loc[i, 'Superpopulation'],
                         'min_dist': row['min_dist'], 'kl_min': row['kl_min'],
                         'entropy_mean': row['entropy_mean']})
    sim_input = pd.DataFrame(rows)
    sim_input.to_csv(SIMULATION_TSV, sep='\t', index=False)
    print(f'wrote {len(sim_input)} rows to {SIMULATION_TSV}')

    DONORS_DIR.mkdir(parents=True, exist_ok=True)
    pool_names = []
    for _, row in canonical.iterrows():
        pool_name = f"{row['universe']}__{row['strategy']}"
        pidx = row['pool']
        donor_file = DONORS_DIR / f'{pool_name}_genomewide.txt'
        donor_file.write_text('\n'.join(sample_ids[i] for i in pidx) + '\n')
        pool_names.append(pool_name)

    POOL_LIST_TXT.write_text('\n'.join(pool_names) + '\n')
    print(f'wrote {len(pool_names)} pool names to {POOL_LIST_TXT}')
    print(f'\n{canonical.groupby(["universe", "strategy"]).size().shape[0]} '
          f'(universe, strategy) combinations -- greedy_maxkl not included, see module docstring')


if __name__ == '__main__':
    main()
