#!/usr/bin/env python3
# Shared nomination + scoring functions, factored out of
# ipynb/vcf_metrics/local/02b_pool_nomination_n8.ipynb (cells 3-5) so that
# scripts (not just the interactive notebook) can generate additional pool
# replicates. Kept verbatim/behavior-identical to the notebook -- this is
# not a rewrite, just a shared import point.

from pathlib import Path

import numpy as np
import pandas as pd

DESIGN_DIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design/csv/designs')

N_POOL = 8
SEED = 42
EPS = 1e-6


def load_design_data():
    npz = np.load(DESIGN_DIR / '1kg_chr22_dist_matrix.npz', allow_pickle=True)
    D = npz['D']
    X = npz['X']
    sample_ids = list(npz['sample_ids'])
    meta = pd.read_csv(DESIGN_DIR / '1kg_sample_meta.tsv', sep='\t')
    assert list(meta['SampleID']) == sample_ids, 'sample order mismatch'
    return D, X, sample_ids, meta


def get_universe_idx(superpops, meta, unrelated_only=True):
    """Return VCF indices for samples from the given superpopulation(s)."""
    mask = meta['Superpopulation'].isin(superpops)
    if unrelated_only:
        mask &= meta['unrelated'] == True
    return meta[mask].index.tolist()


def universe_idx_for(universe_name, meta):
    """Reconstruct the eligible-donor index list for one of the 9 named universes."""
    parts = {
        'EUR_only': ['EUR'], 'AFR_only': ['AFR'], 'EAS_only': ['EAS'],
        'SAS_only': ['SAS'], 'AMR_only': ['AMR'],
        'EUR_AFR': ['EUR', 'AFR'], 'EUR_EAS': ['EUR', 'EAS'],
        'all_3_major': ['EUR', 'AFR', 'EAS'],
        'all_5': ['EUR', 'AFR', 'EAS', 'SAS', 'AMR'],
    }
    return get_universe_idx(parts[universe_name], meta)


### Pool scoring (verbatim from notebook cell 4)

def pool_dist_metrics(pool_idx, D):
    sub = D[np.ix_(pool_idx, pool_idx)]
    ij = np.triu_indices(len(pool_idx), k=1)
    dists = sub[ij]
    return {'min_dist': dists.min(), 'mean_dist': dists.mean(), 'max_dist': dists.max()}


def pool_entropy_score(pool_idx, X):
    g = X[pool_idx]
    n = g.shape[0]
    counts = np.stack([(g == k).sum(axis=0) for k in range(3)]).astype(np.float32)
    probs = counts / n
    H = -np.sum(np.where(probs > 0, probs * np.log(probs), 0), axis=0)
    frac_discrim = (H > 0).mean()
    return {'entropy_mean': H.mean(), 'entropy_p25': np.percentile(H, 25),
            'frac_discrim': frac_discrim}


def pool_kl_divergence(pool_idx, X):
    g = X[pool_idx].astype(np.float32)
    p_alt = g / 2.0
    p_ref = 1.0 - p_alt

    n_pool, n_sites = g.shape
    kl_vals = []
    for i in range(n_pool):
        for j in range(i + 1, n_pool):
            pi_ref, pj_ref = p_ref[i] + EPS, p_ref[j] + EPS
            pi_alt, pj_alt = p_alt[i] + EPS, p_alt[j] + EPS
            kl_ij = (pi_ref * np.log(pi_ref / pj_ref) +
                     pi_alt * np.log(pi_alt / pj_alt))
            pi2, pj2 = p_ref[j] + EPS, p_ref[i] + EPS
            kl_ji = (pj2 * np.log(pj2 / (pi_ref)) +
                     (p_alt[j] + EPS) * np.log((p_alt[j] + EPS) / pi_alt))
            kl_sym = 0.5 * (kl_ij + kl_ji)
            kl_vals.append(kl_sym.mean())

    kl_vals = np.array(kl_vals)
    return {'kl_min': kl_vals.min(), 'kl_mean': kl_vals.mean()}


def pool_det_score(pool_idx, X):
    g = X[pool_idx].astype(np.float64)
    g = g - g.mean(axis=1, keepdims=True)
    cov = g @ g.T / g.shape[1]
    sign, logdet = np.linalg.slogdet(cov)
    return {'log_det': logdet if sign > 0 else -np.inf}


def score_pool(pool_idx, D, X):
    m = {}
    m.update(pool_dist_metrics(pool_idx, D))
    m.update(pool_entropy_score(pool_idx, X))
    m.update(pool_kl_divergence(pool_idx, X))
    m.update(pool_det_score(pool_idx, X))
    return m


### Nomination strategies (verbatim from notebook cell 5)

def nominate_random(universe_idx, n=N_POOL, seed=SEED):
    rng = np.random.default_rng(seed)
    return list(rng.choice(universe_idx, n, replace=False))


def nominate_greedy_maxmin(universe_idx, D, n=N_POOL, seed=SEED):
    rng = np.random.default_rng(seed)
    pool = [int(rng.choice(universe_idx))]
    remaining = set(universe_idx) - {pool[0]}

    while len(pool) < n:
        pool_arr = np.array(pool)
        cands = np.array(list(remaining))
        min_dists_to_pool = D[np.ix_(cands, pool_arr)].min(axis=1)
        best = cands[min_dists_to_pool.argmax()]
        pool.append(int(best))
        remaining.remove(int(best))

    return pool


def nominate_greedy_maxmean(universe_idx, D, n=N_POOL, seed=SEED):
    rng = np.random.default_rng(seed)
    pool = [int(rng.choice(universe_idx))]
    remaining = set(universe_idx) - {pool[0]}

    while len(pool) < n:
        pool_arr = np.array(pool)
        cands = np.array(list(remaining))
        mean_dists = D[np.ix_(cands, pool_arr)].mean(axis=1)
        best = cands[mean_dists.argmax()]
        pool.append(int(best))
        remaining.remove(int(best))

    return pool


def nominate_greedy_maxkl(universe_idx, X, n=N_POOL, seed=SEED):
    rng = np.random.default_rng(seed)
    pool = [int(rng.choice(universe_idx))]
    remaining = set(universe_idx) - {pool[0]}

    def kl_to_pool(cand_idx, pool):
        # NOTE: EPS must be added *after* taking the complement (1 - p), not
        # added to p itself before complementing -- otherwise a homozygous-alt
        # site (dosage=2 -> p=1.0) gives (1 - (1+EPS)) < 0 and log() of a
        # negative number (NaN), which silently corrupts the argmax selection
        # below for every site where either donor is homozygous-alt. This bit
        # pool_kl_divergence() above (which computes p_ref = 1 - p_alt before
        # adding EPS) and is why greedy_maxkl pools scored at/below random.
        p_cand_alt = X[cand_idx].astype(np.float32) / 2.0
        p_cand_ref = 1.0 - p_cand_alt
        p_cand_alt, p_cand_ref = p_cand_alt + EPS, p_cand_ref + EPS
        kls = []
        for pi in pool:
            p_pi_alt = X[pi].astype(np.float32) / 2.0
            p_pi_ref = 1.0 - p_pi_alt
            p_pi_alt, p_pi_ref = p_pi_alt + EPS, p_pi_ref + EPS
            kl_fwd = (p_cand_alt * np.log(p_cand_alt / p_pi_alt) +
                      p_cand_ref * np.log(p_cand_ref / p_pi_ref))
            kl_rev = (p_pi_alt * np.log(p_pi_alt / p_cand_alt) +
                      p_pi_ref * np.log(p_pi_ref / p_cand_ref))
            kls.append(0.5 * (kl_fwd + kl_rev).mean())
        return min(kls)

    while len(pool) < n:
        cands = list(remaining)
        min_kls = [kl_to_pool(c, pool) for c in cands]
        best = cands[np.argmax(min_kls)]
        pool.append(best)
        remaining.remove(best)

    return pool


def nominate_ancestry_balanced(universe_idx, meta, n=N_POOL, seed=SEED):
    rng = np.random.default_rng(seed)
    sp_in_universe = meta.loc[universe_idx, 'Superpopulation'].unique()
    per_sp = int(np.ceil(n / len(sp_in_universe)))
    pool = []
    for sp in sp_in_universe:
        sp_idx = meta.loc[
            (meta.index.isin(universe_idx)) & (meta['Superpopulation'] == sp)
        ].index.tolist()
        k = min(per_sp, len(sp_idx))
        pool += list(rng.choice(sp_idx, k, replace=False))
    rng.shuffle(pool)
    if len(pool) > n:
        pool = pool[:n]
    elif len(pool) < n:
        leftover = [i for i in universe_idx if i not in pool]
        pool += list(rng.choice(leftover, n - len(pool), replace=False))
    return pool


### Adversarial nomination (factored out of
### ipynb/vcf_metrics/local/02d_adversarial_pools.ipynb cells 11/13, so
### 02e_adversarial_eur_afr.ipynb can reuse them for a multi-ancestry universe)

def nominate_greedy_mindist(universe_idx, D, n=N_POOL, seed=SEED):
    """Mirror image of nominate_greedy_maxmin: at each step add the candidate
    that minimizes the max distance to the current pool (complete-linkage-style
    tightening). This is the adversarial floor, not an average-case control
    like 'random'."""
    rng = np.random.default_rng(seed)
    pool = [int(rng.choice(universe_idx))]
    remaining = set(universe_idx) - {pool[0]}

    while len(pool) < n:
        pool_arr = np.array(pool)
        cands = np.array(list(remaining))
        max_dists_to_pool = D[np.ix_(cands, pool_arr)].max(axis=1)
        best = cands[max_dists_to_pool.argmin()]
        pool.append(int(best))
        remaining.remove(int(best))

    return pool


def nominate_adversarial_family(fam_id, ped, sample_ids, fill_universe_idx, n=N_POOL, seed=SEED):
    """A real 1000G family/pedigree (~50% IBD), padded to n with unrelated
    fill drawn from fill_universe_idx. Passing the family's own ancestry
    universe reproduces a same-ancestry adversarial pool; passing a
    *different* ancestry's universe produces a mixed-ancestry adversarial
    pool (related cluster embedded in an otherwise cross-ancestry batch)."""
    fam_sample_ids = ped.loc[ped['FamilyID'] == fam_id, 'SampleID'].tolist()
    fam_idx = [sample_ids.index(s) for s in fam_sample_ids][:n]
    rng = np.random.default_rng(seed)
    leftover = [i for i in fill_universe_idx if i not in fam_idx]
    n_fill = n - len(fam_idx)
    fill = list(rng.choice(leftover, n_fill, replace=False)) if n_fill > 0 else []
    return fam_idx + fill


def nominate(strategy, universe_idx, D, X, meta, n=N_POOL, seed=SEED):
    """Dispatch to the right nominate_* function by strategy name."""
    if strategy == 'random':
        return nominate_random(universe_idx, n=n, seed=seed)
    if strategy == 'greedy_maxmin':
        return nominate_greedy_maxmin(universe_idx, D, n=n, seed=seed)
    if strategy == 'greedy_maxmean':
        return nominate_greedy_maxmean(universe_idx, D, n=n, seed=seed)
    if strategy == 'greedy_maxkl':
        return nominate_greedy_maxkl(universe_idx, X, n=n, seed=seed)
    if strategy == 'ancestry_balanced':
        return nominate_ancestry_balanced(universe_idx, meta, n=n, seed=seed)
    raise ValueError(f'unknown strategy: {strategy}')
