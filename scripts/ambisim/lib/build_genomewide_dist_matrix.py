#!/usr/bin/env python3
# Genome-wide replacement for csv/designs/1kg_chr22_dist_matrix.npz, which
# every min_dist/mean_dist/kl_min/kl_mean value in this whole project (the
# n=8 grid, the n=16 experiment, the orthogonal-sampling pools) was built
# from -- and which turned out to be chr22 only, 20,000 subsampled SNPs
# (see ipynb/vcf_metrics/local/02a_distance_matrix.ipynb). Population
# ancestry structure is consistent across chromosomes, so universe-level
# findings likely replicate, but donor-level selection (exactly who gets
# picked as "most extreme") and cryptic-relatedness safety (see
# github/demux_pool_design NOTES.md) both benefit from genome-wide data.
#
# Source: the SAME already-filtered, already-genome-wide VCF used to build
# every pool's actual simulated genotypes (1000G.merged.common_biallelic.vcf.gz,
# 24 contigs, 6,869,484 common (MAF>=0.05) biallelic SNPs, from
# B01_filter_common_biallelic.sh) -- so the distance matrix and the actual
# simulated genotypes are now built from the same underlying variant set,
# which chr22-only was not.
#
# Two-pass extraction (not a single reservoir-sampling pass over all 6.87M
# sites' full genotypes, which would be extremely slow in Python): first
# pull just CHROM/POS for all sites (cheap, 2 columns), randomly sample
# N_SITES of them in Python, then use bcftools -R to pull only those sites'
# full genotypes (the expensive part, but now ~34x smaller).
#
# Output: csv/designs/1kg_genomewide_dist_matrix.npz (D, X, sample_ids, pcs)
#         csv/designs/1kg_genomewide_sample_meta.tsv

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design')
DESIGN_DIR = PROJDIR / 'csv' / 'designs'
VCF = Path('/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/1000G.merged.common_biallelic.vcf.gz')
PED = Path('/u/project/cluo/terencew/demux_benchmark/pool_design/csv/1000G/meta/20130606_g1k_3202_samples_ped_population.txt')
BCFTOOLS = '/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools'

OUTMAT = DESIGN_DIR / '1kg_genomewide_dist_matrix.npz'
OUTMETA = DESIGN_DIR / '1kg_genomewide_sample_meta.tsv'
SAMPLED_SITES = DESIGN_DIR / '_genomewide_sampled_sites.tsv'  # scratch, kept for reproducibility

N_SITES = 200_000  # ~10x the original chr22-only N_SITES=20,000, ~2.9% of the 6,869,484 available
SEED = 42


def main():
    DESIGN_DIR.mkdir(parents=True, exist_ok=True)

    print(f'{sys.argv[0]}: Step 1/5 -- listing sample IDs', flush=True)
    sample_ids = subprocess.run(
        [BCFTOOLS, 'query', '-l', str(VCF)], capture_output=True, text=True, check=True
    ).stdout.strip().split('\n')
    n_samples = len(sample_ids)
    print(f'{n_samples} samples in VCF', flush=True)

    print(f'{sys.argv[0]}: Step 2/5 -- listing all site positions (CHROM/POS only, cheap)', flush=True)
    proc = subprocess.run(
        [BCFTOOLS, 'query', '-f', '%CHROM\t%POS\n', str(VCF)],
        capture_output=True, text=True, check=True,
    )
    all_sites = proc.stdout.strip().split('\n')
    print(f'{len(all_sites)} total sites available', flush=True)

    print(f'{sys.argv[0]}: Step 3/5 -- randomly sampling {N_SITES} sites (seed={SEED})', flush=True)
    rng = np.random.default_rng(SEED)
    idx = rng.choice(len(all_sites), size=N_SITES, replace=False)
    idx.sort()  # bcftools -R wants sorted regions
    sampled = [all_sites[i] for i in idx]
    with open(SAMPLED_SITES, 'w') as f:
        for line in sampled:
            f.write(line + '\n')
    print(f'wrote {SAMPLED_SITES}', flush=True)

    print(f'{sys.argv[0]}: Step 4/5 -- extracting genotypes for the {N_SITES} sampled sites', flush=True)
    proc = subprocess.run(
        [BCFTOOLS, 'view', '-R', str(SAMPLED_SITES), str(VCF)],
        stdout=subprocess.PIPE, check=True,
    )
    proc2 = subprocess.run(
        [BCFTOOLS, 'query', '-f', '[%GT\t]\n'],
        input=proc.stdout, stdout=subprocess.PIPE, check=True,
    )
    lines = proc2.stdout.decode().strip().split('\n')
    print(f'got genotypes for {len(lines)} sites (expected ~{N_SITES})', flush=True)

    print(f'{sys.argv[0]}: Step 5/5 -- building dosage matrix, distance matrix, PCA', flush=True)
    n_sites_actual = len(lines)
    dosage = np.zeros((n_sites_actual, n_samples), dtype=np.int8)
    for i, line in enumerate(lines):
        gts = line.strip().split('\t')
        for j, gt in enumerate(gts[:n_samples]):
            gt = gt.replace('|', '/')
            alleles = gt.split('/')
            try:
                dosage[i, j] = sum(int(a) for a in alleles if a != '.')
            except ValueError:
                dosage[i, j] = 0  # missing -> 0, matches the original chr22 notebook's convention

    # Samples x sites (transpose of how it was read), matching pool_nomination.py's expectation
    X = dosage.T.astype(np.float32)
    print(f'X shape: {X.shape}', flush=True)

    print('computing pairwise Euclidean distance matrix (XXT trick)...', flush=True)
    XXT = X @ X.T
    sq = (X ** 2).sum(axis=1)
    D_sq = sq[:, None] + sq[None, :] - 2 * XXT
    D_sq = np.maximum(D_sq, 0)
    D = np.sqrt(D_sq).astype(np.float32)
    print(f'D shape: {D.shape}  min={D[D>0].min():.1f}  max={D.max():.1f}  mean={D[D>0].mean():.1f}', flush=True)

    print('PCA...', flush=True)
    pca = PCA(n_components=10, random_state=SEED)
    PCS = pca.fit_transform(X)
    print(f'PC1-2 explained variance: {pca.explained_variance_ratio_[:2]}', flush=True)

    np.savez_compressed(OUTMAT, D=D, X=X, sample_ids=np.array(sample_ids), pcs=PCS)
    print(f'saved {OUTMAT}', flush=True)

    meta = pd.read_csv(PED, sep=' ')
    meta['unrelated'] = (meta['FatherID'] == '0') & (meta['MotherID'] == '0')
    meta_vcf = pd.DataFrame({'SampleID': sample_ids, 'vcf_idx': np.arange(n_samples)})
    meta_vcf = meta_vcf.merge(meta[['SampleID', 'Population', 'Superpopulation', 'unrelated']],
                               on='SampleID', how='left')
    meta_vcf.to_csv(OUTMETA, sep='\t', index=False)
    print(f'saved {OUTMETA}', flush=True)
    print(f'unrelated samples: {meta_vcf["unrelated"].sum()} / {len(meta_vcf)}', flush=True)


if __name__ == '__main__':
    main()
