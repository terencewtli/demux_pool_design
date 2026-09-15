#!/usr/bin/env python3
# Genome-wide, LD-pruned replacement for csv/designs/1kg_chr22_dist_matrix.npz
# using PLINK (compiled, purpose-built for exactly this) instead of a custom
# numpy/bcftools approach. Supersedes build_genomewide_dist_matrix.py, which
# used a raw 200,000-SNP random subsample -- fine for avoiding chr22-only
# bias, but a naive full-variant or arbitrary-subsample approach both miss
# the standard fix for pairwise-distance/PCA/relatedness work: SNPs in LD
# (inherited in blocks) aren't independent evidence, so an un-pruned SNP set
# over-weights whatever genomic regions happen to have long-range LD. The
# standard practice (and what this script does) is to LD-prune first, then
# compute distance on the pruned, approximately-independent set.
#
# Source VCF: 1000G.merged.common_biallelic.vcf.gz (24 contigs, 6,869,484
# common (MAF>=0.05) biallelic SNPs) -- the SAME file A00_prep_pool_vcfs.sh
# subsets for actual pool simulation, so the distance matrix and the real
# simulated genotypes share one variant set.
#
# Pipeline (each step is a fast, compiled PLINK operation, not a Python loop
# over millions of genotype lines):
#   1. VCF -> PLINK binary (.bed/.bim/.fam)
#   2. --indep-pairwise 50 5 0.2 (standard window/step/r^2 LD-pruning params)
#   3. --distance square 1-ibs flat-missing on the pruned set (1-IBS is the
#      standard population-genetics pairwise distance; the original chr22
#      notebook's own comment already noted Euclidean/Manhattan/IBS-based
#      distances are rank-equivalent for genotype data, r>0.99 correlated)
#   4. --recode A on the same pruned set, for the dosage matrix X (still
#      needed by nominate_greedy_maxkl's KL divergence)
#   5. flashpca on a pruned .bed/.bim/.fam (--bfile, not --extract, so a
#      separate pruned bfile is made first) -- replaces an in-Python sklearn
#      PCA on X, which OOM-killed this job (see step 5 comment below)
#   6. This script: parse PLINK + flashpca's outputs into the .npz format
#      pool_nomination.py expects.
#
# NOTE: 1-IBS is on a different numeric scale than the old chr22 matrix's
# raw Euclidean distances (roughly 0-1 here vs. ~0-150 before) -- this is
# expected and doesn't matter, since these values were never meant to carry
# absolute physical meaning, only relative/ranking meaning within one
# analysis (see github/demux_pool_design NOTES.md).
#
# Output: csv/designs/1kg_genomewide_dist_matrix.npz (D, X, sample_ids, pcs)
#         csv/designs/1kg_genomewide_sample_meta.tsv

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJDIR = Path('/u/project/cluo/terencew/claude/project_ideas/pool_design')
DESIGN_DIR = PROJDIR / 'csv' / 'designs'
SCRATCH = DESIGN_DIR / '_plink_genomewide'
VCF = Path('/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/1000G.merged.common_biallelic.vcf.gz')
PED = Path('/u/project/cluo/terencew/demux_benchmark/pool_design/csv/1000G/meta/20130606_g1k_3202_samples_ped_population.txt')
PLINK = 'plink'
FLASHPCA = '/u/project/cluo/terencew/programs/flashpca/flashpca_x86-64'

OUTMAT = DESIGN_DIR / '1kg_genomewide_dist_matrix.npz'
OUTMETA = DESIGN_DIR / '1kg_genomewide_sample_meta.tsv'

SEED = 42


def run(cmd, skip_if_exists=None, **kw):
    if skip_if_exists and Path(skip_if_exists).exists():
        print(f'{skip_if_exists} exists, skipping: ' + ' '.join(str(c) for c in cmd), flush=True)
        return
    print('$ ' + ' '.join(str(c) for c in cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def main():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    prefix = str(SCRATCH / 'genomewide')

    # Every PLINK step below is idempotent (skips if its output already exists) --
    # a prior run got all the way through PLINK and OOM-killed during the Python
    # parsing step afterward (see the parsing code below), so re-running the
    # (expensive, ~1.5h) PLINK steps would be pure waste. Delete the relevant
    # _plink_genomewide/ file yourself if you actually want a step to redo.
    print('Step 1/5: VCF -> PLINK binary', flush=True)
    run([PLINK, '--vcf', str(VCF), '--double-id', '--allow-extra-chr',
         '--vcf-half-call', 'missing', '--make-bed', '--out', prefix],
        skip_if_exists=prefix + '.bed')

    print('Step 2/5: LD pruning (--indep-pairwise 50 5 0.2, standard params)', flush=True)
    run([PLINK, '--bfile', prefix, '--indep-pairwise', '50', '5', '0.2',
         '--out', prefix + '_pruned'], skip_if_exists=prefix + '_pruned.prune.in')
    n_pruned = sum(1 for _ in open(prefix + '_pruned.prune.in'))
    print(f'{n_pruned} SNPs survive LD pruning', flush=True)

    print('Step 3/5: pairwise 1-IBS distance on the pruned set', flush=True)
    run([PLINK, '--bfile', prefix, '--extract', prefix + '_pruned.prune.in',
         '--distance', 'square', '1-ibs', 'flat-missing', '--out', prefix + '_dist'],
        skip_if_exists=prefix + '_dist.mdist')

    # --output-missing-genotype 0 makes every cell numeric (no 'NA' strings) so X can be parsed
    # without pandas at all below -- pandas' wide-dataframe path turned out to be the ACTUAL
    # OOM culprit, not sklearn's PCA as first suspected: swapping in flashpca for the PCA step
    # (see step 5) fixed that step (3m37s, vs. a 3+ hour OOM-kill before) but this job still
    # climbed past the previous run's fatal maxvmem (252G and rising after 34 min, on a job that
    # should need ~6G for X) at the exact same `pd.read_csv(..., usecols=range(6, n_cols))` line
    # below -- caught by watching `qstat -j` on the resubmitted job, not assumed fixed just
    # because the PCA half was. New output name (not overwriting `_dosage.raw`) since this flag
    # changes the file's content, not just how it's read.
    print('Step 4/6: additive dosage recode on the same pruned set (for X / KL)', flush=True)
    run([PLINK, '--bfile', prefix, '--extract', prefix + '_pruned.prune.in',
         '--recode', 'A', '--output-missing-genotype', '0', '--out', prefix + '_dosage_num'],
        skip_if_exists=prefix + '_dosage_num.raw')

    # PCA used to be sklearn's PCA(svd_solver='randomized').fit_transform(X) on the
    # 3202x460k dosage array in step 6 below -- that silently OOM-killed this job
    # (no traceback, log just stops after the D-shape print; maxvmem hit 338G against
    # a 128G request). flashpca is a compiled, purpose-built genotype-PCA tool that
    # streams the bed file in blocks instead of materializing a dense feature matrix
    # for the solver, so it does the same randomized-PCA computation in a fraction of
    # the memory/time. It needs its own pruned .bed/.bim/.fam (it takes --bfile, not
    # --extract), built once here and reused on every rerun.
    print('Step 5/6: PCA via flashpca on the pruned bfile', flush=True)
    run([PLINK, '--bfile', prefix, '--extract', prefix + '_pruned.prune.in',
         '--make-bed', '--out', prefix + '_pruned_bfile'],
        skip_if_exists=prefix + '_pruned_bfile.bed')
    # flashpca writes eigenvalues.txt/eigenvectors.txt/pve.txt to the CWD by
    # default regardless of --outpc's path -- pin all of them under SCRATCH
    # explicitly so a run doesn't scatter files into wherever it happened to
    # be launched from (caught by a stray-file test run during development).
    run([FLASHPCA, '--bfile', prefix + '_pruned_bfile', '-d', '10',
         '--seed', str(SEED), '--numthreads', os.environ.get('NSLOTS', '4'),
         '--outpc', prefix + '_pca.txt', '--outval', prefix + '_pca_eigenvalues.txt',
         '--outvec', prefix + '_pca_eigenvectors.txt', '--outpve', prefix + '_pca_pve.txt'],
        skip_if_exists=prefix + '_pca.txt')

    print('Step 6/6: parsing PLINK + flashpca outputs into .npz', flush=True)
    # 1-ibs + square writes <prefix>.mdist (matrix) and <prefix>.mdist.id (sample order)
    D = np.loadtxt(prefix + '_dist.mdist', dtype=np.float32)
    dist_ids = pd.read_csv(prefix + '_dist.mdist.id', sep='\t', header=None, names=['FID', 'IID'])
    sample_ids = dist_ids['IID'].tolist()
    print(f'D shape: {D.shape}  min={D[D>0].min():.4f}  max={D.max():.4f}  mean={D[D>0].mean():.4f}', flush=True)

    # .raw is sample-major: FID IID PAT MAT SEX PHENOTYPE <snp1> <snp2> ...
    # pandas.read_csv on a 3202-row x ~460,000-column frame is the ACTUAL OOM cause found this
    # session (maxvmem 252G+ and still climbing after 34 min on a resubmitted job, at this exact
    # line) -- pandas' wide-dataframe path has per-column overhead that doesn't show up in a
    # normal (many-rows, few-columns) table, no matter how tight the requested dtype is. Avoided
    # entirely here: read the file as plain text, one row at a time, straight into a
    # preallocated float32 array via np.fromstring (a single C-level call per row, not a Python
    # loop over 460k cells) -- fixed, predictable ~n_samples*n_snps*4 byte footprint (~6GB here),
    # not something to be reconfirmed by watching the resubmitted job again.
    with open(prefix + '_dosage_num.raw') as f:
        header = f.readline().split()
    n_cols = len(header)
    n_snps = n_cols - 6
    n_samples = len(sample_ids)
    raw_ids = []
    X = np.empty((n_samples, n_snps), dtype=np.float32)
    with open(prefix + '_dosage_num.raw') as f:
        next(f)  # header, already read above
        for i, line in enumerate(f):
            fields = line.split(' ', 6)
            raw_ids.append(fields[1])
            X[i] = np.fromstring(fields[6], dtype=np.float32, sep=' ')
    assert raw_ids == sample_ids, 'sample order mismatch between .mdist.id and .raw -- do not proceed'
    print(f'X shape: {X.shape}', flush=True)

    # flashpca's --outpc is FID IID PC1..PC10, in the pruned bfile's .fam order --
    # not guaranteed to match sample_ids, so align explicitly by IID rather than
    # assuming row order survived the plink --make-bed round-trip.
    pcs_df = pd.read_csv(prefix + '_pca.txt', sep='\t')
    pcs_df = pcs_df.set_index('IID').loc[sample_ids]
    PCS = pcs_df[[c for c in pcs_df.columns if c.startswith('PC')]].to_numpy(dtype=np.float32)
    print(f'PCS shape: {PCS.shape}', flush=True)

    np.savez_compressed(OUTMAT, D=D, X=X, sample_ids=np.array(sample_ids), pcs=PCS)
    print(f'saved {OUTMAT}', flush=True)

    meta = pd.read_csv(PED, sep=' ')
    meta['unrelated'] = (meta['FatherID'] == '0') & (meta['MotherID'] == '0')
    meta_vcf = pd.DataFrame({'SampleID': sample_ids, 'vcf_idx': np.arange(len(sample_ids))})
    meta_vcf = meta_vcf.merge(meta[['SampleID', 'Population', 'Superpopulation', 'unrelated']],
                               on='SampleID', how='left')
    meta_vcf.to_csv(OUTMETA, sep='\t', index=False)
    print(f'saved {OUTMETA}', flush=True)
    print(f'unrelated samples: {meta_vcf["unrelated"].sum()} / {len(meta_vcf)}', flush=True)


if __name__ == '__main__':
    main()
