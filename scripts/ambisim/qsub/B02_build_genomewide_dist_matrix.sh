#!/bin/bash
#$ -N B02_build_genomewide_dist_matrix
#$ -cwd
#$ -l h_data=32G,h_rt=8:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim/logs/B02_build_genomewide_dist_matrix.$JOB_ID
#$ -j y

# Build the genome-wide, LD-pruned replacement for
# csv/designs/1kg_chr22_dist_matrix.npz using PLINK -- see
# scripts/ambisim/lib/build_genomewide_dist_matrix_plink.py for full
# rationale (chr22-only, 20k-SNP distance matrix found to underlie every
# min_dist/mean_dist value in this project; a first attempt used a raw
# 200,000-SNP subsample without LD-pruning, corrected here per user
# feedback -- see github/demux_pool_design JOURNAL.md 2026-09-14 entries).
#
# h_rt=8:00:00 is a first-pass estimate -- PLINK's VCF import and distance
# computation are compiled and fast even at millions of variants, but this
# hasn't been measured on this exact file; check
# `grep -a real logs/B02_build_genomewide_dist_matrix.*` after this finishes.
#
# Submit with: qsub scripts/ambisim/qsub/B02_build_genomewide_dist_matrix.sh

source ~/.bashrc
conda activate allcools

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
export PATH=/u/local/apps/plink/1.90b624:$PATH

echo "Start: $(date)"
time python3 $PROJDIR/scripts/ambisim/lib/build_genomewide_dist_matrix_plink.py
echo "End: $(date)"
