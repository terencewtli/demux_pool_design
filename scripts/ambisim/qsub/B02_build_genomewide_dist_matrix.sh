#!/bin/bash
#$ -N B02_build_genomewide_dist_matrix
#$ -cwd
#$ -l h_data=16G,h_rt=8:00:00
#$ -pe shared 4
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim/logs/B02_build_genomewide_dist_matrix.$JOB_ID
#$ -j y

# Build the genome-wide replacement for csv/designs/1kg_chr22_dist_matrix.npz
# -- see scripts/ambisim/lib/build_genomewide_dist_matrix.py for full
# rationale (chr22-only, 20k-SNP distance matrix found to underlie every
# min_dist/mean_dist value in this project; github/demux_pool_design
# NOTES.md documents why this matters).
#
# h_rt=8:00:00 is a first-pass estimate (not measured on this data) --
# the heaviest step is extracting genotypes for 200,000 sites x 3202
# samples via bcftools query, similar order of magnitude to a single pool's
# VCF subsetting (A00_prep_pool_vcfs.sh) but over ~30x more variants;
# check `grep -a real logs/B02_build_genomewide_dist_matrix.*` after this
# finishes and adjust for any future rerun.
#
# Submit with: qsub scripts/ambisim/qsub/B02_build_genomewide_dist_matrix.sh

source ~/.bashrc
conda activate allcools

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design

echo "Start: $(date)"
time python3 $PROJDIR/scripts/ambisim/lib/build_genomewide_dist_matrix.py
echo "End: $(date)"
