#!/bin/bash
#$ -N B04_regenerate_pools_genomewide
#$ -cwd
#$ -l h_data=8G,h_rt=2:00:00
#$ -pe shared 2
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim/logs/B04_regenerate_pools_genomewide.$JOB_ID
#$ -j y

# Regenerates the full 132-pool grid's donor selections + design metrics against the
# genome-wide, LD-pruned distance matrix (1kg_genomewide_dist_matrix.npz) -- see
# scripts/ambisim/lib/regenerate_pools_genomewide.py's own docstring for exactly what is/isn't
# re-derived (adversarial_family/family_mixed donor lists are copied unchanged; their metrics
# are still recomputed).
#
# Single sequential process, not an array -- includes greedy_maxkl for 5 universes x 3 reps,
# estimated ~45-60 min total (measured 116s/pool on the 353-donor AMR universe this session;
# h_rt=2:00:00 is a safety margin, not a measured ceiling -- check
# `grep -a real logs/B04_regenerate_pools_genomewide.*` after this finishes).
#
# Writes to csv/designs/nominated_pools_n8_genomewide.tsv -- the SAME path
# create_pools_genomewide.py (a narrower, redundant script written before this one was found)
# also wrote to; this script's output supersedes that one.

export PATH="/u/home/t/terencew/project-cluo/miniconda3/envs/allcools/bin:$PATH"
set -euo pipefail

cd /u/project/cluo/terencew/claude/project_ideas/pool_design/scripts/ambisim/lib
echo "Start: $(date)"
time python3 regenerate_pools_genomewide.py
echo "End: $(date)"
