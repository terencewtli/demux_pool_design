#!/bin/bash
# SUPERSEDED (2026-09-15): scripts/ambisim/lib/regenerate_pools_genomewide.py (written in an
# undocumented 2026-09-14/15 session, found while investigating min_dist/mean_dist entanglement
# on the genome-wide matrix) already regenerates ALL 132 pools -- including greedy_maxkl for the
# 5 single-ancestry universes -- against 1kg_genomewide_dist_matrix.npz, writing to the SAME
# csv/designs/nominated_pools_n8_genomewide.tsv this script's caller (create_pools_genomewide.py)
# also wrote to. Running both would race/overwrite each other. Job 14752676 (this script) was
# qdel'ed once this was found; use scripts/ambisim/qsub/B04_regenerate_pools_genomewide.sh
# instead. Left in place rather than deleted -- rm it yourself if you want it gone.
#$ -N B03_greedy_maxkl_genomewide
#$ -cwd
#$ -l h_data=8G,h_rt=3:00:00
#$ -pe shared 2
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim/logs/B03_greedy_maxkl_genomewide.$JOB_ID
#$ -j y

# Genome-wide greedy_maxkl nomination for the 5 single-ancestry universes -- see
# scripts/ambisim/lib/run_greedy_maxkl_genomewide.py for scope/rationale. Timed at
# ~116s per 8-donor pool on the smallest universe (AMR, n=353); estimated ~70-90 min
# total across 5 universes x 5 reps. h_rt=3:00:00 is a safety margin over that estimate,
# not a measured ceiling -- check `grep -a real logs/B03_greedy_maxkl_genomewide.*` after
# this finishes and correct the comment above if it's off.
#
# Requires create_pools_genomewide.py to have already been run (appends to the files
# it created).

source ~/.bashrc
conda activate allcools

set -euo pipefail

cd /u/project/cluo/terencew/claude/project_ideas/pool_design/scripts/ambisim/lib
echo "Start: $(date)"
time python3 run_greedy_maxkl_genomewide.py
echo "End: $(date)"
