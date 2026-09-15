#!/bin/bash
#$ -N final_A01a_run_simulation
#$ -cwd
#$ -l h_data=4G,h_rt=12:00:00
#$ -pe shared 8
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim_final/logs/A01a_run_simulation.$JOB_ID.$TASK_ID
#$ -j y

# See A00_prep_pool_vcfs.sh in this directory for the DO NOT SUBMIT YET gating -- unchanged
# here. Mirrors scripts/ambisim/qsub/A01a_run_simulation.sh / ambisim_new's variant, pointed at
# ambisim_final/ and txt/donors_genomewide/.

source ~/.bashrc
conda activate mapping

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
PROP=20   # ambient contamination % -- matched to the n=8 grid

LIST=$PROJDIR/ambisim/txt/pool_experiments.txt
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")
DONORS_FILE=$PROJDIR/txt/donors_genomewide/${POOL}.txt
VCF=$PROJDIR/ambisim_final/vcf/${POOL}.vcf.gz
OUTDIR=$PROJDIR/ambisim_final/$POOL

mkdir -p "$OUTDIR"

echo "Start: $(date)  pool=$POOL"
DROP=$PROJDIR/scripts/ambisim/lib/setup_drop.sh
AMBISIM=$PROJDIR/scripts/ambisim/lib/run_ambisim.sh

time bash $DROP "$OUTDIR" "$SAMPLE" "$DONORS_FILE" "$PROP"
time bash $AMBISIM "$OUTDIR" "$SAMPLE" "$DONORS_FILE" "$VCF"

echo "End: $(date)  pool=$POOL"
