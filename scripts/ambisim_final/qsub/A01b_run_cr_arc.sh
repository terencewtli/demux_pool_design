#!/bin/bash
#$ -N final_A01b_run_cr_arc
#$ -cwd
#$ -l h_data=4G,h_rt=24:00:00
#$ -pe shared 8
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim_final/logs/A01b_run_cr_arc.$JOB_ID.$TASK_ID
#$ -j y

# See A00_prep_pool_vcfs.sh in this directory for the DO NOT SUBMIT YET gating -- unchanged
# here. Mirrors scripts/ambisim_new/qsub/A01b_run_cr_arc.sh (24h h_rt from the start, per that
# script's own note about the n=8 grid's 14h-cap casualty wave), pointed at ambisim_final/.

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0

CRARC=/u/project/cluo/terencew/programs/cellranger-arc-2.0.1/bin/cellranger-arc
REF=/u/project/cluo/terencew/reference/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

LIST=$PROJDIR/ambisim/txt/pool_experiments.txt
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")
OUTDIR=$PROJDIR/ambisim_final/$POOL
SAMPLEDIR=$OUTDIR/$SAMPLE
CR_OUTDIR=$OUTDIR/cr_arc

echo "Start: $(date)  pool=$POOL"

if [ -f "$CR_OUTDIR/$SAMPLE/outs/gex_possorted_bam.bam" ]; then
    echo "$(date): CR-ARC output exists, skipping"
    exit 0
fi

RNA_FASTQ=$SAMPLEDIR/fastq/RNA
ATAC_FASTQ=$SAMPLEDIR/fastq/ATAC

LIB_CSV=$OUTDIR/setup_crarc.csv
cat > "$LIB_CSV" << EOF
fastqs,sample,library_type
$RNA_FASTQ,sim,Gene Expression
$ATAC_FASTQ,sim,Chromatin Accessibility
EOF

mkdir -p "$CR_OUTDIR"
cd "$CR_OUTDIR"

time $CRARC count \
    --id "$SAMPLE" \
    --reference "$REF" \
    --libraries "$LIB_CSV" \
    --localcores=8 \
    --localmem=30

echo "End: $(date)  pool=$POOL"
