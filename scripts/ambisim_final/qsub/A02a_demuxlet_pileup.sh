#!/bin/bash
#$ -N final_A02a_demuxlet_pileup
#$ -cwd
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim_final/logs/A02a_demuxlet_pileup.$JOB_ID.$TASK_ID
#$ -j y

# See A00_prep_pool_vcfs.sh in this directory for the DO NOT SUBMIT YET gating -- unchanged
# here. Mirrors scripts/ambisim_new/qsub/A02a_demuxlet_pileup.sh, pointed at ambisim_final/.
# Fans out two independent qsub jobs per pool (GEX + ATAC popscle pileup) -- unpredictable job
# IDs, so A02b below doesn't hold_jid on them by name pattern automatically; check
# `qstat -u $USER` (or hold_jid on the job names A01a_popscle_gex,A01a_popscle_atac) before
# submitting A02b.

source ~/.bashrc

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
TEMPLATE=/u/project/cluo/terencew/demux_benchmark/template_demux

LIST=$PROJDIR/ambisim/txt/pool_experiments.txt
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")

OUTDIR=$PROJDIR/ambisim_final/$POOL
CR_OUTS=$OUTDIR/cr_arc/$SAMPLE/outs
VCF=$PROJDIR/ambisim_final/vcf/${POOL}.vcf.gz

GEX_BAM=$CR_OUTS/gex_possorted_bam.bam
ATAC_BAM=$CR_OUTS/atac_possorted_bam.bam
BARCODES=$CR_OUTS/filtered_feature_bc_matrix/barcodes.tsv.gz

DEMUX_DIR=$OUTDIR/demux/demuxlet
POPSCLE_GEX=$DEMUX_DIR/gex
POPSCLE_ATAC=$DEMUX_DIR/atac

echo "Start: $(date)  pool=$POOL"

if [ ! -f "$GEX_BAM" ] || [ ! -f "$ATAC_BAM" ] || [ ! -f "$BARCODES" ]; then
    echo "$(date): cr_arc output not found for $POOL, skipping (run A01b_run_cr_arc.sh first)"
    exit 0
fi

if [ ! -f "$VCF" ]; then
    echo "ERROR: pool VCF not found: $VCF" >&2
    exit 1
fi

mkdir -p "$POPSCLE_GEX" "$POPSCLE_ATAC"

if [ -f "$POPSCLE_GEX/$SAMPLE.pileup.var.gz" ]; then
    echo "$(date): GEX pileup exists, skipping"
else
    qsub "$TEMPLATE/popscle_gex.sh" "$GEX_BAM" "$BARCODES" "$SAMPLE" "$VCF" "$POPSCLE_GEX"
fi

if [ -f "$POPSCLE_ATAC/$SAMPLE.pileup.var.gz" ]; then
    echo "$(date): ATAC pileup exists, skipping"
else
    qsub "$TEMPLATE/popscle_atac.sh" "$ATAC_BAM" "$BARCODES" "$SAMPLE" "$VCF" "$POPSCLE_ATAC"
fi

echo "End: $(date)  pool=$POOL"
