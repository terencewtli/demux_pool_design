#!/bin/bash
#$ -N n16_A02b_demuxlet_call
#$ -cwd
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-4:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim_n16/logs/A02b_demuxlet_call.$JOB_ID.$TASK_ID
#$ -j y

# Stage 2 of the demuxlet-only pipeline for the n=16 mini-experiment pools -- adapted from
# scripts/ambisim_new/qsub/A02b_demuxlet_call.sh (same structure, different pool list/paths;
# no scripts/ambisim_n16/qsub/ existed yet even though ambisim_n16/'s cr_arc + ATAC pileup
# output already exists, so whatever ran those steps wasn't checked into scripts/ here).
# Same 32GB call-memory override for the same reason (the shared demuxlet.sh template's 8GB
# default silently OOM-killed ~27/132 calls in the n=8 grid, see NOTES.md).
#
# Run after A02a's pileup jobs have finished for the pool (skips gracefully, not an error, if
# GEX/ATAC pileup isn't there yet -- safe to submit before pileups finish).
#
# Submit with: qsub scripts/ambisim_n16/qsub/A02b_demuxlet_call.sh

source ~/.bashrc

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
TEMPLATE=/u/project/cluo/terencew/demux_benchmark/template_demux
CALL_MEM=32G

LIST=$PROJDIR/txt/n16_pool_experiments.txt
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")

OUTDIR=$PROJDIR/ambisim_n16/$POOL
CR_OUTS=$OUTDIR/cr_arc/$SAMPLE/outs
VCF=$PROJDIR/ambisim_n16/vcf/${POOL}.vcf.gz
BARCODES=$CR_OUTS/filtered_feature_bc_matrix/barcodes.tsv.gz

DEMUX_DIR=$OUTDIR/demux/demuxlet
POPSCLE_GEX=$DEMUX_DIR/gex
POPSCLE_ATAC=$DEMUX_DIR/atac

GEX_PLP=$POPSCLE_GEX/$SAMPLE.pileup
ATAC_PLP=$POPSCLE_ATAC/$SAMPLE.pileup

echo "Start: $(date)  pool=$POOL"

if [ ! -f "${GEX_PLP}.var.gz" ]; then
    echo "$(date): GEX pileup not found for $POOL, skipping (run A02a_demuxlet_pileup.sh first)"
    echo "  expected: ${GEX_PLP}.var.gz"
    exit 0
fi
if [ ! -f "${ATAC_PLP}.var.gz" ]; then
    echo "$(date): ATAC pileup not found for $POOL, skipping (run A02a_demuxlet_pileup.sh first)"
    echo "  expected: ${ATAC_PLP}.var.gz"
    exit 0
fi

if [ -f "$POPSCLE_GEX/$SAMPLE.best" ]; then
    echo "$(date): GEX demuxlet output exists, skipping"
else
    qsub -l h_data=$CALL_MEM "$TEMPLATE/demuxlet.sh" "$BARCODES" "$VCF" "$GEX_PLP" "$POPSCLE_GEX/$SAMPLE"
fi

if [ -f "$POPSCLE_ATAC/$SAMPLE.best" ]; then
    echo "$(date): ATAC demuxlet output exists, skipping"
else
    qsub -l h_data=$CALL_MEM "$TEMPLATE/demuxlet.sh" "$BARCODES" "$VCF" "$ATAC_PLP" "$POPSCLE_ATAC/$SAMPLE"
fi

echo "End: $(date)  pool=$POOL"
