#!/bin/bash
#$ -N A02b_demuxlet_call
#$ -cwd
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim/logs/A02b_demuxlet_call.$JOB_ID.$TASK_ID
#$ -j y

# Stage 2 of the demuxlet-only pipeline: run popscle demuxlet on the GEX and
# ATAC pileups from A02a_demuxlet_pileup.sh, separately, against the pool's
# own ground-truth VCF (its 8 donors) -- same demuxlet.sh template used by
# A01d_call_demux.sh in the IGVF pilot pipeline.
#
# Run after A02a's pileup jobs have finished for the pool (this script does
# not hold_jid on them, since A02a fans out two independent qsub jobs per
# pool with unpredictable job IDs -- check `qstat -u $USER` before running).
#
# Submit with: qsub scripts/ambisim/qsub/A02b_demuxlet_call.sh
# Or for a single pool: qsub -t <N> scripts/ambisim/qsub/A02b_demuxlet_call.sh

source ~/.bashrc

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
TEMPLATE=/u/project/cluo/terencew/demux_benchmark/template_demux

LIST=$PROJDIR/ambisim/txt/pool_experiments.txt
# ID=1
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")

OUTDIR=$PROJDIR/ambisim/$POOL
CR_OUTS=$OUTDIR/cr_arc/$SAMPLE/outs
VCF=$PROJDIR/ambisim/vcf/${POOL}.vcf.gz
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
    qsub "$TEMPLATE/demuxlet.sh" "$BARCODES" "$VCF" "$GEX_PLP" "$POPSCLE_GEX/$SAMPLE"
fi

if [ -f "$POPSCLE_ATAC/$SAMPLE.best" ]; then
    echo "$(date): ATAC demuxlet output exists, skipping"
else
    qsub "$TEMPLATE/demuxlet.sh" "$BARCODES" "$VCF" "$ATAC_PLP" "$POPSCLE_ATAC/$SAMPLE"
fi

echo "End: $(date)  pool=$POOL"
