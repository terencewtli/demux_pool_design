#!/bin/bash
#$ -N A02a_demuxlet_pileup
#$ -cwd
#$ -l h_data=1G,h_rt=1:00:00
#$ -pe shared 1
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/logs/A02a_demuxlet_pileup.$JOB_ID.$TASK_ID
#$ -j y

# Stage 1 of the demuxlet-only pipeline: filter each pool's GEX/ATAC BAM for
# popscle, sort the pool's own ground-truth VCF to match the BAM contig order,
# and run popscle dsc-pileup. RNA and ATAC are queued as separate downstream
# jobs (popscle_gex.sh / popscle_atac.sh), same as A01c_plp.sh does in the
# IGVF pilot pipeline, but restricted to demuxlet -- no cellsnp/varcon/scsplit/
# scavengers, since only demuxlet was asked for here.
#
# Requires Cell Ranger ARC output (A01b_run_cr_arc.sh) to already exist for
# the pool; a pool without it is reported and skipped, not treated as an error,
# since cr_arc has not finished for most pools yet.
#
# Submit with: qsub scripts/ambisim/qsub/A02a_demuxlet_pileup.sh
# Or for a single pool (by its line number in pool_experiments.txt):
#   qsub -t <N> scripts/ambisim/qsub/A02a_demuxlet_pileup.sh

source ~/.bashrc

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
TEMPLATE=/u/project/cluo/terencew/demux_benchmark/template_demux

LIST=$PROJDIR/ambisim/txt/pool_experiments.txt
# ID=10
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")

OUTDIR=$PROJDIR/ambisim/$POOL
CR_OUTS=$OUTDIR/cr_arc/$SAMPLE/outs
VCF=$PROJDIR/ambisim/vcf/${POOL}.vcf.gz

GEX_BAM=$CR_OUTS/gex_possorted_bam.bam
ATAC_BAM=$CR_OUTS/atac_possorted_bam.bam
BARCODES=$CR_OUTS/filtered_feature_bc_matrix/barcodes.tsv.gz

DEMUX_DIR=$OUTDIR/demux/demuxlet
POPSCLE_GEX=$DEMUX_DIR/gex
POPSCLE_ATAC=$DEMUX_DIR/atac

echo "Start: $(date)  pool=$POOL"

if [ ! -f "$GEX_BAM" ] || [ ! -f "$ATAC_BAM" ] || [ ! -f "$BARCODES" ]; then
    echo "$(date): cr_arc output not found for $POOL, skipping (run A01b_run_cr_arc.sh first)"
    echo "  expected: $GEX_BAM"
    echo "  expected: $ATAC_BAM"
    echo "  expected: $BARCODES"
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
