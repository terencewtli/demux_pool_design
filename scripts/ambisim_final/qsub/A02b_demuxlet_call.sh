#!/bin/bash
#$ -N final_A02b_demuxlet_call
#$ -cwd
#$ -l h_data=4G,h_rt=8:00:00
#$ -pe shared 10
#$ -t 1-126:1
#$ -tc 20
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim_final/logs/A02b_demuxlet_call.$JOB_ID.$TASK_ID
#$ -j y

# 126 pools (not 132) -- see A00_prep_pool_vcfs.sh header. Mirrors scripts/ambisim_new/qsub/
# A02b_demuxlet_call.sh, pointed at ambisim_final/.
#
# 2026-09-18 rewrite: demuxlet now runs *inline* in the array task instead of being handed to
# two child `qsub`s per pool. Rationale:
#   - one 126-task array instead of ~252 unpredictable child jobs (was impossible to hold_jid
#     on, impossible to count, and cluttered `qstat`);
#   - kills the "no submit host" failure mode (child qsub from a compute node silently refused,
#     task still exits 0) that lost 4 pools at the A02a stage -- see md/PROGRESS.md;
#   - the log for a pool's GEX and ATAC calls is now one file named after this script, instead
#     of landing in scripts/ambisim_final/qsub/logs/ under the template's own job name.
# The popscle invocation itself is byte-identical to template_demux/demuxlet.sh; that template
# is untouched and still used by the other trees.
#
# Memory: the demuxlet call needs >=32G (the template's 8G default silently OOM-killed ~27/132
# calls in the original n=8 grid). h_data=4G x 10 slots = 40G total, which keeps the per-slot
# request small enough to queue quickly. demuxlet is single-threaded -- the 10 slots are bought
# for the memory, not for parallelism. Check for `Killed` in the log rather than assuming 40G
# is always enough.
#
# Runtime: each call took 3-14 min against the same ~1.3G pool VCFs in the chr22 grid, so ~30
# min/task for GEX+ATAC sequentially. h_rt=8:00:00 is headroom, not an estimate.

source ~/.bashrc

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
SAMPLE=20220928-IGVF-D0
POPSCLE=/u/project/cluo/terencew/programs/popscle/bin/popscle

LIST=$PROJDIR/ambisim/txt/pool_experiments_final126.txt
# ID=1
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")

OUTDIR=$PROJDIR/ambisim_final/$POOL
CR_OUTS=$OUTDIR/cr_arc/$SAMPLE/outs
VCF=$PROJDIR/ambisim_final/vcf/${POOL}.vcf.gz
BARCODES=$CR_OUTS/filtered_feature_bc_matrix/barcodes.tsv.gz

DEMUX_DIR=$OUTDIR/demux/demuxlet
POPSCLE_GEX=$DEMUX_DIR/gex
POPSCLE_ATAC=$DEMUX_DIR/atac

GEX_PLP=$POPSCLE_GEX/$SAMPLE.pileup
ATAC_PLP=$POPSCLE_ATAC/$SAMPLE.pileup

echo "Start: $(date)  pool=$POOL  host=$(hostname -s)"

if [ ! -f "${GEX_PLP}.var.gz" ]; then
    echo "$(date): GEX pileup not found for $POOL, skipping (run A02a_demuxlet_pileup.sh first)"
    exit 0
fi
if [ ! -f "${ATAC_PLP}.var.gz" ]; then
    echo "$(date): ATAC pileup not found for $POOL, skipping (run A02a_demuxlet_pileup.sh first)"
    exit 0
fi

# both calls use the same RNA-space barcode list on purpose -- ipynb/ambisim/01a depends on the
# ATAC .best BARCODE column already being in RNA_BC space for the join against drop_data_rand.txt
if [ -f "$POPSCLE_GEX/$SAMPLE.best" ]; then
    echo "$(date): GEX demuxlet output exists, skipping"
else
    echo "$(date): GEX demuxlet start"
    time $POPSCLE demuxlet --plp "$GEX_PLP" \
      --group-list "$BARCODES" \
      --field GT \
      --vcf "$VCF" --out "$POPSCLE_GEX/$SAMPLE"
    echo "$(date): GEX demuxlet done"
fi

if [ -f "$POPSCLE_ATAC/$SAMPLE.best" ]; then
    echo "$(date): ATAC demuxlet output exists, skipping"
else
    echo "$(date): ATAC demuxlet start"
    time $POPSCLE demuxlet --plp "$ATAC_PLP" \
      --group-list "$BARCODES" \
      --field GT \
      --vcf "$VCF" --out "$POPSCLE_ATAC/$SAMPLE"
    echo "$(date): ATAC demuxlet done"
fi

echo "End: $(date)  pool=$POOL"
