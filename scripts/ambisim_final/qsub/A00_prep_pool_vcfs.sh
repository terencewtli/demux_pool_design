#!/bin/bash
#$ -N final_A00_prep_pool_vcfs
#$ -cwd
#$ -l h_data=2G,h_rt=4:00:00
#$ -pe shared 2
#$ -t 1-132:1
#$ -o /u/project/cluo/terencew/claude/project_ideas/pool_design/ambisim_final/logs/A00_prep_pool_vcfs.$JOB_ID.$TASK_ID
#$ -j y

# "Final" tree: the original 132-pool grid, re-simulated with donor selections regenerated
# against the genome-wide LD-pruned distance matrix (scripts/ambisim/lib/
# regenerate_pools_genomewide.py, job 14752950) instead of the chr22-only one every existing
# ambisim/ pool was built from. Mirrors scripts/ambisim/qsub/A00_prep_pool_vcfs.sh, with two
# differences: donor lists come from txt/donors_genomewide/ (not txt/donors/), and output goes
# to ambisim_final/ (not ambisim/), so this can't collide with or silently get mixed into the
# existing chr22-based grid's results.
#
# DO NOT SUBMIT YET as of 2026-09-15 -- gated behind (1) job 14752950 finishing and its
# donor-selection changes being reviewed (csv/designs/genomewide_vs_chr22_donor_overlap.tsv),
# and (2) the EUR_only orthogonal-sampling pilot (ambisim_new/) actually showing whether
# min_dist/mean_dist separation matters for demux accuracy at all -- see JOURNAL.md 2026-09-15
# entry and NOTES.md's "Recommended path forward" section. A full 132-pool resimulation is a
# multi-week, heavy-disk undertaking (/u/project/cluo was at 95%+ capacity as of 2026-09-14);
# don't launch it on the strength of "the matrix changed" alone.
#
# Output: ambisim_final/vcf/{pool_id}.vcf.gz (+.tbi)

set -euo pipefail

PROJDIR=/u/project/cluo/terencew/claude/project_ideas/pool_design
VCF_MERGED=/u/project/cluo/terencew/demux_benchmark/pool_design/vcf/1000G/1000G.merged.common_biallelic.vcf.gz
BCFTOOLS=/u/local/apps/bcftools/1.11/gcc-4.8.5/bin/bcftools
TABIX=/u/local/apps/htslib/1.12/gcc-4.8.5/bin/tabix
BGZIP=/u/local/apps/htslib/1.12/gcc-4.8.5/bin/bgzip

LIST=$PROJDIR/ambisim/txt/pool_experiments.txt
ID=$SGE_TASK_ID
POOL=$(sed -n "${ID}p" "$LIST")

DONORS_FILE=$PROJDIR/txt/donors_genomewide/${POOL}.txt
OUT=$PROJDIR/ambisim_final/vcf/${POOL}.vcf.gz

echo "Start: $(date)  pool=$POOL"

if [ -f "$OUT" ] && [ -f "${OUT}.tbi" ]; then
    echo "$(date): $OUT exists, skipping"
    exit 0
fi

DONORS_CSV=$(paste -sd',' "$DONORS_FILE")
echo "$(date): subsetting to donors: $DONORS_CSV"

mkdir -p "$PROJDIR/ambisim_final/vcf"
time $BCFTOOLS view --threads 4 -s "$DONORS_CSV" "$VCF_MERGED" \
    | $BCFTOOLS view -g ^miss \
    | $BGZIP -@ 4 > "$OUT"
$TABIX -p vcf "$OUT"

N=$($BCFTOOLS view -H "$OUT" | wc -l)
echo "$(date): $N variants in $OUT  ($(du -sh $OUT | cut -f1))"
echo "End: $(date)"
