#!/bin/bash
#PBS -N RE_RUN_ALL3_v2
#PBS -l select=1:ncpus=20:mem=1024gb:scratch_ssd=500gb:singularity=True
#PBS -l walltime=168:00:00
#PBS -q large_mem@pbs-m1.metacentrum.cz
#PBS -M blanka.jendriskova0000@seznam.cz
#PBS -m bea
#PBS -o /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/RE_RUN_ALL3_v2_Output/stdout.txt
#PBS -e /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/RE_RUN_ALL3_v2_Output/stderr.txt

# ============================================================
# RepeatExplorer2 - Three-sample comparative analysis
# AAZ0 + AAZW + PBrw combined
# Abraxas sylvata - Blanka Jendrisková
# Version 2 (21.02.2026)
#
# Changes from v1:
#   - --mincl changed from 0.001 to 0.01
#   - --assembly_min changed from 3 to 5
#   - Output saved to RE_RUN_ALL3_v2_Output
# ============================================================

trap 'clean_scratch' TERM EXIT

INPUT_FILE="AAZ0_AAZW_PBrw.fa"
INPUT_DIR="/storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs"
OUTPUT_DIR="${INPUT_DIR}/RE_RUN_ALL3_v2_Output"
REPO_IMAGE="library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e"

mkdir -p "$OUTPUT_DIR"

# 1. Copy input data to scratch SSD
echo "=== STEP 1: Copying data to scratch ==="
echo "Scratch: $SCRATCHDIR"
df -h "$SCRATCHDIR"
cp "${INPUT_DIR}/${INPUT_FILE}" "$SCRATCHDIR" || exit 1
cd "$SCRATCHDIR"

# 2. Pull Singularity image
echo "=== STEP 2: Pulling Singularity image ==="
export SINGULARITY_CACHEDIR="$SCRATCHDIR"
export SINGULARITY_LOCALCACHEDIR="$SCRATCHDIR"
export SINGULARITY_TMPDIR="$SCRATCHDIR"
export TMPDIR="$SCRATCHDIR"
singularity pull --arch amd64 "$REPO_IMAGE"
echo "Disk status after image pull:"
df -h "$SCRATCHDIR"

# 3. Run RepeatExplorer2
echo "=== STEP 3: Running RepeatExplorer2 ==="
echo "Input: ${INPUT_FILE}, prefix_length=5"
echo "Samples: AAZ0_ + AAZW_ + PBrw_ (3 samples, 2,676,720 sequences)"
echo "Start: $(date)"

singularity exec --bind "$SCRATCHDIR:/data" "$SCRATCHDIR/repex_tarean_0.3.12-7a7dc9e.sif" \
    seqclust --paired \
    --automatic_filtering \
    --mincl 0.01 \
    --cpu 18 \
    --prefix_length 5 \
    --output_dir=/data/output/ \
    --max_memory 900000000 \
    --cleanup \
    --assembly_min 5 \
    --taxon METAZOA3.0 \
    "/data/${INPUT_FILE}"

SEQCLUST_EXIT=$?
echo "Seqclust finished with exit code: $SEQCLUST_EXIT"
echo "End: $(date)"
df -h "$SCRATCHDIR"

# 4. Copy results back to storage
echo "=== STEP 4: Copying results ==="
if [ -d "$SCRATCHDIR/output" ]; then
    cp -r "$SCRATCHDIR/output"/* "$OUTPUT_DIR" || export CLEAN_SCRATCH=false
    echo "Results successfully copied to: $OUTPUT_DIR"
else
    echo "ERROR: output folder not found!"
    export CLEAN_SCRATCH=false
fi
echo "=== DONE ==="
