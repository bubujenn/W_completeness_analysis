#!/bin/bash
# =============================================================================
# RE.sh - RepeatExplorer High Memory 
# Testing first dataset: AAZ0 vs PBrw (Male vs PacBio)
# =============================================================================

# Define input file and tag length
INPUT_FILE="AAZ0_vs_PBrw.fa"
TAG_LENGTH=5

# Check if file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file '$INPUT_FILE' not found in current directory!"
    exit 1
fi

echo "Preparing PBS job for: $INPUT_FILE"

# Create output directory for this analysis
JOB_ID="RE_RUN_AAZ0"
OUT_DIR="$PWD/${JOB_ID}_Output"
mkdir -p "$OUT_DIR"

PBS_SCRIPT="${OUT_DIR}/${JOB_ID}.pbs"

# Generate PBS script dynamically (based on 'Vadovick' high-mem template)
cat <<EOF > "$PBS_SCRIPT"
#!/bin/bash
#PBS -N $JOB_ID
#PBS -l select=1:ncpus=20:mem=1024gb:scratch_local=1400gb:singularity=True
#PBS -l walltime=144:00:00
#PBS -q large_mem@pbs-m1.metacentrum.cz
#PBS -M blanka.jendriskova0000@seznam.cz
#PBS -m bea
#PBS -o ${OUT_DIR}/stdout.txt
#PBS -e ${OUT_DIR}/stderr.txt

# Clean scratch on exit
trap 'clean_scratch' TERM EXIT

# 1. Copy input data to scratch (fast local disk)
echo "Copying data to scratch..."
cp "$PWD/$INPUT_FILE" "\$SCRATCHDIR" || exit 1
cd "\$SCRATCHDIR"

# 2. Setup Singularity (RE image)
export SINGULARITY_CACHEDIR="\$SCRATCHDIR"
export SINGULARITY_LOCALCACHEDIR="\$SCRATCHDIR"
export SINGULARITY_TMPDIR="\$SCRATCHDIR"
export TMPDIR="\$SCRATCHDIR"

REPO_IMAGE="library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e"
echo "Pulling Singularity image..."
singularity pull --arch amd64 "\$REPO_IMAGE"

# 3. Run RepeatExplorer Analysis (High Memory: 1024GB)
echo "Running seqclust on $INPUT_FILE with prefix length $TAG_LENGTH..."

singularity exec --bind "\$SCRATCHDIR:/data" "\$SCRATCHDIR/repex_tarean_0.3.12-7a7dc9e.sif" \\
seqclust --paired \\
--automatic_filtering --mincl 0.001 --cpu 18 --prefix_length $TAG_LENGTH \\
--output_dir=/data/output/ --max_memory 900000000 --cleanup \\
--assembly_min 3 --taxon METAZOA3.0 "/data/$INPUT_FILE"

# 4. Copy results back to home directory
echo "Analysis complete. Copying results back..."
cp -r "\$SCRATCHDIR/output"/* "$OUT_DIR" || export CLEAN_SCRATCH=false
echo "Done."
EOF

# Submit the job
echo "Submitting job to Metacentrum..."
qsub "$PBS_SCRIPT"

echo "---------------------------------------------------"
echo "Job '$JOB_ID' submitted successfully!"
echo "Check status with: qstat -u jendrb00"
echo "Results will be in: $OUT_DIR"
