#!/bin/bash
#PBS -N RE_RUN_ALL3
#PBS -l select=1:ncpus=20:mem=1024gb:scratch_ssd=500gb:singularity=True
#PBS -l walltime=168:00:00
#PBS -q large_mem@pbs-m1.metacentrum.cz
#PBS -M blanka.jendriskova0000@seznam.cz
#PBS -m bea
#PBS -o /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/RE_RUN_ALL3_Output/stdout.txt
#PBS -e /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/RE_RUN_ALL3_Output/stderr.txt

# ============================================================
# RepeatExplorer2 - Three-sample comparative analysis
# AAZ0 + AAZW + PBrw combined
# Abraxas sylvata - Blanka Jendrisková
# Version 1 (16.02.2026)
#
# BACKGROUND:
# Prior to running this script, three datasets were prepared
# and tagged with unique prefixes so that RepeatExplorer2
#
# 1. AAZ0 (male genome, no W chromosome):
#    - Illumina-like paired-end reads simulated from “male“
#      genome assembly (AA_Z_chr.fasta) using randomreads.sh
#      (BBTools v39.01, read length 100 bp, 20M pairs,
#      seed=42, no error simulation)
#    - Subsampled to 438,646 pairs (0.4x coverage of 463 Mb)
#      using seqtk sample (seed=42)
#    - Forward and reverse reads merged using seqtk mergepe
#      and converted to FASTA using seqtk seq -a
#    - Tagged with prefix AAZ0_ using sed
#
# 2. AAZW (female genome, with W chromosome):
#    - Same simulation approach as AAZ0, applied to female
#      genome assembly (AA_Z_Wasmbl.fasta)
#    - Subsampled to 449,857 pairs (0.4x coverage of 475 Mb)
#    - Tagged with prefix AAZW_ using sed
#
# 3. PBrw (raw female PacBio HiFi reads):
#    - Original PacBio reads (ERR12370385.fastq.gz)
#      subsampled to 2x genome coverage (~14M reads)
#      using seqtk sample (seed=42)
#    - Illumina-like reads simulated from this subset
#      using randomreads.sh (same parameters as above)
#    - Subsampled to 449,857 pairs (0.4x coverage of 475 Mb)
#    - Tagged with prefix PBrw_ using sed
#
# All three tagged FASTA files were concatenated:
#    cat AAZ0_tagged.fa AAZW_tagged.fa PBrw_tagged.fa \
#        > AAZ0_AAZW_PBrw.fa
#    Total input: 2,676,720 sequences
# ============================================================

trap 'clean_scratch' TERM EXIT

INPUT_FILE="AAZ0_AAZW_PBrw.fa"
INPUT_DIR="/storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs"
OUTPUT_DIR="${INPUT_DIR}/RE_RUN_ALL3_Output"
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

# 3. Run RepeatExplorer2 in comparative mode
# --paired:               input is interleaved paired-end reads
# --automatic_filtering:  filter low-complexity sequences
# --mincl 0.001:          minimum cluster size (0.001% of reads)
# --cpu 18:               number of CPU cores
# --prefix_length 5:      length of sample prefix (AAZ0_, AAZW_, PBrw_)
# --max_memory 900000000: max memory in kB (~900 GB)
# --assembly_min 3:       minimum reads for contig assembly
# --taxon METAZOA3.0:     reference database for annotation
echo "=== STEP 3: Running RepeatExplorer2 ==="
echo "Input: ${INPUT_FILE}, prefix_length=5"
echo "Samples: AAZ0_ + AAZW_ + PBrw_ (3 samples, 2,676,720 sequences)"
echo "Start: $(date)"

singularity exec --bind "$SCRATCHDIR:/data" "$SCRATCHDIR/repex_tarean_0.3.12-7a7dc9e.sif" \
    seqclust --paired \
    --automatic_filtering \
    --mincl 0.001 \
    --cpu 18 \
    --prefix_length 5 \
    --output_dir=/data/output/ \
    --max_memory 900000000 \
    --cleanup \
    --assembly_min 3 \
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






OPRAVA SKRIPTU : 
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
