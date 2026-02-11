#!/bin/bash
# =============================================================================
# RepeatExplorer Pipeline 
# =============================================================================
# 1. Simulation from reference genomes (AAZ, AAZW) and PacBio raw data
# 2. Subsample to normalized coverage 0.4x
# 3. Interleave and rename reads
# =============================================================================

# Working directory and modules
WORK_DIR="/storage/brno2/home/jendrb00/20_10_consultation"
cd "$WORK_DIR" || { echo "Directory does not exist"; exit 1; }

module load bbmap
module load seqtk

# =============================================================================
# 1. SIMULATION (randomreads.sh) - PAIRED=t, ADDERRORS=f
# We want perfect reads without errors
# =============================================================================

# --- AAZ (male genome, NO W chromosome) ---
echo "1. Simulating AAZ..."
randomreads.sh \
    ref=AA_Z_chr.fasta \
    out1=AAZ_R1.fastq.gz \
    out2=AAZ_R2.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f || { echo "Error in AAZ simulation"; exit 2; }

# --- AAZW (female genome, WITH W chromosome) ---
echo "2. Simulating AAZW..."
randomreads.sh \
    ref=AA_Z_Wasmbl.fasta \
    out1=AAZW_R1.fastq.gz \
    out2=AAZW_R2.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f || { echo "Error in AAZW simulation"; exit 3; }

# --- PBrw (PacBio female data) ---
echo "3. Preparing PacBio..."
RAW_PACBIO_FILE="/storage/brno2/home/jendrb00/PacBiok100/ERR12370385.fastq.gz"

# Random sample from raw data first (approx 1.4 Gb = 2x genome)
# Estimated average read length ~7kb -> 200k reads
seqtk sample -s42 "$RAW_PACBIO_FILE" 200000 > PB_2xgenome_raw.fastq
seqtk seq -a PB_2xgenome_raw.fastq > PB_2xgenome_raw.fasta || { echo "Error in PacBio sample"; exit 4; }

# Simulation from these long reads
echo "4. Simulating PBrw..."
randomreads.sh \
    ref=PB_2xgenome_raw.fasta \
    out1=PBrw_R1_new.fastq.gz \
    out2=PBrw_R2_new.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f || { echo "Error in PBrw simulation"; exit 5; }

# =============================================================================
# 2. SUBSAMPLE (to coverage 0.4x)
# =============================================================================
echo "--- Subsampling ---"

# Create folder for final inputs
INPUT_DIR="$WORK_DIR/RE_final_inputs"
mkdir -p "$INPUT_DIR"
cd "$INPUT_DIR" || exit 6

# AAZ (smaller genome -> fewer reads: 438646 pairs)
seqtk sample -s42 ../AAZ_R1.fastq.gz 438646 > AAZ_R1.sampled.fq
seqtk sample -s42 ../AAZ_R2.fastq.gz 438646 > AAZ_R2.sampled.fq

# AAZW and PBrw (larger genome with W -> more reads: 449857 pairs)
seqtk sample -s42 ../AAZW_R1.fastq.gz 449857 > AAZW_R1.sampled.fq
seqtk sample -s42 ../AAZW_R2.fastq.gz 449857 > AAZW_R2.sampled.fq

seqtk sample -s42 ../PBrw_R1_new.fastq.gz 449857 > PBrw_R1.sampled.fq
seqtk sample -s42 ../PBrw_R2_new.fastq.gz 449857 > PBrw_R2.sampled.fq

# Check if files exist
if [ ! -s PBrw_R2.sampled.fq ]; then echo "Error in subsampling"; exit 7; fi

# =============================================================================
# 3. INTERLEAVE + TAGGING
# =============================================================================
echo "--- Interleave and tagging ---"

# AAZ
seqtk mergepe AAZ_R1.sampled.fq AAZ_R2.sampled.fq > AAZ_interleaved.fastq
seqtk seq -a AAZ_interleaved.fastq > AAZ_interleaved_new.fa
# Shorten names to first word and add tag AAZ0_
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZ_interleaved_new.fa | sed 's/>/>AAZ0_/' > AAZ0_tagged.fa

# AAZW
seqtk mergepe AAZW_R1.sampled.fq AAZW_R2.sampled.fq > AAZW_interleaved.fastq
seqtk seq -a AAZW_interleaved.fastq > AAZW_interleaved_new.fa
# Tag AAZW_
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZW_interleaved_new.fa | sed 's/>/>AAZW_/' > AAZW_tagged.fa

# PBrw
seqtk mergepe PBrw_R1.sampled.fq PBrw_R2.sampled.fq > PBrw_interleaved.fastq
seqtk seq -a PBrw_interleaved.fastq > PBrw_interleaved_new.fa
# Tag PBrw_
awk '/^>/{split($0,a," "); print a[1]; next} {print}' PBrw_interleaved_new.fa | sed 's/>/>PBrw_/' > PBrw_tagged.fa

# =============================================================================
# 4. CONCATENATION AND RUNNING REPEAT EXPLORER
# =============================================================================
echo "--- Final concatenation ---"

# Dataset 1: W repeats composition (AAZ vs PBrw)
cat AAZ0_tagged.fa PBrw_tagged.fa > AAZ0_vs_PBrw.fa

# Dataset 2: W completeness (AAZW vs PBrw)
cat AAZW_tagged.fa PBrw_tagged.fa > AAZW_vs_PBrw.fa

# Check counts
echo "Checking sequence counts (AAZ0 expects ~1.77M, AAZW expects ~1.80M):"
c1=$(grep -c '^>' AAZ0_vs_PBrw.fa)
c2=$(grep -c '^>' AAZW_vs_PBrw.fa)
echo "AAZ0_vs_PBrw: $c1"
echo "AAZW_vs_PBrw: $c2"

if [ "$c1" -eq 0 ] || [ "$c2" -eq 0 ]; then
    echo "Error: Output files are empty!"; exit 11;
fi

# Run RE
echo "--- Starting RepeatExplorer ---"
RE_SCRIPT="../RE_clustering.sh"  # Script is in parent dir

if [ ! -f "$RE_SCRIPT" ]; then
    echo "Error: Script $RE_SCRIPT not found!"; exit 12;
fi

# Run both in background
bash "$RE_SCRIPT" AAZ0_vs_PBrw.fa 5 &
echo "Started AAZ0 vs PBrw (job in background)"

bash "$RE_SCRIPT" AAZW_vs_PBrw.fa 5 &
echo "Started AAZW vs PBrw (job in background)"

echo "Done. You can exit now."
exit 0
