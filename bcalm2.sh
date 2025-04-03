#!/bin/bash
#PBS -N bcalm_missing_seq
#PBS -l select=1:ncpus=8:mem=60gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -j oe

# Clean scratch space on exit
trap 'clean_scratch' TERM EXIT

# Set base variables
Path=${PWD}
Base="bcalm_missing_seq"
TimeStamp=$(date +"%Y%m%d.%H%M")
Outdir=${Path}/${Base}_${TimeStamp}
SCRATCHDIR=${SCRATCHDIR:-"/scratch"}

# Input files and tools
MERFIN_STATS="${1:-/storage/brno2/home/jendrb00/DIPLOMOVA_PRACE_VYSLEDKY/my_analysis/Abraxas_sylvata_merfin/AA_Z_W/merfin_output4_*/merfin_stats.tsv}"
KMER_SIZE=19
BCALM_BIN="/storage/brno2/home/jendrb00/bcalm/build/bin/bcalm"

# Prepare scratch and output
mkdir -p "$SCRATCHDIR/bcalm_run"
mkdir -p "$Outdir"

# Resolve Merfin stats file (in case of wildcard)
MERFIN_STATS_FILE=$(ls $MERFIN_STATS | tail -n1)
if [[ ! -f "$MERFIN_STATS_FILE" ]]; then
  echo ":x: Error: Merfin stats file not found!"
  exit 1
fi

# Copy Merfin stats to scratch
cp "$MERFIN_STATS_FILE" "$SCRATCHDIR/bcalm_run/"
cd "$SCRATCHDIR/bcalm_run" || exit 2

# Output files
PROBLEM_KMERS_TXT="problem_kmers.txt"
KMER_FASTA="problem_kmers.fasta"
UNITIG_PREFIX="unitigs"

# Filter k-mers: missing + collapsed
echo ":small_blue_diamond: Filtering missing and collapsed k-mers..."
awk -F'\t' '$0 ~ /missing|collapsed/ {print $1}' "$(basename "$MERFIN_STATS_FILE")" > "$PROBLEM_KMERS_TXT"
if [[ ! -s "$PROBLEM_KMERS_TXT" ]]; then
  echo ":warning: No problematic k-mers found. Exiting."
  exit 0
fi

# Convert k-mers to FASTA
echo ":small_blue_diamond: Converting k-mers to FASTA format..."
awk '{print ">kmer" NR "\n" $1}' "$PROBLEM_KMERS_TXT" > "$KMER_FASTA"

# Run BCALM2 to build unitigs
echo ":small_blue_diamond: Running BCALM2..."
"$BCALM_BIN" \
  -in "$KMER_FASTA" \
  -kmer-size "$KMER_SIZE" \
  -abundance-min 1 \
  -out "$UNITIG_PREFIX"

UNITIG_FASTA="${UNITIG_PREFIX}.unitigs.fa"
if [[ ! -f "$UNITIG_FASTA" ]]; then
  echo ":x: Error: BCALM2 failed to produce output!"
  exit 1
fi

# Summarize output
echo ":small_blue_diamond: Summarizing unitigs..."
TOTAL_BP=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "$UNITIG_FASTA")
NUM_UNITIGS=$(grep -c "^>" "$UNITIG_FASTA")

# Copy results back
cp -r "$SCRATCHDIR/bcalm_run"/* "$Outdir" || { export CLEAN_SCRATCH=false; exit 3; }

# Final output info
echo ":white_tick: Done!"
echo "- Output directory: $Outdir"
echo "- Unitigs: $NUM_UNITIGS"
echo "- Total bp in unitigs: $TOTAL_BP"

