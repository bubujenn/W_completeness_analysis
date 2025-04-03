#!/bin/bash
set -euo pipefail

# Configurable parameters
MERFIN_STATS="merfin_stats.tsv"       # Output from merfin -dump with -peak
KMER_SIZE=19
BCALM_BIN="/path/to/bcalm"            # <-- Replace with full path to bcalm binary
OUTDIR="bcalm_output"

# Check input
if [[ ! -f "$MERFIN_STATS" ]]; then
  echo ":x: Error: file $MERFIN_STATS does not exist."
  exit 1
fi

# Setup output folders
mkdir -p "$OUTDIR"
PROBLEM_KMERS_TXT="$OUTDIR/problem_kmers.txt"
KMER_FASTA="$OUTDIR/problem_kmers.fasta"

# Filter k-mers: only missing or collapsed
echo ":small_blue_diamond: Filtering missing and collapsed k-mers from Merfin stats..."
awk -F'\t' '$0 ~ /missing|collapsed/ {print $1}' "$MERFIN_STATS" > "$PROBLEM_KMERS_TXT"
if [[ ! -s "$PROBLEM_KMERS_TXT" ]]; then
  echo ":warning: No missing or collapsed k-mers found. Exiting."
  exit 0
fi

# Convert to FASTA
echo ":small_blue_diamond: Converting k-mer list to FASTA format..."
awk '{print ">kmer" NR "\n" $1}' "$PROBLEM_KMERS_TXT" > "$KMER_FASTA"

# Run BCALM2 to merge overlapping k-mers
echo ":small_blue_diamond: Running BCALM2..."
"$BCALM_BIN" \
  -in "$KMER_FASTA" \
  -kmer-size "$KMER_SIZE" \
  -abundance-min 1 \
  -out "$OUTDIR/unitigs"
UNITIG_FASTA="$OUTDIR/unitigs.unitigs.fa"
if [[ ! -f "$UNITIG_FASTA" ]]; then
  echo ":x: Error: BCALM2 failed to produce $UNITIG_FASTA"
  exit 1
fi

# Sum total missing bases
echo ":small_blue_diamond: Summarizing unitig lengths..."
TOTAL_BP=$(awk '/^>/ {next} {sum += length($0)} END {print sum}' "$UNITIG_FASTA")
NUM_UNITIGS=$(grep -c "^>" "$UNITIG_FASTA")
echo ":white_tick: Estimated missing sequence:"
echo "- Unitigs: $NUM_UNITIGS"
echo "- Total bp: $TOTAL_BP"
