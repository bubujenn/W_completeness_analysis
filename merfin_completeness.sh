#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -N MerfinCompleteness
#PBS -j oe

trap 'clean_scratch' TERM EXIT

READ_COVERAGE=45
KMER_SIZE=19
THREADS=8

WORKDIR="/storage/brno2/home/jendrb00/DIPLOMOVA_PRACE_VYSLEDKY/my_analysis/Abraxas_sylvata_merfin/AA_Z_W/try2"
MERFIN_BIN="/storage/brno2/home/jendrb00/merfin/build/bin/merfin"
MERYL_PATH="/storage/plzen1/home/jendrb00/meryl-1.4.1/bin"
READ_KMER_DB="without_errors.meryl"
ASM_KMER_DB="AA_Z_W.meryl"
READS_FASTQ="ERR12370385.fastq.gz"
ASSEMBLY_FASTA="AA_Z_Wasmbl.fasta"
TIMESTAMP=$(date +"%Y%m%d_%H%M")
OUTDIR="$WORKDIR/merfin_output_completeness_$TIMESTAMP"

mkdir -p "$SCRATCHDIR/merfin_run"
cp -r "$WORKDIR/$READ_KMER_DB" "$SCRATCHDIR/merfin_run/"
cp -r "$WORKDIR/$ASM_KMER_DB" "$SCRATCHDIR/merfin_run/"
cp "$WORKDIR/$READS_FASTQ" "$SCRATCHDIR/merfin_run/"
cp "$WORKDIR/$ASSEMBLY_FASTA" "$SCRATCHDIR/merfin_run/"
cd "$SCRATCHDIR/merfin_run"

"$MERFIN_BIN" -completeness \
  -sequence "$ASSEMBLY_FASTA" \
  -readmers "$READ_KMER_DB" \
  -seqmers "$ASM_KMER_DB" \
  -peak "$READ_COVERAGE" \
  -output merfin_completeness.tsv

  (/storage/brno2/home/jendrb00/merfin/build/bin/merfin \
  -completeness \
  -sequence AA_Z_Wasmbl.fasta \
  -readmers without_errors.meryl \
  -seqmers AA_Z_W.meryl \
  -peak 45 \
  > merfin_completeness.txt 2>&1)

mkdir -p "$OUTDIR"
cp -r "$SCRATCHDIR/merfin_run/"* "$OUTDIR/" || { export CLEAN_SCRATCH=false; echo "Chyba při kopírování!"; exit 9; }

echo "HOTOVO: Výsledky byly uloženy do složky: $OUTDIR"
