#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -N Merfin_job
#PBS -j oe

# Clean scratch on exit
trap 'clean_scratch' TERM EXIT

# VARIABLES
READ_COVERAGE=45       
KMER_SIZE=19           

# Paths 
WORKDIR="/storage/brno2/home/jendrb00/DIPLOMOVA_PRACE_VYSLEDKY/my_analysis/Abraxas_sylvata_merfin/AA_Z_W"
MERFIN_BIN="/storage/brno2/home/jendrb00/merfin/build/bin/merfin"
MERYL_PATH="/storage/plzen1/home/jendrb00/meryl-1.4.1/bin"
READ_KMER_DB="without_errors.meryl"
ASM_KMER_DB="AA_Z_W.meryl"
READS_FASTQ="ERR12370385.fastq.gz"
ASSEMBLY_FASTA="AA_Z_Wasmbl.fasta"
THREADS=8
TIMESTAMP=$(date +"%Y%m%d_%H%M")
OUTDIR="$WORKDIR/merfin_output_$TIMESTAMP"

# Prepare scratch
mkdir -p "$SCRATCHDIR/merfin_run"
cp -r "$WORKDIR/$READ_KMER_DB" "$SCRATCHDIR/merfin_run/"
cp -r "$WORKDIR/$ASM_KMER_DB" "$SCRATCHDIR/merfin_run/"
cp "$WORKDIR/$READS_FASTQ" "$SCRATCHDIR/merfin_run/"
cp "$WORKDIR/$ASSEMBLY_FASTA" "$SCRATCHDIR/merfin_run/"
cd "$SCRATCHDIR/merfin_run"

# Run Merfin
"$MERFIN_BIN" -dump \
  -sequence "$ASSEMBLY_FASTA" \
  -readmers "$READ_KMER_DB" \
  -seqmers "$ASM_KMER_DB" \
  -peak "$READ_COVERAGE" \
  -output merfin_stats.tsv

# Copy result back
cp merfin_stats.tsv "$WORKDIR/"

# Annotate copy number
awk -v cov=$READ_COVERAGE -F'\t' '
BEGIN {OFS="\t"}
NR==1 {print $0, "expected_cn", "deviation", "cn_class"}
NR>1 {
  exp_cn = ($2 / cov);
  dev = $3 - exp_cn;
  cn_class = (exp_cn < 1.5) ? "1x" : ((exp_cn < 2.5) ? "2x" : ((exp_cn < 5) ? "repeat" : "high-copy"));
  print $0, exp_cn, dev, cn_class
}' "$WORKDIR/merfin_stats.tsv" > "$WORKDIR/merfin_stats_annotated.tsv"

# Extract kmers
awk -F'\t' 'NR>1 && $3 == 0 && $2 >= 3 {print $1}' "$WORKDIR/merfin_stats_annotated.tsv" > "$WORKDIR/missing_kmers.txt"
awk -F'\t' -v cov=$READ_COVERAGE 'NR>1 && $3 != 0 && ($3 - ($2 / cov)) < -1 {print $1}' "$WORKDIR/merfin_stats_annotated.tsv" > "$WORKDIR/collapsed_kmers.txt"
cat "$WORKDIR/missing_kmers.txt" "$WORKDIR/collapsed_kmers.txt" | sort -T "$WORKDIR" | uniq > "$WORKDIR/problem_kmers.txt"

# Activate meryl
export PATH="$MERYL_PATH:$PATH"

# Extract FASTA headers
grep "^>" "$WORKDIR/$ASSEMBLY_FASTA" | cut -d' ' -f1 | sed 's/^>//' > "$WORKDIR/all_headers.txt"

# Match headers
grep -Ff "$WORKDIR/problem_kmers.txt" "$WORKDIR/all_headers.txt" > "$WORKDIR/problem_headers.txt"

# Convert kmers list to FASTA
seqtk subseq "$WORKDIR/$ASSEMBLY_FASTA" "$WORKDIR/problem_headers.txt" > "$WORKDIR/problem_scaffolds.fa"

# Create meryl db from scaffolds
meryl count k=$KMER_SIZE "$WORKDIR/problem_scaffolds.fa" output "$WORKDIR/problem_kmers2.meryl"

# Recover reads using meryl-lookup
meryl-lookup -include \
  -sequence "$WORKDIR/$READS_FASTQ" \
  -mers "$WORKDIR/problem_kmers2.meryl" \
  -output "$WORKDIR/problem_reads2.fa"

# Load modules
module purge
module load python/3.9.12-gcc-10.2.1
module load samtools/1.13-gcc-10.2.1
module load minimap2/2.22-gcc-10.2.1
module load bedtools

# Map recovered reads and process
minimap2 -t "$THREADS" -a -x map-ont "$WORKDIR/$ASSEMBLY_FASTA" "$WORKDIR/problem_reads2.fa" | \
  samtools sort -@ "$THREADS" -o "$WORKDIR/problem_reads.sorted.bam"
samtools index "$WORKDIR/problem_reads.sorted.bam"
bedtools bamtobed -i "$WORKDIR/problem_reads.sorted.bam" > "$WORKDIR/problem_reads.bed"

# Copy back
mkdir -p "$OUTDIR"
cp -r "$SCRATCHDIR/merfin_run/"* "$OUTDIR/" || { export CLEAN_SCRATCH=false; echo "Chyba při kopírování"; exit 9; }

# Hotovo
echo "HOTOVO. Výsledky byly uloženy do: $OUTDIR"


