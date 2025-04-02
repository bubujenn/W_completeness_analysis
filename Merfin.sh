#Download Merfin
#Cílová složka 
INSTALL_DIR="/storage/brno12-cerit/home/jendrb00/merfin"
#Přesun do složky
cd "$(dirname "$INSTALL_DIR")"
# Klonuj repozitář
if [ -d "$INSTALL_DIR" ]; then
  echo "❗ Starý Merfin se maže..."
  rm -rf "$INSTALL_DIR"
fi
git clone --recursive https://github.com/arangrhie/merfin.git
# Přesun do složky se zdroji
cd "$INSTALL_DIR/src"
# Pridat modul
module add gcc
#Kompilace
cd ~/merfin/src
make
_________________________________________________________________
#!/bin/bash
#PBS -l select=1:ncpus=8:mem=64gb:scratch_local=200gb
#PBS -l walltime=12:00:00
#PBS -N Merfin_job
#PBS -j oe

# Clean scratch on exit
trap 'clean_scratch' TERM EXIT

# ─── PROMĚNNÉ
READ_COVERAGE=45
KMER_SIZE=19
THREADS=8

# Cesty
WORKDIR="/storage/brno2/home/jendrb00/DIPLOMOVA_PRACE_VYSLEDKY/my_analysis/Abraxas_sylvata_merfin/AA_Z_W"
MERFIN_BIN="/storage/brno2/home/jendrb00/merfin/build/bin/merfin"
MERYL_PATH="/storage/plzen1/home/jendrb00/meryl-1.4.1/bin"
READ_KMER_DB="without_errors.meryl"
ASM_KMER_DB="AA_Z_W.meryl"
READS_FASTQ="ERR12370385.fastq.gz"
ASSEMBLY_FASTA="AA_Z_Wasmbl.fasta"
TIMESTAMP=$(date +"%Y%m%d_%H%M")
OUTDIR="$WORKDIR/merfin_output_sh_$TIMESTAMP"

# ─── PŘÍPRAVA SCRATCHE 
mkdir -p "$SCRATCHDIR/merfin_run"
cp -r "$WORKDIR/$READ_KMER_DB" "$SCRATCHDIR/merfin_run/"
cp -r "$WORKDIR/$ASM_KMER_DB" "$SCRATCHDIR/merfin_run/"
cp "$WORKDIR/$READS_FASTQ" "$SCRATCHDIR/merfin_run/"
cp "$WORKDIR/$ASSEMBLY_FASTA" "$SCRATCHDIR/merfin_run/"
cd "$SCRATCHDIR/merfin_run"

# ─── SPUŠTĚNÍ MERFINU
"$MERFIN_BIN" -dump \
  -sequence "$ASSEMBLY_FASTA" \
  -readmers "$READ_KMER_DB" \
  -seqmers "$ASM_KMER_DB" \
  -peak "$READ_COVERAGE" \
  -output merfin_stats.tsv

# ─── ANOTACE COPY NUMBER DEVIACÍ 
awk -v cov=$READ_COVERAGE -F'\t' '
BEGIN {OFS="\t"}
NR==1 {print $0, "expected_cn", "deviation", "cn_class"}
NR>1 {
  exp_cn = ($2 / cov);
  dev = $3 - exp_cn;
  cn_class = (exp_cn < 1.5) ? "1x" : ((exp_cn < 2.5) ? "2x" : ((exp_cn < 5) ? "repeat" : "high-copy"));
  print $0, exp_cn, dev, cn_class
}' merfin_stats.tsv > merfin_stats_annotated.tsv

# ─── EXTRAKCE K-MERŮ 
awk -F'\t' 'NR>1 && $3 == 0 && $2 >= 3 {print $1}' merfin_stats_annotated.tsv > missing_kmers.txt
awk -F'\t' -v cov=$READ_COVERAGE 'NR>1 && $3 != 0 && ($3 - ($2 / cov)) < -1 {print $1}' merfin_stats_annotated.tsv > collapsed_kmers.txt
sort --temporary-directory="$SCRATCHDIR" collapsed_kmers.txt missing_kmers.txt | uniq > problem_kmers.txt

# ─── VYTVOŘENÍ MERYL DB Z TXT 
"$MERYL_PATH/meryl" count k=$KMER_SIZE output problem_kmers.meryl problem_kmers.txt

# ─── NALEZENÍ READŮ S PROBLÉMOVÝMI K-MERY 
"$MERYL_PATH/meryl-lookup" -include \
  -sequence "$READS_FASTQ" \
  -mers problem_kmers.meryl \
  -output problem_reads.fa

# ─── NAČTENÍ MODULŮ
module purge
module load python/3.9.12-gcc-10.2.1
module load samtools/1.13-gcc-10.2.1
module load minimap2/2.22-gcc-10.2.1

# ─── MAPOVÁNÍ NA ASSEMBLY
minimap2 -t "$THREADS" -a -x map-ont "$ASSEMBLY_FASTA" problem_reads.fa | \
  samtools sort -@ "$THREADS" -o problem_reads.sorted.bam
samtools index problem_reads.sorted.bam

# ─── ZÁVĚR A KOPÍROVÁNÍ 
mkdir -p "$OUTDIR"
cp -r "$SCRATCHDIR/merfin_run/"* "$OUTDIR/" || { export CLEAN_SCRATCH=false; echo "Chyba pri kopirovani"; exit 9; }

echo "HOTOVO: Výsledky byly uloženy do složky: $OUTDIR"



