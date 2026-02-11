#!/bin/bash
# =============================================================================
# RepeatExplorer Pipeline 
# =============================================================================
# 1. Simulace z referenci (genomy AAZ, AAZW a PacBio raw data)
# 2. Subsample na normalizovany coverage 0.4x
# 3. Interleave a rename readu
# =============================================================================

# Pracovni cesta a moduly
WORK_DIR="/storage/brno2/home/jendrb00/20_10_consultation"
cd "$WORK_DIR" || exit 1

module load bbmap
module load seqtk

# =============================================================================
# CAST 1: SIMULACE (randomreads.sh)
# Vsechny simulace s paired=t, adderrors=f (bez chyb)
# =============================================================================

# --- AAZ (genom samce, bez W) ---
# ref: AA_Z_chr.fasta
randomreads.sh \
    ref=AA_Z_chr.fasta \
    out1=AAZ_R1.fastq.gz \
    out2=AAZ_R2.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f

# --- AAZW (genom samice, s W) ---
# ref: AA_Z_Wasmbl.fasta
randomreads.sh \
    ref=AA_Z_Wasmbl.fasta \
    out1=AAZW_R1.fastq.gz \
    out2=AAZW_R2.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f

# --- PBrw (PacBio data samice) ---
# Krok 1: 2x subsample z RAW PacBio dat (1.4 Gb)
# Pouzivam RAW data, ne predpripravene 100bp fragmenty.

RAW_PACBIO_FILE="/storage/brno2/home/jendrb00/PacBiok100/ERR12370385.fastq.gz"

# Nahodny vyber z raw dat (cca 200k readu odpovida 1.4 Gb pri prumerne delce 7kb)
seqtk sample -s42 "$RAW_PACBIO_FILE" 200000 > PB_2xgenome_raw.fastq
seqtk seq -a PB_2xgenome_raw.fastq > PB_2xgenome_raw.fasta

# Krok 2: Simulace z dlouhych readu
# randomreads.sh naseka dlouhe ready na kratke parove ready
randomreads.sh \
    ref=PB_2xgenome_raw.fasta \
    out1=PBrw_R1_new.fastq.gz \
    out2=PBrw_R2_new.fastq.gz \
    len=100 \
    reads=20000000 \
    seed=42 \
    paired=t \
    adderrors=f

# =============================================================================
# CAST 2: SUBSAMPLE (na coverage 0.4x)
# =============================================================================
# AAZ: 438646 paru (mensi genom)
# AAZW a PBrw: 449857 paru (vetsi genom s W)

INPUT_DIR="$WORK_DIR/RE_final_inputs"
mkdir -p "$INPUT_DIR"
cd "$INPUT_DIR" || exit 1

# AAZ
seqtk sample -s42 ../AAZ_R1.fastq.gz 438646 > AAZ_R1.sampled.fq
seqtk sample -s42 ../AAZ_R2.fastq.gz 438646 > AAZ_R2.sampled.fq

# AAZW
seqtk sample -s42 ../AAZW_R1.fastq.gz 449857 > AAZW_R1.sampled.fq
seqtk sample -s42 ../AAZW_R2.fastq.gz 449857 > AAZW_R2.sampled.fq

# PBrw
seqtk sample -s42 ../PBrw_R1_new.fastq.gz 449857 > PBrw_R1.sampled.fq
seqtk sample -s42 ../PBrw_R2_new.fastq.gz 449857 > PBrw_R2.sampled.fq

# =============================================================================
# CAST 3: INTERLEAVE + FASTA
# =============================================================================

# AAZ
seqtk mergepe AAZ_R1.sampled.fq AAZ_R2.sampled.fq > AAZ_interleaved.fastq
seqtk seq -a AAZ_interleaved.fastq > AAZ_interleaved_new.fa

# AAZW
seqtk mergepe AAZW_R1.sampled.fq AAZW_R2.sampled.fq > AAZW_interleaved.fastq
seqtk seq -a AAZW_interleaved.fastq > AAZW_interleaved_new.fa

# PBrw
seqtk mergepe PBrw_R1.sampled.fq PBrw_R2.sampled.fq > PBrw_interleaved.fastq
seqtk seq -a PBrw_interleaved.fastq > PBrw_interleaved_new.fa

# =============================================================================
# CAST 4: ZKRACENI NAZVU A TAGOVANI
# =============================================================================
# Zkraceni nazvu na prvni slovo pred mezerou a pridani prefixu

awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZ_interleaved_new.fa | sed 's/>/>AAZ0_/' > AAZ0_tagged.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZW_interleaved_new.fa | sed 's/>/>AAZW_/' > AAZW_tagged.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' PBrw_interleaved_new.fa | sed 's/>/>PBrw_/' > PBrw_tagged.fa

# =============================================================================
# CAST 5: SPOJENI (CONCATENACE)
# =============================================================================

cat AAZ0_tagged.fa PBrw_tagged.fa > AAZ0_vs_PBrw.fa   # Slozeni repetic W
cat AAZW_tagged.fa PBrw_tagged.fa > AAZW_vs_PBrw.fa   # Kompletnost W

# Kontrola poctu sekvenci
# Ocekavano: AAZ0_vs_PBrw ~1.77M, AAZW_vs_PBrw ~1.80M
echo "Pocet sekvenci:"
grep -c '^>' AAZ0_vs_PBrw.fa
grep -c '^>' AAZW_vs_PBrw.fa

# =============================================================================
# CAST 6: SPUSTENI RepeatExplorer
# =============================================================================

RE_SCRIPT="/storage/plzen1/home/p817n421/Shared/blues/scripts/RE_clustering.sh"

echo "Spoustim AAZ0 vs PBrw..."
bash "$RE_SCRIPT" AAZ0_vs_PBrw.fa 5 &

echo "Spoustim AAZW vs PBrw..."
bash "$RE_SCRIPT" AAZW_vs_PBrw.fa 5 &

wait
echo "Hotovo."
