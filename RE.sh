#!/bin/bash
# RepeatExplorer Pipeline for W Chromosome Completeness Analysis

# STEP 1: Zkrácení názvů sekvencí (awk)
# STEP 2: Tagování (sed)  
# STEP 3: Concatenace (cat)
# STEP 4: Spuštění RepeatExplorer (bash RE_clustering.sh)
============================================================================
#  RE_clustering.sh
============================================================================

Fasta=$(basename ${1})
TimeStamp=$(date +"%Y%m%d.%H%M")
Outdir=$PWD/RE_${TimeStamp} 
Outfile=${Outdir}/RE_clustering.pbs

mkdir ${Outdir}
cd ${Outdir}

(
cat << EndOfJobScript
#!/bin/bash
#PBS -N RE_clustering
#PBS -l select=1:ncpus=20:mem=1024gb:scratch_local=1400gb:singularity=True
#PBS -l walltime=144:00:00
#PBS -q large_mem@pbs-m1.metacentrum.cz
#PBS -M blanka.jendriskova0000@seznam.cz
#PBS -m bea
#PBS -o ${Outdir}/stdout.txt
#PBS -e ${Outdir}/stderr.txt

trap 'clean_scratch' TERM EXIT
export TMPDIR=\${SCRATCHDIR}
cp ${1} \${SCRATCHDIR}
cd \${SCRATCHDIR}

test -n "\$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

singularity pull --arch amd64 library://repeatexplorer/default/repex_tarean:0.3.12-7a7dc9e 

export SINGULARITY_CACHEDIR=\${SCRATCHDIR}
export SINGULARITY_LOCALCACHEDIR=\${SCRATCHDIR}
export SINGULARITY_TMPDIR=\${SCRATCHDIR}

singularity exec --bind \${SCRATCHDIR}:/data/ --bind \$SCRATCHDIR:/tmp/ \${SCRATCHDIR}/repex_tarean_0.3.12-7a7dc9e.sif seqclust --paired --automatic_filtering --mincl 0.001 --cpu 18 --prefix_length ${2} \
--output_dir=/data/output/ --max_memory 900000000 --cleanup --assembly_min 3  --taxon METAZOA3.0 /data/${Fasta}

rm -r \${SCRATCHDIR}/repex_tarean_0.3.12-7a7dc9e.sif
cp -r \${SCRATCHDIR}/output ${Outdir} || export CLEAN_SCRATCH=false
EndOfJobScript
) > ${Outfile}

chmod +x ${Outfile}
qsub ${Outfile}
cd ..

============================================================================
# PART 2: Data Preparation
============================================================================

cd /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs

# Zkrácení názvů sekvencí
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZ_interleaved.fa > AAZ_short.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZW_interleaved.fa > AAZW_short.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' PBrw_interleaved.fa > PBrw_short.fa

# Tagování
sed 's/>/&AAZ0_/' AAZ_short.fa > AAZ0_tagged.fa
sed 's/>/&AAZW_/' AAZW_short.fa > AAZW_tagged.fa
sed 's/>/&PBrw_/' PBrw_short.fa > PBrw_tagged.fa

# Concatenace
cat AAZ0_tagged.fa PBrw_tagged.fa > AAZ0_vs_PBrw.fa
cat AAZW_tagged.fa PBrw_tagged.fa > AAZW_vs_PBrw.fa

============================================================================
# PART 3: Spuštění RepeatExplorer
============================================================================

# Dataset 1: AAZ0 vs PBrw
bash /storage/plzen1/home/p817n421/Shared/blues/scripts/RE_clustering.sh /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/AAZ0_vs_PBrw.fa 5

# Dataset 2: AAZW vs PBrw
bash /storage/plzen1/home/p817n421/Shared/blues/scripts/RE_clustering.sh /storage/brno2/home/jendrb00/20_10_consultation/RE_final_inputs/AAZW_vs_PBrw.fa 5



______________________________________________________________________________________________________________________________

 =============================================================================
# PART 1: SIMULATION (randomreads.sh)
 =============================================================================
# AAZ (genom samce, bez W) ---
# ref: AA_Z_chr.fasta
randomreads.sh ref=AA_Z_chr.fasta out=AAZ_100bp.fastq.gz \
    len=100 reads=20000000 seed=42 paired=f adderrors=f

# AAZW (genom samice, s W) ---
# ref: AA_Z_Wasmbl.fasta
randomreads.sh ref=AA_Z_Wasmbl.fasta out=AAZW_100bp.fastq.gz \
    len=100 reads=20000000 seed=42 paired=f adderrors=f

# PBrw (PacBio 2× vzorek) ---
# Krok 1: 2× vzorek z PacBio (14M readů = 1.4 Gb)
seqtk sample -s42 pb_100_se.fastq.gz 14000000 > PB_2xgenome.fastq
seqtk seq -a PB_2xgenome.fastq > PB_2xgenome.fasta

# Krok 2: Simulace
randomreads.sh ref=PB_2xgenome.fasta out1=PBrw_R1.fastq.gz out2=PBrw_R2.fastq.gz \
    len=100 reads=20000000 seed=42 paired=t adderrors=f

 =============================================================================
# PART 2: SUBSAMPLE (normalizace na 0.4× coverage)
 =============================================================================
# Počty: AAZ 438646 párů, AAZW/PBrw 449857 párů

# AAZ (438646 párů)
seqtk sample -s42 ../AAZ_R1.fastq.gz 438646 > AAZ_R1.sampled.fq
seqtk sample -s42 ../AAZ_R2.fastq.gz 438646 > AAZ_R2.sampled.fq

# AAZW (449857 párů)
seqtk sample -s42 ../AAZW_R1.fastq.gz 449857 > AAZW_R1.sampled.fq
seqtk sample -s42 ../AAZW_R2.fastq.gz 449857 > AAZW_R2.sampled.fq

# PBrw (449857 párů)
seqtk sample -s42 ../PBrw_2xgenome/PBrw_R1.fastq.gz 449857 > PBrw_R1.sampled.fq
seqtk sample -s42 ../PBrw_2xgenome/PBrw_R2.fastq.gz 449857 > PBrw_R2.sampled.fq

=============================================================================
# PART 3: INTERLEAVE + FASTA
=============================================================================

seqtk mergepe AAZ_R1.sampled.fq AAZ_R2.sampled.fq > AAZ_interleaved.fastq
seqtk seq -a AAZ_interleaved.fastq > AAZ_interleaved.fa

seqtk mergepe AAZW_R1.sampled.fq AAZW_R2.sampled.fq > AAZW_interleaved.fastq
seqtk seq -a AAZW_interleaved.fastq > AAZW_interleaved.fa

seqtk mergepe PBrw_R1.sampled.fq PBrw_R2.sampled.fq > PBrw_interleaved.fastq
seqtk seq -a PBrw_interleaved.fastq > PBrw_interleaved.fa

=============================================================================
# PART 4: ZKRÁCENÍ NÁZVŮ
=============================================================================

awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZ_interleaved.fa > AAZ_short.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' AAZW_interleaved.fa > AAZW_short.fa
awk '/^>/{split($0,a," "); print a[1]; next} {print}' PBrw_interleaved.fa > PBrw_short.fa

=============================================================================
# PART 5: TAGOVÁNÍ (prefix_length = 5)
=============================================================================

sed 's/>/>AAZ0_/' AAZ_short.fa > AAZ0_tagged.fa
sed 's/>/>AAZW_/' AAZW_short.fa > AAZW_tagged.fa
sed 's/>/>PBrw_/' PBrw_short.fa > PBrw_tagged.fa

=============================================================================
# PART 6: CONCATENACE
=============================================================================

cat AAZ0_tagged.fa PBrw_tagged.fa > AAZ0_vs_PBrw.fa   # Složení repetic W
cat AAZW_tagged.fa PBrw_tagged.fa > AAZW_vs_PBrw.fa   # Kompletnost W

# Ověření
grep -c '^>' AAZ0_vs_PBrw.fa   # 1777006
grep -c '^>' AAZW_vs_PBrw.fa   # 1799428

