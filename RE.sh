#!/bin/bash

# script to comparative clustering by RepeatExplorer
# expected name of fasta files /PATH/PrefixSampleNumber*.fasta 
# usage: bash /PATH/script.sh /PATH/inputs/file.fasta prefix_length

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





Priprava dat - 
#!/bin/bash
# Data preparation for RE 

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

Spusteni -
# Dataset 1: AAZ0 vs PBrw (W chromosome)
bash RE_clustering.sh AAZ0_vs_PBrw_v3.fa 5

# Dataset 2: AAZW vs PBrw (W chromosome completeness)  
bash RE_clustering.sh AAZW_vs_PBrw_v3.fa 5
