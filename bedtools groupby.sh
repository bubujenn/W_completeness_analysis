#check - meryl print is working right
qsub -I -l select=1:ncpus=16:mem=128gb:scratch_local=200gb -l walltime=24:00:00 
export PATH=/storage/plzen1/home/jendrb00/meryl-1.4.1/bin:$PATH
meryl print expected_W_size.meryl > expected_W_size.tsv
meryl print W_missing_in_assembly.meryl > W_missing_in_assembly.tsv
module load bedtools 
bedtools groupby -i expected_W_size.tsv -c 2 -o sum > expected_W_size_bedtools_sum.txt
bedtools groupby -i W_missing_in_assembly.tsv -c 2 -o sum > W_missing_in_assembly_bedtools_sum.txt
