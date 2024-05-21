###nucleotide diversity analysis 
module load miniconda3
source activate bowtie2
bowtie2-build allinput.fa allinput.fa
for i in *_1.fastq; do
bowtie2 -p 32 -x allgenomes.fa -1 ${i} -2 ${i/_1.fastq}_2.fastq > allinput.fa-vs-${i/_1.fastq}.sam
done
conda deactiavte
source activate samtools
for i in *sam; do
samtools view -@ 32 -bS $i | samtools sort -@ 32 > ${i/.sam/}.sorted.bam
done
conda deactivate
module unload miniconda3

#instrain
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate instrain
for i in *bam; do
inStrain profile $i allgenomes.fa -o ${i/.bam}.gene.IS -p 32 -s genomes.stb -g genes.fna
done
inStrain compare -i *gene.IS -s genomes.stb -p 32 -o allinput.gene.IS.COMPARE
