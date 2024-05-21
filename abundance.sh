### abundance analysis
#mask reference genome
genome_file="allinput.fa"
masked_genomefn="allinput_masked.fasta"
bash /lomi_home/wenxiu/pipeline/bowtie_align_material/masked.sh $genome_file

#bwa
module load miniconda3
source activate bwa
bwa index $masked_genomefn
for i in *_1.fastq; do
bwa mem -5SP -t 32 ${1/.fa}_masked.fasta $i ${i/_1.fastq/_2.fastq} | samblaster  > ${i/_1.fastq}.sam
done
conda deactivate
source activate samtools
for i in *sam; do
samtools view -@ 32 -bS -h -b $i > ${i/.sam/}.bam
samtools sort -@ 32 ${i/.sam/}.bam -o ${i/.sam/}.sorted.bam
done
conda deactivate

#filtering
source activate bamm
for i in *sorted.bam; do
bamm filter -b $i --percentage_id 0.99 --percentage_aln 0.9
done
conda deactivate

#bbmap calculate reads
source activate bbmap
for i in *sorted.bam; do
pileup.sh in=$i out=${i/.bam/}.cov rpkm=${i/.bam/}.rpkm overwrite=true
done
