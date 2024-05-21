###01 remove the 3-end of the ligation site in each HiC read
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate hicstuff
for fq in *clean_1.fastq; do
fwd_read=$fq
rev_read=${fwd_read/_1.fastq/_2.fastq}
FILESTEM=${fwd_read/_1.fastq}
hicstuff cutsite -1 $fwd_read -2 $rev_read -p ${FILESTEM}_cutsite -e Sau3AI -s 75 -t 32
done

###02 hic reads mapped to reference genome
##mask reference genome
genome_file="allinput.fa"
masked_genomefn="allinput_masked.fasta"
bash /lomi_home/wenxiu/pipeline/bowtie_align_material/masked.sh $genome_file

##map HiC reads
export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/bamm/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/bwa/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/samtools/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/samblaster/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/sambamba/bin
bwa index $masked_genomefn
for fq in *cutsite_1.fastq; do
fwd_read=$fq
rev_read=${fwd_read/_1.fastq/_2.fastq}
FILESTEM=${fwd_read/_1.fastq}
bwa mem -5SP -t 32 ${masked_genomefn} $fwd_read $rev_read | samblaster | samtools view -S -h -b -q 20 -F 2316 > ${FILESTEM}_${masked_genomefn}.aligned.bam
sambamba view -t 32 -h -f bam -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" ${FILESTEM}_${masked_genomefn}.aligned.bam -o ${FILESTEM}_${masked_genomefn}.aligned.uniq.bam
samtools sort ${FILESTEM}_${masked_genomefn}.aligned.uniq.bam > ${FILESTEM}_${masked_genomefn}.aligned.sorted.bam
bamm filter -b ${FILESTEM}_${masked_genomefn}.aligned.sorted.bam --percentage_id 0.99 --percentage_aln 0.9
done

##get read pair mapped to different contigs
export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/bedtools/bin
for i in *filtered.bam; do bamToBed -i $i > ${i/.bam/}.bed; done
for i in *filtered.bed; do
cat $i | awk -v FS='\t' -v OFS='\t' '{print $1, $4}' | awk -F'/' '{print $1}' | awk -v FS='\t' -v OFS='\t' '{print $2, $1}' > ${i}.aligned.contig
python /lomi_home/wenxiu/script/hic/hiclink_generate_difflink.py ${i}.aligned.contig ${i}.difflink.txt
done
cat *difflink.txt > all_difflink.txt

###03 remove the read pairs with overlap alignment with greater than 80% identity and 40 bp overlap length
##do pairwise blastn in each contig pair region (for example: fwdread1_mapped_contig_region.fa, revread1_mapped_contig_region.fa)
blastn -subject fwdread1_mapped_contig_region.fa -query revread1_mapped_contig_region.fa -out read1_mapped_contig_region.out -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 1 #filter read1 if its mapped contig regions has greater than 80 identity and 40 bp alignemnt length

###04 remove duplicate-like reads in a contig-pair region using CDHIT 
##do cdhit in each contig pair has all mapped regions (for example: contig1-contig2_allmappedregions.fa)
/lomi_home/wenxiu/software/cdhit-4.8.1/cd-hit-est -i contig1-contig2_allmappedregions.fa -o contig1-contig2_allmappedregions.99.fa -T 1 -c 0.99 -aS 0.95 -M 0 -d 0 #Use the results of cdhit as an aid, manually remove duplicate-like reads with same/near-same mapping location in a contig-pair.

##05 get final links after getting drop.list from blastn and cdhit results of all contig-pair.
python /lomi_home/wenxiu/script/hic/alllink_filterReads.py all_difflink.txt drop.list all_difflink_filtered.txt
python /lomi_home/wenxiu/script/hic/hiclink_from_alllink.py all_difflink_filtered.txt difflink_filtered.xlsx
