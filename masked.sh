#Cover highly conserved and repetitive regions of the genome:
#Find rRNA region
module load miniconda3
source activate metaxa-2.2
metaxa2 --plus --mode m --cpu 32 --multi_thread T --table T -g ssu --not_found T  -i $1 -o ${1/.fa/}.metaxa2_ssu
metaxa2 --plus --mode m --cpu 32 --multi_thread T --table T -g lsu --not_found T  -i $1 -o ${1/.fa/}.metaxa2_lsu
 
cut -f 1,9,10 ${1/.fa/}.metaxa2_ssu.extraction.results >> masked_metaxa.bed
cut -f 1,9,10 ${1/.fa/}.metaxa2_lsu.metaxa2_lsu.extraction.results >> masked_metaxa.bed

conda deactivate
source activate barrnap
barrnap --kingdom bac --threads 32 --reject 0.3 $1 > out_bac.gff
barrnap --kingdom arc --threads 32 --reject 0.3 $1 > out_arc.gff
barrnap --kingdom euk --threads 32 --reject 0.3 $1 > out_euk.gff
 
cut -f 1,4,5 out_bac.gff >> masked_barrnap.bed
cut -f 1,4,5 out_arc.gff >> masked_barrnap.bed
cut -f 1,4,5 out_euk.gff >> masked_barrnap.bed

cat masked_barrnap.bed masked_metaxa.bed > masked_rrna.bed
conda deactivate

#Find tRNA region
source activate bioperl
perl /share/pipelines/tools/SEQMODIFY_FASTA_aliquot.pl $1 1000

source activate trnascan-se
export TMPDIR="/lomi_home/wenxiu/tmp"
export PATH=$PATH:/lomi_home/wenxiu/software/GNU/parallel/bin

ls ${1}_trunc* |  parallel -j 32 "tRNAscan-SE -A -b {}.ar.bedformat --thread 1 {} "
cat *.ar.bedformat >> ${1/.fa}.bedformat
find . -name "${1}_trunc*bedformat" | xargs rm

ls ${1}_trunc* |  parallel -j 32 "tRNAscan-SE -B -b {}.bac.bedformat --thread 1 {} "
cat *.bac.bedformat >> ${1/.fa}.bedformat
cut -f 1,2,3 ${1/.fa}.bedformat > masked_trnascan.bed
conda deactivate

find . -name "${1}_trunc*" | xargs rm

#Find duplicate areas
module load ncbi-cxx-toolkit
dustmasker -in $1 -out ${1/.fa/}.lowcom.out -outfmt acclist
cut -d '>' -f 2 ${1/.fa/}.lowcom.out > masked_dustmasker.bed 
module unload ncbi-cxx-toolkit 

#Cover the above area with NNNN
cat masked_rrna.bed masked_trnascan.bed masked_dustmasker.bed > ${1/.fa}_mask.bed
module load miniconda3
source activate bedtools
awk -F'\t' '$2>=0' ${1/.fa}_mask.bed > ${1/.fa}_mask_0.bed
bedtools maskfasta -fi $1 -bed ${1/.fa}_mask_0.bed -fo ${1/.fa}_masked.fasta
conda deactivate
