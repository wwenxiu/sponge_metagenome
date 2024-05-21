###provirus link
##using phaster_api get allbin_provirus.fa
viralset="/lomi_home/wenxiu/HIC-meta/CCs_analyze/virus_inden/fnl_viralset/vp_mrna.min5k.cir.0.99-0.85.fasta"
prophageset="/lomi_home/wenxiu/HIC-meta/CCs_analyze/virus_inden/MAG_prophage_pip/phaster_api/allbin_provirus.fa"

#minimap2
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/samtools/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate minimap2
minimap2 -d catalogue.mmi ${prophageset} # make index
minimap2 -t 32 -N 5 -a catalogue.mmi ${viralset} | samtools view -F 3584 -b --threads 32 > provirus2viralset.bam

#bedtools
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/bedtools/bin
bamToBed -i provirus2viralset.bam > provirus2viralset.bed.txt

###cripsr link
mag_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/binning/promag_set"
viralset="/lomi_home/wenxiu/HIC-meta/CCs_analyze/virus_iden/fnl_viralset/vp_mrna.min5k.cir.0.99-0.85.fasta"

##CRISPRdect
CRISPRdect_file=${workdir}'/CRISPRdect'
mkdir -p ${CRISPRdect_file}
export PATH=$PATH:/lomi_home/wenxiu/software/hmmer/bin
export PATH=$PATH:/lomi_home/wenxiu/software/genemark/genemark_suite_linux_64/gmsuite
export PATH=$PATH:/lomi_home/wenxiu/software/CRISPRDetect_3.0
export PATH=$PATH:/lomi_home/wenxiu/software/cdhit-4.8.1
export PATH=$PATH:/lomi_home/wenxiu/software/EMBOSS/EMBOSS/bin
export PATH=$PATH:/lomi_home/wenxiu/software/blast/ncbi-blast-2.2.28+/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/ViennaRNA/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate crispr_env
for i in `ls ${mag_file}/*fa`; do
fn=$(basename $i)
/lomi_home/wenxiu/software/CRISPRDetect_3.0/CRISPRDetect3 -f $i -o ${CRISPRdect_file}/${fn/.fa/}_CRISPRDetect3 -T 0 -wgs 1 -tmp_dir /lomi_home/wenxiu/tmp -annotate_cas_genes 1 -check_plasmid 1 > ${CRISPRdect_file}/${fn/.fa/}.log
done
conda deactivate

##spacepharer
spacepharer_file=${workdir}'/spacepharer'
mkdir -p ${spacepharer_file}
source activate spacepharer
#create a database of the phage genomes targetSetDB and control sequences targetSetDB_rev
spacepharer createsetdb ${viralset} ${spacepharer_file}/targetSetDB/viralset /lomi_home/wenxiu/tmp
spacepharer createsetdb ${viralset} ${spacepharer_file}/targetSetDB/viralset_rev /lomi_home/wenxiu/tmp --reverse-fragments 1
for i in `ls ${CRISPRdect_file}/bin*CRISPRDetect3.spacers.fa`; do
fn=$(basename $i)
spacepharer easy-predict $i ${spacepharer_file}/targetSetDB/viralset ${spacepharer_file}/${fn}_spacepharer_predictions.tsv /lomi_home/wenxiu/tmp
done
conda deactivate
