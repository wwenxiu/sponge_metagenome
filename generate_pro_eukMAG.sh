###binning to get archaeal/bacterial/eukaryotic MAG
##BASALT binning
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3/bin
source activate BASALT
BASALT -a CC1M_spades.fasta, CC2M_spades.fasta, CC3M_spades.fasta, CC1T_spades.fasta, CCS_illumina_contigs -s CC1T-HiC_citsite_1.fastq,CC1T-HiC_citsite_2.fastq/CC1M-HiC_citsite_1.fastq,CC1M-HiC_citsite_2.fastq/CC2M-HiC_citsite_1.fastq,CC2M-HiC_citsite_2.fastq/CC3M-HiC_citsite_1.fastq,CC3M-HiC_citsite_2.fastq -l PC_CCS_clean.fq -t 32 -m 500 --autopara more-sensitive --refinepara deep --max-ctn 100 --min-cpn 0 --mode continue

conda deactivate

##VAMB binning
#minimap2
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/samtools/bin
source activate minimap2
minimap2 -d catalogue_ont.mmi assembly_ont.fasta # make index
minimap2 -t 32 -N 5 -ax map-ont catalogue_ont.mmi CCS_clean.fq | samtools view -F 3584 -b --threads 32 > CCS_ont.bam

conda deactivate
source activate vamb
vamb --outdir vamb_out --fasta assembly_ont.fasta --bamfiles CCS_ont.bam

##eukaryotic mag classification
##tiara
source activate tiara-env
for i in raw_binset/*.fa; do
base=$(basename $i)
tiara -i $i -o tiara/${base}_out -t 32
done
conda deactivate

##busco 
source activate busco
for i in possible_euk_binset/*fa; do
BASE=`basename $i`
busco -m genome -i $i -o ${BASE/.fa/}.busco -c 1 -e 1e-03
done

##MAG checkm and drep
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate checkm2
checkm2 predict --threads 56 --input raw_binset --output-directory checkm2_out -x fa --database_path /lomi_home/wenxiu/database/CheckM2_database/uniref100.KO.1.dmnd
conda deactivate
source activate checkm-genome
checkm lineage_wf -t 32 -x fa raw_binset checkm_out

#get high quality MAG
module load miniconda3
source activate drep-2.6.2
dRep dereplicate derep_file -g HQ_MAG/*.fa -pa 0.99 -sa 0.99 -p 32 -comp 0 -con 100
conda deactivate 
module unload miniconda3

##archaeal/bacterial classification
export PATH=$PATH:/share/software/miniconda3_py39/bin
source activate gtdbtk-1.7.0
gtdbtk classify_wf --cpus 32 --genome_dir final_probinset --out_dir gtdbtk_out -x fa

##archaeal/bacterial phylogeny
#marker find
export PATH=$PATH:/share/software/hmmer/hmmer-3.2.1/bin
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/envs/clustalomega/bin
export PATH=$PATH:/tiagor_home/jianchang/software/miniconda3_4.10/envs/python3.10/bin
python3 /lomi_home/wenxiu/software/markerfinder/markerfinder.py -i final_binset pro_tree -t 32 -c 

#trimal
export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/trimal/bin
trimal -in pro_tree.concat.aln -out pro_tree.concat.aut1.aln -automated1

#iqtree
export PATH=$PATH:/share/software/iqtree/iqtree-2.0/bin
iqtree -s pro_tree.concat.aut1.aln -m MFP -B 1000 --bnni -T 56
