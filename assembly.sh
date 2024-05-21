###prokaryotic/eukaryotic sample reads assembly
##spades
export PATH=$PATH:/lomi_home/wenxiu/software/SPAdes-3.15.3-Linux/bin
cleanReads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/PC_DNA_cleanReads"
for fq in ${cleanReads_file}/*_1.fastq; do
filename=$(basename $fq)
python /lomi_home/wenxiu/software/SPAdes-3.15.3-Linux/bin/metaspades.py --meta -s1 $fq -s2 ${fq/_1.fq/_2.fq} -o ${filename/_1.fastq}_spades -t 32 -m 500
done

##megahit
export PATH=$PATH:/share/software/megahit
megahit --presets meta-sensitive -t 32 -m 0.9 -r ${cleanReads_file}/CC?M*cutsite*.fastq -o CC_meghit_asm

##opera-ms, after combine the contigs assembled from short reads and 454 contigs using CDHIT
ass_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/assembly/PC_DNA_illumina_ass"
long_reads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/ST_DNA_ont_Reads"

export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate perl
perl /lomi_home/wenxiu/software/OPERA-MS/OPERA-MS.pl \
    --contig-file ${ass_file}/pooled_assembly.fa \
    --short-read1 ../combine_1.fastq \
    --short-read2 ../combine_2.fastq \
    --long-read ${long_reads_file}/PC_CCS.fq \
    --no-ref-clustering \
    --num-processors 32 \
    --out-dir PC_CCS_res > log.err
conda deactivate

##flye
longReads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/ST_DNA_ont_Reads"
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate flye
for fq in `ls ${longReads_file}/*fq`; do
filename=$(basename $fq)
flye --nano-raw $fq --out-dir ${filename/.fq} --threads 32 --meta
done
conda deactivate

##nextpolish using HiC cutsite reads, after combine the hybrid assembled contigs and contigs assembled from ONT reads of sponges DNA using CDHIT
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate NextPolish
nextPolish ./run.cfg
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 16
multithread_jobs = 32
genome = ./PC_assembly_min1k.fa
genome_size = auto
workdir = ./outdir
polish_options = -p {multithread_jobs}

[sgs_option] #optional
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa --unpaired

###viral sample reads assembly
##spades 
shortReads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/VP_DNA_illumina_cleanReads"
export PATH=$PATH:/lomi_home/wenxiu/software/SPAdes-3.15.3-Linux/bin

for fq in `ls ${shortReads_file}/*_1.fastq`; do
filename=$(basename $fq)
python /lomi_home/wenxiu/software/SPAdes-3.15.3-Linux/bin/spades.py --metaviral -1 ${fq} -2 ${fq/_1/_2} -o ${filename/_1.fastq}_metaspades -t 32 -m 500
done

##opera-ms, after combine the contigs assembled from short reads using CDHIT
ass_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/assembly/VP_DNA_illumina_ass"
longReads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/VP_DNA_ont_Reads"
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate perl
perl /lomi_home/wenxiu/software/OPERA-MS/OPERA-MS.pl \
    --contig-file ${ass_file}/pooled_assembly.fa \
    --short-read1 ../combine_1.fastq \
    --short-read2 ../combine_2.fastq \
    --long-read ${long_reads_file}/VP_CCS_pooled.fq \
    --no-ref-clustering \
    --num-processors 32 \
    --out-dir VP_CCS_res > log.err
conda deactivate

##flye
longReads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/VP_DNA_ont_Reads"
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate flye
for fq in `ls ${longReads_file}/*fq`; do
filename=$(basename $fq)
flye --nano-raw $fq --out-dir ${filename/.fq} --threads 32 --meta
done
conda deactivate

##nextpolish using short illumina reads, after combine the hybrid assembled contigs and contigs assembled from ONT reads of  viral DNA using CDHIT
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate NextPolish

nextPolish ./run.cfg
##run_paired.cfg
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 16
multithread_jobs = 32
genome = ./VP_assembly_min1k.fa
genome_size = auto
workdir = ./outdir
polish_options = -p {multithread_jobs}

[sgs_option] #optional
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa 

###mRNA reads assembly
cleanReads_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/mRNA_illumina_cleanReads"
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate trinity
for fq in `ls ${cleanReads_file}/*_1.fastq`; do
filename=$(basebame $fq)
Trinity --seqType fq --left $fq --right ${fq/_1.fastq/_2.fastq} --CPU 32 --max_memory 500G --output ${filename/_1.fastq}_trinity
done

##rnaspades
export PATH=$PATH:/lomi_home/wenxiu/software/SPAdes-3.15.3-Linux/bin
for fq in `ls ${cleanReads_file}/*_1.fastq`; do
filename=$(basename $fq)
python /lomi_home/wenxiu/software/SPAdes-3.15.3-Linux/bin/spades.py --rna -1 ${fq} -2 ${fq/_1/_2} -o ${filename/_1.fastq}_rnaspades -t 32 -m 500
done