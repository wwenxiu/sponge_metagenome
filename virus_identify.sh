### identify vOTUs using virsorter,virsorter2, virfinder,deepvirfinder, and CAT
Assembly_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/assembly/VP_mrna_ass"
viral_iden_file="/lomi_home/wenxiu/HIC-meta/CCs_analyze/virus_iden/VP_virus_pip"
virsorter_file=${viral_iden_file}"/virsorter"
virfinder_file=${viral_iden_file}"/virfinder"
virsorter2_file=${viral_iden_file}"/virsorter2"
deepvirfinder_file=${viral_iden_file}"/deepvirfinder"
threads=32

#1.5kb
export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/seqmagick/bin
for fa in `ls ${Assembly_file}/*fa`; do
seqmagick convert --min-length 1500 $fa ${fa/.fa}.min1.5k.fa
done

#virsorter
export PERL5LIB=/lomi_home/wenxiu/software/miniconda3_4.10/envs/virsorter/lib/perl5/site_perl/5.22.0/
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate virsorter
virsorter_file=${viral_iden_file}"/virsorter"
mkdir -p ${virsorter_file} 
for fa in `ls ${Assembly_file}/*min1.5k.fa`; do
assembly_bn=`basename $fa`
assembly_fs=${assembly_bn/.*/}
wrapper_phage_contigs_sorter_iPlant.pl -f $fa --db 1 --wdir ${virsorter_file}/${assembly_fs} --ncpu $threads --data-dir /lomi_home/wenxiu/software/virsorter-data
done
conda deactivate 

#virfinder
virfinder_file=${viral_iden_file}"/virfinder"
source activate virfinder
mkdir -p ${virfinder_file} 
for fa in `ls ${Assembly_file}/*min1.5k.fa`; do
assembly_bn=`basename $fa`
assembly_fs=${assembly_bn/.*/}
Rscript /lomi_home/wenxiu/software/miniconda3_4.10/envs/virfinder/bin/virfinder.R $fa ${virfinder_file}/${assembly_fs}.csv
done
conda deactivate 

#virsorter2
virsorter2_file=${viral_iden_file}"/virsorter2"
source activate vs2
mkdir -p ${virsorter2_file} 
for fa in `ls ${Assembly_file}/*min1.5k.fa`; do
assembly_bn=`basename $fa`
assembly_fs=${assembly_bn/.*/}
virsorter run -w ${virsorter2_file}/${assembly_fs} -i $fa --min-length 2000 -j $threads all
done
conda deactivate 

#deepvirfinder
deepvirfinder_file=${viral_iden_file}"/deepvirfinder"
source activate dvf
mkdir -p ${deepvirfinder_file} 
cd ${deepvirfinder_file}
for fa in `ls ${Assembly_file}/*min1.5k.fa`; do
assembly_bn=`basename $fa`
assembly_fs=${assembly_bn/.*/}
python /lomi_home/wenxiu/software/DeepVirFinder/dvf.py -i $fa -o ${deepvirfinder_file}/${assembly_fs} -l 300 -c 4
done
conda deactivate

#CAT
CAT_file=${viral_iden_file}"/CAT"
source activate cat
for fa in `ls ${Assembly_file}/*min1.5k.fa`; do
assembly_bn=`basename $fa`
assembly_fs=${assembly_bn/.*/}
CAT contigs -c $fa -d /lomi_home/wenxiu/database/CAT_db/CAT_prepare_20210107/2021-01-07_CAT_database -t /lomi_home/wenxiu/database/CAT_db/CAT_prepare_20210107/2021-01-07_taxonomy  -o ${CAT_file}/${assembly_fs}
CAT add_names -i ${CAT_file}/${assembly_fs}.contig2classification.txt -o ${CAT_file}/${assembly_fs}_CAT_tax.txt -t /lomi_home/wenxiu/database/CAT_db/CAT_prepare_20210107/2021-01-07_taxonomy --only_official
done
conda deactivate

##cluster virus sequnce after summarizing outputs and get vp_mrna.min5k.cir.fa
module load blast/2.8.1+
virus_file="vp_mrna.min5k.cir.fa"
#ani cluster
makeblastdb -in ${virus_file} -dbtype nucl -out ${virus_file}
blastn -query ${virus_file} -db ${virus_file} -outfmt '6 std qlen slen' -max_target_seqs 10000 -out ${virus_file}.blast.tsv -num_threads 32 
python /lomi_home/wenxiu/software/checkv/scripts/anicalc.py -i ${virus_file}.blast.tsv -o ${virus_file}.blast.ani.tsv
python /lomi_home/wenxiu/software/checkv/scripts/aniclust.py --fna ${virus_file} --ani ${virus_file}.blast.ani.tsv --out ${virus_file}.blast.ani.cluster.tsv --min_ani 99 --min_tcov 85 --min_qcov 0
cat ${virus_file}.blast.ani.cluster.tsv | awk -F'\t' '{print $1}' | sort | uniq > cluster.list

export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/seqkit/bin
seqkit grep -f cluster.list ${virus_file} > ${virus_file/.fa/}.0.95-0.85.fasta

##vOTU taxonmony classification
#PhaGCN2
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
source activate phagcn2

python /lomi_home/wenxiu/software/PhaGCN2.0/run_Speed_up.py --contigs viral.contigs.min5k.cir.0.99-0.85.fasta --len 10000

