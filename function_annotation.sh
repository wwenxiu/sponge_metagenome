### function annotation for virus
viralset="/lomi_home/wenxiu/HIC-meta/CCs_analyze/virus_iden/fnl_viralset/vp_mrna.min5k.cir.0.99-0.85.fasta"
viralset_protein="/lomi_home/wenxiu/HIC-meta/CCs_analyze/virus_iden/fnl_viralset/vp_mrna.min5k.cir.0.99-0.85.protein.faa"

#orf predict
export PATH=$PATH:/lomi_home/wenxiu/software/miniconda3_4.10/bin
export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/seqkit/bin
source activate phanotate
phanotate.py -o ${viralset}.gene.fa -f fasta $viralset
seqkit translate ${viralset/.fasta}.gene.fa > ${viralset/.fasta}.protein.faa 
conda deactivate

#pfam
export PATH=$PATH:/lomi_home/wenxiu/software/hmmer/bin
hmmscan --cpu 32 --domtblout viralset.pfam35.out.dm /lomi_home/wenxiu/database/Pfam-A/Pfam35.0/Pfam-A.hmm ${viralset/.fasta}.protein.faa 
sh /lomi_home/wenxiu/database/CAzyme/db_V10/hmmscan-parser_pfam.sh viralset.pfam35.out.dm > viralset.pfam34.out.ps
cat viralset.pfam35.out.ps | awk '$5<1e-15&&$10>0.35' > viralset.pfam35.out.stringent.txt

#kofam
hmmscan --cpu 32 --domtblout viralset.kofam.out.dm /lomi_home/wenxiu/database/kofam/kofam_profiles.hmm ${viralset/.fasta}.protein.faa 
sh /lomi_home/wenxiu/database/CAzyme/db_V10/hmmscan-parser.sh viralset.kofam.out.dm > viralset.kofam.out.ps
cat viralset.kofam.out.ps | awk '$5<1e-15&&$10>0.35' > viralset.kofam.out.stringent.txt

#VOG
hmmscan --cpu 32 --domtblout viralset.vog.out.dm /lomi_home/wenxiu/database/VOGDB/VOG_v214/All_vog.hmm ${viralset/.fasta}.protein.faa
sh /lomi_home/wenxiu/database/CAzyme/db_V10/hmmscan-parser.sh viralset.vog.out.dm > viralset.vog.out.ps
cat viralset.vog.out.ps | awk '$5<1e-15&&$10>0.35' > viralset.vog.out.stringent.txt


### function annotation for archaeal/bacterial genomes
bin_dir="/lomi_home/wenxiu/HIC-meta/CCs_analyze/bin_annotation/bins"
inputset="pro_genomes"

#orf predict
export PATH=$PATH:/share/software/miniconda/miniconda3-latest/envs/prodigal/bin
for i in ${bin_dir}/*fa
do prodigal -a ${i/.fa/}.protein.faa -c -d ${i/.fa/}.gene.fna -i $i -m -f gff -o ${i/.fa/}.gff -p meta
done

#pfam and kofam 
export PATH=$PATH:/lomi_home/wenxiu/software/GNU/parallel/bin
export PATH=$PATH:/lomi_home/wenxiu/software/hmmer/bin
ln -s ${bin_dir}/*faa .

ls *faa | parallel -j 15 "hmmscan --cpu 2 --domtblout {}.single.pfam34.out.dm /lomi_home/wenxiu/database/Pfam-A/Pfam34.0/Pfam-A.hmm {}"
cat *.single.pfam34.out.dm > ${inputset}.pfam34.out.dm
sh /lomi_home/wenxiu/database/CAzyme/db_V10/hmmscan-parser_pfam.sh ${inputset}.pfam34.out.dm > ${inputset}.pfam34.out.ps
cat ${inputset}.pfam34.out.ps | awk '$5<1e-15&&$10>0.35' > ${inputset}.pfam34.out.stringent.txt

ls *faa | parallel -j 15 "hmmscan --cpu 2 --domtblout {}.single.kofam.out.dm /lomi_home/wenxiu/database/kofam/kofam_profiles.hmm {}"
cat *.single.kofam.out.dm > ${inputset}.kofam.out.dm
sh /lomi_home/wenxiu/database/CAzyme/db_V10/hmmscan-parser.sh ${inputset}.kofam.out.dm > ${inputset}.kofam.out.ps
cat ${inputset}.kofam.out.ps | awk '$5<1e-15&&$10>0.35' > ${inputset}.kofam.out.stringent.txt