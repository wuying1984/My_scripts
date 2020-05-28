##### Convert gff to gtf #####
gffread my.gff3 -T -o my.gtf

##### get sequence from fasta file#####
for i in `cat Nb_LRR-RLK.list`; do cat Niben101_annotation.proteins.fa| grep -A1 $i; done >Nb_LRR-RLK.fa

##### do blastp #####
makeblastdb -in protein.fa -input_type fasta -dbtype prot -titleb TITLE -out TITLE
blastp -db NB101 -query AtPTI_genes.fasta  -evalue 0.0001 -seg no  -num_threads 8 -outfmt "7 qacc sacc qlen slen qstart qend sstart send evalue score nident mismatch length pident qcov qcovhsp"  -out AtPTI_genes.blastNB101

----------------------------------------------------------------------------------------------------------------------------------------
##### split gff into non-seq part (repeatmask, augustus,genemark,snap,est2genome,protein2genome) #####
cat Bfra_R1V1.fa.all.gff | grep ^tig | grep -v "match_part" |awk '{print $_>"part."$2".gff"}'
mv part...gff part.genome.gff

##### merge bed file #####
cat part.est2genome.gff| awk '{print $1"\t"$4"\t"$5}' | sort -k1,1 -k2,2n -k3,3n >est2genome.bed
bedtools merge -i est2genome.bed | awk '{print $1"\test2genome\texpressed_sequence_match\t"$2"\t"$3"\t111\t+\t.\t"$1}' >part.est2genome.gff
cat part.protein2genome.gff|awk '{print $1"\t"$4"\t"$5}' | sort -k1,1 -k2,2n -k3,3n >protein2genome.bed
bedtools merge -i protein2genome.bed | awk '{print $1"\tprotein2genome\tprotein_match\t"$2"\t"$3"\t111\t+\t.\t"$1}' >part.protein2genome.gff
cat Bfra_R1V1.fa.all.gff | awk '$2=="maker"' | awk '$3=="gene"' | grep ^tig | awk '{print $_>"gene."$2".gff"}'
cat Bfra_R1V1.fa.all.gff | awk '$2=="maker"' | awk '$3=="CDS"' | grep ^tig | awk '{print $_>"CDS."$2".gff"}'

##### get intersection between CDS/gene and other features
##### 第一列是overlap的百分比
for i in `echo repeatmasker est2genome protein2genome augustus_masked`;
	bedtools intersect -wo -a gene.maker.gff -b part.$i.gff | awk '{print $19/($5-$4+1)"\t"$_}' >intersect.$i.txt
done
for i in `echo repeatmasker est2genome protein2genome augustus_masked`;do bedtools intersect -wo -a CDS.maker.gff -b part.$i.gff | awk '{print $19/($5-$4+1)"\t"$_}' >intersect.$i.txt;done

##### use get_geneID_for_checkup.pl to get the overlapping status
perl get_geneID_for_checkup.pl
---------------------------------------------------------------------------------------------------------------------------------------

##########get intron bed##################################
perl get_all_intron_bed.pl annotation/gencode.v28.annotation.gtf | bedtools getfasta -s -name -fi annotation/GRCh38.p12.genome.fa -bed - -fo annotation/gencode.v28.intron.fa
perl get_all_intron_bed.pl Bfra_R1V1.fa.all.gtf >part.intron.bed

##### gtf to CDS fasta################################
 perl ../../get_all_CDS_bed.pl Bfra_R1V1.fa.all.gtf |  bedtools getfasta -s -name -fi ~/Xiaolab/Botrytis/Maker/Bfra_R1V1.fa -bed - -fo Bfra_individual_CDS.fa
 perl ../../connect_individual_CDS.pl Bfra_add_individual_CDS.fa Bfra_add_CDS.fa
 
 ##### CDS to protein##################################
 ~/miniconda3/bin/python convert_cds_to_pep_for_fasta.py Bfra_CDS.fa
 
 ##### get how many *(stop codon) in the sequence#################
 for i in `cat ../../Bfra_CDS.pep.fa | grep -v ">" | grep "*$"` ; do echo -ne $i"\t"; echo $i| tr -dc "*"|length;done| awk '$2=="2"'
 
 ##### check middle *##########################
 for i in `cat Bfra_add_CDS.pep.fa | grep -v ">" | grep "*$"` ; do echo -ne $i"\t"; echo $i| tr -dc "*"|length;done| awk '$2=="2"'
 
 ##### for jbrowse ########
 cp ${i/_trinity/}_6th.gff ~/mnt/jbrowse/;cd ../;
 maker2jbrowse UCSC1.final.gff -o UCSC1
 
 ##### if the jbrowse is not available on cluster, just download the desttop style of jbrowse and selecte the gff containing the probmatic gene for desplay
 ##### gene needs to be checked
 ##### 1. overlap with repeat >30%
 ##### 2. CDS length <10bp
 ##### 3. intron <10bp (not allowed when upload to NCBI)
 ##### 4. gene overlap with each other
 ##### 5. gene too long >10000
 ##### 6. not start with M or end with *
 ##### 7. overlap with more than one times of augustus/genemark/snap (not include yet)
c Bfra_CDS.pep.fa| grep ">"  |sed "s#>##" >all_gene_list.txt
c all_gene_list.txt | sort >A.txt
c qualified_gene_list.txt | sort >B.txt
comm A.txt B.txt -2 -3
comm A.txt B.txt -2 -3| wc -l
comm A.txt B.txt -2 -3 >gene_to_be_checked_list.txt
c gene_list_for_checkup_intron.tab >>gene_to_be_checked_list.txt
for i in `cat gene_to_be_checked_list.txt`; do cat gene.maker.gff| grep $i; done >gene_to_be_checked_list.gff
c ../Bfra_R1V1.fa.all.gff| grep ^tig >Bfra_R1V1.fa.all_no_seq.gff
bedtools intersect -wa -a Bfra_R1V1.fa.all_no_seq.gff -b gene_to_be_checked_list.gff |uniq >gene_for_jbrowse.gff
c gene_for_jbrowse.gff| awk '$2~/maker/' >gene_for_jbrowse_gene.gff
c gene_for_jbrowse.gff| awk '$2~/protein2genome/' >gene_for_jbrowse_protein2genome.gff
c gene_for_jbrowse.gff| awk '$2~/repeatmasker/' >gene_for_jbrowse_repeatmasker.gff
c gene_for_jbrowse.gff | awk '$2~/cdna2genome/||$2~/est2genome/' >gene_for_jbrowse_cna.gff
c gene_for_jbrowse.gff | awk '$2~/augustus_masked/||$2~/snap_masked/ ||$2~/genemark/' >gene_for_jbrowse_software.gff
 
 https://jbrowse.org/blog/2020/02/04/jbrowse-1-16-8.html
for i in `head test1.txt`; do cat gene.maker.gff| grep $i| awk '{print "'$i'""\t"$1"\t"$4".."$5}';done
~/miniconda3/bin/python replace_add_rm.py >output/Bfra_R1V1.fa.all_no_seq.modified.gff
~/miniconda3/bin/python sort_gff.py >output/Bfra_R1V1.fa.all_no_seq.modified.sorted.gff

#############get overlapped gene##########################
c  Bfra_R1V1.fa.all_no_seq.modified.sorted.gff| awk '$3~/gene/'>gene.gff
intersectBed -wo -a gene.gff -b gene.gff | awk '$9!=$18' | awk '{print $9}' | awk -F\; '{print $1}' | sed 's@ID=@@' | us >genes_with_intersect.txt
intersectBed -wo -a gene.gff -b gene.gff | awk '$9!=$18' | awk '{print $9}' | awk -F\; '{print $1}' | sed 's@ID=@@' | us

for i in `cat genes_with_intersect.txt`;do cat gene.gff| grep $i; done >overlap_gene.gff
bedtools intersect -wa -a ../../Bfra_R1V1.fa.all_no_seq.gff -b overlap_gene.gff |uniq >overlap_gene_for_jbrowse.gff
bedtools intersect -wa -a ../Bfra_R1V1.fa.all_no_seq.modified.sorted.gff -b overlap_gene.gff |uniq| awk '$2~/maker/' >overlap_gene_for_jbrowse_gene.gff
c overlap_gene_for_jbrowse.gff| awk '$2~/maker/' >overlap_gene_for_jbrowse_gene.gff
c overlap_gene_for_jbrowse.gff| awk '$2~/protein2genome/' >overlap_gene_for_jbrowse_protein2genome.gff
c overlap_gene_for_jbrowse.gff| awk '$2~/repeatmasker/' >overlap_gene_for_jbrowse_repeatmasker.gff
c overlap_gene_for_jbrowse.gff | awk '$2~/cdna2genome/||$2~/est2genome/' >overlap_gene_for_jbrowse_cna.gff
for i in `cat genes_with_intersect.txt`; do cat overlap_gene_for_jbrowse_gene.gff |grep $i| awk '$3~/gene/{print "'$i'""\t"$1":"$4".."$5}';done >overlap_gene_for_checkup.txt
c overlap_gene_for_jbrowse.gff | awk '$2~/augustus_masked/||$2~/snap_masked/ ||$2~/genemark/' >overlap_gene_for_jbrowse_software.gff

 
##### contig assembly stat
~/miniconda3/bin/assembly-stats ../Bfra_R1V1.fa
#########scaffolds contig assembly stat
seqtk cutN 1 ../Bfra_R1V1.fa >contig.fasta

##### http://tiramisutes.github.io/2016/08/30/NCBI-submission.html

#add header to fasta
c Bfra_R1V1.fa| sed 'n;d' >A
c Bfra_R1V1.fa| sed -n 'n;p' >B
c A | wc
c B | wc
cat A | awk '{print $1" [organism=Botrytis fragariae] [isolate=BVB16]"}' >C
paste C B | tr '\t' '\n' >D
./linux64.tbl2asn -p Bfra -t template_Bfra_2020.sbt -M n -Z discrepency.report.txt  
##### template_Bfra_2020.sbt fill the website online and save it
##### ./linux64.tbl2asn is downloaded from NCBI
-a r50k -l
 
 
 
