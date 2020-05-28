### data folder containing some files used as input for the scripts.

### 1. combineTABwith.pl 
#### conbine two files based on the specified column
```
combineTABwith.pl A 1 B 1
#1 is the column to be combined
#A column 1 can have duplicate, B column 1 is unique
#A, B is tab split file
#  A: a 1    B: a aaa   output:  a 1 a aaa
      b 2       b bbb            b 2 b bbb
      c 3       c ccc            c 3 c ccc
      c 1       d ddd            c 1 c ccc
```

### 2. connect_individual_CDS.pl
#### conect individual CDS from one gene into whole CDS
```
# use get_all_CDS_bed.pl to get individual CDS file
 perl connect_individual_CDS.pl individual_CDS.fa Whole_CDS.fa
```

### 3. convert_cds_to_pep_for_fasta.py
#### convert CDS sequences into peptides
```
 ~/miniconda3/bin/python convert_cds_to_pep_for_fasta.py CDS.fa
# outfile is CDS.pep.fa
```

### 4.get_all_CDS_bed.pl
```
## first convert gff to gtf
gffread ../Bfra_R1V1.fa.all.gff -T -o Bfra_R1V1.fa.all.gtf
## get individual CDS sequence
perl get_all_CDS_bed.pl input.gtf |  bedtools getfasta -s -name -fi genome.fa -bed - -fo individual_CDS.fa
```
 
### 5. get_all_exon_bed.pl
#### similar to get_all_CDS_bed.pl, but get the exon bed or exon fasta
 
### 6. get_all_intron_bed.pl
#### similar to get_all_CDS_bed.pl, but get the intron bed or intron fasta
 
### 7. get_geneID_for_checkup.pl
#### get the gene ID that overlap with repeatmask (may need a bit change to exclude simple repeat)
####                                   est2genome
####                                   protein2genome
####                                   AUGUSTUS,GeneMark,SNAP (may need a bit change to output how many AUGUSTUS id overlap with the ID
                                       
 ```
 perl get_geneID_for_checkup.pl
 ```

### 8. get_partial.pl
#### check if a gene starts with M and ends with *
#### input a fasta file, change the name in the script
```
perl get_partial.pl
```

### 9.get_seq_need.pl
#### input cds file, pep file and the id need, output the seq needed

### 10.get_tbl_new_Bfra.pl & 11. get_tbl_new.pl
#### convert gff file into tbl file
#### change the input file and check the script before use it

### 12. replace_add_rm.py
#### replace the changed gff for specific gene
#### add gff for specific genes
#### remove gff for specific genes
#### file format is in `data folder`
```
# replace_dct = get_add("data/bfra_replace.txt")
# add_dct = get_add("data/bfra_add.txt")
# rm_dct = get_rm("data/bfra_rm.txt")
g1 = GFF("data/Bfra_R1V1.fa.all_no_seq.maker.gff")
g2 = GFF("data/Bfra_R1V1.fa.all_no_seq.other.gff")
~/miniconda3/bin/python replace_add_rm.py >output/Bfra_R1V1.fa.all_no_seq.modified.gff
```

### 13.sort_gff.py
#### sort gff file after replace_add_rm
```
~/miniconda3/bin/python sort_gff.py >output/Bfra_R1V1.fa.all_no_seq.modified.sorted.gff
```
