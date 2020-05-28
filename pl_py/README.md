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

`
 ~/miniconda3/bin/python convert_cds_to_pep_for_fasta.py CDS.fa
`

### 4.get_all_CDS_bed.pl
```
## first convert gff to gtf
gffread ../Bfra_R1V1.fa.all.gff -T -o Bfra_R1V1.fa.all.gtf
## get individual CDS sequence
perl get_all_CDS_bed.pl input.gtf |  bedtools getfasta -s -name -fi genome.fa -bed - -fo individual_CDS.fa
```
 
 ### 5. 2333 May 11 14:57 get_all_exon_bed.pl
-rw-r--r-- 1 ywu domain users  3018 May 12 08:13 get_all_intron_bed.pl
-rw-r--r-- 1 ywu domain users  4841 May 12 08:13 get_geneID_for_checkup.pl
-rw-r--r-- 1 ywu domain users   725 May 27 18:07 get_partial.pl
-rw-r--r-- 1 ywu domain users  1838 May 12 08:13 get_seq_need.pl
-rwxr-xr-x 1 ywu domain users 17886 May 27 18:06 get_tbl_new_Bfra.pl
-rwxr-xr-x 1 ywu domain users 16632 May 27 18:07 get_tbl_new.pl
-rw-r--r-- 1 ywu domain users  5412 May 27 18:16 replace_add_rm.py
-rw-r--r-- 1 ywu domain users  2024 May 27 18:16 sort_gff.py
