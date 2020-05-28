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
```
 perl connect_individual_CDS.pl individual_CDS.fa Whole_CDS.fa
```
-rw-r--r-- 1 ywu domain users   893 May 12 08:13 convert_cds_to_pep_for_fasta.py
-rw-r--r-- 1 ywu domain users  2331 May 12 08:13 get_all_CDS_bed.pl
-rw-r--r-- 1 ywu domain users  2333 May 11 14:57 get_all_exon_bed.pl
-rw-r--r-- 1 ywu domain users  3018 May 12 08:13 get_all_intron_bed.pl
-rw-r--r-- 1 ywu domain users  4841 May 12 08:13 get_geneID_for_checkup.pl
-rw-r--r-- 1 ywu domain users   725 May 27 18:07 get_partial.pl
-rw-r--r-- 1 ywu domain users  1838 May 12 08:13 get_seq_need.pl
-rwxr-xr-x 1 ywu domain users 17886 May 27 18:06 get_tbl_new_Bfra.pl
-rwxr-xr-x 1 ywu domain users 16632 May 27 18:07 get_tbl_new.pl
-rw-r--r-- 1 ywu domain users  5412 May 27 18:16 replace_add_rm.py
-rw-r--r-- 1 ywu domain users  2024 May 27 18:16 sort_gff.py
