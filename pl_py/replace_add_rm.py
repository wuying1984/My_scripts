#!/usr/bin/env python

# cat Bfra_R1V1.fa.all_no_seq.gff | awk '$2=="maker"' >Bfra_R1V1.fa.all_no_seq.maker.gff
# cat Bfra_R1V1.fa.all_no_seq.gff | awk '$2!="maker"' >Bfra_R1V1.fa.all_no_seq.other.gff
from __future__ import print_function, division
from  more_itertools import unique_everseen
import sys

def get_replace(filename):
    """ replace """
    replace_dct = dict()
    with open(filename, "r") as f:
       for line in f.readlines()[1:]:
           action, oldid, newid, scaffold = line.rstrip().split("\t")
           replace_dct[oldid] = newid
    return(replace_dct)

def get_add(filename):
    """ add """
    add_dct = dict()
    with open(filename, "r") as f:
       for line in f.readlines()[1:]:
           action, newid, oldid, scaffold = line.rstrip().split("\t")
           add_dct[newid] = oldid
    return(add_dct)

def get_rm(filename):
    """ rm """
    rm_dct = dict()
    with open(filename, "r") as f:
       for line in f.readlines()[1:]:
           action, oldid, scaffold = line.rstrip().split("\t")
           rm_dct[oldid] = 1
    return(rm_dct)

def print_replaced_lines(lines2, id, id2):
    """ replace lines2 (non-maker gff lines) """
    dct_iCDS_lines = dict()
    for line in lines2:
        scaffold, feature, structure, start, end, dot1, strand, dot2, other = line.rstrip().split("\t")
        if structure == "match":
            other = "ID=" + id + ";Name=" + id + ";Alias=" + id2 + ";"
            # print gene line
            print("\t".join([scaffold, "maker", "gene", start, end, dot1, strand, dot2, other]))
            # print mRNA line

            other = "ID=" + id + "-RA;Parent=" + id + ";Name=" + id + "-RA;Alias=" + id2 + ";"
            print("\t".join([scaffold, "maker", "mRNA", start, end, dot1, strand, dot2, other]))
        elif structure == "match_part":
            other = "ID=" + id + "-RA:cds;Parent=" + id + "-RA;"
            if not strand in dct_iCDS_lines.keys():
                dct_iCDS_lines[strand] = dict()
            dct_iCDS_lines[strand][int(start)] = "\t".join([scaffold, "maker", "CDS", start, end, dot1, strand, dot2, other])

    phase = 0
    current_iCDSs_len = 0
    for strand in dct_iCDS_lines.keys():
        if strand == "+":
            sorted_starts = sorted(dct_iCDS_lines[strand].keys())
        elif strand == "-":
            sorted_starts = sorted(dct_iCDS_lines[strand].keys(), reverse=True)

        for start in sorted_starts:
            # print("###", start, phase)
            scaffold, feature, structure, start, end, dot1, strand, dot2, other = dct_iCDS_lines[strand][start].split("\t")
            # exon
            print("\t".join([scaffold, feature, "exon", start, end, dot1, strand, dot2, other.replace("cds","exon")]))


        for start in sorted_starts:
            scaffold, feature, structure, start, end, dot1, strand, dot2, other = dct_iCDS_lines[strand][start].split("\t")
            # CDS
            print("\t".join([scaffold, feature, structure, start, end, dot1, strand, str(phase), other]))

            iCDS_len = int(end) - int(start) + 1
            current_iCDSs_len += iCDS_len
            phase = 3 - current_iCDSs_len % 3
            if phase == 3: phase = 0

def print_lines(lines):
    for line in lines:
        print(line.rstrip())







class GFF:
    """ GFF object (now only use for maker) """
    def __init__(self, filename):
        self.filename = filename
        self.get_id2lines()
        
    def get_id2lines(self):
        self.id2lines = dict()
        with open(self.filename, "r") as f:
            for line in f.readlines():
                #print(line)
                if not line.startswith("#"):
                    scaffold, feature, structure, start, end, dot1, strand, dot2, other = line.rstrip().split("\t")

                    id = other.rstrip().split("ID=")[1].split(";")[0].split("-RA")[0]
                    if not id in self.id2lines.keys():
                        self.id2lines[id] = list()
                    self.id2lines[id].append(line.rstrip())

                    if "Name=" in other:
                        id2 = other.rstrip().split("Name=")[1].split(";")[0].split("-RA")[0]
                        if not id2 in self.id2lines.keys():
                            self.id2lines[id2] = list()
                        self.id2lines[id2].append(line.rstrip())

                    if "Target=" in other:
                        id3 = other.rstrip().split("Target=")[1].split(" ")[0].split("-RA")[0]
                        if not id3 in self.id2lines.keys():
                            self.id2lines[id3] = list()
                        self.id2lines[id3].append(line.rstrip())


if __name__ == '__main__':
    replace_dct = get_add("data/bfra_replace.txt")
    add_dct = get_add("data/bfra_add.txt")
    rm_dct = get_rm("data/bfra_rm.txt")

    g1 = GFF("data/Bfra_R1V1.fa.all_no_seq.maker.gff")
    g2 = GFF("data/Bfra_R1V1.fa.all_no_seq.other.gff")

    for id in g1.id2lines.keys():
        if id in replace_dct:
            id2 = replace_dct[id]
            print_lines(list(unique_everseen(g2.id2lines[id2])))
        elif id in rm_dct:
            pass
        else:
            print_lines(list(unique_everseen(g1.id2lines[id])))
            
    for id2 in add_dct.keys():
        id = add_dct[id2]
        lines2 = g2.id2lines[id2]
        print_replaced_lines(list(unique_everseen(lines2)), id, id2)

         

