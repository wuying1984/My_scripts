#!/usr/bin/env python

from __future__ import print_function, division
from  more_itertools import unique_everseen
import sys

def print_lines(lines):
    for line in lines:
        print(line.rstrip())




class GFF:
    """ GFF object (now only use for maker) """
    def __init__(self, filename):
        self.filename = filename
        self.dct = {} # record gene location
        self.get_id2lines()

    def get_id2lines(self):
        self.id2lines = dict()
        with open(self.filename, "r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    scaffold, feature, structure, start, end, dot1, strand, dot2, other = line.rstrip().split("\t")
                    start = int(start)
                    end = int(end)
                    id = other.rstrip().split("ID=")[1].split(";")[0].split("-RA")[0]
                    self.dct[id] = [scaffold, start, end, strand]
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


#if __name__ == '__main__':
g3 = GFF("output/Bfra_R1V1.fa.all_no_seq.modified.gff")

sorted_ids = [_[0] for _ in sorted(g3.dct.items(),  key=lambda x: (x[1][0], x[1][1], x[1][2]))]

for id in sorted_ids:
    print_lines(list(unique_everseen(g3.id2lines[id])))
            
         

