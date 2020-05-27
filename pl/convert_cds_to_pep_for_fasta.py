from pyfaidx import Fasta
from Bio.Seq import Seq
import sys

# Usage: python -i convert_cds_to_protein_for_fasta.py /projects/b1080/hy/ribo/ribo/yeast/annotation/candidateORF.fa

def cds2pep(cds):
    return(str(Seq(cds).translate()))


def cds_fasta_file_to_pep_fasta_file(cdsFastaPath):
    cdsFasta = Fasta(cdsFastaPath,  as_raw = True)

    suffix = "fasta"
    if cdsFastaPath.endswith("fa"):
        suffix = "fa"
    elif cdsFastaPath.endswith("fas"):
        suffix = "fas"

    with open(cdsFastaPath.rstrip("." + suffix) + ".pep." + suffix, "w") as w:
        for name in cdsFasta.keys():
            cds = str(cdsFasta[name])
            pep = cds2pep(cds)
            w.write(">" + name + "\n" + pep + "\n") 

if __name__ == "__main__":
    cdsFastaPath = sys.argv[1]
    print("converting this cds fasta file: ", cdsFastaPath)
    cds_fasta_file_to_pep_fasta_file(cdsFastaPath)

