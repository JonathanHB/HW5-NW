# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch
import numpy as np

def main():
    """
    This function should
    (1) Align all species to humans and print species in order of most similar to human BRD
    (2) Print all alignment scores between each species BRD2 and human BRD2
    """
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/tursiops_truncatus_BRD2.fa")

    # TODO Align all species to humans and print species in order of most similar to human BRD
    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix

    #nonhuman species, sequences, and alignments to the human sequence
    species = ["Gallus gallus", "Mus musculus", "Balaeniceps_rex", "tursiops_truncatus"]
    seqs_nh = [gg_seq, mm_seq, br_seq, tt_seq]
    alignments = []
    alnscores = []

    # using gap opening penalty of -10 and a gap extension penalty of -1 and BLOSUM62 matrix
    #align species to human
    for i in seqs_nh:
        alignments.append(NeedlemanWunsch(f"./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1).align(hs_seq, i))
        alnscores.append(alignments[-1][0])
        #print(alignments[-1][1])
        #print(alignments[-1][2])

    #sort by similarity and print results
    asort = np.flip(np.argsort(alnscores))

    print("Species by BRD2 similarity to human BRD2")
    for i in asort:
        print(f"{species[i]}: {alnscores[i]}")



if __name__ == "__main__":
    main()
