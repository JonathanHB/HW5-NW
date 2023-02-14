# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

lt = ""

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """

    gA_mat = np.array([[0., 0., 0., 0.],[0., -11., -1., -1.],[0., -11., -11., -1.],[0., -11., -11., -11.],[0., -11., -11., -11.]])
    gB_mat = np.array([[0., 0., 0., 0.],[0., -11., -11., -11.],[0., -1., -11., -11.],[0., -1., -11., -11.],[0., -1., -11., -11.]])
    aln_mat = np.array([[0., -11., -12., -13.],[-11., 5., -6., -7.],[-12., -6., 4., -7.],[-13., -7., -1., 5.],[-14., -8., -6., 4.]])

    seq1, _ = read_fasta(f"{lt}./data/test_seq1.fa")
    seq2, _ = read_fasta(f"{lt}./data/test_seq2.fa")

    nw = NeedlemanWunsch(f"{lt}./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    nw.align(seq1, seq2)

    assert (nw._gapA_matrix == gA_mat).all(), "mismatched gap_A matrix"
    assert (nw._gapB_matrix == gB_mat).all(), "mismatched gap_B matrix"

    assert (nw._align_matrix == aln_mat).all(), "mismatched alignment matrix"


def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """

    aln_score = 17.0
    aln_seqA = 'MAVHQLIRRP'
    aln_seqB = 'M---QLIRHP'

    seq3, _ = read_fasta(f"{lt}./data/test_seq3.fa")
    seq4, _ = read_fasta(f"{lt}./data/test_seq4.fa")

    nw = NeedlemanWunsch(f"{lt}./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    alignment = nw.align(seq3, seq4)

    assert alignment[0] == aln_score, "wrong alignment score"
    assert alignment[1] == aln_seqA, "wrong sequence B alignment"
    assert alignment[2] == aln_seqB, "wrong sequence B alignment"


