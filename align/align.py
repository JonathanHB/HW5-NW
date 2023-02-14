# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub

    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA
        self._seqB = seqB
        
        # Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing

        n = len(seqA)+1
        m = len(seqB)+1

        self._align_matrix = np.zeros([n, m])
        # 0 if gap does not start, self.gap_open if it does (except on the first row and column, which already include the gap opening penalty)
        self._gapA_matrix = np.ones([n, m])*(self.gap_open + self.gap_extend)
        self._gapB_matrix = np.ones([n, m])*(self.gap_open + self.gap_extend)
        self._back_A = np.zeros([n,m])
        self._back_B = np.zeros([n,m])

        #tree = [[[0,0] for i in range(m)] for i in range(n)] #backwards pointers for debugging

        #global alignment

        #initialize corner values
        self._gapA_matrix[0, 0] = 0
        self._gapB_matrix[0, 0] = 0

        #initialize edge values
        for i in range(1, n):
            self._align_matrix[i, 0] = self.gap_extend*(i) + self.gap_open
            self._gapA_matrix[i, 0] = 0
            self._gapB_matrix[i, 0] = 0
            self._back_A[i, 0] = 1

        for i in range(1, m):
            self._align_matrix[0, i] = self.gap_extend*(i) + self.gap_open
            self._gapA_matrix[0, i] = 0
            self._gapB_matrix[0, i] = 0
            self._back_B[0, i] = 1

        #fill in matrix
        for i in range(1, n):
            for j in range(1, m):
                #fill in the score in the current cell based on the best of the above, left, and above left cells
                diagscore = self._align_matrix[i-1, j-1] + self.sub_dict[(self._seqA[i-1], self._seqB[j-1])]
                vertscore = self._align_matrix[i-1, j] + self._gapB_matrix[i-1, j]
                horiscore = self._align_matrix[i, j-1] + self._gapA_matrix[i, j-1]

                self._align_matrix[i, j] = max(diagscore, vertscore, horiscore)

                #debugging code---------------------------
                #print(dvhscores)
                #print(np.argsort(dvhscores))
                # if j == 1:
                #     print(f"{i}, {j}")
                #     print("diagonal: " + str(diagscore))
                #     print("vertical: " + str(vertscore))
                #     print("horizontal: " + str(horiscore))
                #-----------------------------------------

                #fill in gap and backtrace matrices

                #get the direction of the best upstream square to tell us which backtrace matrix square to fill
                dvhscores = [diagscore, vertscore, horiscore]
                best_direction = np.argsort(dvhscores)[-1]

                #a and b may be [consistently] reversed here
                if best_direction == 1:
                    self._gapB_matrix[i, j] = self.gap_extend
                    self._back_B[i, j] = 1
                    #tree[i][j] = [-1, 0]

                elif best_direction == 2:
                    self._gapA_matrix[i, j] = self.gap_extend
                    self._back_A[i, j] = 1
                    #tree[i][j] = [0, -1]

                #else:
                    #tree[i][j] = [-1, -1]

        #debugging
        #for k in tree:
        #    print(k)

        return self._backtrace()


    def _backtrace(self) -> Tuple[float, str, str]:

        self.alignment_score = self._align_matrix[-1,-1]

        n = len(self._seqA)+1
        m = len(self._seqB)+1

        self.seqA_align = []
        self.seqB_align = []

        #initial position
        x = n-1 #A
        y = m-1 #B

        for i in range(max(m,n)-1):

            #print(f"{x}, {y}")
            #if this pair is a gap in B, move up the column
            if self._back_B[x,y] == 1:
                #print(f"A {x}, {y}")
                self.seqA_align.append(self._seqA[x-1])
                self.seqB_align.append("-")
                x-=1

            #if this pair is a gap in A, move left along the row
            elif self._back_A[x,y] == 1:
                #print(f"B {x}, {y}")
                self.seqB_align.append(self._seqB[y-1])
                self.seqA_align.append("-")
                y-=1

            #if this pair is a (mis)match, move diagonally up and left
            else:
                self.seqA_align.append(self._seqA[x-1])
                self.seqB_align.append(self._seqB[y-1])
                x-=1
                y-=1

        #join the lists and fix [flip] the order
        self.seqA_align = "".join(self.seqA_align)[::-1]
        self.seqB_align = "".join(self.seqB_align)[::-1]

        # print(self._seqA)
        # print(self._seqB)
        #
        # print(self.seqA_align)
        # print(self.seqB_align)
        #
        # print(self._back_A)
        # print(self._back_B)
        #
        # print(list(self._gapA_matrix))
        # print(list(self._gapB_matrix))
        # #
        # print(list(self._align_matrix))

        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header

#
# #testing code
# seq1, _ = read_fasta("../data/test_seq1.fa")
# seq2, _ = read_fasta("../data/test_seq2.fa")
# nw = NeedlemanWunsch("../substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
# c=nw.align(seq1,seq2)
# print(nw._align_matrix)
# print(c)

# # print(seq1)
#
# #seq3xseq4
# #M[3-gap]QLIR[H/R]P
# print(sum([5,-13, 5, 4, 4, 5, 0, 7]))