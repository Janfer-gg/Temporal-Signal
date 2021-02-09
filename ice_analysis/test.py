# from Bio import pairwise2, Seq
# seq1 = Seq.Seq('AAATAA')
# seq2 = Seq.Seq('AAACAA')
#
# score_matrix = {("A", "A"): 2,("A", "T"): -1,("A", "C"): -1,("A", "G"): -1,
#                 ("T", "A"): -1,("T", "T"): 2,("T", "C"): -1,("T", "G"): -1,
#                 ("C", "A"): -1,("C", "T"): -1,("C", "C"): 2,("C", "G"): -1,
#                 ("G", "A"): -1,("G", "T"): -1,("G", "C"): -1,("G", "G"): 2}
#
# alignments = pairwise2.align.localms(seq1, seq2, 5, -4, -3, -1)
# print(alignments)
# for alignment in alignments:
#     print(pairwise2.format_alignment(*alignment))
from Bio.SubsMat import MatrixInfo as matlist

