from Bio.Seq import MutableSeq
#原始序列
dna_seq = MutableSeq("TAGCTACAGTGAAATCTCGA")
print(dna_seq)

print()
#反向
dna_seq_rev = dna_seq[::-1]
print(dna_seq_rev)

#互补
dna_seq.complement()
print(dna_seq)

#反向互补
dna_seq_rev.complement()
print(dna_seq_rev)





#翻译
def dna_trans_protein(dnaSeq):
    codonTable = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    proteinSeq = ""
    for codonStart in range(0, len(dnaSeq), 3):
        codon = dnaSeq[codonStart:codonStart + 3]
        if codon in codonTable:
            proteinSeq += codonTable[codon]
    return proteinSeq

codon=dna_trans_protein("ATGGTTGTCTATTAACTTGTTCAAAAAAGTATCAGGAGTTGTCAAGGCAGAGAAGAGAGTGTTTGC")
print(codon)
