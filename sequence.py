# 检查gRNA的碱基
def check_nucleotide(sequence):
    sequence_nucleotide = ['A', 'C', 'G', 'T', 'U', 'a', 'c', 'g', 't', 'u']
    return all(nuc in sequence_nucleotide for nuc in sequence)


# gRNA序列标准化
def gRNA_normalize(gRNA):
    gRNA_normal =[]
    for g in gRNA:
        gRNA_normal.append(g.upper().replace('U','T'))
    return gRNA_normal

def DNA_rev_complement(DNA_sequence):
    reverse = DNA_sequence[::-1]
    complement=''
    for nt in reverse:
        if nt =='A':
            complement += 'T'
        if nt =='T':
            complement += 'A'
        if nt =='G':
            complement += 'C'
        if nt =='C':
            complement += 'G'
    return complement