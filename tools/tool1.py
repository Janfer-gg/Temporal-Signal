import rpy2.robjects as robjects
robjects.r("source('dot_plot.R')")

filepath ='C://Users//41518//Desktop'
seq1 = 'TTACTTGGGTGTGGAGGCATGTGCCTATAATCCAAGGTACTCAGGAGGCTGAGGCATGAGAATTGCTTGAACCCGGGAGGTAGAGGTTGCAGTGAGCTGAGATTATGCCACTGCACTCCAGCCTGGATGACAGAGTGAGGCTCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAGAAGTAGAATCATATGGTAATTCTGTGTTTAACTTTTTGAGGAACGCCAGACTCTTCCAAACCACTGCCCTGTGTTACATTCCCAGTAGCAATGTGTTCAGTTGCTTCACATCTTTCCCAACGCTTTTGATTTTCCTCTTTATTTTTCATTGAATCTAATCCTAGTGGGTATGAAGTATCTTACTGTGTTTTGATTTGCATTTCTCTAATGAATAATGATGTTGAGCATATTTTCATGTGCTTATGGGTCATTTGTAAACTTTTAGAGAAATGTCTGTTCCAATCTTTTATCCGTTTTAAAAATTGGATTGTCTTTTTGTTGTTGAGGTATAAGATTTCTTTCTATGTTCTAGATATTAGATGCTTAACAGATAAGCAGGAGAATCGCCTGAACCCGGGAGGCAGAGGCTGCAGTGAGCGGAGATTGTACCGCTGCCCTCCAGCCTGGGCGACAGAGGAGGGAGACTGTCTGAAACAAACAACAACAACAACAACAAGATGAAAATAATACGGCAAATACTGTTCTTGACTTTTACTGTATCGAGTGAGAACTGTAAGGCTTAGAACAGCATGTGGTGGACTGTAATGATTCTATTTGTTTGCTTCTTATCCCTTTGCTTCACCCAGAAAAGGATTATAGCAGTCTTTGTGACAAGCAACCGATAGGAAGACGTCTCTTCAGGCAGTTCTGTGATACCAAACCCACTCTAAAGAGGCACATTGAATTCTTGGATGCAGTGGTGAGCAGTTTATCTCCATATTGAGCAACCACCCAATCTTATGCTTTTGAAAATGTAAAAACTTGGCCAGGCGCAGTAGCTCATGCTGTAATCCCAGCATTTTGGAAGGCTGAGGCGGGCGGATCATGAGGTCAGGAGATCGAGACCATCCTGACCAACACAGTGAAACCCGTCTCTACTAAAATACAAAAAGTTAGCCGGGCGTGGTGGCGCGTGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGGGAATCACTTGAACCAGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCACCACTGCCTTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAGAAAATGTAAAAACTAGACCTGGGCAATGATGTAGGAGTGGAGGGACTGGTTGTCTTGGTGTCATATTATCTTATTAGGACAGAAATCCCTTTAGTCTGAAGTGGTATTTTGTGTAGAATTACGTCTCAAGTGTTGGAGAATCACATGTAGTCATTGAATGACTTTGAAACTTGAGGCTTGAATTGTATGAAGATGACAAACTAATTAAGATGATAATGCTTTTTTTATTTTTTATTTTTTTGAAATGGAGTCTCTCTCTGTCTCCCAGGCTGGACTGCAGTGGCATGATCTTGGCTTACTGCAACCTCCGCCTCCTGGGTTCAAGCAATTCTTCTGCCTCAGCCTTGTGATCCACCTGCCTTGGCCCCCCAAAGTGCTGGGATTACAGGCATAAGCCACTGTGCCCGGCTCTAAGATGATAGTGCTTGTTCGGCTACCAGAAAAGTGTC'
seq2 = 'TTACTTGGGTGTGGAGGCATGTGCCTATAATCCAAGGTACTCAGGAGGCTGAGGCATGAGAATTGCTTGAACCCGGGAGGTAGAGGTTGCAGTGAGCTGAGATTATGCCACTGCACTCCAGCCTGGATGACAGAGTGAGGCTCTGTCTCAAAAAAAAAAAAAAAAAAAAAAAAGAAGTAGAATCATATGGTAATTCTGTGTTTAACTTTTTGAGGAACGCCAGACTCTTCCAAACCACTGCCCTGTGTTACATTCCCAGTAGCAATGTGTTCAGTTGCTTCACATCTTTCCCAACGCTTTTGATTTTCCTCTTTATTTTTCATTGAATCTAATCCTAGTGGGTATGAAGTATCTTACTGTGTTTTGATTTGCATTTCTCTAATGAATAATGATGTTGAGCATATTTTCATGTGCTTATGGGTCATTTGTAAACTTTTAGAGAAATGTCTGTTCCAATCTTTTATCCGTTTTAAAAATTGGATTGTCTTTTTGTTGTTGAGGTATAAGATTTCTTTCTATGTTCTAGATATTAGATGCTTAACAGATAAGCAGGAGAATCGCCTGAACCCGGGAGGCAGAGGCTGCAGTGAGCGGAGATTGTACCGCTGCCCTCCAGCCTGGGCGACAGAGGAGGGAGACTGTCTGAAACAAACAACAACAACAACAACAAGATGAAAATAATACGGCAAATACTGTTCTTGACTTTTACTGTATCGAGTGAGAACTGTAAGGCTTAGAACAGCATGTGGTGGACTGTAATGATTCTATTTGTTTGCTTCTTATCCCTTTGCTTCACCCAGAAAAGGATTATAGCAGTCTTTGTGACAAGCAACCGATAGGAAGACGTCTCTTCAGGCAGTTCTGTGATACCAAACCCACTCTAAAGAGGCACATTGAATTCTTGGATGCAGTGGTGAGCAGTTTATCTCCATATTGAGCAACCACCCAATCTTATGCTTTTGAAAATGTAAAAACTTGGCCAGGCGCAGTAGCTCATGCTGTAATCCCAGCATTTTGGAAGGCTGAGGCGGGCGGATCATGAGGTCAGGAGATCGAGACCATCCTGACCAACACAGTGAAACCCGTCTCTACTAAAATACAAAAAGTTAGCCGGGCGTGGTGGCGCGTGCCTGTAGTCCCAGCTACTTGGGAGGCTGAGGCAGGGGAATCACTTGAACCAGGAGGCGGAGGTTGCAGTGAGCCGAGATCGCACCACTGCCTTCCAGCCTGGGCGACAGAGCAAGACTCCGTCTCAAAAAAAAAGAAAATGTAAAAACTAGACCTGGGCAATGATGTAGGAGTGGAGGGACTGGTTGTCTTGGTGTCATATTATCTTATTAGGACAGAAATCCCTTTAGTCTGAAGTGGTATTTTGTGTAGAATTACGTCTCAAGTGTTGGAGAATCACATGTAGTCATTGAATGACTTTGAAACTTGAGGCTTGAATTGTATGAAGATGACAAACTAATTAAGATGATAATGCTTTTTTTATTTTTTATTTTTTTGAAATGGAGTCTCTCTCTGTCTCCCAGGCTGGACTGCAGTGGCATGATCTTGGCTTACTGCAACCTCCGCCTCCTGGGTTCAAGCAATTCTTCTGCCTCAGCCTTGTGATCCACCTGCCTTGGCCCCCCAAAGTGCTGGGATTACAGGCATAAGCCACTGTGCCCGGCTCTAAGATGATAGTGCTTGTTCGGCTACCAGAAAAGTGTC'
size =10
mis = 0
robjects.r['dot_plot'](seq1,seq2,size,mis,filepath)