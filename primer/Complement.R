complement_analysis<-function(table){
  for(i in 1:nrow(table)){
    analysis_seq<-DNAString(table[i,]$PCR)
    seq_com<-reverseComplement(analysis_seq)
    for (j in 1:(length(seq_com)-5)) {
      n<-matchPattern(substring(seq_com,j,j+5),analysis_seq,max.mismatch = 1,min.mismatch = 0)
      if(length(n)!=0){
        print(i)
        print(j)
      }
    }
  }
}

s1<-DNAString(s1)

hairpin<-findPalindromes(s1,min.armlength = 4,max.looplength = 25)

self_annealing
if(length(hairpin)!=0){
  
}

