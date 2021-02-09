WT.calculation<-function(result,Gene,Gene_rev){
  primer1<-result[1,][1]
  primer2<-result[2,][2]
  primer3<-result[2,][1]
  primer4<-result[1,][2]
  
  {
    #反向基因
    if(Gene_rev) {
      
      primer1.pos <-
        matchPattern(reverse(DNAString(primer1)), subject = Gene)
      primer2.pos <-
        matchPattern(complement(DNAString(primer2)), subject = Gene)
      primer3.pos <-
        matchPattern(reverse(DNAString(primer3)), subject = Gene)
      primer4.pos <-
        matchPattern(complement(DNAString(primer4)), subject = Gene)
      
      analysis_seq1<-substring(Gene,min(start(primer1.pos), start(primer2.pos)),max(end(primer1.pos), end(primer2.pos)))
      analysis_seq2<-substring(Gene,min(start(primer1.pos), start(primer4.pos)),max(end(primer1.pos), end(primer4.pos)))
      analysis_seq3<-substring(Gene,min(start(primer2.pos), start(primer3.pos)),max(end(primer2.pos), end(primer3.pos)))
      GC1<-analysis_GC(analysis_seq1)
      GC2<-analysis_GC(analysis_seq2)
      GC3<-analysis_GC(analysis_seq3)
      
      WT1.distance <- 
        max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) +1
      WT2.distance <-
        max(end(primer1.pos), end(primer4.pos)) - min(start(primer1.pos), start(primer4.pos)) + 1
      WT3.distance <-
        max(end(primer2.pos), end(primer3.pos)) - min(start(primer2.pos), start(primer3.pos)) + 1
    }
    
    #正向基因
    else{
      primer1.pos <-
        matchPattern(DNAString(primer1), subject = Gene)
      primer2.pos <-
        matchPattern(reverseComplement(DNAString(primer2)), subject = Gene)
      primer3.pos <-
        matchPattern(DNAString(primer3), subject = Gene)
      primer4.pos <-
        matchPattern(reverseComplement(DNAString(primer4)), subject = Gene)
      
      analysis_seq1<-substring(Gene,min(start(primer1.pos), start(primer2.pos)),max(end(primer1.pos), end(primer2.pos)))
      analysis_seq2<-substring(Gene,min(start(primer1.pos), start(primer4.pos)),max(end(primer1.pos), end(primer4.pos)))
      analysis_seq3<-substring(Gene,min(start(primer2.pos), start(primer3.pos)),max(end(primer2.pos), end(primer3.pos)))
      GC1<-analysis_GC(analysis_seq1)
      GC2<-analysis_GC(analysis_seq2)
      GC3<-analysis_GC(analysis_seq3)
      
      
      WT1.distance <-
        max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) + 1
      WT2.distance <-
        max(end(primer1.pos), end(primer4.pos)) - min(start(primer1.pos), start(primer4.pos)) + 1
      WT3.distance <-
        max(end(primer2.pos), end(primer3.pos)) - min(start(primer2.pos), start(primer3.pos)) + 1
    }
  }
  primer.data<-data.frame(F1=character(),R1=character(),R2=character(),F2=character(),WT_F1R1=numeric(),WT_F1R2=numeric(),WT_F2R1=numeric(),GC1=character(),GC2=character(),GC3=character())
  primer.data[1, 1] <- primer1
  primer.data[1, 2] <- primer2
  primer.data[1, 3] <- primer4
  primer.data[1, 4] <- primer3
  primer.data[1, 5] <- WT1.distance
  primer.data[1, 6] <- WT2.distance
  primer.data[1, 7] <- WT3.distance
  primer.data[1, 8] <- GC1
  primer.data[1, 9] <- GC2
  primer.data[1, 10] <- GC3
  
  return(primer.data)
}



