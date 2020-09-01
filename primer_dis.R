#计算片段敲除的靶位点距离
pri1<-read.csv("primer_result1.csv",header = TRUE)
pri2<-read.csv("primer_result2.csv",header = TRUE)
primer3<-pri1[1,2]
primer4<-pri2[1,1]
{
  #反向基因
  if(Gene_rev) {
    if (min(gRNA$start) < 750 | max(gRNA$end) + 750 > nchar(Gene)) {
      primer3.pos <-
        matchPattern(complement(DNAString(primer3)), subject = Gene2)
      primer4.pos <-
        matchPattern(reverse(DNAString(primer4)), subject = Gene2)
    }
    
    else{
      primer3.pos <-
        matchPattern(complement(DNAString(primer3)), subject = Gene)
      primer4.pos <-
        matchPattern(reverse(DNAString(primer4)), subject = Gene)
    }
    WT2.distance <-
      max(end(primer1.pos), end(primer3.pos)) - min(start(primer1.pos), start(primer3.pos)) + 1
    WT3.distance <-
      max(end(primer2.pos), end(primer4.pos)) - min(start(primer2.pos), start(primer4.pos)) + 1
  }
  
  #正向基因
  else{
    if (min(gRNA$start) < 750 | max(gRNA$end) + 750 > nchar(Gene)) {
      primer3.pos <-
        matchPattern(reverseComplement(DNAString(primer3)), subject = Gene2)
      primer4.pos <-
        matchPattern(DNAString(primer4), subject = Gene2)
    }
    
    else{
      primer3.pos <-
        matchPattern(reverseComplement(DNAString(primer3)), subject = Gene)
      primer4.pos <-
        matchPattern(DNAString(primer4), subject = Gene)
    }
    WT2.distance <-
      max(end(primer1.pos), end(primer3.pos)) - min(start(primer1.pos), start(primer3.pos)) + 1
    WT3.distance <-
      max(end(primer2.pos), end(primer4.pos)) - min(start(primer2.pos), start(primer4.pos)) + 1
  }
}

primer.data<-data.frame(F1=character(),R1=character(),R2=character(),F2=character(),F1R1=numeric(),F1R2=numeric(),F2R1=numeric())
primer.data[1,1]<-primer1
primer.data[1,2]<-primer2
primer.data[1,5]<-WT.distance
if(exists("primer3")){
  primer.data[1,3]<-primer3
  primer.data[1,6]<-WT2.distance
}
if(exists("primer4")){
  primer.data[1,4]<-primer4
  primer.data[1,7]<-WT3.distance
}
write.csv(primer.data,paste0("C://Users//41518//Desktop//",term,"靶位点.csv"),row.names = FALSE)



