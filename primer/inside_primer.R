
inside_primer<-function(Gene, PCR_seq, Gene_rev,result1,filepath){
  if(nchar(PCR_seq)>=30 & nchar(PCR_seq)<=200){
    result_in <- primer_design3_in(Gene, PCR_seq, Gene_rev, result1,filepath)
  }
  if(nchar(PCR_seq)>200){
    result_in <- primer_design_in(Gene, PCR_seq, Gene_rev, result1,filepath)
  }
  return(result_in)
}


inside_primer2<-function(Gene,Gene_rev,result1,Exon_region,filepath){
  {
    if (Gene_rev) {
      F1.pos <-matchPattern(reverse(DNAString(result1$primer1)), subject = Gene)
      R1.pos <-matchPattern(complement(DNAString(result1$primer2)), subject = Gene)
      {
        if (end(F1.pos) - 999 < Exon_region$start) {
          PCR_seq1 <- substring(Gene, end(F1.pos) - 999, Exon_region$end)
          Exon_pos<-data.frame(start=Exon_region$end-Exon_region$end+1,
                               end=Exon_region$end-Exon_region$start+1)
        }
        else{
          PCR_seq1 <-substring(Gene, Exon_region$start, Exon_region$end)
          Exon_pos<-data.frame(start=1,end=1000-(end(F1.pos)-Exon_region$end))
        }
      }
      
      {
        if (start(R1.pos) + 999 > Exon_region$end) {
          PCR_seq2 <-substring(Gene, Exon_region$start, start(R1.pos) +999)
          Exon_pos2<-data.frame(start=(start(R1.pos) + 999)-Exon_region$end+1,
                                end=(start(R1.pos) + 999)-Exon_region$start+1)
        }
        else{
          PCR_seq2 <-substring(Gene, Exon_region$start, Exon_region$end)
          Exon_pos2<-data.frame(start=(Exon_region$end-start(R1.pos))-999,end=2000)
        }
      }
    }
    else{
      F1.pos <-matchPattern(DNAString(result1$primer1), subject = Gene)
      R1.pos <-matchPattern(reverseComplement(DNAString(result1$primer2)), subject = Gene)
      {
        if(start(F1.pos)+999>Exon_region$end){
          PCR_seq1 <- substring(Gene,Exon_region$start,start(F1.pos)+999)
          Exon_pos<-data.frame(start=1,end=Exon_region$end-Exon_region$start+1)
        }
        else{
          PCR_seq1 <- substring(Gene,Exon_region$start,Exon_region$end)
          Exon_pos<-data.frame(start=1,end=1000-(Exon_region$start-start(F1.pos)))
        }
      }
      
      {
        if(end(R1.pos)-999<Exon_region$start){
          PCR_seq2 <- substring(Gene,end(R1.pos)-999,Exon_region$end)
          Exon_pos2<-data.frame(start=Exon_region$start-(end(R1.pos)-999)+1,
                                end=Exon_region$end-(end(R1.pos)-999)+1)
        }
        else{
          PCR_seq2 <- substring(Gene,Exon_region$start,Exon_region$end)
          Exon_pos2<-data.frame(start=(end(R1.pos)-Exon_region$start)-999,end=2000)
        }
      }
    }
  }
  result_in<-primer_design2_in(Gene,PCR_seq1,PCR_seq2,Gene_rev,result1,Exon_pos,Exon_pos2,filepath)
  return(result_in)
}



inside_primer3<-function(Gene,Gene_rev,result1,Exon_region,filepath){
  {
    if (Gene_rev) {
      F1.pos <-matchPattern(reverse(DNAString(result1$primer1)), subject = Gene)
      R1.pos <-matchPattern(complement(DNAString(result1$primer2)), subject = Gene)
      {
        if (end(F1.pos) - 999 < Exon_region[1, ]$start) {
          PCR_seq1 <- substring(Gene, end(F1.pos) - 999, Exon_region[1, ]$end)
          Exon_pos<-data.frame(start=Exon_region[1, ]$end-Exon_region[which(Exon_region$end>=(end(F1.pos) - 999)),]$end+1,
                               end=Exon_region[1, ]$end-Exon_region[which(Exon_region$end>=(end(F1.pos) - 999)),]$start+1)
        }
        else{
          PCR_seq1 <-substring(Gene, Exon_region[1, ]$start, Exon_region[1, ]$end)
          Exon_pos<-data.frame(start=1,end=1000-(end(F1.pos)-Exon_region[1, ]$end))
        }
      }
      
      {
        if (start(R1.pos) + 999 > Exon_region[nrow(Exon_region), ]$end) {
          PCR_seq2 <-substring(Gene, Exon_region[nrow(Exon_region), ]$start, start(R1.pos) +999)
          Exon_pos2<-data.frame(start=(start(R1.pos) + 999)-Exon_region[which(Exon_region$start<=(start(R1.pos) + 999)),]$end+1,
                                end=(start(R1.pos) + 999)-Exon_region[which(Exon_region$start<=(start(R1.pos) + 999)),]$start+1)
        }
        else{
          PCR_seq2 <-substring(Gene, Exon_region[nrow(Exon_region), ]$start, Exon_region[nrow(Exon_region), ]$end)
          Exon_pos2<-data.frame(start=(Exon_region[nrow(Exon_region), ]$end-start(R1.pos))-999,end=2000)
        }
      }
    }
    else{
      F1.pos <-matchPattern(DNAString(result1$primer1), subject = Gene)
      R1.pos <-matchPattern(reverseComplement(DNAString(result1$primer2)), subject = Gene)
      {
        if(start(F1.pos)+999>Exon_region[1,]$end){
          PCR_seq1 <- substring(Gene,Exon_region[1,]$start,start(F1.pos)+999)
          Exon_pos<-data.frame(start=Exon_region[which(Exon_region$start<=(start(F1.pos)+999)),]$start-Exon_region[1,]$start+1,
                               end=Exon_region[which(Exon_region$start<=(start(F1.pos)+999)),]$end-Exon_region[1,]$start+1)
        }
        else{
          PCR_seq1 <- substring(Gene,Exon_region[1,]$start,Exon_region[1,]$end)
          Exon_pos<-data.frame(start=1,end=1000-(Exon_region[1,]$start-start(F1.pos)))
        }
      }
      {
        if(end(R1.pos)-999<Exon_region[nrow(Exon_region),]$start){
          PCR_seq2 <- substring(Gene,end(R1.pos)-999,Exon_region[nrow(Exon_region),]$end)
          Exon_pos2<-data.frame(start=Exon_region[which(Exon_region$end>=(end(R1.pos)-999)),]$start-(end(R1.pos)-999)+1,
                                end=Exon_region[which(Exon_region$end>=(end(R1.pos)-999)),]$end-(end(R1.pos)-999)+1)
        }
        else{
          PCR_seq2 <- substring(Gene,Exon_region[nrow(Exon_region),]$start,Exon_region[nrow(Exon_region),]$end)
          Exon_pos2<-data.frame(start=(end(R1.pos)-Exon_region[nrow(Exon_region),]$start)-999,end=2000)
        }
      }
    }
  }
  result_in<-primer_design2_in(Gene,PCR_seq1,PCR_seq2,Gene_rev,result1,Exon_pos,Exon_pos2,filepath)
  return(result_in)
}
