analysis_GC<-function(seq){
  seq<-DNAString(seq)
  GC<-sum(letterFrequency(seq,c("G","C")))/nchar(seq)
  {  
    if(GC >= 0.35 &GC <= 0.65) {
      return("Normal")
    }
    else if (GC >= 0.25 & GC < 0.35) {
      return("Low")
    }
    else if (GC > 0.65 & GC <= 0.80) {
      return("High")
    }
    else{
      return("FALSE")
    }
  }
}
