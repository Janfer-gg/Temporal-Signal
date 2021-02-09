Get_dot_region6 <- function(KO_region) {
  if(KO_region$start<300 | KO_region$end+300>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$end+500, KO_region$end+500 + 400)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$end, KO_region$end+400)
  }
  #出现连续10个及以上的A碱基重复，或连续15个及以上G\C\T单碱基重复
  # a1<-matchPattern("AAAAAAAAAA",subject = analysis_seq)
  a2<-matchPattern("GGGGGGGGGGGGGGG",subject = analysis_seq)
  a3<-matchPattern("TTTTTTTTTTTTTTT",subject = analysis_seq)
  a4<-matchPattern("ccccccccccccccc",subject = analysis_seq)
  if(length(a2)!=0 | length(a3)!=0 | length(a4)!=0){
    print("下游连续单碱基重复")
    return(TRUE)
    break
  }
  
  #70bp以内含有3个及以上的10-20bp连续重复
  # analysis_pos <- numeric()
  # len <- nchar(analysis_seq) - 9
  # for (i in 1:len) {
  #   pattern <- substring(analysis_seq, i, i + 9)
  #   pos <- gregexpr(pattern, analysis_seq)[[1]]
  #   if(length(pos)>=3){
  #     if(pos[3]-pos[1]<=60){
  #       print("70bp内3次重复")
  #       return(TRUE)
  #       break
  #     }
  #   }
  # }
  
  # analysis_seq<-DNAString(analysis_seq)
  # pos<-findPalindromes(analysis_seq, min.armlength = 20, max.looplength = 50, min.looplength = 0, max.mismatch = 0)
  # if(length(pos)>0){
  #   return(TRUE)
  # }
  
  #间隔≤50bp的20bp及以上的反向互补序列
  # len <- nchar(analysis_seq) - 19
  # for (i in 1:len) {
  #   pattern <- substring(analysis_seq, i, i + 19)
  #   pattern <- reverseComplement(DNAString(pattern))
  #   pos <- matchPattern(pattern, analysis_seq, max.mismatch = 4)
  #   if (length(pos) > 0) {
  #     for (j in 1:length(pos)) {
  #       if (end(pos[j]) - i <= 90 & end(pos[j]) - i >= 39) {
  #         print("回文")
  #         return(TRUE)
  #         break
  #       }
  #     }
  #   }
  # }
  
  #≥7bp的全GC序列重复
  # len2 <- nchar(analysis_seq) - 6
  # for(i in 1:len2){
  #   pattern2<-substring(analysis_seq,i,i+6)
  #   #7bp的全GC序列
  #   if(grepl("T",pattern2)==FALSE & grepl("A",pattern2)==FALSE){
  #     pos<-gregexpr(pattern2, analysis_seq)[[1]]
  #     if(length(pos)>=2){
  #       print("GC重复")
  #       return(TRUE)
  #     }
  #   }
  # }
  return(FALSE)
}

