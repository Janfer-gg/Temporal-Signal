#反向重复
#含有反向重复>100 bp的序列，此区域弃用
Get_dot_region2 <- function(KO_region) {
  if(KO_region$start<800 | KO_region$end+800>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 800, KO_region$end+500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 99
  for (i in 1:len) {
    pattern <- substring(analysis_seq_rev, i, i + 99)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    if (length(pos) > 1) {
      return(TRUE)
      break
    }
  }
  return(FALSE)

}




