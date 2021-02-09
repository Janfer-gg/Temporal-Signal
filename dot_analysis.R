#有>=25 bp的反向重复序列，间隔<50 bp，发卡结构
dot_analysis1 <- function(KO_region) {
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
  len <- nchar(analysis_seq) - 24
  for (i in 1:len) {
    pattern <- substring(analysis_seq_rev, i, i + 24)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != i) {
        if(pos[j] -i <=75 &pos[j] -i >=24 ){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}


#70bp以内含有3个及以上的10-20bp连续重复
dot_analysis2 <- function(KO_region) {
  if(KO_region$start<800 | KO_region$end+800>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 800, KO_region$end+500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- matchPattern(pattern, analysis_seq)
    if(length(pos)>=3){
      for(j in 1:(length(pos)-2)){
        if(end(pos)[j+2]-end(pos)[j]<=60){
          return(TRUE)
        }
      }
    }
  }
  return(FALSE)
}








