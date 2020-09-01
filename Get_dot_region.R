#正向重复
Get_dot_region1 <- function(KO_region) {
  if(KO_region$start<900 | KO_region$end+900>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 900, KO_region$end+500 + 900)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 900, KO_region$end + 900)
  }
  
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != i) {
        analysis_pos <- append(analysis_pos, pos[j])
      }
    }
  }
  if(length(analysis_pos)==0){
    return(FALSE)
  }
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  analysis_frame <- data.frame()
  for (i in 1:length(analysis_pos)) {
    analysis_frame[i, 1] <- analysis_pos[i] 
    analysis_frame[i, 2] <- analysis_pos[i] + 9
  }
  times <- 0
  i <- 1
  j <- 1
  dot.frame <- data.frame()
  repeat {
    if(nrow(analysis_frame)<2){
      return(FALSE)
    }
    if (analysis_frame[i + 1, 1] == analysis_frame[i, 1]+1 ) {
      times <- times + 1
      if (times == 1) {
        dot.frame[j, 1] <- analysis_frame[i, 1]
      }
      i <- i + 1
      if (i == length(analysis_frame[, 1])) {
        dot.frame[j, 2] <- analysis_frame[i, 2]
        break
      }
    }
    else{
      if (times == 0) {
        dot.frame[j, 1] <- analysis_frame[i, 1]
      }
      dot.frame[j, 2] <- analysis_frame[i , 2]
      i <- i + 1
      if (i == length(analysis_frame[, 1])) {
        j <- j + 1
        dot.frame[j, 1] <- analysis_frame[i, 1]
        dot.frame[j, 2] <- analysis_frame[i, 2]
        break
      }

      times <- 0
      j <- j + 1
    }
  }
  dot.frame[, 1] <-dot.frame[, 1] + KO_region$start - 900- 1
  dot.frame[, 2] <-dot.frame[, 2] + KO_region$start - 900- 1
  names(dot.frame) <- c("start", "end")
  length<-dot.frame$end-dot.frame$start+1
  dot.frame<-cbind(dot.frame,length)
  #print(dot.frame)
  #含有直接重复>100 bp的序列，此区域弃用
  dot.frame_100<-dot.frame[which(dot.frame$length>100),]
  if(nrow(dot.frame_100)!=0){
    print("直接重复>100 bp的序列")
    return(TRUE)
    break
  }
  #含有连续5个以上，间隔小于50 bp的20 bp以上的重复序列，此区域弃用
  dot.frame_20<-dot.frame[which(dot.frame$length>20),]
  times <- 0
  i <- 1
  if(nrow(dot.frame_20 != 0)) {
    repeat {
      if (i == nrow(dot.frame_20)) {
        break
      }
      if (dot.frame_20[i + 1, 1] < dot.frame_20[i, 1] + 50) {
        times <- times + 1
        i <- i + 1
      }
      else{
        if (times >= 4) {
          print("连续5个以上，间隔小于50 bp的20 bp以上的重复序列")
          #连续5个以上，间隔小于50 bp的20 bp以上的重复序列
          return(TRUE)
          break
        }
        times <- 0
        i <- i + 1
      }
    }
  }
  return(FALSE)
}







