#反向重复
#含有反向重复>100 bp的序列，此区域弃用
Get_dot_region2 <- function(KO_region) {
  if(KO_region$start<900 | KO_region$end+900>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 900, KO_region$end+500 + 900)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 900, KO_region$end + 900)
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
    # for (j in 1:length(pos)) {
    #   if (pos[j] != -1) {
    #     analysis_pos <- append(analysis_pos,pos[j])
    #   }
    # }
  }
  return(FALSE)
  # analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  # analysis_frame <- data.frame()
  # for (i in 1:length(analysis_pos)) {
  #   analysis_frame[i, 1] <- analysis_pos[i] 
  #   analysis_frame[i, 2] <- analysis_pos[i] + 9
  # }
  # times <- 0
  # i <- 1
  # j <- 1
  # dot.rev_frame <- data.frame()
  # repeat {
  #   if (analysis_frame[i + 1, 1] == analysis_frame[i, 1]+1 ) {
  #     times <- times + 1
  #     if (times == 1) {
  #       dot.rev_frame[j, 1] <- analysis_frame[i, 1]
  #     }
  #     i <- i + 1
  #     if (i == length(analysis_frame[, 1])) {
  #       dot.rev_frame[j, 2] <- analysis_frame[i, 2]
  #       break
  #     }
  #   }
  #   else{
  #     if (times == 0) {
  #       dot.rev_frame[j, 1] <- analysis_frame[i, 1]
  #     }
  #     dot.rev_frame[j, 2] <- analysis_frame[i , 2]
  #     i <- i + 1
  #     if (i == length(analysis_frame[, 1])) {
  #       j <- j + 1
  #       dot.rev_frame[j, 1] <- analysis_frame[i, 1]
  #       dot.rev_frame[j, 2] <- analysis_frame[i, 2]
  #       break
  #     }
  #     times <- 0
  #     j <- j + 1
  #   }
  # }
  # dot.rev_frame[, 1] <-dot.rev_frame[, 1] + KO_region$start - 900 - 1
  # dot.rev_frame[, 2] <-dot.rev_frame[, 2] + KO_region$start - 900 - 1
  # names(dot.rev_frame) <- c("start", "end")
  # length<-dot.rev_frame$end-dot.rev_frame$start+1
  # dot.rev_frame<-cbind(dot.rev_frame,length)
  # #print(dot.rev_frame)
  # dot.rev_frame_100<-dot.rev_frame[which(dot.rev_frame$length>100),]
  # if(nrow(dot.rev_frame_100)!=0){
  #   return(TRUE)
  #   break
  # }
  # return(FALSE)
}
#待定
# dot.rev_frame_100[which(dot.rev_frame_100$start),]
# 
# #含有反向重复>100 bp的序列，且重复区域在靶位点邻近500 bp内，此区域弃用
# dot.rev_frame_100<-dot.rev_frame[which(dot.rev_frame$length>100),]
# if(nrow(dot.frame_100)!=0){
#   for(i in 1:nrow(dot.rev_frame_100)){
#     if(dot.rev_frame_100[i,]$start<KO_region$end & dot.rev_frame_100[i,]$end>KO_region$start-500){
#       return(TRUE)
#       break
#     }
#     else if(dot.rev_frame_100[i,]$end>KO_region$start & dot.rev_frame_100[i,]$start<KO_region$end+500){
#       return(TRUE)
#       break
#     }
#   }
# }

#有>=25 bp的反向重复序列，间隔<50 bp，区域弃用
Get_dot_region3 <- function(KO_region) {
  if(KO_region$start<900 | KO_region$end+900>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 900, KO_region$end+500 + 900)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 900, KO_region$end + 900)
  }
  
  analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 24
  for (i in 1:len) {
    pattern <- substring(analysis_seq_rev, i, i + 24)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != i) {
        if(pos[j] -i <=70 &pos[j] -i >=20 ){
          return(TRUE)
        }
      }
    }
    # if (length(pos) > 1) {
    #   for (j in 1:length(pos)) {
    #     if (end(pos[j]) - i <= 90 & end(pos[j]) - i >= 39) {
    #       print("间隔<50 bp,>=25 bp的反向")
    #       return(TRUE)
    #       break
    #     }
    #   }
    # }
  }
  return(FALSE)
}



