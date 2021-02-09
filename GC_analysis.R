#平均GC含量
GC_analysis1<-function(KO_region){
  if(KO_region$start<800 | KO_region$end+800>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 800, KO_region$end+500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_seq<-DNAString(analysis_seq)
  GC<-sum(letterFrequency(analysis_seq,c("G","C")))/nchar(analysis_seq)
  if(GC > 0.7| GC < 0.25){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}


#GC含量分析,连续50bp大于80%或小于25%的区域
GC_analysis2 <- function(KO_region) {
  k <- 1
  GC_avoid_region <- data.frame(start = numeric(), end = numeric())
  if(KO_region$start<800 | KO_region$end+800>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 800, KO_region$end+500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30
  for (i in 1:(length(GC) - 49)){
    j <- i + 49
    GC_50 <- GC[i:j]
    if (all(GC_50 >= 0.80) | all(GC_50 <= 0.25)) {
      GC_avoid_region[k, ]$start <- i
      GC_avoid_region[k, ]$end <- j
      k <- k + 1
    }
  }
  if (nrow(GC_avoid_region) != 0) {
    times <- 0
    i <- 1
    j <- 1
    avoid_district <- data.frame()
    repeat {
      #待修改
      if(nrow(GC_avoid_region)==1){
        avoid_district[j, 1] <- GC_avoid_region[i, 1]
        avoid_district[j, 2] <- GC_avoid_region[i, 2]
        break
      }
      if (GC_avoid_region[i + 1, 1] >= GC_avoid_region[i, 1] &
          GC_avoid_region[i + 1, 1] <= GC_avoid_region[i, 2]) {
        times <- times + 1
        if (times == 1) {
          avoid_district[j, 1] <- GC_avoid_region[i, 1]
        }
        i <- i + 1
        if (i == length(GC_avoid_region[, 1])) {
          avoid_district[j, 2] <- GC_avoid_region[i, 2]
          break
        }
      }
      else{
        if (times == 0) {
          avoid_district[j, 1] <- GC_avoid_region[i, 1]
        }
        avoid_district[j, 2] <- GC_avoid_region[i , 2]
        i <- i + 1
        if (i == length(GC_avoid_region[, 1])) {
          avoid_district[j, 2] <- GC_avoid_region[i-1, 2]
          j<-j+1
          avoid_district[j,1] <-GC_avoid_region[i,1]
          avoid_district[j,2]<-GC_avoid_region[i,2]
          break
        }
        times <- 0
        j <- j + 1
      }
    }
    avoid_district[, 1] <-
      avoid_district[, 1] + KO_region$start - 800 - 1
    avoid_district[, 2] <-
      avoid_district[, 2] + KO_region$start - 800 - 1
    names(avoid_district) <- c("start", "end")
    return(avoid_district)
  }
  else{
    return(FALSE)
  }
}

#平均GC含量
GC_analysis3<-function(KO_region){
  if(KO_region$start<800 | KO_region$end+800>nchar(Gene)){
    analysis_seq <-
      substring(Gene2, KO_region$start+500 - 800, KO_region$end+500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_seq<-DNAString(analysis_seq)
  GC<-sum(letterFrequency(analysis_seq,c("G","C")))/nchar(analysis_seq)
  return(GC)
}
