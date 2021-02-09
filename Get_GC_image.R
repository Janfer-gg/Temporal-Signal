Get_GC_image <- function(KO_region) {
  if (KO_region$start < 800 | KO_region$end + 800 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$start + 500 - 800, KO_region$end + 500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30
  GC.data <- data.frame(x = c(1:(nchar(analysis_seq) - 29)), y = GC)
  p4 <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())
  p4<- p4 + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p4)
}

#大片段敲除的上游GC分析
Get_GC_image2 <- function(KO_region) {
  if (KO_region$start < 500 | KO_region$end + 500 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$start + 500 - 800, KO_region$start+500)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$start)
  }
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30 
  GC.data <- data.frame(x = c(1:(nchar(analysis_seq) - 29)), y = GC)
  p4 <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())
  p4<- p4 + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p4)
}

#大片段敲除的下游GC分析
Get_GC_image3 <- function(KO_region) {
  if (KO_region$start < 500 | KO_region$end + 500 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$end+500, KO_region$end+500+800)
  }
  else{
    analysis_seq <- substring(Gene, KO_region$end, KO_region$end+800)
  }
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30 
  GC.data <- data.frame(x = c(1:(nchar(analysis_seq) - 29)), y = GC)
  p4 <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())
  p4<- p4 + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p4)
}

#大片段敲除的上游GC分析
Get_GC_image4 <- function(KO_region) {
  analysis_seq <-substring(Gene3, KO_region$start - 800, KO_region$start)
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30 
  GC.data <- data.frame(x = c(1:(nchar(analysis_seq) - 29)), y = GC)
  p4 <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())
  p4<- p4 + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p4)
}

#大片段敲除的下游GC分析
Get_GC_image5 <- function(KO_region) {
  analysis_seq <- substring(Gene3, KO_region$end, KO_region$end+800)
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30 
  GC.data <- data.frame(x = c(1:(nchar(analysis_seq) - 29)), y = GC)
  p4 <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())
  p4<- p4 + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p4)
}
