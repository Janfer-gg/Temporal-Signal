Get_dot_plot<-function(KO_region){
  if (KO_region$start < 800 | KO_region$end + 800 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$start + 500 - 800, KO_region$end + 500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  plot.f<-data.frame(x=numeric(),y=numeric())
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j])
      plot.f<-rbind(plot.f,f)
    }
  }
  
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  
  plot.rev.f<-data.frame(x=numeric(),y=numeric())
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  #analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pattern <- DNAString(pattern)
    pattern <- reverseComplement(pattern)
    pattern <- as.character(pattern)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != -1) {
        f<-data.frame(x=i,y=pos[j])
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  p3 <- ggplot(data = plot.f, aes(x = x, y = y)) + geom_point(size=0.6,color="red") + theme_bw() +
    theme(panel.grid = element_blank())+ geom_point(data = plot.rev.f, aes(x = x, y = y),size=0.6,color="green4")
  p3<- p3 + xlab("Sequence1") + ylab("Sequence2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p3)
}


#大片段敲除的上游点阵分析
Get_dot_plot2<-function(KO_region){
  if (KO_region$start < 500 | KO_region$end + 500 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$start + 500 - 800, KO_region$start+500)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$start)
  }
  
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  plot.f<-data.frame(x=numeric(),y=numeric())
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j])
      plot.f<-rbind(plot.f,f)
    }
  }
  
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  
  plot.rev.f<-data.frame(x=numeric(),y=numeric())
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  #analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pattern <- DNAString(pattern)
    pattern <- reverseComplement(pattern)
    pattern <- as.character(pattern)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != -1) {
        f<-data.frame(x=i,y=pos[j])
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  p3 <- ggplot(data = plot.f, aes(x = x, y = y)) + geom_point(size=0.6,color="red") + theme_bw() +
    theme(panel.grid = element_blank())+ geom_point(data = plot.rev.f, aes(x = x, y = y),size=0.6,color="green4")
  p3<- p3 + xlab("Sequence1") + ylab("Sequence2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p3)
}


#大片段敲除的下游点阵分析
Get_dot_plot3<-function(KO_region){
  if (KO_region$start < 500 | KO_region$end + 500 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$end+500, KO_region$end+500+800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$end, KO_region$end+800)
  }
  
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  plot.f<-data.frame(x=numeric(),y=numeric())
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j])
      plot.f<-rbind(plot.f,f)
    }
  }
  
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  
  plot.rev.f<-data.frame(x=numeric(),y=numeric())
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  #analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pattern <- DNAString(pattern)
    pattern <- reverseComplement(pattern)
    pattern <- as.character(pattern)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != -1) {
        f<-data.frame(x=i,y=pos[j])
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  p3 <- ggplot(data = plot.f, aes(x = x, y = y)) + geom_point(size=0.6,color="red") + theme_bw() +
    theme(panel.grid = element_blank())+ geom_point(data = plot.rev.f, aes(x = x, y = y),size=0.6,color="green4")
  p3<- p3 + xlab("Sequence1") + ylab("Sequence2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p3)
}


#大片段敲除的上游点阵分析
Get_dot_plot4<-function(KO_region){
  analysis_seq <-substring(Gene3, KO_region$start - 800, KO_region$start)
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  plot.f<-data.frame(x=numeric(),y=numeric())
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j])
      plot.f<-rbind(plot.f,f)
    }
  }
  
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  
  plot.rev.f<-data.frame(x=numeric(),y=numeric())
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  #analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pattern <- DNAString(pattern)
    pattern <- reverseComplement(pattern)
    pattern <- as.character(pattern)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != -1) {
        f<-data.frame(x=i,y=pos[j])
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  p3 <- ggplot(data = plot.f, aes(x = x, y = y)) + geom_point(size=0.6,color="red") + theme_bw() +
    theme(panel.grid = element_blank())+ geom_point(data = plot.rev.f, aes(x = x, y = y),size=0.6,color="green4")
  p3<- p3 + xlab("Sequence1") + ylab("Sequence2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p3)
}


#大片段敲除的下游点阵分析
Get_dot_plot5<-function(KO_region){
  analysis_seq <-substring(Gene3, KO_region$end, KO_region$end+800)
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  plot.f<-data.frame(x=numeric(),y=numeric())
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j])
      plot.f<-rbind(plot.f,f)
    }
  }
  
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  
  plot.rev.f<-data.frame(x=numeric(),y=numeric())
  analysis_pos <- numeric()
  len <- nchar(analysis_seq) - 9
  #analysis_seq_rev<-as.character(reverseComplement(DNAString(analysis_seq)))
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pattern <- DNAString(pattern)
    pattern <- reverseComplement(pattern)
    pattern <- as.character(pattern)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      if (pos[j] != -1) {
        f<-data.frame(x=i,y=pos[j])
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  p3 <- ggplot(data = plot.f, aes(x = x, y = y)) + geom_point(size=0.6,color="red") + theme_bw() +
    theme(panel.grid = element_blank())+ geom_point(data = plot.rev.f, aes(x = x, y = y),size=0.6,color="green4")
  p3<- p3 + xlab("Sequence1") + ylab("Sequence2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  return(p3)
}




