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
  plot.f<-data.frame()
  for (i in 1:len) {
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- gregexpr(pattern, analysis_seq)[[1]]
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j],z="forward")
      plot.f<-rbind(plot.f,f)
    }
  }
  
  analysis_pos <- sort(analysis_pos[!duplicated(analysis_pos)])
  
  plot.rev.f<-data.frame()
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
        f<-data.frame(x=i,y=pos[j],z="reverseComplement")
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  dot.table<-rbind(plot.f,plot.rev.f)
  p3 <-
    ggplot(data = dot.table, aes(x = x, y = y,colour=z)) + geom_point(size=0.6) + theme_bw() +
    theme(panel.grid = element_blank())+theme(legend.position="top")+theme(legend.title=element_blank())
  p3
  return(p3)
}


