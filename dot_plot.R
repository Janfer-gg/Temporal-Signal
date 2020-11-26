dot_plot<-function(seq1,seq2,filepath){
  analysis_pos <- numeric()
  len <- nchar(seq1) - 9
  plot.f<-data.frame(x=numeric(),y=numeric())
  for (i in 1:len) {
    pattern <- substring(seq1, i, i + 9)
    pos <- start(matchPattern(pattern, seq2,max.mismatch = 0))
    for (j in 1:length(pos)) {
      f<-data.frame(x=i,y=pos[j])
      plot.f<-rbind(plot.f,f)
    }
  }
  
  plot.rev.f<-data.frame(x=numeric(),y=numeric())
  analysis_pos <- numeric()
  for (i in 1:len) {
    pattern <- substring(seq1, i, i + 9)
    pattern <- DNAString(pattern)
    pattern <- reverseComplement(pattern)
    pattern <- as.character(pattern)
    pos <- start(matchPattern(pattern, seq2,max.mismatch = 0))
    if(length(pos)!=0){
      for (j in 1:length(pos)) {
        f<-data.frame(x=i,y=pos[j])
        plot.rev.f<-rbind(plot.rev.f,f)
      }
    }
  }
  p1 <- ggplot(data = plot.f, aes(x = x, y = y)) + geom_point(size=0.5,color="red") + theme_bw() +
    theme(panel.grid = element_blank())+ geom_point(data = plot.rev.f, aes(x = x, y = y),size=0.5,color="green4")+
    scale_x_continuous(n.breaks=8)+scale_y_continuous(n.breaks=8)
  p1<- p1 + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))
  png(file = paste0(filepath,"//","dot.png"),width = 480*2,height = 480*2,res = 72*2)
  print(p1)
  dev.off()

}


