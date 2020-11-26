GC_plot <- function(seq,size,filepath) {
  seq <- DNAString(seq)

  GC <-rowSums(letterFrequencyInSlidingView(seq, size, c("G","C"))) / size
  GC.data <- data.frame(x = c(1:(nchar(seq) - (size-1))), y = GC)
  p <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())+theme(plot.margin=unit(rep(3,4),'lines'))
  p<- p + xlab("Sequence 1") + ylab("Sequence 2")+ theme(axis.text.x = element_text(size=12,color = "black"),axis.text.y = element_text(size=12,color = "black"))+
    scale_x_continuous(n.breaks=8)
  png(file = paste0(filepath,"//","GC.png"),width = 480*3,height = 480*2,res = 72*2)
  print(p)
  dev.off()
}

