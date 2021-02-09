get_large<-function(region1,region2,gRNA,label){
  ff <-data.frame(x = gRNA$start, y = 0.59)
  large_start <- gRNA$start-300
  large_end <- gRNA$end+300
  y <- 0.5
  f <- data.frame(x = c(large_start:large_end), y = y)
  p1_large <-
    ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "black", size = 2) + theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_blank()) +
    scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
  {
    if(nrow(region1)==0) {
      start <- large_start
      end <- large_end
      p1_large <- p1_large + annotate(
        "rect",
        xmin = start,
        xmax = end,
        ymin = y - 0.04,
        ymax = y + 0.04,
        colour = "#D01027",
        alpha = .0
      )
    }
    else{
      start <- as.numeric(region1$Exon_start)
      if (start < large_start) {
        start <- large_start
      }
      end <- as.numeric(region1$Exon_end)
      if (end > large_end) {
        end <- large_end
      }
      p1_large <- p1_large + annotate(
        "rect",
        xmin = start,
        xmax = end,
        ymin = y - 0.04,
        ymax = y + 0.04,
        colour = "#D01027",
        alpha = .0
      )
    }
  }
  
  start <- as.numeric(region2$start)
  if(start<large_start){
    start<-large_start
  }
  end <- as.numeric(region2$end)
  if(end>large_end){
    end<-large_end
  }
  p1_large <- p1_large + annotate(
    "rect",
    xmin = start,
    xmax = end,
    ymin = y - 0.04,
    ymax = y + 0.04,
    fill = "#D01027",
  )+ annotate(
    "text",
    label = sub("xon ", "", region2$Exon),
    x = (start+end)/2,
    y = y - 0.07,
    size = 4,
    hjust = 0,
    color = "black"
  )
  img<-'cut.png'
  p1_large <-
    p1_large + geom_image(data = ff,aes(x = x, y = y),image=img)+annotate("text",label=label,x = gRNA$start,y = 0.67,size=6)
  return(p1_large)
}

