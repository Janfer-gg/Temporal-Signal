#各转录本及其外显子示意图
KO_region_image2 <- function(Gene,allinfo) {
  start <- allinfo$start         #基因组起始位置
  y <- 0
  f <- data.frame(x = c(-(nchar(Gene)/50):nchar(Gene)*1.02), y = y)
  p2 <-
    ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "white") + theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_blank())+scale_y_continuous(breaks=NULL,limits = c(-19,0))+scale_x_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
  transcript.pos<-data.frame(name=character(),start=numeric(),end=numeric(),y=numeric())
  
  for (i in 1:length(allinfo$Transcript)){
    Transcript <- allinfo$Transcript[[i]]
    transcript.pos[i,]$name<-Transcript$display_name
    Transcript_start <-
      Transcript$start - start + 1         #转录本起始位置
    Transcript_end <-
      Transcript$end - start + 1             #转录本终止位置
    transcript.pos[i,]$start<-Transcript_start
    transcript.pos[i,]$end<-Transcript_end
    transcript.pos2<-transcript.pos
  }
  for(i in 1:nrow(transcript.pos2)){
    if(is.na(transcript.pos2[i,]$y)){
      y<-y-0.7
      transcript.pos2[i,]$y<-y
      for(j in i:nrow(transcript.pos2)){
        if(is.na(transcript.pos2[j,]$y)){
          if(transcript.pos2[j,]$start>=transcript.pos2[i,]$end+(nchar(Gene)/100)){
            transcript.pos2[j,]$y<-y
            transcript.pos2[i,]$end<-transcript.pos2[j,]$end
          }
        }
      }
    }
  }
  transcript.pos$y<-transcript.pos2$y
  
  
  for (i in 1:length(allinfo$Transcript))
  {
    Transcript <- allinfo$Transcript[[i]]                #转录本
    {
      if (Transcript$display_name %in% color_blue) {
        color = "blue"
      }
      else if(Transcript$source=="ensembl_havana"){
        color ="#f8b957"
      }
      else{
        color ="#D01027"
      }
    }
    CDS_start <-
      Transcript$Translation$start - start + 1    #CDS起始位置
    CDS_end <-
      Transcript$Translation$end - start + 1        #CDS终止位置
    
    y<-transcript.pos[i,]$y
    ff <- data.frame(x = c(transcript.pos[i,]$start:transcript.pos[i,]$end),
                     y = y)
    #转录本线
    p2 <-
      p2 + geom_line(data = ff,
                     aes(x = x, y = y),
                     color = "#bdc4ca",
      ) + annotate(
        "text",
        label = paste0(Transcript$display_name),
        x = transcript.pos[i,]$start,
        y = y - 0.3,
        size = 2,color=color,hjust = 0
      )
    
    #外显子框
    for (j in 1:length(Transcript$Exon))
    {
      Exon_start <- Transcript$Exon[[j]]$start - start + 1   #外显子起始位置
      Exon_end <-
        Transcript$Exon[[j]]$end - start + 1       #外显子终止位置
      p2 <-
        p2 + annotate(
          "rect",
          xmin = Exon_start,
          xmax = Exon_end,
          ymin = y - 0.1,
          ymax = y + 0.1,
          col = color,
          alpha = .0
        )
      #完全在编码区
      if (length(CDS_start) != 0 & length(CDS_end != 0)) {
        if (Exon_start >= CDS_start & Exon_end <= CDS_end) {
          p2 <-
            p2 + annotate(
              "rect",
              xmin = Exon_start,
              xmax = Exon_end,
              ymin = y - 0.1,
              ymax = y + 0.1,
              fill = color
            )
        }
        #部分在编码区(上游)
        else if (Exon_start <= CDS_start &
                 Exon_end >= CDS_start & Exon_end <= CDS_end) {
          p2 <-
            p2 + annotate(
              "rect",
              xmin = CDS_start,
              xmax = Exon_end,
              ymin = y - 0.1,
              ymax = y + 0.1,
              fill = color
            )
        }
        #部分在编码区(下游)
        else if (Exon_start >= CDS_start &
                 Exon_start <= CDS_end & Exon_end >= CDS_end) {
          p2 <-
            p2 + annotate(
              "rect",
              xmin = Exon_start,
              xmax = CDS_end,
              ymin = y - 0.1,
              ymax = y + 0.1,
              fill = color
            )
        }
      }
    }
  }
  #调用KO_region
  if (exists("KO_region_20")) {
    p2 <-
      p2 + annotate(
        "rect",
        xmin = KO_region_20[t,]$start - 300,
        xmax = KO_region_20[t,]$end + 300,
        ymin = min(transcript.pos$y) - 0.7,
        ymax = 0,
        alpha = .0,
        color = "red"
      )
  }
  #裁剪
  png(file=paste0(filepath,"//","transcript_region2.png"),width = 480*4,height = 480*3,res = 72*3)
  print(p2)
  dev.off()
  cut_large<-min(transcript.pos$y)
  {
    if (cut_large >= -2.1) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x280")
    }
    else if (cut_large >= -4.2) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x440")
    }
    else if (cut_large >= -6.3) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x600")
    }
    else if (cut_large >= -8.4) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x720")
    }
    else if (cut_large >= -10.5) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x960")
    }
    else if (cut_large >= -12.6) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x980")
    }
    else if (cut_large >= -15) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region2.png")), "x1140")
    }
  }
  return(img)
}
