#各转录本及其外显子示意图
KO_region_image3 <- function(Gene, allinfo) {
  start <- allinfo$start         #基因组起始位置
  y <- 0
  f <- data.frame(x = c(-(nchar(Gene)/50):nchar(Gene)*1.02), y = y)
  p2 <-
    ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "white") + theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_blank())+scale_y_continuous(breaks=NULL,limits = c(-21.7,0))+scale_x_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
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
          if(transcript.pos2[j,]$start>=transcript.pos2[i,]$end+(nchar(Gene)/100)& transcript.pos[j,]$end-transcript.pos[j,]$start>=(nchar(Gene)/10) | transcript.pos2[j,]$start>=transcript.pos2[i,]$end+(nchar(Gene)/20)){
            transcript.pos2[j,]$y<-y
            transcript.pos2[i,]$end<-transcript.pos2[j,]$end
          }
          else if(transcript.pos2[j,]$end<=transcript.pos2[i,]$start-(nchar(Gene)/100)& transcript.pos[j,]$end-transcript.pos[j,]$start>=(nchar(Gene)/10) | transcript.pos2[j,]$end<=transcript.pos2[i,]$start-(nchar(Gene)/20)){
            transcript.pos2[j,]$y<-y
            transcript.pos2[i,]$start<-transcript.pos2[j,]$start
          }
        }
      }
    }
  }
  transcript.pos$y<-transcript.pos2$y
  
  
  Exon_data<-data.frame(start=numeric(),end=numeric(),y=numeric(),fill=character(),color=character())
  p<-1
  for (i in 1:length(allinfo$Transcript)){
    Transcript <- allinfo$Transcript[[i]]                #转录本
    {
      if (Transcript$display_name %in% color_blue) {
        color = "blue"
      }
      else if(Transcript$source=="ensembl_havana"){
        color ="orange2"
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
    if(y<(-21)){
      next
    }
    ff <- data.frame(x = c(transcript.pos[i,]$start,transcript.pos[i,]$end),y = y)
    {
      if (Gene_rev) {
        label <- paste0(Transcript$display_name, "<")
      }
      else{
        label <- paste0(Transcript$display_name, ">")
      }
    }
    #转录本线
    p2 <-
      p2 + geom_line(data = ff,
                     aes(x = x, y = y),
                     color = "#bdc4ca", ) + annotate(
                       "text",
                       label = label,
                       x = transcript.pos[i, ]$start,
                       y = y - 0.3,
                       size = 2,
                       color = color,
                       hjust = 0
                     )
    
    #外显子框
    for (j in 1:length(Transcript$Exon)){
      Exon_start <- Transcript$Exon[[j]]$start - start + 1   #外显子起始位置
      Exon_end <-
        Transcript$Exon[[j]]$end - start + 1       #外显子终止位置
      
      Exon_data[p,]$start<-Exon_start
      Exon_data[p,]$end<-Exon_end
      Exon_data[p,]$y<-y
      Exon_data[p,]$fill<-0
      Exon_data[p,]$color<-color
      p<-p+1
      
      #完全在编码区
      if (length(CDS_start) != 0 & length(CDS_end != 0)) {
        if (Exon_start >= CDS_start & Exon_end <= CDS_end) {
          Exon_data[p,]$start<-Exon_start
          Exon_data[p,]$end<-Exon_end
          Exon_data[p,]$y<-y
          Exon_data[p,]$fill<-color
          p<-p+1
        }
        #部分在编码区(上游)
        else if (Exon_start <= CDS_start &
                 Exon_end >= CDS_start & Exon_end <= CDS_end) {
          Exon_data[p,]$start<-CDS_start
          Exon_data[p,]$end<-Exon_end
          Exon_data[p,]$y<-y
          Exon_data[p,]$fill<-color
          p<-p+1
        }
        
        #部分在编码区(下游)
        else if (Exon_start >= CDS_start &
                 Exon_start <= CDS_end & Exon_end >= CDS_end) {
          Exon_data[p,]$start<-Exon_start
          Exon_data[p,]$end<-CDS_end
          Exon_data[p,]$y<-y
          Exon_data[p,]$fill<-color
          p<-p+1
          
        }
        #部分在编码区(中间)
        else if (Exon_start<CDS_start & Exon_end>CDS_start &Exon_end>CDS_end){
          Exon_data[p,]$start<-CDS_start
          Exon_data[p,]$end<-CDS_end
          Exon_data[p,]$y<-y
          Exon_data[p,]$fill<-color
          p<-p+1
        }
        
      }
    }
  }
  
  Exon_data1<-Exon_data[which(Exon_data$fill==0),]
  Exon_data2<-Exon_data[which(Exon_data$fill!=0),]
  p2 <-
    p2 + annotate(
      "rect",
      xmin = Exon_data1$start,
      xmax = Exon_data1$end,
      ymin = Exon_data1$y - 0.1,
      ymax = Exon_data1$y + 0.1,
      color=Exon_data1$color,
      alpha=.0
    )
  p2 <-
    p2 + annotate(
      "rect",
      xmin = Exon_data2$start,
      xmax = Exon_data2$end,
      ymin = Exon_data2$y - 0.1,
      ymax = Exon_data2$y + 0.1,
      fill = Exon_data2$fill,
    )
  #调用KO_region
  {
    if (exists("Gene3")) {
      {
        if(min(transcript.pos$y)<(-21)){
          p2 <-
            p2 + annotate(
              "rect",
              xmin = min(gRNA_planC$start)-1200,
              xmax = max(gRNA_planC$end)-1200,
              ymin = -21.7,
              ymax = 0,
              alpha = .0,
              color = "grey"
            )
        }
        else{
          p2 <-
            p2 + annotate(
              "rect",
              xmin = min(gRNA_planC$start)-1200,
              xmax = max(gRNA_planC$end)-1200,
              ymin = min(transcript.pos$y) - 0.7,
              ymax = 0,
              alpha = .0,
              color = "grey"
            )
        }
      }
    }
    else{
      {
        if(min(transcript.pos$y)<(-21)){
          p2 <-
            p2 + annotate(
              "rect",
              xmin = min(gRNA_planC$start),
              xmax = max(gRNA_planC$end),
              ymin = -21.7,
              ymax = 0,
              alpha = .0,
              color = "grey"
            )
        }
        else{
          p2 <-
            p2 + annotate(
              "rect",
              xmin = min(gRNA_planC$start),
              xmax = max(gRNA_planC$end),
              ymin = min(transcript.pos$y) - 0.7,
              ymax = 0,
              alpha = .0,
              color = "grey"
            )
        }
      }
    }
  }
  #裁剪
  png(file = paste0(filepath,"//","transcript_region3.png"),width = 480*4,height = 480*3,res = 72*3)
  print(p2)
  dev.off()
  cut_large<-min(transcript.pos$y)
  {
    if (cut_large >= -2.1) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x260")
    }
    else if (cut_large >= -4.2) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x420")
    }
    else if (cut_large >= -6.3) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x580")
    }
    else if (cut_large >= -8.4) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x700")
    }
    else if (cut_large >= -10.5) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x840")
    }
    else if (cut_large >= -12.6) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x960")
    }
    else if (cut_large >= -15) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x1120")
    }
    else if (cut_large >= -18) {
      img <-
        image_crop(image_read(paste0(filepath,"//","transcript_region3.png")), "x1320")
    }
    else{
      img <- image_read(paste0(filepath,"//","transcript_region3.png"))
    }
  }
  return(img)
}
