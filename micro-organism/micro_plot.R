micro_plot<-function(gene,length,seq,filepath,pos1,pos2=NULL){
  library(Biostrings)
  library(ggplot2)
  library(ggimage)
  # source("dot_plot.R")
  # source("GC_plot.R")
  f1<-data.frame(start=numeric(),end=numeric(),color=character(),text=character())
  f1[1,]<-c(-500,1,"#666666","HA-F")
  f1[2,]<-c(1,length,"#D01027",gene)
  f1[3,]<-c(length,length+500,"#666666","HA-R")
  
  if(length<500){
    f2<-data.frame(x=c(-600,length+600),y=1)
  }
  if(length>=500){
    f2<-data.frame(x=c(-800,length+800),y=1)
  }
  
  p <-
    ggplot(data = f2, aes(x = x, y = y)) + geom_path(color = "#202020", size = 2) + theme_bw() +
    theme(panel.grid = element_blank(), panel.border = element_blank()) +
    scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)+
    annotate("rect",xmin = as.numeric(f1$start),xmax = as.numeric(f1$end),ymin = 0.96,ymax=1.04,fill=f1$color,color="#dddddd")+
    annotate("text",x=-250,y=1,label="HA-F",color="white",fontface="bold",size=6)+annotate("text",x=length+250,y=1,label="HA-R",color="white",fontface="bold",size=6)+
    annotate("text",x=length/2,y=1,label=gene,color="white",fontface="bold",size=6)
  
  img<-"cut.png"
  f3<-data.frame(x=pos1,y=1.09)
  
  p<-p+geom_image(data = f3,aes(x = x, y = y),image=img)+
    annotate("text",x=pos1,y=1.17,label="g1",color="black",size=5)
  
  if(!is.null(pos2)){
    img2<-"cut2.png"
    f4<-data.frame(x=pos2,y=0.91)
    p<-p+geom_image(data = f4,aes(x = x, y = y),image=img2)+
      annotate("text",x=pos2,y=0.85,label="g2",color="black",size=5)
  }
  
  png(file = paste0(filepath,"//","plot1.png"),width = 480*3,height = 480*2,res = 72*2)
  print(p)
  dev.off()
  
  pp<-ggplot()+annotate("segment",x=1,xend = 1400,y=1,yend = 1,color="#202020",size=2)+
    annotate("rect",xmin = 200,xmax = 700,ymin = 0.96,ymax=1.04,fill="#666666",color="#dddddd")+
    annotate("rect",xmin = 700,xmax = 1200,ymin = 0.96,ymax=1.04,fill="#666666",color="#dddddd")+
    annotate("text",x=450,y=1,label="HA-F",color="white",fontface="bold",size=6)+annotate("text",x=950,y=1,label="HA-R",color="white",fontface="bold",size=6)+
    theme_bw()+theme(panel.grid = element_blank(), panel.border = element_blank()) +
    scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
  
  png(file = paste0(filepath,"//","plot2.png"),width = 480*3,height = 480*2,res = 72*2)
  print(pp)
  dev.off()
  
  # dot_plot(seq,seq,filepath)
  # GC_plot(seq,30,filepath)
}


sample<-read.csv("C://Users//41518//Desktop//微生物//数据安全/micro-全.csv")
sample<-sample[which(sample$species=="K-12 substr. MG1655"),]

sam<-sample[which(sample$gene=="envC"),]

gene<-sam$gene
length<-sam$end-sam$start+1
seq<-sam$analysis_seq
pos1<-sam$pos1
pos2<-sam$pos2
filepath<-"C://Users//41518//Desktop"
micro_plot(gene,length,seq,filepath,pos1,pos2)

