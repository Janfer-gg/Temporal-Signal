#内侧靶位点和外侧靶位点都有
# primer.table<-read.csv(paste0(filepath2,"//resultA.csv"),header = TRUE)

image1<-function(output,Exon_region,Gene,Gene_rev){
  primer.table<-output
  {
    #反向基因
    if(Gene_rev){
      F1.pos<-matchPattern(reverse(DNAString(primer.table$F1)), subject = Gene)
      F2.pos<-matchPattern(reverse(DNAString(primer.table$F2)), subject = Gene)
      R1.pos<-matchPattern(complement(DNAString(primer.table$R1)),subject = Gene)
      R2.pos<-matchPattern(complement(DNAString(primer.table$R2)),subject = Gene)
      
      f <- data.frame(x=c((start(R1.pos)-200):(end(F1.pos)+200)),y=1)
      p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
      
      f1<-data.frame(x=c(end(F1.pos):start(F1.pos)),y=1.10)
      f2<-data.frame(x=c(start(R1.pos):end(R1.pos)),y=1.10)
      f3<-data.frame(x=c(end(F2.pos):start(F2.pos)),y=1.10)
      f4<-data.frame(x=c(start(R2.pos):end(R2.pos)),y=1.07)
      
      p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        geom_line(data = f4,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        annotate("text",label="F1",x=start(F1.pos),y=1.13,size=4)+
        annotate("text",label="R1",x=start(R1.pos),y=1.13,size=4)+
        annotate("text",label="F2",x=start(F2.pos),y=1.13,size=4)+
        annotate("text",label="R2",x=start(R2.pos),y=1.05,size=3)
      
    }
    
    #正向基因
    else{
      F1.pos<-matchPattern(DNAString(primer.table$F1), subject = Gene)
      F2.pos<-matchPattern(DNAString(primer.table$F2), subject = Gene)
      R1.pos<-matchPattern(reverseComplement(DNAString(primer.table$R1)),subject = Gene)
      R2.pos<-matchPattern(reverseComplement(DNAString(primer.table$R2)),subject = Gene)
      f <- data.frame(x=c((start(F1.pos)-200):(end(R1.pos)+200)),y=1)
      p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
      
      f1<-data.frame(x=c(start(F1.pos):end(F1.pos)),y=1.10)
      f2<-data.frame(x=c(end(R1.pos):start(R1.pos)),y=1.10)
      f3<-data.frame(x=c(start(F2.pos):end(F2.pos)),y=1.10)
      f4<-data.frame(x=c(end(R2.pos):start(R2.pos)),y=1.07)
      
      p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        geom_line(data = f4,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        annotate("text",label="F1",x=start(F1.pos),y=1.13,size=4)+
        annotate("text",label="R1",x=start(R1.pos),y=1.13,size=4)+
        annotate("text",label="F2",x=start(F2.pos),y=1.13,size=4)+
        annotate("text",label="R2",x=start(R2.pos),y=1.05,size=3)
    }
  }
  
  for(i in 1:nrow(Exon_region)){
    p2 <- p2 + annotate("rect",xmin = Exon_region[i,]$start,xmax = Exon_region[i,]$end,ymin = 0.96,ymax = 1.04,fill = "#D01027")+
      annotate("text",label=sub("xon ", "", Exon_region[i,]$Exon),x=(Exon_region[i,]$start+Exon_region[i,]$end)/2,y=0.94,size=4)
  }
  return(p2)
}




#只有外侧靶位点
image2<-function(output,Exon_region,Gene,Gene_rev){
  primer.table<-output
  {
    #反向基因
    if(Gene_rev){
      F1.pos<-matchPattern(reverse(DNAString(primer.table$primer1)), subject = Gene)
      R1.pos<-matchPattern(complement(DNAString(primer.table$primer2)),subject = Gene)
      f <- data.frame(x=c((start(R1.pos)-200):(end(F1.pos)+200)),y=1)
      p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
      
      f1<-data.frame(x=c(end(F1.pos):start(F1.pos)),y=1.10)
      f2<-data.frame(x=c(start(R1.pos):end(R1.pos)),y=1.10)
      
      p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        annotate("text",label="F1",x=start(F1.pos),y=1.13,size=4)+
        geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        annotate("text",label="R1",x=start(R1.pos),y=1.13,size=4)
      
    }
    
    #正向基因
    else{
      F1.pos<-matchPattern(DNAString(primer.table$primer1), subject = Gene)
      R1.pos<-matchPattern(reverseComplement(DNAString(primer.table$primer2)),subject = Gene)
      f <- data.frame(x=c((start(F1.pos)-200):(end(R1.pos)+200)),y=1)
      p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
      
      f1<-data.frame(x=c(start(F1.pos):end(F1.pos)),y=1.10)
      f2<-data.frame(x=c(end(R1.pos):start(R1.pos)),y=1.10)
      
      p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        annotate("text",label="F1",x=start(F1.pos),y=1.13,size=4)+
        geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        annotate("text",label="R1",x=start(R1.pos),y=1.13,size=4)
 
    }
  }
  
  for (i in 1:nrow(Exon_region)) {
    p2 <-
      p2 + annotate(
        "rect",
        xmin = Exon_region[i, ]$start,
        xmax = Exon_region[i, ]$end,
        ymin = 0.96,
        ymax = 1.04,
        fill = "#D01027"
      )+ annotate("text",label=sub("xon ", "", Exon_region[i,]$Exon),x=(Exon_region[i,]$start+Exon_region[i,]$end)/2,y=0.94,size=4)
  }
  return(p2)
}


#有一个内侧靶位点
image3<-function(output,Exon_region,Gene,Gene_rev){
  primer.table<-output
  {
    #反向基因
    if(Gene_rev){
      F1.pos<-matchPattern(reverse(DNAString(primer.table$F1)), subject = Gene)
      R1.pos<-matchPattern(complement(DNAString(primer.table$R1)),subject = Gene)
      
      f <- data.frame(x=c((start(R1.pos)-200):(end(F1.pos)+200)),y=1)
      p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
      
      f1<-data.frame(x=c(end(F1.pos):start(F1.pos)),y=1.10)
      f2<-data.frame(x=c(start(R1.pos):end(R1.pos)),y=1.10)
      p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        annotate("text",label="F1",x=start(F1.pos),y=1.13,size=4)+
        geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        annotate("text",label="R1",x=start(R1.pos),y=1.13,size=4)
      
      if(!is.na(primer.table$F2)){
        F2.pos<-matchPattern(reverse(DNAString(primer.table$F2)), subject = Gene)
        f3<-data.frame(x=c(end(F2.pos):start(F2.pos)),y=1.10)
        p2<-p2+geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
          annotate("text",label="F2",x=start(F2.pos),y=1.13,size=4)
      }
      else if(!is.na(primer.table$R2)){
        R2.pos<-matchPattern(complement(DNAString(primer.table$R2)),subject = Gene)
        f4<-data.frame(x=c(start(R2.pos):end(R2.pos)),y=1.10)
        p2<-p2+geom_line(data = f4,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
          annotate("text",label="R2",x=start(R2.pos),y=1.13,size=4)
      }
      
    }
    
    #正向基因
    else{
      F1.pos<-matchPattern(DNAString(primer.table$F1), subject = Gene)
      R1.pos<-matchPattern(reverseComplement(DNAString(primer.table$R1)),subject = Gene)
      f <- data.frame(x=c((start(F1.pos)-200):(end(R1.pos)+200)),y=1)
      p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
      
      f1<-data.frame(x=c(start(F1.pos):end(F1.pos)),y=1.10)
      f2<-data.frame(x=c(end(R1.pos):start(R1.pos)),y=1.10)
      p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
        annotate("text",label="F1",x=start(F1.pos),y=1.13,size=4)+
        geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
        annotate("text",label="R1",x=start(R1.pos),y=1.13,size=4)
      
      if(!is.na(primer.table$F2)){
        F2.pos<-matchPattern(DNAString(primer.table$F2), subject = Gene)
        f3<-data.frame(x=c(start(F2.pos):end(F2.pos)),y=1.10)
        p2<-p2+geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
          annotate("text",label="F2",x=start(F2.pos),y=1.13,size=4)
      }
      if(!is.na(primer.table$R2)){
        R2.pos<-matchPattern(reverseComplement(DNAString(primer.table$R2)),subject = Gene)
        f4<-data.frame(x=c(end(R2.pos):start(R2.pos)),y=1.10)        
        p2<-p2+geom_line(data = f4,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
          annotate("text",label="R2",x=start(R2.pos),y=1.13,size=4)
      }
      
    }
  }
  
  
  for(i in 1:nrow(Exon_region)){
    p2 <- p2 + annotate("rect",xmin = Exon_region[i,]$start,xmax = Exon_region[i,]$end,ymin = 0.96,ymax = 1.04,fill = "#D01027")+
      annotate("text",label=sub("xon ", "", Exon_region[i,]$Exon),x=(Exon_region[i,]$start+Exon_region[i,]$end)/2,y=0.94,size=4)
  }
  return(p2)
}


#移码方案画图
# image4<-function(output){
#   primer.table<-output
#   {
#     #反向基因
#     if(Gene_rev){
#       F1.pos<-matchPattern(reverse(DNAString(primer.table$F1)), subject = Gene)
#       R1.pos<-matchPattern(complement(DNAString(primer.table$R1)),subject = Gene)
#       f <- data.frame(x=c((start(R1.pos)-200):(end(F1.pos)+200)),y=1)
#       p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
#         theme(panel.grid = element_blank(), panel.border = element_blank()) +
#         scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
#       
#       f1<-data.frame(x=c(end(F1.pos):start(F1.pos)),y=1.10)
#       f2<-data.frame(x=c(start(R1.pos):end(R1.pos)),y=1.10)
#       
#       p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
#         geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")
#       
#     }
#     
#     #正向基因
#     else{
#       F1.pos<-matchPattern(DNAString(primer.table$F1), subject = Gene)
#       R1.pos<-matchPattern(reverseComplement(DNAString(primer.table$R1)),subject = Gene)
#       f <- data.frame(x=c((start(F1.pos)-200):(end(R1.pos)+200)),y=1)
#       p1 <-ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
#         theme(panel.grid = element_blank(), panel.border = element_blank()) +
#         scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
#       
#       f1<-data.frame(x=c(start(F1.pos):end(F1.pos)),y=1.10)
#       f2<-data.frame(x=c(end(R1.pos):start(R1.pos)),y=1.10)
#       
#       p2<-p1+geom_line(data = f1,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "last",type = "closed"),color="red")+
#         geom_line(data = f2,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")
#       
#     }
#   }
#   
#   # for (i in 1:nrow(Exon_region)) {
#   #   p2 <-
#   #     p2 + annotate(
#   #       "rect",
#   #       xmin = Exon_region[i, ]$start,
#   #       xmax = Exon_region[i, ]$end,
#   #       ymin = 0.96,
#   #       ymax = 1.04,
#   #       fill = "#D01027"
#   #     )
#   # }
#   return(p2)
# }
