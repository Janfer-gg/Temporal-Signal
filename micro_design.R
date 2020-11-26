micro_design<-function(species,gene,seq1,seq2,filepath){
  
  setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
  library(reticulate)
  library(ggplot2)
  library(Biostrings)
  source("table_select.R")
  source("dot_plot.R")
  source("GC_plot.R")
  source_python("crispor_table_download3.py")
  
  {
    if(nchar(seq1)<=3000 &nchar(seq1)>=70){
      seq<-substring(seq1,ceiling(nchar(seq1)/3),floor(nchar(seq1)*2/3))
      py$run(seq,species,filepath)
      gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
      gRNA.table<-table_select(gRNA.table)
      if(length(gRNA.table)!=0){
        
        w<-1 
        repeat{
          s<-which(abs(gRNA.table$V1 - gRNA.table[w,]$V1) <= 23 & abs(gRNA.table$V1 - gRNA.table[w,]$V1)!=0)
          if(length(s)!=0){
            gRNA.table<-gRNA.table[-s,]
          }
          w<-w+1
          if(w >= nrow(gRNA.table)){
            break
          }
        }
        
        if(nrow(gRNA.table)>=2){
          gRNA<-data.frame(strand=gRNA.table[1:2,]$V4,gRNA=gRNA.table[1:2,]$V2,pos=gRNA.table[1:2,]$V1+ceiling(nchar(seq1)/3)-1,score=gRNA.table[1:2,]$V3)
        }
        
        if(nrow(gRNA.table)==1){
          gRNA<-data.frame(strand=gRNA.table$V4,gRNA=gRNA.table$V2,pos=gRNA.table$V1+ceiling(nchar(seq1)/3)-1,score=gRNA.table$V3)
        }
        
      }
    }
    else if(nchar(seq1)>3000){
      
      seq_up<-substring(seq1,500,1000)
      py$run(seq_up,species,filepath)
      gRNA.table1<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
      gRNA.table1<-table_select(gRNA.table1)
      
      seq_down<-substring(seq1,nchar(seq1)-1000,nchar(seq1)-500)
      py$run(seq_down,species,filepath)
      gRNA.table2<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
      gRNA.table2<-table_select(gRNA.table2)
      
      if(nrow(gRNA.table1)!=0 & nrow(gRNA.table2)!=0){
        gRNA<-data.frame(strand=c(gRNA.table1[1,]$V4,gRNA.table2[1,]$V4),gRNA=c(gRNA.table1[1,]$V2,gRNA.table2[1,]$V2),pos=c(gRNA.table1[1,]$V1+500-1,gRNA.table2[1,]$V1+(nchar(seq1)-1000-1)),score=c(gRNA.table1[1,]$V3,gRNA.table2[1,]$V3))
      }
    }
  }
  
  {
    if(exists("gRNA")){
      
      f1<-data.frame(start=numeric(),end=numeric(),color=character(),text=character())
      f1[1,]<-c(-500,1,"grey","HA-F")
      f1[2,]<-c(1,nchar(seq1),"#D01027",gene)
      f1[3,]<-c(nchar(seq1),nchar(seq1)+500,"grey","HA-R")
      
      f2<-data.frame(x=c(-600,nchar(seq1)+600),y=1)
      p <-
        ggplot(data = f2, aes(x = x, y = y)) + geom_path(color = "black", size = 4) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)+
        annotate("rect",xmin = as.numeric(f1$start),xmax = as.numeric(f1$end),ymin = 0.97,ymax=1.03,fill=f1$color,color="black")+
        annotate("text",x=-250,y=1,label="HA-F",color="yellow")+annotate("text",x=nchar(seq1)+250,y=1,label="HA-R",color="yellow")+
        annotate("text",x=nchar(seq1)/2,y=1,label=gene,color="yellow",size=5)
      
      if(nrow(gRNA)==2){
        f3<-data.frame(x=gRNA[1,]$pos,y=c(1.06,1.09))
        f4<-data.frame(x=gRNA[2,]$pos,y=c(1.06,1.09))
        p<-p+geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
          geom_line(data = f4,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
          annotate("text",x=gRNA[1,]$pos,y=1.12,label="g1",color="red",size=5)+
          annotate("text",x=gRNA[2,]$pos,y=1.12,label="g2",color="red",size=5)
        
        output<-data.frame(Gene=gene,gRNA1=gRNA[1,]$gRNA,strand1=gRNA[1,]$strand,score1=gRNA[1,]$score,
                           gRNA2=gRNA[2,]$gRNA,strand2=gRNA[2,]$strand,score2=gRNA[2,]$score)
      }
      if(nrow(gRNA)==1){
        f3<-data.frame(x=gRNA[1,]$pos,y=c(1.06,1.09))
        p<-p+geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
          annotate("text",x=gRNA[1,]$pos,y=1.12,label="g1",color="red",size=5)
        
        output<-data.frame(Gene=gene,gRNA1=gRNA[1,]$gRNA,strand1=gRNA[1,]$strand,score1=gRNA[1,]$score,
                           gRNA2=NA,strand2=NA,score2=NA)
      }
      
      dot_plot(seq2,seq2,10,0,filepath)
      GC<-GC_plot(seq2,30,filepath)
      
      output$GC<-GC
      
      write.csv(output,paste0(filepath,"//result.csv"),row.names = FALSE)
    }
    
    else{
      write.csv("fail",paste0(filepath,"//fail.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
  return(output)
}