setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
library(reticulate)
# library(ggplot2)
library(Biostrings)
source("table_select.R")
# source("dot_plot.R")
# source("GC_plot.R")
source_python("crispor_table_download3.py")

# 
filepath<-"C://Users//41518//Desktop//微生物//run1"
# dir.create(filepath)

gene.table<-read.csv("C://Users//41518//Desktop//微生物//micro-全.csv")


for(i in 3371:18191){
  
  if(nchar(gene.table[i,]$seq)<=500 &nchar(gene.table[i,]$seq)>=70){
    
    seq<-substring(gene.table[i,]$seq,21,nchar(gene.table[i,]$seq)-20)
    py$run(seq,gene.table[i,]$species,filepath)
    gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
    gRNA.table<-table_select(gRNA.table)
    if(length(gRNA.table)!=0){
      #删除重叠的gRNA序列
      w<-1 
      repeat{
        s<-which(abs(gRNA.table$V1 - gRNA.table[w,]$V1) < 23 & abs(gRNA.table$V1 - gRNA.table[w,]$V1)!=0)
        if(length(s)!=0){
          gRNA.table<-gRNA.table[-s,]
        }
        w<-w+1
        if(w >= nrow(gRNA.table)){
          break
        }
      }
      
      if(nrow(gRNA.table)>=2){
        gRNA<-data.frame(strand=gRNA.table[1:2,]$V4,gRNA=gRNA.table[1:2,]$V2,pos=gRNA.table[1:2,]$V1+20,score=gRNA.table[1:2,]$V3,lindel=gRNA.table[1:2,]$V8)
      }
      
      if(nrow(gRNA.table)==1){
        gRNA<-data.frame(strand=gRNA.table$V4,gRNA=gRNA.table$V2,pos=gRNA.table$V1+20,score=gRNA.table$V3,lindel=gRNA.table$V8)
      }
    }
  }
  
  #如果序列小于等于3000bp
  # if(nchar(gene.table[i,]$seq)<=3000 &nchar(gene.table[i,]$seq)>=70){
  #   seq<-substring(gene.table[i,]$seq,ceiling(nchar(gene.table[i,]$seq)/3),floor(nchar(gene.table[i,]$seq)*2/3))
  #   py$run(seq,gene.table[i,]$species,filepath)
  #   gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  #   gRNA.table<-table_select(gRNA.table)
  #   if(length(gRNA.table)!=0){
  #     #删除重叠的gRNA序列
  #     w<-1 
  #     repeat{
  #       s<-which(abs(gRNA.table$V1 - gRNA.table[w,]$V1) <= 23 & abs(gRNA.table$V1 - gRNA.table[w,]$V1)!=0)
  #       if(length(s)!=0){
  #         gRNA.table<-gRNA.table[-s,]
  #       }
  #       w<-w+1
  #       if(w >= nrow(gRNA.table)){
  #         break
  #       }
  #     }
  #     
  #     if(nrow(gRNA.table)>=2){
  #       gRNA<-data.frame(strand=gRNA.table[1:2,]$V4,gRNA=gRNA.table[1:2,]$V2,pos=gRNA.table[1:2,]$V1+ceiling(nchar(gene.table[i,]$seq)/3)-1,score=gRNA.table[1:2,]$V3)
  #     }
  #     
  #     if(nrow(gRNA.table)==1){
  #       gRNA<-data.frame(strand=gRNA.table$V4,gRNA=gRNA.table$V2,pos=gRNA.table$V1+ceiling(nchar(gene.table[i,]$seq)/3)-1,score=gRNA.table$V3)
  #     }
  #   }
  # }
  # else if(nchar(gene.table[i,]$seq)>3000){
  #   #上游
  #   seq1<-substring(gene.table[i,]$seq,500,1000)
  #   py$run(seq1,gene.table[i,]$species,filepath)
  #   gRNA.table1<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  #   gRNA.table1<-table_select(gRNA.table1)
  #   #下游
  #   seq2<-substring(gene.table[i,]$seq,nchar(gene.table[i,]$seq)-1000,nchar(gene.table[i,]$seq)-500)
  #   py$run(seq2,gene.table[i,]$species,filepath)
  #   gRNA.table2<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  #   gRNA.table2<-table_select(gRNA.table2)
  #   
  #   if(nrow(gRNA.table1)!=0 & nrow(gRNA.table2)!=0){
  #     gRNA<-data.frame(strand=c(gRNA.table1[1,]$V4,gRNA.table2[1,]$V4),gRNA=c(gRNA.table1[1,]$V2,gRNA.table2[1,]$V2),pos=c(gRNA.table1[1,]$V1+500-1,gRNA.table2[1,]$V1+(nchar(gene.table[i,]$seq)-1000-1)),score=c(gRNA.table1[1,]$V3,gRNA.table2[1,]$V3))
  #   }
  # }
  
  
  if(exists("gRNA")){
    if(nrow(gRNA)==2){
      gene.table[i,]$gRNA1<-gRNA[1,]$gRNA
      gene.table[i,]$strand1<-gRNA[1,]$strand
      gene.table[i,]$pos1<-gRNA[1,]$pos
      gene.table[i,]$lindel1<-gRNA[1,]$lindel
      gene.table[i,]$score1<-100
      gene.table[i,]$gRNA2<-gRNA[2,]$gRNA
      gene.table[i,]$strand2<-gRNA[2,]$strand
      gene.table[i,]$pos2<-gRNA[2,]$pos
      gene.table[i,]$lindel2<-gRNA[2,]$lindel
      gene.table[i,]$score2<-100
      
    }
    if(nrow(gRNA)==1){
      gene.table[i,]$gRNA1<-gRNA$gRNA
      gene.table[i,]$strand1<-gRNA$strand
      gene.table[i,]$pos1<-gRNA$pos
      gene.table[i,]$lindel1<-gRNA$lindel
      gene.table[i,]$score1<-100
    }
    print(paste(i,"success"))
    write.csv(gene.table,"C://Users//41518//Desktop//微生物//micro-全.csv",row.names = FALSE)
    rm(gRNA)
  }
}


# con <- dbConnect(MySQL(), host="47.92.215.36", dbname="YuanJingData", user="ubigene_data", password="@ubigene2020@..")







# {
#   #输出
#   if(exists("gRNA")){
#     f1<-data.frame(start=numeric(),end=numeric(),color=character(),text=character())
#     f1[1,]<-c(-500,1,"grey","HA-F")
#     f1[2,]<-c(1,nchar(gene.table$seq),"#D01027",gene.table$gene)
#     f1[3,]<-c(nchar(gene.table$seq),nchar(gene.table$seq)+500,"grey","HA-R")
#     
#     f2<-data.frame(x=c(-600,nchar(gene.table$seq)+600),y=1)
#     p <-
#       ggplot(data = f2, aes(x = x, y = y)) + geom_path(color = "black", size = 4) + theme_bw() +
#       theme(panel.grid = element_blank(), panel.border = element_blank()) +
#       scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)+
#       annotate("rect",xmin = as.numeric(f1$start),xmax = as.numeric(f1$end),ymin = 0.97,ymax=1.03,fill=f1$color,color="black")+
#       annotate("text",x=-250,y=1,label="HA-F",color="yellow")+annotate("text",x=nchar(gene.table$seq)+250,y=1,label="HA-R",color="yellow")+
#       annotate("text",x=nchar(gene.table$seq)/2,y=1,label=gene.table$gene,color="yellow",size=5)
#     
#     if(nrow(gRNA)==2){
#       f3<-data.frame(x=gRNA[1,]$pos,y=c(1.06,1.09))
#       f4<-data.frame(x=gRNA[2,]$pos,y=c(1.06,1.09))
#       p<-p+geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
#         geom_line(data = f4,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
#         annotate("text",x=gRNA[1,]$pos,y=1.12,label="g1",color="red",size=5)+
#         annotate("text",x=gRNA[2,]$pos,y=1.12,label="g2",color="red",size=5)
#       
#       output<-data.frame(Gene=gene.table$gene,gRNA1=gRNA[1,]$gRNA,strand1=gRNA[1,]$strand,score1=gRNA[1,]$score,
#                          gRNA2=gRNA[2,]$gRNA,strand2=gRNA[2,]$strand,score2=gRNA[2,]$score)
#     }
#     if(nrow(gRNA)==1){
#       f3<-data.frame(x=gRNA[1,]$pos,y=c(1.06,1.09))
#       p<-p+geom_line(data = f3,aes(x=x,y=y),arrow=arrow(length = unit(0.15, "cm"),ends = "first",type = "closed"),color="red")+
#         annotate("text",x=gRNA[1,]$pos,y=1.12,label="g1",color="red",size=5)
#       
#       output<-data.frame(Gene=gene.table$gene,gRNA1=gRNA[1,]$gRNA,strand1=gRNA[1,]$strand,score1=gRNA[1,]$score,
#                          gRNA2=NA,strand2=NA,score2=NA)
#     }
#     
#     dot_plot(gene.table$analysis_seq,gene.table$analysis_seq,10,0,filepath)
#     GC<-GC_plot(gene.table$analysis_seq,30,"GC",filepath)
#     output$GC<-GC
#     
#     write.csv(output,paste0(filepath,"//result.csv"),row.names = FALSE)
#   }
#   else{
#     write.csv("fail",paste0(filepath,"//fail.txt"),row.names = FALSE,col.names = FALSE)
#   }
# }








