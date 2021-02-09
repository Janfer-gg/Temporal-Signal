setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
library(reticulate)
# library(ggplot2)
# library(Biostrings)
source("table_select.R")
# source("dot_plot.R")
# source("GC_plot.R")
source_python("crispor_table_download3.py")

# 
filepath<-"C://Users//41518//Desktop//微生物//run1"
# dir.create(filepath)

gene.table<-read.csv("C://Users//41518//Desktop//微生物//micro-补充.csv")

for(i in 8953:9030){
  
  if(is.na(gene.table[i,]$pos1)){
    #如果序列小于等于3000bp
    if(nchar(gene.table[i,]$seq)<=3000 &nchar(gene.table[i,]$seq)>=70){
      seq<-substring(gene.table[i,]$seq,ceiling(nchar(gene.table[i,]$seq)/3),floor(nchar(gene.table[i,]$seq)*2/3))
      py$run(seq,gene.table[i,]$species,filepath)
      gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
      gRNA.table<-table_select(gRNA.table)
      if(length(gRNA.table)!=0){
        #删除重叠的gRNA序列
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
          gRNA<-data.frame(strand=gRNA.table[1:2,]$V4,gRNA=gRNA.table[1:2,]$V2,pos=gRNA.table[1:2,]$V1+ceiling(nchar(gene.table[i,]$seq)/3)-1,score=gRNA.table[1:2,]$V3)
        }
        
        if(nrow(gRNA.table)==1){
          gRNA<-data.frame(strand=gRNA.table$V4,gRNA=gRNA.table$V2,pos=gRNA.table$V1+ceiling(nchar(gene.table[i,]$seq)/3)-1,score=gRNA.table$V3)
        }
      }
    }
    else if(nchar(gene.table[i,]$seq)>3000){
      #上游
      seq1<-substring(gene.table[i,]$seq,500,1000)
      py$run(seq1,gene.table[i,]$species,filepath)
      gRNA.table1<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
      gRNA.table1<-table_select(gRNA.table1)
      if(length(gRNA.table1)==0){
        next
      }
      #下游
      seq2<-substring(gene.table[i,]$seq,nchar(gene.table[i,]$seq)-1000,nchar(gene.table[i,]$seq)-500)
      py$run(seq2,gene.table[i,]$species,filepath)
      gRNA.table2<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
      gRNA.table2<-table_select(gRNA.table2)
      if(length(gRNA.table2)==0){
        next
      }
      
      if(nrow(gRNA.table1)!=0 & nrow(gRNA.table2)!=0){
        gRNA<-data.frame(strand=c(gRNA.table1[1,]$V4,gRNA.table2[1,]$V4),gRNA=c(gRNA.table1[1,]$V2,gRNA.table2[1,]$V2),pos=c(gRNA.table1[1,]$V1+500-1,gRNA.table2[1,]$V1+(nchar(gene.table[i,]$seq)-1000-1)),score=c(gRNA.table1[1,]$V3,gRNA.table2[1,]$V3))
      }
    }
    if(exists("gRNA")){
      if(nrow(gRNA)==2){
        gene.table[i,]$gRNA1<-gRNA[1,]$gRNA
        gene.table[i,]$strand1<-gRNA[1,]$strand
        gene.table[i,]$pos1<-gRNA[1,]$pos
        gene.table[i,]$gRNA2<-gRNA[2,]$gRNA
        gene.table[i,]$strand2<-gRNA[2,]$strand
        gene.table[i,]$pos2<-gRNA[2,]$pos
        
      }
      if(nrow(gRNA)==1){
        gene.table[i,]$gRNA1<-gRNA$gRNA
        gene.table[i,]$strand1<-gRNA$strand
        gene.table[i,]$pos1<-gRNA$pos
      }
      print(paste(i,"success"))
      rm(gRNA)
    }
  }
 
}









