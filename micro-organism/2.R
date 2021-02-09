setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
library(reticulate)
# library(ggplot2)
# library(Biostrings)
source("table_select.R")
# source("dot_plot.R")
# source("GC_plot.R")
source_python("crispor_table_download3.py")

# 
filepath<-"C://Users//41518//Desktop//微生物//run2"
# dir.create(filepath)

gene.table<-read.csv("C://Users//41518//Desktop//微生物//micro2701-3500.csv")

# species<-"K-12 substr. MG1655"
# gene.table<-gene.table[which(gene.table$species==species),]
# 
# species<-"Escherichia coli BL21(DE3)"
# gene.table<-gene.table[which(gene.table$species==species),]
# 
# species<-"Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S"
# gene.table<-gene.table[which(gene.table$species==species),]
# 
# species<-"shigella flexneri 2a str. 301"
# gene.table<-gene.table[which(gene.table$species==species),]
# 
# gene.table<-gene.table[which(gene.table$gene=="carB"),]
# 
# if(nrow(gene.table)==0| nchar(gene.table$seq)<70){
#   stop("no gene")
# }
# gene.table<-transform(gene.table,gRNA1=character(nrow(gene.table)),strand1=character(nrow(gene.table)),pos1=numeric(nrow(gene.table)),gRNA2=character(nrow(gene.table)),strand2=character(nrow(gene.table)),pos2=numeric(nrow(gene.table)))
for(i in 2701:3500){
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
      #排除完全重叠的情况
      gRNA.table<-gRNA.table[!duplicated(gRNA.table$V1),]
      
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
    #下游
    seq2<-substring(gene.table[i,]$seq,nchar(gene.table[i,]$seq)-1000,nchar(gene.table[i,]$seq)-500)
    py$run(seq2,gene.table[i,]$species,filepath)
    gRNA.table2<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
    gRNA.table2<-table_select(gRNA.table2)
    
    if(length(gRNA.table1)!=0 & length(gRNA.table2)!=0){
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
    write.csv(gene.table,"C://Users//41518//Desktop//微生物//micro2701-3500.csv",row.names = FALSE)
    rm(gRNA)
  }
}









