setwd("C://Users//41518//Desktop//work//ubigene//micro-organism/")
filepath<-"C://Users//41518//Desktop//微生物//run"
library(reticulate)
source("table_select2.R")
source_python("crispor_table_download3.py")
species<-"Saccharomyces cerevisiae S288c"

gene.table<-read.csv("C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv")

for(i in a){
  
  if(nchar(gene.table[i,]$seq)>=400){
    
    #上游
    seq1<-substring(gene.table[i,]$seq,1,200)
    py$run(seq1,species,filepath)
    gRNA.table1<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
    gRNA.table1<-table_select1(gRNA.table1)
    #下游
    seq2<-substring(gene.table[i,]$seq,nchar(gene.table[i,]$seq)-200,nchar(gene.table[i,]$seq))
    py$run(seq2,species,filepath)
    gRNA.table2<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
    gRNA.table2<-table_select2(gRNA.table2)
    
    if(!is.null(gRNA.table1) & !is.null(gRNA.table2)){
      if(nrow(gRNA.table1)!=0 & nrow(gRNA.table2)!=0){
        gRNA<-data.frame(strand=c(gRNA.table1[1,]$V4,gRNA.table2[1,]$V4),gRNA=c(gRNA.table1[1,]$V2,gRNA.table2[1,]$V2),pos=c(gRNA.table1[1,]$V1,gRNA.table2[1,]$V1+(nchar(gene.table[i,]$seq)-200-1)),score=c(gRNA.table1[1,]$V3,gRNA.table2[1,]$V3))
      }
    }

    
    if(exists("gRNA")){
      if(nrow(gRNA)==2){
        gene.table[i,]$gRNA1<-gRNA[1,]$gRNA
        gene.table[i,]$strand1<-gRNA[1,]$strand
        gene.table[i,]$pos1<-gRNA[1,]$pos
        gene.table[i,]$score1<-gRNA[1,]$score
        gene.table[i,]$gRNA2<-gRNA[2,]$gRNA
        gene.table[i,]$strand2<-gRNA[2,]$strand
        gene.table[i,]$pos2<-gRNA[2,]$pos
        gene.table[i,]$score2<-gRNA[2,]$score
        
      }
      print(paste(i,"success"))
      write.csv(gene.table,"C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv",row.names = FALSE)
      rm(gRNA)
    }
    else{
      print(paste(i,"fail"))
    }
  }
}



