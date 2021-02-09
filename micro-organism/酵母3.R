setwd("C://Users//41518//Desktop//work//ubigene//micro-organism/")
filepath<-"C://Users//41518//Desktop//微生物//run2"
library(reticulate)
source("table_select2.R")
source_python("crispor_table_download3.py")
species<-"Saccharomyces cerevisiae S288c"

gene.table<-read.csv("C://Users//41518//Desktop//微生物//Sac-single.csv")

for(i in 6233:6400){
  if(nchar(gene.table[i,]$seq)<=2000){
    seq<-substring(gene.table[i,]$seq,4,nchar(gene.table[i,]$seq)-3)
  }
  if(nchar(gene.table[i,]$seq)>2000){
    seq<-substring(gene.table[i,]$seq,4,2000)
  }
  
  py$run(seq,species,filepath)
  gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  
  gRNA.table<-gRNA.table[which(!grepl("Inefficient", gRNA.table$V2)),]
  
  gRNA.table<-gRNA.table[which(as.numeric(gRNA.table$V3)>=70 & as.numeric(gRNA.table$V8)>=60),]
  if(nrow(gRNA.table)==0){
    print(paste(i,"fail"))
    next
  }
  
  gRNA.table$V2 <-gsub(" ", "", substring(gRNA.table$V2, 1, 24))
  strand <- sub("[^a-zA-Z]+", "", gRNA.table$V1)
  gRNA.table$V4<-strand
  pos <- as.numeric(sub("[^0-9]+", "", gRNA.table$V1))
  for(j in 1:length(strand)){
    if(strand[j]=="fw"){
      pos[j]=pos[j]+2
    }
    else if(strand[j]=="rev"){
      pos[j]=pos[j]+22
    }
  }
  gRNA.table$V1<-pos
  
  gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=90 & as.numeric(gRNA.table$V8)>=80 & gRNA.table$V10=="0-0-0-0-0"),]
  if(nrow(gRNA.table1)==0){
    gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
    if(nrow(gRNA.table1)==0){
      gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70),]
      if(nrow(gRNA.table1)==0){
        gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=70 & as.numeric(gRNA.table$V8)>=60),]
      }
    }
  }
  
  if(!is.null(gRNA.table1)){
    if(nrow(gRNA.table1)!=0){
      gRNA<-data.frame(strand=gRNA.table1[1,]$V4,gRNA=gRNA.table1[1,]$V2,pos=gRNA.table1[1,]$V1,score=gRNA.table1[1,]$V3)
    }
  }
    
  {
    if(exists("gRNA")){
      gene.table[i, ]$gRNA1 <- gRNA$gRNA
      gene.table[i,]$strand1 <- gRNA$strand
      gene.table[i,]$pos1 <- gRNA$pos
      gene.table[i,]$score1 <- gRNA$score

      print(paste(i,"success"))
      write.csv(gene.table,"C://Users//41518//Desktop//微生物//Sac-single.csv",row.names = FALSE)
      rm(gRNA)
    }
    else{
      print(paste(i,"fail"))
    }
  }
}
