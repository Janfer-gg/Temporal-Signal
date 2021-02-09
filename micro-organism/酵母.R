gene.table<-read.csv("C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv")

setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
library(reticulate)
source_python("crispor_table_download3.py")
source("table_select2.R")
# 
filepath<-"C://Users//41518//Desktop//微生物//run"

species<-"Saccharomyces cerevisiae S288c"


# gene.table<-read.csv("C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv")
for(k in d){
  # if(is.na(gene.table[k,]$gene)){
  # 
  #     py$run(gene.table[k,]$seq,species,filepath)
  #     gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  # 
  #     gRNA.table<-gRNA.table[which(!grepl("Inefficient", gRNA.table$V2)),]
  # 
  #     gRNA.table<-gRNA.table[which(as.numeric(gRNA.table$V3)>70 & as.numeric(gRNA.table$V8)>60),]
  #     if(nrow(gRNA.table)<2){
  #       print(paste(k,"fail"))
  #       next
  #     }
  # 
  #     gRNA.table$V2 <-gsub(" ", "", substring(gRNA.table$V2, 1, 24))
  #     strand <- sub("[^a-zA-Z]+", "", gRNA.table$V1)
  #     gRNA.table$V4<-strand
  #     pos <- as.numeric(sub("[^0-9]+", "", gRNA.table$V1))
  #     for(j in 1:length(strand)){
  #       if(strand[j]=="fw"){
  #         pos[j]=pos[j]+2
  #       }
  #       else if(strand[j]=="rev"){
  #         pos[j]=pos[j]+22
  #       }
  #     }
  #     gRNA.table$V1<-pos
  # 
  # 
  #     table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
  # 
  #     {
  #       if(nrow(table1)>=2){
  #         for(i in 1:(nrow(table1)-1)){
  #           for(j in (i+1):nrow(table1)){
  #             if(abs(table1[i,]$V1 - table1[j,]$V1) >= 23){
  #               gRNA<-data.frame(strand=c(table1[i,]$V4,table1[j,]$V4),gRNA=c(table1[i,]$V2,table1[j,]$V2),pos=c(table1[i,]$V1,table1[j,]$V1),score=c(table1[i,]$V3,table1[j,]$V3))
  #             }
  #           }
  #         }
  # 
  #         if(exists("gRNA")==FALSE){
  #           for(i in 1:(nrow(gRNA.table)-1)){
  #             for(j in (i+1):nrow(gRNA.table)){
  #               if(abs(gRNA.table[i,]$V1 - gRNA.table[j,]$V1) >= 23){
  #                 gRNA<-data.frame(strand=c(gRNA.table[i,]$V4,gRNA.table[j,]$V4),gRNA=c(gRNA.table[i,]$V2,gRNA.table[j,]$V2),pos=c(gRNA.table[i,]$V1,gRNA.table[j,]$V1),score=c(gRNA.table[i,]$V3,gRNA.table[j,]$V3))
  #               }
  #             }
  #           }
  #         }
  #       }
  # 
  #       else{
  #         for(i in 1:(nrow(gRNA.table)-1)){
  #           for(j in (i+1):nrow(gRNA.table)){
  #             if(abs(gRNA.table[i,]$V1 - gRNA.table[j,]$V1) >= 23){
  #               gRNA<-data.frame(strand=c(gRNA.table[i,]$V4,gRNA.table[j,]$V4),gRNA=c(gRNA.table[i,]$V2,gRNA.table[j,]$V2),pos=c(gRNA.table[i,]$V1,gRNA.table[j,]$V1),score=c(gRNA.table[i,]$V3,gRNA.table[j,]$V3))
  #             }
  #           }
  #         }
  #       }
  #     }
  # 
  #     {
  #       if(exists("gRNA")){
  #         gene.table[k, ]$gRNA1 <- gRNA[1, ]$gRNA
  #         gene.table[k, ]$strand1 <- gRNA[1, ]$strand
  #         gene.table[k, ]$pos1 <- gRNA[1, ]$pos
  #         gene.table[k, ]$score1 <- gRNA[1, ]$score
  #         gene.table[k, ]$gRNA2 <- gRNA[2, ]$gRNA
  #         gene.table[k, ]$strand2 <- gRNA[2, ]$strand
  #         gene.table[k, ]$pos2 <- gRNA[2, ]$pos
  #         gene.table[k, ]$score2 <- gRNA[2, ]$score
  #         print(paste(k,"success"))
  #         write.csv(gene.table,"C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv",row.names = FALSE)
  #         rm(gRNA)
  #       }
  #       else{
  #         print(paste(k,"fail"))
  #       }
  #     }
  # }
  #     
  
  if(nchar(gene.table[k,]$seq)>=400){

    #上游
    seq1<-substring(gene.table[k,]$seq,1,200)
    py$run(seq1,species,filepath)
    gRNA.table1<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
    gRNA.table1<-table_select1(gRNA.table1)
    #下游
    seq2<-substring(gene.table[k,]$seq,nchar(gene.table[k,]$seq)-200,nchar(gene.table[k,]$seq))
    py$run(seq2,species,filepath)
    gRNA.table2<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
    gRNA.table2<-table_select2(gRNA.table2)

    # gRNA.table1<-gRNA.table1[which(gRNA.table1$V1>abs(gene.table[k,]$up_dis)),]


    if(!is.null(gRNA.table1) & !is.null(gRNA.table2)){
      if(nrow(gRNA.table1)!=0 & nrow(gRNA.table2)!=0){
        gRNA<-data.frame(strand=c(gRNA.table1[1,]$V4,gRNA.table2[1,]$V4),gRNA=c(gRNA.table1[1,]$V2,gRNA.table2[1,]$V2),pos=c(gRNA.table1[1,]$V1,gRNA.table2[1,]$V1+(nchar(gene.table[k,]$seq)-200-1)),score=c(gRNA.table1[1,]$V3,gRNA.table2[1,]$V3))
      }
    }


    if(exists("gRNA")){
      if(nrow(gRNA)==2){
        gene.table[k,]$gRNA1<-gRNA[1,]$gRNA
        gene.table[k,]$strand1<-gRNA[1,]$strand
        gene.table[k,]$pos1<-gRNA[1,]$pos
        gene.table[k,]$score1<-gRNA[1,]$score
        gene.table[k,]$gRNA2<-gRNA[2,]$gRNA
        gene.table[k,]$strand2<-gRNA[2,]$strand
        gene.table[k,]$pos2<-gRNA[2,]$pos
        gene.table[k,]$score2<-gRNA[2,]$score

      }
      print(paste(k,"success"))
      write.csv(gene.table,"C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv",row.names = FALSE)
      rm(gRNA)
    }
    else{
      print(paste(k,"fail"))
    }
  }

  
  # if(nchar(gene.table[k,]$seq)>=200 & nchar(gene.table[k,]$seq)<400){
  #   
  #   py$run(gene.table[k,]$seq,species,filepath)
  #   gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  #   
  #   gRNA.table<-gRNA.table[which(!grepl("Inefficient", gRNA.table$V2)),]
  #   
  #   gRNA.table<-gRNA.table[which(as.numeric(gRNA.table$V3)>70 & as.numeric(gRNA.table$V8)>60),]
  #   if(nrow(gRNA.table)==0){
  #     print(paste(k,"fail"))
  #     next
  #   }
  #   
  #   gRNA.table$V2 <-gsub(" ", "", substring(gRNA.table$V2, 1, 24))
  #   strand <- sub("[^a-zA-Z]+", "", gRNA.table$V1)
  #   gRNA.table$V4<-strand
  #   pos <- as.numeric(sub("[^0-9]+", "", gRNA.table$V1))
  #   for(j in 1:length(strand)){
  #     if(strand[j]=="fw"){
  #       pos[j]=pos[j]+2
  #     }
  #     else if(strand[j]=="rev"){
  #       pos[j]=pos[j]+22
  #     }
  #   }
  #   gRNA.table$V1<-pos
  #   
  #   gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V1)<=100 & as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
  #   gRNA.table1<-gRNA.table1[which((nchar(gene.table[k,]$seq)-gRNA.table1$V1)>abs(gene.table[k,]$down_dis)),]
  #   if(nrow(gRNA.table1)==0){
  #     gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V1)<=200 & as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
  #     gRNA.table1<-gRNA.table1[which((nchar(gene.table[k,]$seq)-gRNA.table1$V1)>abs(gene.table[k,]$down_dis)),]
  #     if(nrow(gRNA.table1)==0){
  #       gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V1)<=100 & as.numeric(gRNA.table$V3)>70 & as.numeric(gRNA.table$V8)>60),]
  #       gRNA.table1<-gRNA.table1[which((nchar(gene.table[k,]$seq)-gRNA.table1$V1)>abs(gene.table[k,]$down_dis)),]
  #       if(nrow(gRNA.table1)==0){
  #         gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V1)<=200 & as.numeric(gRNA.table$V3)>70 & as.numeric(gRNA.table$V8)>60),]
  #         gRNA.table1<-gRNA.table1[which((nchar(gene.table[k,]$seq)-gRNA.table1$V1)>abs(gene.table[k,]$down_dis)),]
  #       }
  #     }
  #   }
  #   
  #   gRNA.table2<-gRNA.table[which(as.numeric(gRNA.table$V1)>=nchar(gene.table[k,]$seq)-100 & as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
  #   gRNA.table2<-gRNA.table2[which((nchar(gene.table[k,]$seq)-gRNA.table2$V1)>abs(gene.table[k,]$down_dis)),]
  #   if(nrow(gRNA.table2)==0){
  #     gRNA.table2<-gRNA.table[which(as.numeric(gRNA.table$V1)>=nchar(gene.table[k,]$seq)-200 & as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
  #     gRNA.table2<-gRNA.table2[which((nchar(gene.table[k,]$seq)-gRNA.table2$V1)>abs(gene.table[k,]$down_dis)),]
  #     if(nrow(gRNA.table2)==0){
  #       gRNA.table2<-gRNA.table[which(as.numeric(gRNA.table$V1)>=nchar(gene.table[k,]$seq)-100 & as.numeric(gRNA.table$V3)>70 & as.numeric(gRNA.table$V8)>60),]
  #       gRNA.table2<-gRNA.table2[which((nchar(gene.table[k,]$seq)-gRNA.table2$V1)>abs(gene.table[k,]$down_dis)),]
  #       if(nrow(gRNA.table2)==0){
  #         gRNA.table2<-gRNA.table[which(as.numeric(gRNA.table$V1)>=nchar(gene.table[k,]$seq)-200 & as.numeric(gRNA.table$V3)>70 & as.numeric(gRNA.table$V8)>60),]
  #         gRNA.table2<-gRNA.table2[which((nchar(gene.table[k,]$seq)-gRNA.table2$V1)>abs(gene.table[k,]$down_dis)),]
  #       }
  #     }
  #   }
  #   
  #   if(nrow(gRNA.table1)!=0 & nrow(gRNA.table2)!=0){
  #     gRNA<-data.frame(strand=c(gRNA.table1[1,]$V4,gRNA.table2[1,]$V4),gRNA=c(gRNA.table1[1,]$V2,gRNA.table2[1,]$V2),pos=c(gRNA.table1[1,]$V1,gRNA.table2[1,]$V1),score=c(gRNA.table1[1,]$V3,gRNA.table2[1,]$V3))
  #   }
  #   
  #   if(exists("gRNA")){
  #     if(nrow(gRNA)==2){
  #       gene.table[k,]$gRNA1<-gRNA[1,]$gRNA
  #       gene.table[k,]$strand1<-gRNA[1,]$strand
  #       gene.table[k,]$pos1<-gRNA[1,]$pos
  #       gene.table[k,]$score1<-gRNA[1,]$score
  #       gene.table[k,]$gRNA2<-gRNA[2,]$gRNA
  #       gene.table[k,]$strand2<-gRNA[2,]$strand
  #       gene.table[k,]$pos2<-gRNA[2,]$pos
  #       gene.table[k,]$score2<-gRNA[2,]$score
  #       
  #     }
  #     print(paste(k,"success"))
  #     write.csv(gene.table,"C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv",row.names = FALSE)
  #     rm(gRNA)
  #   }
  #   else{
  #     print(paste(k,"fail"))
  #   }
  # }
  

}


