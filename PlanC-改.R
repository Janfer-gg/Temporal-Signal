setwd("C://Users//41518//Desktop//work/ubigene")
#创建文件夹
filepath<-"C://Users//41518//Desktop//SPATA5"
dir.create(filepath)
library(ggplot2)
library(ggimage)
library(magick)
library(ggpubr)
library(httr)
library(jsonlite)
library(xml2)
library(rvest)
library(Biostrings)
library(dplyr)
library(reticulate)
library(stringr)
source("Get_ID.R")
source("Get_allinfo.R")
source("Get_seq.R")
source("Get_transcript_table.R")
source("delete_noprotein.R")
source("Get_max_transcript.R")
source("del_othergene.R")
source("Get_dot_region3.R")
source("Get_dot_region5.R")
source("Get_dot_region6.R")
source("Get_dot_region7.R")
source("Get_dot_region8.R")
source("Get_dot_region9.R")
source("Get_dot_region10.R")
source("Get_result1.R")
source("Get_result2.R")
source("Get_result3.R")
source("GC_analysis.R")
source("gRNA2_C.R")
source("Get_dot_plot.R")
source("Get_GC_image.R")
source("KO_longregion.R")
source("KO_region_300.R")
source("KO_region_image1.R")
source("KO_region_image2.R")
source("KO_region_image3.R")
source("Get_avoid_region.R")
source("Get_mark_region.R")

#进度条10%
write.table("1",paste0(filepath,"//","10%.txt"),row.names = FALSE,col.names = FALSE)

term <- ("SPATA5")
species<-"Human"
ID <- Get_ID(term,species)

# 获取信息 --------------------------------------------------------------------
Gene <- Get_seq(ID)                    #基因序列
Gene2 <- Get_seq2(ID)                 #5'和3'端各增加500bp
Gene3 <- Get_seq3(ID)
allinfo <- Get_allinfo(ID)
start <- allinfo$start
source_python("ensembl_table_download.py")
py$download_csv(ID)
transcript.table <- read.csv("transcript.csv",header = TRUE)

#如果要敲除的基因与其他基因有重叠
othergene<-Get_othergene(species,allinfo$seq_region_name,allinfo$start,allinfo$end)
othergene.table<-othergene[which(othergene$gene_id!=ID),]
avoid_region<-data.frame()
mark_region<-data.frame()
if(nrow(othergene.table)!=0){
  avoid_region<-Get_avoid_region(othergene.table)
  mark_region<-Get_mark_region(othergene.table)
}

#进度条20%
write.table("1",paste0(filepath,"//","20%.txt"),row.names = FALSE,col.names = FALSE)

print("ensembl数据调用成功")
color_blue<-character()
incomplete.transcript<-character()
for (i in 1:nrow(transcript.table)) {
  if (!grepl("Protein coding", transcript.table[i, ]$Biotype)) {
    color_blue <- append(color_blue,transcript.table[i,]$Name)
  }
  else if (grepl("incomplete", transcript.table[i, ]$Flags)) {
    incomplete.transcript <- append(incomplete.transcript, transcript.table[i,]$Name)
  }
}

transcript.table <-
  delete_noprotein(transcript.table)           #删除非编码蛋白和不完整的转录本
transcript.count <- nrow(transcript.table)     #转录本数量
transcript.name <-
  transcript.table[which(transcript.table$bp == max(transcript.table$bp)), ]$Name[1]     #最长的转录本
t_Exon_region <-
  Get_max_transcript(allinfo, transcript.table)         #最长转录本的外显子位置
t_Exon_region$Exon_start <- as.numeric(t_Exon_region$Exon_start)
t_Exon_region$Exon_end <- as.numeric(t_Exon_region$Exon_end)
t_Exon_region_sort <-
  t_Exon_region[order(t_Exon_region$Exon_start), ]

#是否为反向基因
Gene_rev <-
  all(t_Exon_region$Exon_start == sort(t_Exon_region$Exon_start, decreasing = TRUE))
#反向序列，将基因序列颠倒
if (Gene_rev) {
  gene <- DNAString(Gene)
  rev <- reverse(gene)
  Gene <- as.character(rev)
  gene2 <- DNAString(Gene2)
  rev2 <- reverse(gene2)
  Gene2 <- as.character(rev2)
  gene3 <- DNAString(Gene3)
  rev3 <- reverse(gene3)
  Gene3 <- as.character(rev3)
}



for (j in 1:length(allinfo$Transcript)) {
  if (allinfo$Transcript[[j]]$display_name == transcript.name) {
    if (Gene_rev) {
      t_CDS_start <- allinfo$Transcript[[j]]$Translation$end - start + 1
    }
    else{
      t_CDS_start <- allinfo$Transcript[[j]]$Translation$start - start + 1
    }
  }
}#最长转录本的CDS起始位置


#Ensembl上转录本的实心区域
ko.data <- data.frame()
name <- paste(transcript.table$Name, collapse = "")
for (j in 1:length(allinfo$Transcript)) {
  if (grepl(allinfo$Transcript[[j]]$display_name, name)) {
    Transcript <- allinfo$Transcript[[j]]
    if (length(Transcript$Translation) != 0) {
      #如果有编码区
      CDS_start <-
        Transcript$Translation$start - start + 1     #CDS起始位置
      CDS_end <-
        Transcript$Translation$end - start + 1       #CDS终止位置
      Transcript_start <-
        Transcript$start - start + 1           #转录本起始位置
      Transcript_end <-
        Transcript$end - start + 1
      Exon_info_CDS <- data.frame()            #实心外显子的位置
      k <- 1
      for (i in 1:length(Transcript$Exon)) {
        Exon_start <- Transcript$Exon[[i]]$start - start + 1
        Exon_end <- Transcript$Exon[[i]]$end - start + 1
        #完全在编码区的外显子
        if (Exon_start >= CDS_start & Exon_end <= CDS_end) {
          Exon_info_CDS[k, 1] <- paste("Exon", i)
          Exon_info_CDS[k, 2] <- Exon_start
          Exon_info_CDS[k, 3] <- Exon_end
          k <- k + 1
        }
        #部分在编码区(上游)
        else if (Exon_start <= CDS_start &
                 Exon_end >= CDS_start & Exon_end <= CDS_end) {
          Exon_info_CDS[k, 1] <- paste("Exon", i)
          Exon_info_CDS[k, 2] <- CDS_start
          Exon_info_CDS[k, 3] <- Exon_end
          k <- k + 1
        }
        #部分在编码区(下游)
        else if (Exon_start >= CDS_start &
                 Exon_start <= CDS_end & Exon_end >= CDS_end) {
          Exon_info_CDS[k, 1] <- paste("Exon", i)
          Exon_info_CDS[k, 2] <- Exon_start
          Exon_info_CDS[k, 3] <- CDS_end
          k <- k + 1
        }
        else if (Exon_start < CDS_start & Exon_end > CDS_end) {
          Exon_info_CDS[k, 1] <- paste("Exon", i)
          Exon_info_CDS[k, 2] <- CDS_start
          Exon_info_CDS[k, 3] <- CDS_end
          k <- k + 1
        }
      }
      transcript <-
        rep(allinfo$Transcript[[j]]$display_name,
            length(Exon_info_CDS[, 1]))
      names(Exon_info_CDS) = c("Exon", "start", "end")
      ko_info <- cbind(transcript, Exon_info_CDS)
      ko.data <- rbind(ko.data, ko_info)
    }
  }
}
t_Exon_CDS <-ko.data[which(ko.data$transcript == transcript.name),]
ATG_Exon<-t_Exon_CDS[1,]$Exon
stop_Exon<-t_Exon_CDS[nrow(t_Exon_CDS),]$Exon

KO_region2 <- KO_longregion(t_Exon_CDS)
print(KO_region2)

if(nrow(KO_region2)!=0) {
  
  #如果与其他基因重叠
  if(nrow(avoid_region)!=0){
    avoid_ko_region<-data.frame()
    avoid_ko_del<-numeric()
    for (i in 1:nrow(KO_region2)) {
      for(j in 1:nrow(avoid_region)){
        #如果重叠
        if(!(KO_region2[i,]$end<avoid_region[j,]$start | KO_region2[i,]$start>avoid_region[j,]$end)){
          avoid_ko_del<-append(avoid_ko_del,i)
          avoid_ko_region<-rbind(avoid_ko_region,KO_region2[i,])
        }
      }
    }
    if(length(avoid_ko_del)!=0){
      KO_region2 <- KO_region2[-avoid_ko_del, ]
    }
  }
  #GC含量分析:平均GC含量大于70%或小于30%，则删除该区域
  GC_del <- numeric()
  for (i in 1:nrow(KO_region2)) {
    analysis_GC <- GC_analysis1(KO_region2[i,])
    if (analysis_GC == TRUE) {
      GC_del <- append(GC_del, i)
    }
  }
  if (length(GC_del) != 0) {
    KO_region2 <- KO_region2[-GC_del, ]
  }
  
  #点阵图分析
  dot_del <- numeric()
  for (i in 1:nrow(KO_region2)) {
    analysis_dot5 <- Get_dot_region5(KO_region2[i,])
    analysis_dot6 <- Get_dot_region6(KO_region2[i,])
    analysis_dot7 <- Get_dot_region7(KO_region2[i,])
    analysis_dot8 <- Get_dot_region8(KO_region2[i,])
    analysis_dot9 <- Get_dot_region9(KO_region2[i,])
    analysis_dot10 <- Get_dot_region10(KO_region2[i,])
    analysis_dot11 <- Get_dot_region11(KO_region2[i,])
    analysis_dot12 <- Get_dot_region12(KO_region2[i,])
    if (analysis_dot7 == TRUE | analysis_dot8 == TRUE |
        analysis_dot9 == TRUE | analysis_dot5 == TRUE |
        analysis_dot6 == TRUE | analysis_dot10 == TRUE |
        analysis_dot11 == TRUE | analysis_dot12 == TRUE) {
      dot_del <- append(dot_del, i)
    }
  }
  if (length(dot_del) != 0) {
    KO_region2 <- KO_region2[-dot_del,]
  }
}
print(KO_region2)
#进度条50%
write.table("1",paste0(filepath,"//","50%.txt"),row.names = FALSE,col.names = FALSE)

# 整个敲除 --------------------------------------------------------------------
ko_seq<-data.frame(seq1=character(),seq2=character())
if (nrow(KO_region2) != 0) {
  for(t in 1:nrow(KO_region2)){
    # 长度不够 --------------------------------------------------------------------
    {
      if (KO_region2[t, ]$start < 1200 |
          KO_region2[t, ]$end + 1200 > nchar(Gene)) {
        #外显子上游
        {
          if(KO_region2[t, ]$start==min(t_Exon_CDS$start)){
            judge<-"TRUE"
          }
          else if(KO_region2[t, ]$start!=min(t_Exon_CDS$start) & t_Exon_region_sort[which(t_Exon_region_sort$Exon_start==KO_region2[t, ]$start)-1,]$Exon_end+1000<=KO_region2[t, ]$start ){
            judge<-"TRUE"
          }
          else{
            judge<-"FALSE"
          }
        }
        {
          #800bp
          if (judge) {
            ko_start1 <- KO_region2[t, ]$start + 1200 - 800
            ko_end1 <- KO_region2[t, ]$start + 1200
            ko_seq1 <- substring(Gene3, ko_start1, ko_end1)
            if (Gene_rev) {
              seq1 <- DNAString(ko_seq1)
              seq1_rev <- reverse(seq1)
              ko_seq1 <- as.character(seq1_rev)
            }
            ko_seq[t, 1] <- ko_seq1
          }
          #400bp
          else{
            ko_start1 <- KO_region2[t, ]$start + 1200 - 400
            ko_end1 <- KO_region2[t, ]$start + 1200
            ko_seq1 <- substring(Gene3, ko_start1, ko_end1)
            if (Gene_rev) {
              seq1 <- DNAString(ko_seq1)
              seq1_rev <- reverse(seq1)
              ko_seq1 <- as.character(seq1_rev)
            }
            ko_seq[t,1]<-ko_seq1
          }
        }
        rm(judge)
        #外显子下游
        #800bp
        {
          if(KO_region2[t, ]$end==max(t_Exon_CDS$end)){
            judge<-"TRUE"
          }
          else if(KO_region2[t, ]$end!=max(t_Exon_CDS$end) & t_Exon_region_sort[which(t_Exon_region_sort$Exon_end==KO_region2[t, ]$end)+1,]$Exon_start-1000>=KO_region2[t, ]$end ){
            judge<-"TRUE"
          }
          else{
            judge<-"FALSE"
          }
        }
        {
          if (judge) {
            ko_start2 <- KO_region2[t, ]$end + 1200
            ko_end2 <- KO_region2[t, ]$end + 1200 + 800
            ko_seq2 <- substring(Gene3, ko_start2, ko_end2)
            if (Gene_rev) {
              seq2 <- DNAString(ko_seq2)
              seq2_rev <- reverse(seq2)
              ko_seq2 <- as.character(seq2_rev)
            }
            ko_seq[t, 2] <- ko_seq2
          }
          #400bp
          else{
            ko_start2 <- KO_region2[t, ]$end + 1200
            ko_end2 <- KO_region2[t, ]$end + 1200 + 400
            ko_seq2 <- substring(Gene3, ko_start2, ko_end2)
            if (Gene_rev) {
              seq2 <- DNAString(ko_seq2)
              seq2_rev <- reverse(seq2)
              ko_seq2 <- as.character(seq2_rev)
            }
            ko_seq[t,2]<-ko_seq2
          }
        }
      }
      

# 长度够 ---------------------------------------------------------------------
      else{
        #外显子上游
        #800bp
        {
          if(KO_region2[t, ]$start==min(t_Exon_CDS$start)){
            judge<-"TRUE"
          }
          else if(KO_region2[t, ]$start!=min(t_Exon_CDS$start) & t_Exon_region_sort[which(t_Exon_region_sort$Exon_start==KO_region2[t, ]$start)-1,]$Exon_end+1000<=KO_region2[t, ]$start ){
            judge<-"TRUE"
          }
          else{
            judge<-"FALSE"
          }
        }
        {
          if (judge) {
            ko_start1 <- KO_region2[t, ]$start - 800
            ko_end1 <- KO_region2[t, ]$start
            ko_seq1 <- substring(Gene, ko_start1, ko_end1)
            if (Gene_rev) {
              seq1 <- DNAString(ko_seq1)
              seq1_rev <- reverse(seq1)
              ko_seq1 <- as.character(seq1_rev)
            }
            ko_seq[t, 1] <- ko_seq1
          }
          else{
            ko_start1 <- KO_region2[t, ]$start - 400
            ko_end1 <- KO_region2[t, ]$start
            ko_seq1 <- substring(Gene, ko_start1, ko_end1)
            if (Gene_rev) {
              seq1 <- DNAString(ko_seq1)
              seq1_rev <- reverse(seq1)
              ko_seq1 <- as.character(seq1_rev)
            }
            ko_seq[t,1]<-ko_seq1
          }
        }
        
        #外显子下游
        #800bp
        {
          if(KO_region2[t, ]$end==max(t_Exon_CDS$end)){
            judge<-"TRUE"
          }
          else if(KO_region2[t, ]$end!=max(t_Exon_CDS$end) & t_Exon_region_sort[which(t_Exon_region_sort$Exon_end==KO_region2[t, ]$end)+1,]$Exon_start-1000>=KO_region2[t, ]$end ){
            judge<-"TRUE"
          }
          else{
            judge<-"FALSE"
          }
        }
        {
          if (judge) {
            ko_start2 <- KO_region2[t, ]$end
            ko_end2 <- KO_region2[t, ]$end + 800
            ko_seq2 <- substring(Gene, ko_start2, ko_end2)
            if (Gene_rev) {
              seq2 <- DNAString(ko_seq2)
              seq2_rev <- reverse(seq2)
              ko_seq4 <- as.character(seq2_rev)
            }
            ko_seq[t, 2] <- ko_seq2
          }
          else{
            ko_start2 <- KO_region2[t,]$end
            ko_end2 <- KO_region2[t,]$end + 400
            ko_seq2 <- substring(Gene, ko_start2, ko_end2)
            if (Gene_rev) {
              seq2 <- DNAString(ko_seq2)
              seq2_rev <- reverse(seq2)
              ko_seq2 <- as.character(seq2_rev)
            }
            ko_seq[t,2]<-ko_seq2
          }
        }
      }
    }
  }
}

temp<-list.files(paste0(filepath,"//","gRNA"))
if (nrow(KO_region2) != 0) {
  for(t in 1:length(temp)){
    {
      #上游
      gRNA.table1 <- read.csv(paste0(filepath, "//gRNA//", temp[t], "//gRNA1.csv"),header = FALSE,encoding = "UTF-8")
      gRNA.table2 <- read.csv(paste0(filepath, "//gRNA//", temp[t], "//gRNA2.csv"),header = FALSE,encoding = "UTF-8")
      #特异性得分70下，和Inefficient的gRNA排除掉
      gRNA.del <- numeric()
      for (i in 1:length(gRNA.table1[, 1])) {
        if (gRNA.table1[i, 3] < 70) {
          gRNA.del <- append(gRNA.del, i)
        }
        if (gRNA.table1[i, 3] == "No matches") {
          gRNA.del <- append(gRNA.del, i)
        }
        else if (grepl("Inefficient", gRNA.table1[i, 2])) {
          gRNA.del <- append(gRNA.del, i)
        }
      }
      gRNA.table1 <- gRNA.table1[-gRNA.del, ]
      if(nrow(gRNA.table1)==0){
        next
      }
      #0-0-0(优化)
      count_0 <- numeric()
      for (p in 1:nrow(gRNA.table1)) {
        target <-
          str_extract_all(gRNA.table1[p, ]$V9, "\\d+\\s\\-\\s\\d+\\s\\-\\s\\d+")[[1]][2]
        target <- gsub("\\s", "", target)
        if (target == "0-0-0") {
          count_0 <- append(count_0, p)
        }
      }
      gRNA.table_min <- gRNA.table1[-count_0, ]
      gRNA.table1 <- rbind(gRNA.table1[count_0, ], gRNA.table_min)
      if (nrow(gRNA.table1) > 10) {
        gRNA.table1 <- gRNA.table1[1:10,]
      }
      
      #下游
      #特异性得分70下，和Inefficient的gRNA排除掉
      gRNA.del <- numeric()
      for (i in 1:length(gRNA.table2[, 1])) {
        if (gRNA.table2[i, 3] < 70) {
          gRNA.del <- append(gRNA.del, i)
        }
        if (gRNA.table2[i, 3] == "No matches") {
          gRNA.del <- append(gRNA.del, i)
        }
        else if (grepl("Inefficient", gRNA.table2[i, 2])) {
          gRNA.del <- append(gRNA.del, i)
        }
      }
      gRNA.table2 <- gRNA.table2[-gRNA.del, ]
      
      if(nrow(gRNA.table2)==0){
        next
      }
      
      #0-0-0(优化)
      count_0 <- numeric()
      for (p in 1:nrow(gRNA.table2)) {
        target <-
          str_extract_all(gRNA.table2[p, ]$V9, "\\d+\\s\\-\\s\\d+\\s\\-\\s\\d+")[[1]][2]
        target <- gsub("\\s", "", target)
        if (target == "0-0-0") {
          count_0 <- append(count_0, p)
        }
      }
      gRNA.table_min <- gRNA.table2[-count_0, ]
      gRNA.table2 <- rbind(gRNA.table2[count_0, ], gRNA.table_min)
      
      if (nrow(gRNA.table2) > 10) {
        gRNA.table2 <- gRNA.table2[1:10,]
      }
      if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
        #筛选不到gRNA时要及时退出
        next
      }
      
      #长度不够
      if (KO_region2[t, ]$start < 1200 | KO_region2[t, ]$end + 1200 > nchar(Gene)) {
        KO_region2[t, ]$start <- KO_region2[t, ]$start + 1200
        KO_region2[t, ]$end <- KO_region2[t, ]$end + 1200
        #上下游的gRNA合并
        gRNA.table <- rbind(gRNA.table1, gRNA.table2)
        #获取每个gRNA在基因上的位置
        strand <- sub("[^a-zA-Z]+", "", gRNA.table[, 1])
        gRNA_seq <- substring(gRNA.table[, 2], 1, 20)
        Score1 <- gRNA.table[, 3]
        analysis_seq <-
          gsub(" ", "", substring(gRNA.table[, 2], 1, 24))
        gRNA.table <- cbind(strand, gRNA_seq, analysis_seq, Score1)
        gRNA.table <- as.data.frame(gRNA.table)
        gene <- DNAString(Gene3)
        {
          if (Gene_rev) {
            for (i in 1:length(gRNA.table[, 1])) {
              if (gRNA.table[i,]$strand == "fw") {
                gRNA_rev <- DNAString(gRNA.table[i,]$gRNA_seq)
                gRNA_rev <- reverse(gRNA_rev)
                rev <-
                  matchPattern(pattern = gRNA_rev, subject = gene)
                rev_start <- start(rev)
                rev_end <- end(rev)
                gRNA.table[i, 5] <- rev_start
                gRNA.table[i, 6] <- rev_end
              }
              else if (gRNA.table[i,]$strand == "rev") {
                gRNA_fw <- DNAString(gRNA.table[i,]$gRNA_seq)
                gRNA_fw <- complement(gRNA_fw)
                fw <-
                  matchPattern(pattern = gRNA_fw, subject = gene)
                fw_start <- start(fw)
                fw_end <- end(fw)
                gRNA.table[i, 5] <- fw_start
                gRNA.table[i, 6] <- fw_end
              }
            }
            names(gRNA.table)[5:6] <- c("start", "end")
            print(gRNA.table)
          }
          
          else{
            for (i in 1:length(gRNA.table[, 1])) {
              if (gRNA.table[i, ]$strand == "rev") {
                gRNA_rev <- DNAString(gRNA.table[i, ]$gRNA_seq)
                gRNA_rev <- reverseComplement(gRNA_rev)
                rev <-
                  matchPattern(pattern = gRNA_rev, subject = gene)
                rev_start <- start(rev)
                rev_end <- end(rev)
                gRNA.table[i, 5] <- rev_start[1]
                gRNA.table[i, 6] <- rev_end[1]
              }
              else if (gRNA.table[i, ]$strand == "fw") {
                fw <-
                  matchPattern(pattern = gRNA.table[i, ]$gRNA_seq,
                               subject = gene)
                fw_start <- start(fw)
                fw_end <- end(fw)
                gRNA.table[i, 5] <- fw_start
                gRNA.table[i, 6] <- fw_end
              }
            }
            names(gRNA.table)[5:6] <- c("start", "end")
            print(gRNA.table)
          }
        }
        #局部GC含量大于80%或小于25%，避免该区域
        GC_avoid_region <- GC_analysis2(KO_region2[t,])
        if (GC_avoid_region != FALSE) {
          GC_del <- numeric()
          for (i in 1:nrow(gRNA.table)) {
            for (j in 1:nrow(GC_avoid_region)) {
              if (gRNA.table[i,]$start %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2]) |
                  gRNA.table[i,]$end %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2])) {
                GC_del <- append(GC_del, i)
              }
            }
          }
          if (length(GC_del) != 0) {
            gRNA.table <- gRNA.table[-GC_del, ]
          }
        }
        
        #符合条件的gRNA进行切割效率预测
        write.csv(gRNA.table, file = "CCTOP-predictor.csv", row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer("CCTOP-predictor.csv",species)
        gRNA.table <- read.csv("CCTOP-predictor.csv", header = TRUE)
        print(gRNA.table)
        #切割效率得分低于0.60的删除
        gRNA.table <-gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
        
        #上下游分开
        gRNA.table1 <-gRNA.table[which(gRNA.table$end <= KO_region2[t,]$start),]
        gRNA.table2 <-gRNA.table[which(gRNA.table$start >= KO_region2[t,]$end),]
        
        #切割得分大于0.65的优先
        gRNA.table1 <-rbind(gRNA.table1[which(gRNA.table1$crispr_score >= 0.65), ], gRNA.table1[which(gRNA.table1$crispr_score <0.65), ])
        gRNA.table2 <-rbind(gRNA.table2[which(gRNA.table2$crispr_score >= 0.65), ], gRNA.table2[which(gRNA.table2$crispr_score <0.65), ])
        
        #重新合并
        gRNA.table <- rbind(gRNA.table1, gRNA.table2)
        
        #table1是外显子上游的gRNA,table2是外显子下游的gRNA
        {
          if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
            rm(Gene3)
            next
          }
          else{
            gRNA_planC <- Get_result1(gRNA.table, KO_region2[t,])        #相差0.05分以内优选
          }
        }
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA_planC) == "NULL") {
          gRNA_planC <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
        }
        
        #第二对gRNA
        if (class(gRNA_planC) != "NULL" &
            nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
          #在gRNA表格中删除第一对gRNA
          n = which(gRNA.table$analysis_seq %in% gRNA_planC$analysis_seq)
          gRNA.table_2 <- gRNA.table[-n,]
          #删除重叠的
          gRNA.table_3 <-
            gRNA.table_2[which(
              abs(gRNA.table_2$start - gRNA_planC[1, ]$start) >= 20 &
                abs(gRNA.table_2$start - gRNA_planC[2, ]$start) >= 20
            ), ]
          #上下游分开
          gRNA.table1_2 <-
            gRNA.table_3[which(gRNA.table_3$end <= KO_region2[t,]$start),]
          gRNA.table2_2 <-
            gRNA.table_3[which(gRNA.table_3$start >= KO_region2[t,]$end),]
          if (nrow(gRNA.table1_2) != 0 & nrow(gRNA.table2_2) != 0) {
            #相差0.05分以内优选
            gRNA2_planC <- Get_result1(gRNA.table_3, KO_region2[t, ])
            #如果没有相差0.05分以内的gRNA
            if (class(gRNA2_planC) == "NULL") {
              gRNA2_planC <- rbind(gRNA.table1_2[1, ], gRNA.table2_2[1, ])
            }
          }
        }
        if(exists("gRNA_planC")==TRUE){
          judge3<-"TRUE"
        }   
      }
      #长度够
      else{
  
      }
    }
  }
}


