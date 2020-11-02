setwd("C://Users//41518//Desktop//work/ubigene")
#创建文件夹
filepath<-"C://Users//41518//Desktop//靶位点测试3//GREM1"
dir.create(filepath)
library(ggplot2)
library(ggimage)
library(magick)
library(ggpubr)
library(httr)
library(jsonlite)
library(xml2)
library(Biostrings)
library(reticulate)
library(stringr)
source("Get_ID.R")
source("Get_allinfo.R")
source("Get_seq.R")
source("Get_transcript_table.R")
source("delete_noprotein.R")
source("Get_max_transcript.R")
source("del_othergene.R")
source("Get_dot_region5.R")
source("Get_dot_region6.R")
source("Get_dot_region7.R")
source("Get_dot_region8.R")
source("Get_dot_region9.R")
source("Get_dot_region10.R")
source("Get_result1.R")
source("Get_result2.R")
source("Get_result3.R")
source("GC_analysis_C.R")
source("gRNA2_C.R")
source("Get_dot_plot.R")
source("Get_GC_image.R")
source("KO_longregion.R")
source("KO_region_image3.R")
source("Get_avoid_region.R")
source("Get_mark_region.R")
source("dot_analysis_C.R")

#进度条10%
write.table("1",paste0(filepath,"//","C10%.txt"),row.names = FALSE,col.names = FALSE)

term <- ("GREM1")
species<-"Human"
ID <- Get_ID(term,species)

# 获取信息 --------------------------------------------------------------------
Gene <- Get_seq(ID)                    #基因序列
Gene2 <- Get_seq2(ID)                 #5'和3'端各增加500bp
allinfo <- Get_allinfo(ID)
start <- allinfo$start
#如果要敲除的基因与其他基因有重叠
othergene<-Get_othergene(species,allinfo$seq_region_name,allinfo$start,allinfo$end)
othergene.table<-othergene[which(othergene$Parent!=ID),]
avoid_region <- data.frame(start = numeric(), end = numeric())
mark_region<-data.frame()
if(nrow(othergene.table)!=0){
  avoid_region<-Get_avoid_region(othergene.table)
  mark_region<-Get_mark_region(othergene.table)
}

source_python("ensembl_table_download.py")
py$download_csv(filepath,ID)
transcript.table <- read.csv(paste0(filepath,"//transcript.csv"),header = TRUE)

#进度条20%
write.table("1",paste0(filepath,"//","C20%.txt"),row.names = FALSE,col.names = FALSE)

print("ensembl success")
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
t_Exon_CDS_sort <- t_Exon_CDS[order(t_Exon_CDS$start),]
ATG_Exon<-t_Exon_CDS[1,]$Exon
stop_Exon<-t_Exon_CDS[nrow(t_Exon_CDS),]$Exon

# 整个敲除 --------------------------------------------------------------------
KO_region2 <- KO_longregion(t_Exon_CDS)

#排序
KO_region2<-rbind(KO_region2[which(KO_region2$end-KO_region2$start<=10000),],KO_region2[which(KO_region2$end-KO_region2$start>10000),])
print(KO_region2)
if(nrow(KO_region2)==0){
  write.table("4",paste0(filepath,"//","result3.txt"),row.names = FALSE,col.names = FALSE)
}

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
    if(nrow(KO_region2)==0){
      write.table("1",paste0(filepath,"//","result3.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}
if(nrow(KO_region2)!=0) {
  #GC含量分析:平均GC含量大于70%或小于30%，则删除该区域
  GC_del <- numeric()
  for (i in 1:nrow(KO_region2)) {
    analysis_GC1 <- GC_analysis1(KO_region2[i,])
    analysis_GC2 <- GC_analysis2(KO_region2[i,])
    if (analysis_GC1 == TRUE | analysis_GC2==TRUE) {
      GC_del <- append(GC_del, i)
    }
  }
  if (length(GC_del) != 0) {
    KO_region2 <- KO_region2[-GC_del, ]
  }
  if(nrow(KO_region2)==0){
    write.table("2",paste0(filepath,"//","result3.txt"),row.names = FALSE,col.names = FALSE)
  }
}

if(nrow(KO_region2)!=0) {
  #点阵图分析
  dot_del <- numeric()
  for (i in 1:nrow(KO_region2)) {
    analysis_dot5 <- Get_dot_region5(KO_region2[i,])
    analysis_dot6 <- Get_dot_region6(KO_region2[i,])
    analysis_dot7 <- Get_dot_region7(KO_region2[i,])
    analysis_dot8 <- Get_dot_region8(KO_region2[i,])
    analysis_dot9 <- Get_dot_region9(KO_region2[i,])
    analysis_dot11 <- Get_dot_region11(KO_region2[i,])
    if (analysis_dot7 == TRUE | analysis_dot8 == TRUE |
        analysis_dot9 == TRUE | analysis_dot5 == TRUE |
        analysis_dot6 == TRUE | analysis_dot11 == TRUE ) {
      dot_del <- append(dot_del, i)
    }
  }
  if (length(dot_del) != 0) {
    KO_region2 <- KO_region2[-dot_del,]
  }
  if(nrow(KO_region2)==0){
    write.table("2",paste0(filepath,"//","result3.txt"),row.names = FALSE,col.names = FALSE)
  }
}
print(KO_region2)
#进度条60%
write.table("1",paste0(filepath,"//","C60%.txt"),row.names = FALSE,col.names = FALSE)

if (nrow(KO_region2) != 0) {
  for(t in 1:nrow(KO_region2)){
    # 长度不够 --------------------------------------------------------------------
    if(KO_region2[t,]$start<1200 | KO_region2[t,]$end+1200>nchar(Gene)){
      Gene3 <- Get_seq3(ID)
      if (Gene_rev) {
        gene3 <- DNAString(Gene3)
        rev3 <- reverse(gene3)
        Gene3 <- as.character(rev3)
      }
      KO_region2[t,]$start<-KO_region2[t,]$start+1200
      KO_region2[t,]$end<-KO_region2[t,]$end+1200
      ko_start1 <- KO_region2[t,]$start -400
      ko_end1 <- KO_region2[t,]$start
      ko_seq1 <- substring(Gene3, ko_start1, ko_end1)
      if (Gene_rev) {
        seq1 <- DNAString(ko_seq1)
        seq1_rev <- reverse(seq1)
        ko_seq1 <- as.character(seq1_rev)
      }
      source_python("crispor_table_download.py")
      py$run(ko_seq1, species,filepath,"3")
      gRNA.table1 <-
        read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
      #特异性得分60下，和Inefficient的gRNA排除掉
      gRNA.del <- numeric()
      gRNA.table1<-gRNA.table1[which(gRNA.table1$V3!="No matches"),]
      for (i in 1:length(gRNA.table1[, 1])) {
        if (as.numeric(gRNA.table1[i, 3]) < 60) {
          gRNA.del <- append(gRNA.del, i)
        }
        else if (grepl("Inefficient", gRNA.table1[i, 2])) {
          gRNA.del <- append(gRNA.del, i)
        }
      }
      if(length(gRNA.del)!=0){
        gRNA.table1 <- gRNA.table1[-gRNA.del, ]
      }
      if (nrow(gRNA.table1) == 0) {
        #往前400bp
        {
          if (KO_region2[t,]$start - 1200 == min(t_Exon_CDS$start)) {
            judge <- "TRUE"
          }
          else if ((KO_region2[t,]$start - 1200) != min(t_Exon_CDS$start) &t_Exon_region_sort[which(t_Exon_region_sort$Exon_start == (KO_region2[t,]$start -1200)) - 1, ]$Exon_end + 1000 <= (KO_region2[t,]$start - 1200)) {
            judge <- "TRUE"
          }
          else{
            judge <- "FALSE"
          }
        }
        if (judge) {
          ko_start1 <- KO_region2[t, ]$start - 800
          ko_end1 <- KO_region2[t, ]$start - 400
          ko_seq1 <- substring(Gene, ko_start1, ko_end1)
          if (Gene_rev) {
            seq1 <- DNAString(ko_seq1)
            seq1_rev <- reverse(seq1)
            ko_seq1 <- as.character(seq1_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq1, species, filepath,"3")
          gRNA.table1 <-read.csv(paste0(filepath, "//gRNA3.csv"),header = FALSE,encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table1 <-gRNA.table1[which(gRNA.table1$V3 != "No matches"), ]
          for (i in 1:length(gRNA.table1[, 1])) {
            if (as.numeric(gRNA.table1[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table1[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if (length(gRNA.del) != 0) {
            gRNA.table1 <- gRNA.table1[-gRNA.del,]
          }
          if (nrow(gRNA.table1) == 0) {
            rm(Gene3)
            next
          }
        }
        else{
          rm(Gene3)
          next
        }
      }
      #0-0-0(优化)
      count_0 <- numeric()
      for (p in 1:nrow(gRNA.table1)) {
        if (gRNA.table1[p,]$V9 == "0-0-0") {
          count_0 <- append(count_0, p)
        }
      }
      if (length(count_0) != 0) {
        gRNA.table_min <- gRNA.table1[-count_0, ]
        gRNA.table1 <- rbind(gRNA.table1[count_0, ], gRNA.table_min)
      }
      if (nrow(gRNA.table1) > 10) {
        gRNA.table1 <- gRNA.table1[1:10,]
      }
      
      #外显子下游
      ko_start2 <- KO_region2[t,]$end 
      ko_end2 <- KO_region2[t,]$end + 400
      ko_seq2 <- substring(Gene3, ko_start2, ko_end2)
      if (Gene_rev) {
        seq2 <- DNAString(ko_seq2)
        seq2_rev <- reverse(seq2)
        ko_seq2 <- as.character(seq2_rev)
      }
      #读取ko_seq的gRNA表格
      source_python("crispor_table_download.py")
      py$run(ko_seq2, species,filepath,"3")
      gRNA.table2 <-
        read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
      #特异性得分60下，和Inefficient的gRNA排除掉
      gRNA.del <- numeric()
      gRNA.table2<-gRNA.table2[which(gRNA.table2$V3!="No matches"),]
      for (i in 1:length(gRNA.table2[, 1])) {
        if (as.numeric(gRNA.table2[i, 3]) < 60) {
          gRNA.del <- append(gRNA.del, i)
        }
        else if (grepl("Inefficient", gRNA.table2[i, 2])) {
          gRNA.del <- append(gRNA.del, i)
        }
      }
      if(length(gRNA.del)!=0){
        gRNA.table2 <- gRNA.table2[-gRNA.del, ]
      }
      
      if(nrow(gRNA.table2)==0){
        {
          if((KO_region2[t, ]$end-1200)==max(t_Exon_CDS$end)){
            judge<-"TRUE"
          }
          else if((KO_region2[t, ]$end-1200)!=max(t_Exon_CDS$end) & t_Exon_region_sort[which(t_Exon_region_sort$Exon_end==(KO_region2[t, ]$end-1200))+1,]$Exon_start-1000>=(KO_region2[t, ]$end-1200) ){
            judge<-"TRUE"
          }
          else{
            judge<-"FALSE"
          }
        }
        if(judge){
          ko_start2 <- KO_region2[t,]$end +400
          ko_end2 <- KO_region2[t,]$end + 800
          ko_seq2 <- substring(Gene, ko_start2, ko_end2)
          if (Gene_rev) {
            seq2 <- DNAString(ko_seq2)
            seq2_rev <- reverse(seq2)
            ko_seq2 <- as.character(seq2_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq2, species,filepath,"3")
          gRNA.table2 <-read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table2<-gRNA.table2[which(gRNA.table2$V3!="No matches"),]
          for (i in 1:length(gRNA.table2[, 1])) {
            if (as.numeric(gRNA.table2[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table2[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table2 <- gRNA.table2[-gRNA.del, ]
          }
          if(nrow(gRNA.table2)==0){
            rm(Gene3)
            next
          }
        }
        else{
          rm(Gene3)
          next
        }
      }
      
      #0-0-0(优化)
      count_0 <- numeric()
      for (p in 1:nrow(gRNA.table2)) {
        if (gRNA.table2[p,]$V9 == "0-0-0") {
          count_0 <- append(count_0, p)
        }
      }
      if (length(count_0) != 0) {
        gRNA.table_min <- gRNA.table2[-count_0, ]
        gRNA.table2 <- rbind(gRNA.table2[count_0, ], gRNA.table_min)
      }
      if (nrow(gRNA.table2) > 10) {
        gRNA.table2 <- gRNA.table2[1:10,]
      }
      if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
        #筛选不到gRNA时要及时退出
        rm(Gene3)
        next
      }
      
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
      GC_avoid_region <- GC_analysis4(KO_region2[t,],Gene3)
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
      write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor3.csv"), row.names = FALSE)
      source_python("crispr_get_score.py")
      py$reader_writer(paste0(filepath,"//CCTOP-predictor3.csv"),species)
      gRNA.table <- read.csv(paste0(filepath,"//CCTOP-predictor3.csv"), header = TRUE)
      print(gRNA.table)
      #切割效率得分低于0.60的删除
      gRNA.table <-
        gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
      
      #上下游分开
      gRNA.table1 <-
        gRNA.table[which(gRNA.table$end <= KO_region2[t,]$start),]
      gRNA.table2 <-
        gRNA.table[which(gRNA.table$start >= KO_region2[t,]$end),]
      
      #切割得分大于0.65的优先
      gRNA.table1 <-
        rbind(gRNA.table1[which(gRNA.table1$crispr_score >= 0.65), ], gRNA.table1[which(gRNA.table1$crispr_score <
                                                                                          0.65), ])
      gRNA.table2 <-
        rbind(gRNA.table2[which(gRNA.table2$crispr_score >= 0.65), ], gRNA.table2[which(gRNA.table2$crispr_score <
                                                                                          0.65), ])
      
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
        {
          if (nrow(gRNA.table1_2) != 0 & nrow(gRNA.table2_2) != 0) {
            #相差0.05分以内优选
            gRNA2_planC <- Get_result1(gRNA.table_3, KO_region2[t,])
            #如果没有相差0.05分以内的gRNA
            if (class(gRNA2_planC) == "NULL") {
              gRNA2_planC <- rbind(gRNA.table1_2[1,], gRNA.table2_2[1,])
            }
          }
        }
      }
    }
    
    # 长度够 ---------------------------------------------------------------------
    else{
      ko_start1 <- KO_region2[t,]$start - 400
      ko_end1 <- KO_region2[t,]$start
      # if (ko_start1 < 0) {
      #   KO_region2[t,]$start <- min(t_Exon_CDS$start)
      # }
      # ko_start1 <- KO_region2[t,]$start - 400
      # ko_end1 <- KO_region2[t,]$start
      
      ko_seq1 <- substring(Gene, ko_start1, ko_end1)
      if (Gene_rev) {
        seq1 <- DNAString(ko_seq1)
        seq1_rev <- reverse(seq1)
        ko_seq1 <- as.character(seq1_rev)
      }
      #读取ko_seq的gRNA表格
      source_python("crispor_table_download.py")
      py$run(ko_seq1, species,filepath,"3")
      gRNA.table1 <-
        read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
      #特异性得分60下，和Inefficient的gRNA排除掉
      gRNA.del <- numeric()
      gRNA.table1<-gRNA.table1[which(gRNA.table1$V3!="No matches"),]
      for (i in 1:length(gRNA.table1[, 1])) {
        if (as.numeric(gRNA.table1[i, 3]) < 60) {
          gRNA.del <- append(gRNA.del, i)
        }
        else if (grepl("Inefficient", gRNA.table1[i, 2])) {
          gRNA.del <- append(gRNA.del, i)
        }
      }
      if(length(gRNA.del)!=0){
        gRNA.table1 <- gRNA.table1[-gRNA.del, ]
      }
      
      if(nrow(gRNA.table1)==0){
        #往前400bp
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
        if(judge){
          ko_start1 <- KO_region2[t,]$start - 800
          ko_end1 <- KO_region2[t,]$start - 400
          ko_seq1 <- substring(Gene, ko_start1, ko_end1)
          if (Gene_rev) {
            seq1 <- DNAString(ko_seq1)
            seq1_rev <- reverse(seq1)
            ko_seq1 <- as.character(seq1_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq1, species,filepath,"3")
          gRNA.table1 <-
            read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table1<-gRNA.table1[which(gRNA.table1$V3!="No matches"),]
          for (i in 1:length(gRNA.table1[, 1])) {
            if (as.numeric(gRNA.table1[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table1[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table1 <- gRNA.table1[-gRNA.del, ]
          }
          if(nrow(gRNA.table1)==0){
            next
          }
        }
        else{
          next
        }
      }
      
      #0-0-0(优化)
      count_0 <- numeric()
      for (p in 1:nrow(gRNA.table1)) {
        if (gRNA.table1[p,]$V9 == "0-0-0") {
          count_0 <- append(count_0, p)
        }
      }
      if (length(count_0) != 0) {
        gRNA.table_min <- gRNA.table1[-count_0, ]
        gRNA.table1 <- rbind(gRNA.table1[count_0, ], gRNA.table_min)
      }
      if (nrow(gRNA.table1) > 10) {
        gRNA.table1 <- gRNA.table1[1:10,]
      }
      
      #外显子下游
      ko_start2 <- KO_region2[t,]$end
      ko_end2 <- KO_region2[t,]$end + 400
      # if (ko_end2 > nchar(Gene)) {
      #   KO_region2[t,]$end <- max(t_Exon_CDS$end)
      # }
      # ko_start2 <- KO_region2[t,]$end
      # ko_end2 <- KO_region2[t,]$end + 400
      ko_seq2 <- substring(Gene, ko_start2, ko_end2)
      if (Gene_rev) {
        seq2 <- DNAString(ko_seq2)
        seq2_rev <- reverse(seq2)
        ko_seq2 <- as.character(seq2_rev)
      }
      #读取ko_seq的gRNA表格
      source_python("crispor_table_download.py")
      py$run(ko_seq2, species,filepath,"3")
      gRNA.table2 <-
        read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
      #特异性得分60下，和Inefficient的gRNA排除掉
      gRNA.del <- numeric()
      gRNA.table2<-gRNA.table2[which(gRNA.table2$V3!="No matches"),]
      for (i in 1:length(gRNA.table2[, 1])) {
        if (as.numeric(gRNA.table2[i, 3]) < 60) {
          gRNA.del <- append(gRNA.del, i)
        }
        else if (grepl("Inefficient", gRNA.table2[i, 2])) {
          gRNA.del <- append(gRNA.del, i)
        }
      }
      if(length(gRNA.del)!=0){
        gRNA.table2 <- gRNA.table2[-gRNA.del, ]
      }
      
      if(nrow(gRNA.table2)==0){
        #往后找400bp
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
        if(judge){
          ko_start2 <- KO_region2[t,]$end +400
          ko_end2 <- KO_region2[t,]$end + 800
          ko_seq2 <- substring(Gene, ko_start2, ko_end2)
          if (Gene_rev) {
            seq2 <- DNAString(ko_seq2)
            seq2_rev <- reverse(seq2)
            ko_seq2 <- as.character(seq2_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq2, species,filepath,"3")
          gRNA.table2 <-
            read.csv(paste0(filepath,"//gRNA3.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table2<-gRNA.table2[which(gRNA.table2$V3!="No matches"),]
          for (i in 1:length(gRNA.table2[, 1])) {
            if (as.numeric(gRNA.table2[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table2[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table2 <- gRNA.table2[-gRNA.del, ]
          }
          if(nrow(gRNA.table2)==0){
            next
          }
        }
        else{
          next
        }
      }
      
      #0-0-0(优化)
      count_0 <- numeric()
      for (p in 1:nrow(gRNA.table2)) {
        if (gRNA.table2[p,]$V9 == "0-0-0") {
          count_0 <- append(count_0, p)
        }
      }
      if (length(count_0) != 0) {
        gRNA.table_min <- gRNA.table2[-count_0, ]
        gRNA.table2 <- rbind(gRNA.table2[count_0, ], gRNA.table_min)
      }
      if (nrow(gRNA.table2) > 10) {
        gRNA.table2 <- gRNA.table2[1:10,]
      }
      if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
        #筛选不到gRNA时要及时退出
        next
      }
      
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
      gene <- DNAString(Gene)
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
      GC_avoid_region <- GC_analysis3(KO_region2[t,])
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
      write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor3.csv"), row.names = FALSE)
      source_python("crispr_get_score.py")
      py$reader_writer(paste0(filepath,"//CCTOP-predictor3.csv"),species)
      gRNA.table <- read.csv(paste0(filepath,"//CCTOP-predictor3.csv"), header = TRUE)
      print(gRNA.table)
      #切割效率得分低于0.60的删除
      gRNA.table <-
        gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
      
      #上下游分开
      gRNA.table1 <-
        gRNA.table[which(gRNA.table$end <= KO_region2[t,]$start),]
      gRNA.table2 <-
        gRNA.table[which(gRNA.table$start >= KO_region2[t,]$end),]
      
      #切割得分大于0.65的优先
      gRNA.table1 <-
        rbind(gRNA.table1[which(gRNA.table1$crispr_score >= 0.65), ], gRNA.table1[which(gRNA.table1$crispr_score <
                                                                                          0.65), ])
      gRNA.table2 <-
        rbind(gRNA.table2[which(gRNA.table2$crispr_score >= 0.65), ], gRNA.table2[which(gRNA.table2$crispr_score <
                                                                                          0.65), ])
      
      #重新合并
      gRNA.table <- rbind(gRNA.table1, gRNA.table2)
      
      #table1是外显子上游的gRNA,table2是外显子下游的gRNA
      {
        if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
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
      else{
        print("change region")
      }
    }
    if (exists("gRNA_planC") == TRUE) {
      break
    }
  }
}

#进度条90%
write.table("8",paste0(filepath,"//","C90%.txt"),row.names = FALSE,col.names = FALSE)

if (exists("gRNA_planC") == TRUE) {
  # 敲了哪些外显子 -----------------------------------------------------------------
  if(exists("Gene3") == TRUE){
    t_Exon_region2 <-
      t_Exon_region[which(
        t_Exon_region$Exon_start +1200 >= min(gRNA_planC$start) &
          t_Exon_region$Exon_end +1200 <= max(gRNA_planC$end)
      ), ]
    t_Exon_CDS2 <-
      t_Exon_CDS[which(t_Exon_CDS$start +1200 >= min(gRNA_planC$start) &
                         t_Exon_CDS$end  +1200 <= max(gRNA_planC$end)), ]
    which_ko <- t_Exon_CDS2$Exon
  }
  else{
    t_Exon_region2 <-
      t_Exon_region[which(
        t_Exon_region$Exon_start >= min(gRNA_planC$start) &
          t_Exon_region$Exon_end <= max(gRNA_planC$end)
      ), ]
    t_Exon_CDS2 <-
      t_Exon_CDS[which(t_Exon_CDS$start >= min(gRNA_planC$start) &
                         t_Exon_CDS$end <= max(gRNA_planC$end)), ]
    which_ko <- t_Exon_CDS2$Exon
  }
}

# 画图 ----------------------------------------------------------------------
{
  if (exists("gRNA_planC") == TRUE) {
    {
      if (exists("Gene3") == TRUE) {
        y <- 0.5
        f <- data.frame(x = c(1:nchar(Gene3)), y = y)
        p1 <-
          ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#bdc4ca", size = 2) + theme_bw() +
          theme(panel.grid = element_blank(), panel.border = element_blank()) +
          scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
          xlab(NULL) + ylab(NULL)
        {
          if (Gene_rev) {
            label <- paste0(transcript.name, "<")
          }
          else{
            label <- paste0(transcript.name, ">")
          }
        }
        for (i in 1:nrow(t_Exon_region)) {
          start <- as.numeric(t_Exon_region[i, ]$Exon_start + 1200)
          end <- as.numeric(t_Exon_region[i, ]$Exon_end + 1200)
          p1 <- p1 + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            colour = "orange2",
            alpha = .0
          )
        }  #最长的转录本的展示
        for (i in 1:nrow(t_Exon_CDS)) {
          start <- as.numeric(t_Exon_CDS[i, ]$start + 1200)
          end <- as.numeric(t_Exon_CDS[i, ]$end + 1200)
          p1 <- p1 + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            fill = "orange2"
          )+annotate(
            "text",
            label = label,
            x = 1,
            y = y - 0.07,
            size = 4,
            hjust = 0,
            color = "orange2"
          )
        }
        
        p1 <- p1 + annotate(
          "rect",
          xmin = min(gRNA_planC$start) ,
          xmax = max(gRNA_planC$end) ,
          ymin = y - 0.07,
          ymax = y + 0.07,
          alpha = .0,
          color = "#D01027"
        )
        
        # 敲除区域放大图 -----------------------------------------------------------------
        ff <-data.frame(x = gRNA_planC[1,]$start,y = 0.59)
        fff <-data.frame(x = gRNA_planC[2,]$start, y = 0.59)
        large_start <- min(gRNA_planC$start) - 300
        large_end <- max(gRNA_planC$end) + 300
        y <- 0.5
        f <- data.frame(x = c(large_start:large_end), y = y)
        p1_large <-
          ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
          theme(panel.grid = element_blank(), panel.border = element_blank()) +
          scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
          xlab(NULL) + ylab(NULL)
        
        for (i in 1:nrow(t_Exon_region2)) {
          start <- as.numeric(t_Exon_region2[i, ]$Exon_start) +1200
          end <- as.numeric(t_Exon_region2[i, ]$Exon_end) +1200
          p1_large <- p1_large + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            colour = "#D01027",
            alpha = .0
          )
        }
        for (i in 1:nrow(t_Exon_CDS2)) {
          start <- as.numeric(t_Exon_CDS2[i, ]$start) +1200
          end <- as.numeric(t_Exon_CDS2[i, ]$end)+1200
          p1_large <- p1_large + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            fill = "#D01027"
          ) + annotate(
            "text",
            label = sub("xon ", "", t_Exon_CDS2[i, ]$Exon),
            x = (start+end)/2,
            y = y - 0.07,
            size = 4,
            hjust = 0,
            color = "black"
          )
        }
        
        img <- 'cut.png'
        p1_large <-
          p1_large + geom_image(data = ff, aes(x = x, y = y), image = img) +
          geom_image(data = fff, aes(x = x, y = y), image = img) +
          annotate(
            "text",
            label = "g1",
            x = gRNA_planC[1,]$start,
            y = 0.67,
            size = 6
          ) +
          annotate(
            "text",
            label = "g2",
            x = gRNA_planC[2,]$start,
            y = 0.67,
            size = 6
          )
      }

# 画图 ----------------------------------------------------------------------
      else{
        y <- 0.5
        f <- data.frame(x = c(1:nchar(Gene)), y = y)
        {
          if (Gene_rev) {
            label <- paste0(transcript.name, "<")
          }
          else{
            label <- paste0(transcript.name, ">")
          }
        }
        p1 <-
          ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#bdc4ca", size = 2) + theme_bw() +
          theme(panel.grid = element_blank(), panel.border = element_blank()) +
          scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) + xlab(NULL) +
          ylab(NULL) + annotate(
            "text",
            label = label,
            x = 1,
            y = y - 0.07,
            size = 4,
            hjust = 0,
            color = "orange2"
          )
        
        for (i in 1:nrow(t_Exon_region)) {
          start <- as.numeric(t_Exon_region[i, ]$Exon_start)
          end <- as.numeric(t_Exon_region[i, ]$Exon_end)
          p1 <- p1 + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            colour = "orange2",
            alpha = .0
          )
        }  #最长的转录本的展示
        for (i in 1:nrow(t_Exon_CDS)) {
          start <- as.numeric(t_Exon_CDS[i, ]$start)
          end <- as.numeric(t_Exon_CDS[i, ]$end)
          p1 <- p1 + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            fill = "orange2",
          )
        }
        
        p1 <- p1 + annotate(
          "rect",
          xmin = min(gRNA_planC$start) ,
          xmax = max(gRNA_planC$end) ,
          ymin = y - 0.07,
          ymax = y + 0.07,
          alpha = .0,
          color = "#D01027"
        )
        
        # 敲除区域放大图 -----------------------------------------------------------------
        ff <-data.frame(x = gRNA_planC[1,]$start,y= 0.59)
        fff <-data.frame(x = gRNA_planC[2,]$start,y = 0.59)
        large_start <- min(gRNA_planC$start) - 300
        large_end <- max(gRNA_planC$end) + 300
        y <- 0.5
        f <- data.frame(x = c(large_start:large_end), y = y)
        p1_large <-
          ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
          theme(panel.grid = element_blank(), panel.border = element_blank()) +
          scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) +
          xlab(NULL) + ylab(NULL)
        for (i in 1:nrow(t_Exon_region2)) {
          start <- as.numeric(t_Exon_region2[i, ]$Exon_start)
          end <- as.numeric(t_Exon_region2[i, ]$Exon_end)
          p1_large <- p1_large + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            colour = "#D01027",
            alpha = .0
          )
        }
        for (i in 1:nrow(t_Exon_CDS2)) {
          start <- as.numeric(t_Exon_CDS2[i, ]$start)
          end <- as.numeric(t_Exon_CDS2[i, ]$end)
          p1_large <- p1_large + annotate(
            "rect",
            xmin = start,
            xmax = end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            fill = "#D01027"
          ) + annotate(
            "text",
            label = sub("xon ", "", t_Exon_CDS2[i, ]$Exon),
            x = (start+end)/2,
            y = y - 0.07,
            size = 4,
            hjust = 0,
            color = "black"
          )
        }
        
        img <- 'cut.png'
        p1_large <-
          p1_large + geom_image(data = ff, aes(x = x, y = y), image = img) +
          geom_image(data = fff, aes(x = x, y = y), image = img) +
          annotate(
            "text",
            label = "g1",
            x = gRNA_planC[1,]$start,
            y = 0.67,
            size = 6
          ) +
          annotate(
            "text",
            label = "g2",
            x = gRNA_planC[2,]$start,
            y = 0.67,
            size = 6
          )
      }
    }
        
    png(file = paste0(filepath, "//", "gRNA_position3.png"),width = 480 * 3,height = 480 * 2,res = 72 * 2)
    print(p1)
    dev.off()
    png(file = paste0(filepath, "//", "gRNA_position_large3.png"),width = 480 * 3,height = 480 * 2,res = 72 * 2)
    print(p1_large)
    dev.off()
  }
  else{
    print("fail")
    if(nrow(KO_region2)!=0){
      write.table("3",paste0(filepath,"//","result3.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}

# 输出 ----------------------------------------------------------------------
#敲除大小
if (exists("gRNA_planC") == TRUE) {
  if (Gene_rev) {
    if (all(gRNA_planC$strand == "rev") |
        all(gRNA_planC$strand == "fw")) {
      KO_length <- abs(gRNA_planC[1,]$end - gRNA_planC[2,]$end) 
    }
    else{
      if (gRNA_planC[1, ]$strand == "rev" & gRNA_planC[2, ]$strand == "fw") {
        if(gRNA_planC[1,]$start>gRNA_planC[2,]$start){
          pos1 <- gRNA_planC[1, ]$end-3
          pos2 <- gRNA_planC[2, ]$start+3
          KO_length <- abs(pos1 - pos2) +1
        }
        else{
          pos1 <- gRNA_planC[1, ]$end-3
          pos2 <- gRNA_planC[2, ]$start+3
          KO_length <- abs(pos1 - pos2) -1
        }
      }
      else{
        if(gRNA_planC[1,]$start>gRNA_planC[2,]$start){
          pos1 <- gRNA_planC[1, ]$start+3
          pos2 <- gRNA_planC[2, ]$end-3
          KO_length <- abs(pos1 - pos2) - 1
        }
        else{
          pos1 <- gRNA_planC[1, ]$start+3
          pos2 <- gRNA_planC[2, ]$end-3
          KO_length <- abs(pos1 - pos2) + 1
        }
      }
    }
  }
  #正向
  else{
    if (all(gRNA_planC$strand == "rev") |
        all(gRNA_planC$strand == "fw")) {
      KO_length <- abs(gRNA_planC[1,]$end - gRNA_planC[2,]$end) 
    }
    else{
      if (gRNA_planC[1, ]$strand == "rev" & gRNA_planC[2, ]$strand == "fw") {
        pos1 <- gRNA_planC[1, ]$start+3
        pos2 <- gRNA_planC[2, ]$end-3
        KO_length <- abs(pos1 - pos2-1)
      }
      else{
        pos1 <- gRNA_planC[1, ]$end-3
        pos2 <- gRNA_planC[2, ]$start+3
        KO_length <- abs(pos1 - pos2+1)
      }
    }
  }
  #敲除的CDS
  KO_length_CDS <- KO_region2[t,]$Exon_length
  
  {
    if (exists("gRNA2_planC") == TRUE) {
      if (Gene_rev) {
        if (all(gRNA2_planC$strand == "rev") |
            all(gRNA2_planC$strand == "fw")) {
          KO_length2 <- abs(gRNA2_planC[1, ]$end - gRNA2_planC[2, ]$end)
        }
        else{
          if (gRNA2_planC[1,]$strand == "rev" & gRNA2_planC[2,]$strand == "fw") {
            if(gRNA2_planC[1,]$start>gRNA2_planC[1,]$start){
              pos1 <- gRNA2_planC[1,]$end-3
              pos2 <- gRNA2_planC[2,]$start+3
              KO_length2 <- abs(pos1 - pos2) +1
            }
            else{
              pos1 <- gRNA2_planC[1,]$end-3
              pos2 <- gRNA2_planC[2,]$start+3
              KO_length2 <- abs(pos1 - pos2) -1
            }
          }
          else{
            if(gRNA2_planC[1,]$start>gRNA2_planC[2,]$start){
              pos1 <- gRNA2_planC[1, ]$start+3
              pos2 <- gRNA2_planC[2, ]$end-3
              KO_length2 <- abs(pos1 - pos2) - 1
            }
            else{
              pos1 <- gRNA2_planC[1, ]$start+3
              pos2 <- gRNA2_planC[2, ]$end-3
              KO_length2 <- abs(pos1 - pos2) + 1
            }
          }
        }
      }
      #正向
      else{
        if (all(gRNA2_planC$strand == "rev") |
            all(gRNA2_planC$strand == "fw")) {
          KO_length2 <- abs(gRNA2_planC[1, ]$end - gRNA2_planC[2, ]$end)
        }
        else{
          if (gRNA2_planC[1,]$strand == "rev" & gRNA2_planC[2,]$strand == "fw") {
            pos1 <- gRNA2_planC[1,]$start+3
            pos2 <- gRNA2_planC[2,]$end-3
            KO_length2 <- abs(pos1 - pos2-1)
          }
          else{
            pos1 <- gRNA2_planC[1,]$end-3
            pos2 <- gRNA2_planC[2,]$start+3
            KO_length2 <- abs(pos1 - pos2+1)
          }
        }
      }
      #敲除的CDS
      KO_length_CDS2 <- KO_region2[t,]$Exon_length
      mark<-"FALSE"
    }
    
    else{
      gRNA2_planC<-Get_gRNA2_planC(KO_region2)
      mark<-"TRUE"
    }
  }
}
if(exists("gRNA2_planC")==TRUE){
  if(class(gRNA2_planC)=="NULL"){
    rm(mark)
    rm(gRNA2_planC)
  }
}


if (exists("gRNA_planC") == TRUE) {
  output <- data.frame(Gene=character(),gRNA1=character(),strand1=character(),score1=numeric(),score2=numeric(),
                       gRNA2=character(),strand2=character(),score2_1=numeric(),score2_2=numeric(),region=numeric(),cds=numeric(),
                       gRNA3=character(),strand3=character(),score3_1=numeric(),score3_2=numeric(),
                       gRNA4=character(),strand4=character(),score4_1=numeric(),score4_2=numeric(),region2=numeric(),cds2=numeric(),
                       incomplete_transcript=character(),transcript=character(),Exon_count=numeric(),
                       start_condon=character(),stop_condon=character(),ko_condon=character(),
                       overlap1=character(),overlap2=character(),tip1=character(),
                       tip2=character(),tip3=character(),tip4=character(),GC1=numeric(),GC2=numeric(),
                       mark=character(),ko_condon2=character(),species=character())
  output[1, 1] <- term
  output[1, 2] <- gRNA_planC[1,]$analysis_seq
  output[1, 3] <- gRNA_planC[1,]$strand
  output[1, 4] <- gRNA_planC[1,]$Score1
  output[1, 5] <- gRNA_planC[1,]$crispr_score
  output[1, 6] <- gRNA_planC[2,]$analysis_seq
  output[1, 7] <- gRNA_planC[2,]$strand
  output[1, 8] <- gRNA_planC[2,]$Score1
  output[1, 9] <- gRNA_planC[2,]$crispr_score
  output[1, 10] <- KO_length
  output[1, 11] <- KO_length_CDS
  output[1, 38] <- species
  #有没有重叠的lncRNA,microRNA...
  if(nrow(mark_region)!=0){
    for(j in 1:nrow(mark_region)){
      if(min(gRNA_planC$start)>mark_region[j,]$end | max(gRNA_planC$end)<mark_region[j,]$start){
        overlap1<-"FALSE"
      }
      else{
        overlap1<-"TRUE"
      }
    }
    output[1, 28] <- overlap1
  }
  if(exists("gRNA2_planC")==TRUE){
    output[1, 12] <- gRNA2_planC[1,]$analysis_seq
    output[1, 13] <- gRNA2_planC[1,]$strand
    output[1, 14] <- gRNA2_planC[1,]$Score1
    output[1, 15] <- gRNA2_planC[1,]$crispr_score
    output[1, 16] <- gRNA2_planC[2,]$analysis_seq
    output[1, 17] <- gRNA2_planC[2,]$strand
    output[1, 18] <- gRNA2_planC[2,]$Score1
    output[1, 19] <- gRNA2_planC[2,]$crispr_score
    if(mark=="TRUE"){
      output[1, 20] <- gRNA2_planC[1,8]
      output[1, 21] <- gRNA2_planC[1,9]
      output[1, 37] <- gRNA2_planC[2,8]
    }
    else if(mark=="FALSE"){
      output[1, 20] <- KO_length2
      output[1, 21] <- KO_length_CDS2
    }
    output[1, 36]<-mark
    #有没有重叠的lncRNA,microRNA...
    if(nrow(mark_region)!=0){
      for(j in 1:nrow(mark_region)){
        if(min(gRNA2_planC$start)>mark_region[j,]$end | max(gRNA2_planC$end)<mark_region[j,]$start){
          overlap2<-"FALSE"
        }
        else{
          overlap2<-"TRUE"
        }
      }
      output[1, 29] <- overlap2
    }
  }
  
  output[1, 22] <- paste(incomplete.transcript,collapse = ",")
  output[1, 23] <- transcript.name
  output[1, 24] <- nrow(t_Exon_region)
  output[1, 25] <- ATG_Exon
  output[1, 26] <- stop_Exon
  output[1, 27] <- paste(which_ko,collapse = ",")
  
  print(gRNA_planC)

  
  #点阵分析和GC含量分析
  {
    if (exists("Gene3")) {
      p3 <- Get_dot_plot4(KO_region2[t, ])
      p4 <- Get_dot_plot5(KO_region2[t, ])
      p5 <- Get_GC_image4(KO_region2[t, ])
      p6 <- Get_GC_image5(KO_region2[t, ])
      #上游发卡结构
      tip1 <- dot_analysis5(KO_region2[t, ])
      #下游发卡结构
      tip2 <- dot_analysis6(KO_region2[t, ])
      #上游片段重复
      tip3 <- dot_analysis7(KO_region2[t, ])
      #下游片段重复
      tip4 <- dot_analysis8(KO_region2[t, ])
      #上游GC含量
      GC1 <- GC_analysis7(KO_region2[t, ])
      #下游GC含量
      GC2 <- GC_analysis8(KO_region2[t, ])
      
    }
    else{
      p3 <- Get_dot_plot2(KO_region2[t, ])
      p4 <- Get_dot_plot3(KO_region2[t, ])
      p5 <- Get_GC_image2(KO_region2[t, ])
      p6 <- Get_GC_image3(KO_region2[t, ])
      #上游发卡结构
      tip1 <- dot_analysis1(KO_region2[t, ])
      #下游发卡结构
      tip2 <- dot_analysis2(KO_region2[t, ])
      #上游片段重复
      tip3 <- dot_analysis3(KO_region2[t, ])
      #下游片段重复
      tip4 <- dot_analysis4(KO_region2[t, ])
      #上游GC含量
      GC1 <- GC_analysis5(KO_region2[t, ])
      #下游GC含量
      GC2 <- GC_analysis6(KO_region2[t, ])
    }
  }
  output[1, 30] <- tip1
  output[1, 31] <- tip2
  output[1, 32] <- tip3
  output[1, 33] <- tip4
  output[1, 34] <- round(GC1*100,2)
  output[1, 35] <- round(GC2*100,2)
  
  
  png(file = paste0(filepath,"//","Lattice and diagram3.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p3)
  dev.off()
  png(file = paste0(filepath,"//","Lattice and diagram4.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p4)
  dev.off()
  
  png(file = paste0(filepath,"//","GC content3.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p5)
  dev.off()
  png(file = paste0(filepath,"//","GC content4.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p6)
  dev.off()
  # 所有转录本的展示图 ---------------------------------------------------------------
  p2<-KO_region_image3(Gene,allinfo)
  image_write(p2,paste0(filepath,"//","transcript_region3.png"))
  
  write.csv(output, file = paste0(filepath,"//","planC.csv"), row.names = FALSE)
  print("success")
}



