Get_gRNA<-function(ID,species){
  setwd("C://Users//41518//Desktop//work/ubigene")
  library(ggplot2)
  library(ggpubr)
  library(httr)
  library(jsonlite)
  library(xml2)
  library(rvest)
  library(Biostrings)
  library(dplyr)
  library(reticulate)
  source("Get_ID.R")
  source("Get_allinfo.R")
  source("Get_seq.R")
  source("Get_transcript_table.R")
  source("delete_noprotein.R")
  source("Get_max_transcript.R")
  source("Get_dot_region.R")
  source("Get_dot_region2.R")
  source("Get_result1.R")
  source("Get_result2.R")
  source("Get_result3.R")
  source("GC_analysis.R")
  source("KO_longregion.R")
  source("KO_region_300.R")
  source_python("C://Users//41518//Desktop//work//Ubigene//crispr_get_score(2)(1).py")
  source_python("C://Users//41518//Desktop//work//Ubigene//crispor_table_download.py")
  
  Gene <- Get_seq(ID)                    #基因序列
  Gene2 <- Get_seq2(ID)                 #5'和3'端各增加500bp
  allinfo <- Get_allinfo(ID)
  start <- allinfo$start
  transcript.table <- Get_transcript_table(ID,species)
  transcript.table <-
    delete_noprotein(transcript.table)           #删除非编码蛋白和不完整的转录本
  transcript.count <- nrow(transcript.table)     #转录本数量
  transcript.name <-
    transcript.table[which(transcript.table$bp == max(transcript.table$bp)), ]$Name     #最长的转录本
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
  t_Exon_CDS <-
    ko.data[which(ko.data$transcript == transcript.name),]
  #计算覆盖情况
  ko.data <-
    transform(ko.data, times = numeric(length(ko.data$transcript)))
  for (i in 1:length(ko.data$transcript)) {
    times <- 0
    for (j in 1:length(ko.data$transcript)) {
      if (max(ko.data[i, ]$start, ko.data[j, ]$start) <= min(ko.data[i, ]$end, ko.data[j, ]$end)) {
        times <- times + 1
      }
    }
    ko.data$times[i] <- times
  }
  
  #起始密码子
  Exon_CDS <- data.frame()
  for (i in 1:length(transcript.table$Name)) {
    Exon_CDS <-
      rbind(Exon_CDS, ko.data[which(ko.data$transcript == transcript.table$Name[i]),][1,])
  }
  
  #参考转录本的KO区域
  KO_region <-
    ko.data[which(ko.data$transcript == transcript.name), ]
  
  #删除不在前30%CDS的
  l <- 0
  CDS_length <- KO_region$end - KO_region$start + 1
  CDS_length_30 <- floor(sum(CDS_length) * 0.3)
  for (q in 1:length(CDS_length)) {
    l <- l + CDS_length[q]
    if (l > CDS_length_30) {
      break
    }
  }
  t_CDS_30 <- KO_region[q,]$start + l - CDS_length_30
  KO_region <- KO_region[1:q,]
  
  #在前30%的起始密码子
  n <- which(Exon_CDS$start > t_CDS_30)
  if (length(n) != 0) {
    if (Gene_rev) {
      Exon_CDS <- Exon_CDS[n,]
    }
    else{
      Exon_CDS <- Exon_CDS[-n,]
    }
  }
  
  Exon_CDS <- Exon_CDS[!duplicated(Exon_CDS$start), ]
  Exon_length <- Exon_CDS$end - Exon_CDS$start + 1
  Exon_CDS <- cbind(Exon_CDS, Exon_length)
  
  #判断该外显子是否为起始密码子所在的外显子
  region_del <- numeric()
  for (i in 1:nrow(KO_region)) {
    if (Gene_rev) {
      if (KO_region[i,]$end == t_CDS_start) {
        region_del <- append(region_del, i)
      }
    }
    else{
      if (KO_region[i,]$start == t_CDS_start) {
        region_del <- append(region_del, i)
      }
    }
  }
  if (length(region_del) != 0) {
    KO_region <- KO_region[-region_del,]
  }
  
  KO_region3 <- KO_region
  #KO_region
  Exon_300 <- data.frame()
  Exon_length_500 <- numeric()
  k <- 1
  Exon_length <- numeric()
  for (i in 1:nrow(KO_region)) {
    #判断该区域前后500bp是否有其他外显子，有则合并
    Exon.name <- KO_region[i, ]$Exon
    j <- which(t_Exon_region_sort$Exon_name == Exon.name)
    Exon_length[i] <- KO_region[i, ]$end - KO_region[i, ]$start + 1
    ko_start <- KO_region[i, ]$start - 500
    #Exon_300是长度大于300bp的外显子
    if (Exon_length[i] >= 300) {
      Exon_length_500[k] <- KO_region[i, ]$end - KO_region[i, ]$start + 1
      Exon_300 <- rbind(Exon_300, KO_region[i,])
      k <- k + 1
    }
    repeat {
      if (ko_start <= t_Exon_region_sort[j - 1, ]$Exon_end) {
        ko_start <- t_Exon_region_sort[j - 1, ]$Exon_start - 500
        Exon_length[i] <-
          Exon_length[i] + t_Exon_region_sort[j - 1, ]$Exon_end - t_Exon_region_sort[j - 1, ]$Exon_start +
          1
        j <- j - 1
        if (j == 1) {
          break
        }
      }
      else{
        break
      }
    }
    j <- which(t_Exon_region_sort$Exon_name == Exon.name)
    ko_end <- KO_region[i, ]$end + 500
    repeat {
      if (ko_end >= t_Exon_region_sort[j + 1, ]$Exon_start) {
        ko_end <- t_Exon_region_sort[j + 1, ]$Exon_end + 500
        Exon_length[i] <-
          Exon_length[i] + t_Exon_region_sort[j + 1, ]$Exon_end - t_Exon_region_sort[j + 1, ]$Exon_start +
          1
        j <- j + 1
        if (j == nrow(t_Exon_region_sort)) {
          break
        }
      }
      else{
        break
      }
    }
    KO_region[i, ]$start <- ko_start + 500
    KO_region[i, ]$end <- ko_end - 500
  }
  #KO区域排序
  KO_region <- cbind(KO_region, Exon_length)
  {
    if (Gene_rev) {
      KO_region <- KO_region[order(-KO_region$times,-KO_region$end),]
    }
    else{
      KO_region <- KO_region[order(-KO_region$times, KO_region$end),]
    }
  }
  
  KO_region <- KO_region[!duplicated(KO_region$Exon_length), ]
  
  KO_region <-
    KO_region[which(KO_region$Exon_length %% 3 != 0), ]              #外显子非3的倍数
  KO_region_500 <-
    KO_region[which(KO_region$Exon_length >= 100),]         #外显子不小于100bp
  
  Exon_300 <- cbind(Exon_300, Exon_length_500)          #单个大于300bp的外显子
  names(Exon_300)[6] <- "Exon_length"
  
  #times-1
  for (i in 1:nrow(KO_region_500)) {
    for (j in 1:nrow(Exon_CDS)) {
      if (Exon_CDS[j,]$start >= KO_region_500[i,]$start &
          Exon_CDS[j,]$end <= KO_region_500[i,]$end) {
        KO_region_500[i,]$times <- KO_region_500[i,]$times - 1
      }
    }
  }
  
  for (i in 1:nrow(Exon_300)) {
    for (j in 1:nrow(Exon_CDS)) {
      if (Exon_CDS[j,]$start >= Exon_300[i,]$start &
          Exon_CDS[j,]$end <= Exon_300[i,]$end) {
        Exon_300[i,]$times <- Exon_300[i,]$times - 1
      }
    }
  }
  
  #KO区域中包含了最后一个外显子时
  {
    if (Gene_rev) {
      KO_region_500 <- KO_region_500[order(KO_region_500$start),]
      if(KO_region_500[1,]$start==t_Exon_region_sort[1,]$Exon_start){
        KO_region_500[1,]$start<-t_Exon_CDS[nrow(t_Exon_CDS),]$start
      }
    }
    else{
      KO_region_500<- KO_region_500[order(-KO_region_500$end),]
      if(KO_region_500[1,]$end==t_Exon_region[nrow(t_Exon_region),]$Exon_end){
        KO_region_500[1,]$end<-t_Exon_CDS[nrow(t_Exon_CDS),]$end
        KO_region_500[1,]$Exon_length<-KO_region_500[1,]$end-KO_region_500[1,]$start
      }
    }
  }
  
  KO_region <- rbind(KO_region_500, Exon_300)
  Exon_CDS <- Exon_CDS[which(Exon_CDS$Exon_length >= 300),]
  KO_region <- rbind(KO_region, Exon_CDS)
  
  #排序
  {
    if (Gene_rev) {
      KO_region <- KO_region[order(-KO_region$times,-KO_region$end),]
    }
    else{
      KO_region <- KO_region[order(-KO_region$times, KO_region$end),]
    }
  }
  print(KO_region)
  print(Exon_300)
  print(KO_region_500)
  print(Exon_CDS)    #含有起始密码子的外显子
  
  #GC含量分析:平均GC含量大于70%或小于30%，则删除该区域
  GC_del <- numeric()
  for (i in 1:nrow(KO_region)) {
    analysis_GC <- GC_analysis1(KO_region[i, ])
    if (analysis_GC == TRUE) {
      GC_del <- append(GC_del, i)
    }
  }
  if (length(GC_del) != 0) {
    KO_region <- KO_region[-GC_del,]
  }
  
  #点阵图分析
  dot_del <- numeric()
  for (i in 1:nrow(KO_region)) {
    analysis_dot1 <- Get_dot_region1(KO_region[i,])
    analysis_dot2 <- Get_dot_region2(KO_region[i,])
    analysis_dot3 <- Get_dot_region3(KO_region[i,])
    #analysis_dot4 <- Get_dot_region4(KO_region[i,])
    if (analysis_dot1 == TRUE | analysis_dot2 == TRUE |
        analysis_dot3 == TRUE | analysis_dot4 == TRUE ) {
      dot_del <- append(dot_del, i)
    }
  }
  if (length(dot_del) != 0) {
    KO_region <- KO_region[-dot_del,]
  }
  
  print(KO_region)
  
  # gRNA设计方案 ----------------------------------------------------------------
  if (nrow(KO_region) != 0) {
    for (t in 1:nrow(KO_region)) {
      # 设计在起始密码子后面 -----------------------------------------------------------------
      for (j in 1:nrow(Exon_CDS)) {
        if (all(KO_region[t, ] %in% Exon_CDS[j, ])) {
          ko_start <- KO_region[t, ]$start + 3
          ko_end <- KO_region[t, ]$end
          ko_seq <- substring(Gene, ko_start, ko_end)
          if (Gene_rev) {
            seq <- DNAString(ko_seq)
            seq_rev <- reverse(seq)
            ko_seq <- as.character(seq_rev)
          }
          #读取ko_seq的gRNA表格
          py$run(ko_seq)
          gRNA.table <- read.csv("gRNA.csv", header = FALSE)
          
          #特异性得分70以下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          for (i in 1:length(gRNA.table[, 1])) {
            if (gRNA.table[i, 3] <= 70) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          gRNA.table <- gRNA.table[-gRNA.del,]
          
          if (nrow(gRNA.table) < 2) {
            #筛选不到gRNA时要及时退出
            next
          }
          
          #获取每个gRNA在基因上的位置
          strand <- sub("[^a-zA-Z]+", "", gRNA.table[, 1])
          gRNA_seq <- substring(gRNA.table[, 2], 1, 20)
          Score1 <- gRNA.table[, 3]
          analysis_seq <-
            gsub(" ", "", substring(gRNA.table[, 2], 1, 24))
          gRNA.table <-
            cbind(strand, gRNA_seq, analysis_seq, Score1)
          gRNA.table <- as.data.frame(gRNA.table)
          gene <- DNAString(Gene)
          {
            if (Gene_rev) {
              for (i in 1:length(gRNA.table[, 1])) {
                if (gRNA.table[i, ]$strand == "fw") {
                  gRNA_rev <- DNAString(gRNA.table[i, ]$gRNA_seq)
                  gRNA_rev <- reverse(gRNA_rev)
                  rev <-
                    matchPattern(pattern = gRNA_rev, subject = gene)
                  rev_start <- start(rev)
                  rev_end <- end(rev)
                  gRNA.table[i, 5] <- rev_start
                  gRNA.table[i, 6] <- rev_end
                }
                else if (gRNA.table[i, ]$strand == "rev") {
                  gRNA_fw <- DNAString(gRNA.table[i, ]$gRNA_seq)
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
                if (gRNA.table[i,]$strand == "rev") {
                  gRNA_rev <- DNAString(gRNA.table[i,]$gRNA_seq)
                  gRNA_rev <- reverseComplement(gRNA_rev)
                  rev <-
                    matchPattern(pattern = gRNA_rev, subject = gene)
                  rev_start <- start(rev)
                  rev_end <- end(rev)
                  gRNA.table[i, 5] <- rev_start
                  gRNA.table[i, 6] <- rev_end
                }
                else if (gRNA.table[i,]$strand == "fw") {
                  fw <-
                    matchPattern(pattern = gRNA.table[i,]$gRNA_seq,
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
          GC_avoid_region <- GC_analysis2(KO_region[t, ])
          if (GC_avoid_region != FALSE) {
            GC_del <- numeric()
            for (i in 1:nrow(gRNA.table)) {
              for (j in 1:nrow(GC_avoid_region)) {
                if (gRNA.table[i, ]$start %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2]) |
                    gRNA.table[i, ]$end %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2])) {
                  GC_del <- append(GC_del, i)
                }
              }
            }
            if (length(GC_del) != 0) {
              gRNA.table <- gRNA.table[-GC_del,]
            }
          }
          
          if (nrow(gRNA.table) > 50) {
            gRNA.table <- gRNA.table[1:50,]
          }
          
          #符合条件的gRNA进行切割效率预测
          write.csv(gRNA.table, file = "CCTOP-predictor.csv", row.names = FALSE)
          source_python("C://Users//41518//Desktop//work//Ubigene//crispr_get_score(2)(1).py")
          py$reader_writer("C:\\Users\\41518\\Desktop\\work\\Ubigene\\CCTOP-predictor.csv")
          gRNA.table <-
            read.csv("CCTOP-predictor.csv", header = TRUE)
          
          #切割效率得分低于0.60的删除
          gRNA.table <-
            gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]
          
          #相差0.05分以内优选
          gRNA <- Get_result2(gRNA.table)
          
          
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA) == "NULL") {
            gRNA <- Get_result3(gRNA.table)
          }
          
          #第二对gRNA
          if (class(gRNA) != "NULL") {
            #在gRNA表格中删除第一对gRNA
            n=which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
            gRNA.table_2<-gRNA.table[-n,]
            #相差0.05分以内优选
            gRNA2<-Get_result2(gRNA.table_2)
            #如果没有相差0.05分以内的gRNA
            if (class(gRNA2) == "NULL") {
              gRNA2 <- Get_result3(gRNA.table_2)
            }
            
            if (class(gRNA2) == "NULL") {
              rm(gRNA2)
              #######换区域
            }
          }
          
          if (class(gRNA) == "NULL") {
            rm(gRNA)
          }
          
          if (exists("gRNA") == TRUE) {
            if (KO_region[t, ]$times != transcript.count) {
              No.transcript <-
                ko.data[which(ko.data$end == KO_region[t,]$end),]$transcript
              which.transcript <-
                transcript.table$Name %in% No.transcript
              Not_KO_transcript <-
                transcript.table$Name[which(which.transcript == "FALSE")]
              sprintf("%s转录本可能无法被影响到", Not_KO_transcript)
            }
          }
        }
      }
      if (exists("gRNA") == TRUE) {
        break
      }
      
      # 设计在内含子上 -----------------------------------------------------------------
      for (j in 1:nrow(KO_region_500)) {
        if (all(KO_region[t, ] %in% KO_region_500[j, ])) {
          #外显子上游
          ko_start1 <- KO_region[t, ]$start - 400
          ko_end1 <- KO_region[t, ]$start
          ko_seq1 <- substring(Gene, ko_start1, ko_end1)
          if (Gene_rev) {
            seq1 <- DNAString(ko_seq1)
            seq1_rev <- reverse(seq1)
            ko_seq1 <- as.character(seq1_rev)
          }
          #读取ko_seq的gRNA表格
          py$run(ko_seq1)
          gRNA.table1 <- read.csv("gRNA.csv", header = FALSE)
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
          gRNA.table1 <- gRNA.table1[-gRNA.del,]
          if (nrow(gRNA.table1) > 7) {
            gRNA.table1 <- gRNA.table1[1:7, ]
          }
          
          #外显子下游
          ko_start2 <- KO_region[t, ]$end
          ko_end2 <- KO_region[t, ]$end + 400
          ko_seq2 <- substring(Gene, ko_start2, ko_end2)
          if (Gene_rev) {
            seq2 <- DNAString(ko_seq2)
            seq2_rev <- reverse(seq2)
            ko_seq2 <- as.character(seq2_rev)
          }
          #读取ko_seq的gRNA表格
          py$run(ko_seq2)
          gRNA.table2 <- read.csv("gRNA.csv", header = FALSE)
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
          gRNA.table2 <- gRNA.table2[-gRNA.del,]
          if (nrow(gRNA.table2) > 7) {
            gRNA.table2 <- gRNA.table2[1:7, ]
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
          gRNA.table <-
            cbind(strand, gRNA_seq, analysis_seq, Score1)
          gRNA.table <- as.data.frame(gRNA.table)
          gene <- DNAString(Gene)
          {
            if (Gene_rev) {
              for (i in 1:length(gRNA.table[, 1])) {
                if (gRNA.table[i, ]$strand == "fw") {
                  gRNA_rev <- DNAString(gRNA.table[i, ]$gRNA_seq)
                  gRNA_rev <- reverse(gRNA_rev)
                  rev <-
                    matchPattern(pattern = gRNA_rev, subject = gene)
                  rev_start <- start(rev)
                  rev_end <- end(rev)
                  gRNA.table[i, 5] <- rev_start
                  gRNA.table[i, 6] <- rev_end
                }
                else if (gRNA.table[i, ]$strand == "rev") {
                  gRNA_fw <- DNAString(gRNA.table[i, ]$gRNA_seq)
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
                if (gRNA.table[i,]$strand == "rev") {
                  gRNA_rev <- DNAString(gRNA.table[i,]$gRNA_seq)
                  gRNA_rev <- reverseComplement(gRNA_rev)
                  rev <-
                    matchPattern(pattern = gRNA_rev, subject = gene)
                  rev_start <- start(rev)
                  rev_end <- end(rev)
                  gRNA.table[i, 5] <- rev_start
                  gRNA.table[i, 6] <- rev_end
                }
                else if (gRNA.table[i,]$strand == "fw") {
                  fw <-
                    matchPattern(pattern = gRNA.table[i,]$gRNA_seq,
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
          GC_avoid_region <- GC_analysis2(KO_region[t, ])
          if (GC_avoid_region != FALSE) {
            GC_del <- numeric()
            for (i in 1:nrow(gRNA.table)) {
              for (j in 1:nrow(GC_avoid_region)) {
                if (gRNA.table[i, ]$start %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2]) |
                    gRNA.table[i, ]$end %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2])) {
                  GC_del <- append(GC_del, i)
                }
              }
            }
            if (length(GC_del) != 0) {
              gRNA.table <- gRNA.table[-GC_del,]
            }
          }
          
          #符合条件的gRNA进行切割效率预测
          write.csv(gRNA.table, file = "CCTOP-predictor.csv", row.names = FALSE)
          source_python("C://Users//41518//Desktop//work//Ubigene//crispr_get_score(2)(1).py")
          py$reader_writer("C:\\Users\\41518\\Desktop\\work\\Ubigene\\CCTOP-predictor.csv")
          gRNA.table <-
            read.csv("CCTOP-predictor.csv", header = TRUE)
          
          #切割效率得分低于0.60的删除
          gRNA.table <-
            gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]
          
          #上下游分开
          gRNA.table1 <-
            gRNA.table[which(gRNA.table$end <= KO_region[t, ]$start), ]
          gRNA.table2 <-
            gRNA.table[which(gRNA.table$start >= KO_region[t, ]$end), ]
          #table1是外显子上游的gRNA,table2是外显子下游的gRNA
          {
            if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
              next
            }
            else{
              gRNA <- Get_result1(gRNA.table, KO_region[t,])        #相差0.05分以内优选
            }
          }
          
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA) == "NULL") {
            gRNA.table1 <-
              gRNA.table1[order(-gRNA.table1$crispr_score, -gRNA.table1$Score1), ]
            gRNA.table2 <-
              gRNA.table2[order(-gRNA.table2$crispr_score, -gRNA.table2$Score1), ]
            gRNA <- rbind(gRNA.table1[1,], gRNA.table2[1,])
          }
          
          #第二对gRNA
          if (class(gRNA) != "NULL" &
              nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
            #在gRNA表格中删除第一对gRNA
            n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
            gRNA.table_2 <- gRNA.table[-n, ]
            #相差0.05分以内优选
            gRNA2 <- Get_result1(gRNA.table_2, KO_region[t, ])
            #如果没有相差0.05分以内的gRNA
            if (class(gRNA2) == "NULL") {
              if (gRNA.table1[1, ]$analysis_seq %in% gRNA$analysis_seq) {
                if (gRNA.table2[1, ]$analysis_seq %in% gRNA$analysis_seq) {
                  gRNA2 <- rbind(gRNA.table1[2, ], gRNA.table2[2, ])
                }
                else{
                  gRNA2 <- rbind(gRNA.table1[2, ], gRNA.table2[1, ])
                }
              }
              else{
                if (gRNA.table2[1, ]$analysis_seq %in% gRNA$analysis_seq) {
                  gRNA2 <- rbind(gRNA.table1[1, ], gRNA.table2[2, ])
                }
                else{
                  gRNA2 <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
                }
              }
            }
            
            if (class(gRNA2) == "NULL") {
              rm(gRNA2)
              #######换区域
            }
          }
          
          if (exists("gRNA") == TRUE) {
            if (KO_region[t, ]$times != transcript.count) {
              No.transcript <-
                ko.data[which(ko.data$start == KO_region[t,]$start),]$transcript
              which.transcript <-
                transcript.table$Name %in% No.transcript
              Not_KO_transcript <-
                transcript.table$Name[which(which.transcript == "FALSE")]
            }
          }
        }
      }
      if (exists("gRNA") == TRUE) {
        judge <- "TRUE"
        break
      }
      
      # 设计在外显子上 -----------------------------------------------------------------
      for (j in 1:nrow(Exon_300)) {
        if (all(KO_region[t, ] %in% Exon_300[j, ])) {
          ko_start <- KO_region[t, ]$start
          ko_end <- KO_region[t, ]$end
          ko_seq <- substring(Gene, ko_start, ko_end)
          if (Gene_rev) {
            seq <- DNAString(ko_seq)
            seq_rev <- reverse(seq)
            ko_seq <- as.character(seq_rev)
          }
          #读取ko_seq的gRNA表格
          py$run(ko_seq)
          gRNA.table <- read.csv("gRNA.csv", header = FALSE)
          
          #特异性得分80以下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          for (i in 1:length(gRNA.table[, 1])) {
            if (gRNA.table[i, 3] <= 80) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          gRNA.table <- gRNA.table[-gRNA.del,]
          
          if (nrow(gRNA.table) < 2) {
            #筛选不到gRNA时要及时退出
            next
          }
          
          #获取每个gRNA在基因上的位置
          strand <- sub("[^a-zA-Z]+", "", gRNA.table[, 1])
          gRNA_seq <- substring(gRNA.table[, 2], 1, 20)
          Score1 <- gRNA.table[, 3]
          analysis_seq <-
            gsub(" ", "", substring(gRNA.table[, 2], 1, 24))
          gRNA.table <-
            cbind(strand, gRNA_seq, analysis_seq, Score1)
          gRNA.table <- as.data.frame(gRNA.table)
          gene <- DNAString(Gene)
          {
            if (Gene_rev) {
              for (i in 1:length(gRNA.table[, 1])) {
                if (gRNA.table[i, ]$strand == "fw") {
                  gRNA_rev <- DNAString(gRNA.table[i, ]$gRNA_seq)
                  gRNA_rev <- reverse(gRNA_rev)
                  rev <-
                    matchPattern(pattern = gRNA_rev, subject = gene)
                  rev_start <- start(rev)
                  rev_end <- end(rev)
                  gRNA.table[i, 5] <- rev_start
                  gRNA.table[i, 6] <- rev_end
                }
                else if (gRNA.table[i, ]$strand == "rev") {
                  gRNA_fw <- DNAString(gRNA.table[i, ]$gRNA_seq)
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
                if (gRNA.table[i,]$strand == "rev") {
                  gRNA_rev <- DNAString(gRNA.table[i,]$gRNA_seq)
                  gRNA_rev <- reverseComplement(gRNA_rev)
                  rev <-
                    matchPattern(pattern = gRNA_rev, subject = gene)
                  rev_start <- start(rev)
                  rev_end <- end(rev)
                  gRNA.table[i, 5] <- rev_start
                  gRNA.table[i, 6] <- rev_end
                }
                else if (gRNA.table[i,]$strand == "fw") {
                  fw <-
                    matchPattern(pattern = gRNA.table[i,]$gRNA_seq,
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
          GC_avoid_region <- GC_analysis2(KO_region[t, ])
          if (GC_avoid_region != FALSE) {
            GC_del <- numeric()
            for (i in 1:nrow(gRNA.table)) {
              for (j in 1:nrow(GC_avoid_region)) {
                if (gRNA.table[i, ]$start %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2]) |
                    gRNA.table[i, ]$end %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2])) {
                  GC_del <- append(GC_del, i)
                }
              }
            }
            if (length(GC_del) != 0) {
              gRNA.table <- gRNA.table[-GC_del,]
            }
          }
          if (nrow(gRNA.table) > 50) {
            gRNA.table <- gRNA.table[1:50,]
          }
          #符合条件的gRNA进行切割效率预测
          write.csv(gRNA.table, file = "CCTOP-predictor.csv", row.names = FALSE)
          source_python("C://Users//41518//Desktop//work//Ubigene//crispr_get_score(2)(1).py")
          py$reader_writer("C:\\Users\\41518\\Desktop\\work\\Ubigene\\CCTOP-predictor.csv")
          gRNA.table <-
            read.csv("CCTOP-predictor.csv", header = TRUE)
          
          #切割效率得分低于0.60的删除
          gRNA.table <-
            gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]
          
          #相差0.05分以内优选
          gRNA <- Get_result2(gRNA.table)
          
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA) == "NULL") {
            gRNA <- Get_result3(gRNA.table)
          }
          
          #第二对gRNA
          if (class(gRNA) != "NULL") {
            #在gRNA表格中删除第一对gRNA
            n=which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
            gRNA.table_2<-gRNA.table[-n,]
            #相差0.05分以内优选
            gRNA2<-Get_result2(gRNA.table_2)
            #如果没有相差0.05分以内的gRNA
            if (class(gRNA2) == "NULL") {
              gRNA2 <- Get_result3(gRNA.table_2)
            }
            
            if (class(gRNA2) == "NULL") {
              rm(gRNA2)
              
              #######换区域
            }
          }
          
          if (class(gRNA) == "NULL") {
            rm(gRNA)
          }
          
          if (exists("gRNA") == TRUE) {
            if (KO_region[t, ]$times != transcript.count) {
              No.transcript <-
                ko.data[which(ko.data$start == KO_region[t,]$start),]$transcript
              which.transcript <-
                transcript.table$Name %in% No.transcript
              Not_KO_transcript <-
                transcript.table$Name[which(which.transcript == "FALSE")]
              sprintf("%s转录本可能无法被影响到", Not_KO_transcript)
            }
          }
        }
      }
      if (exists("gRNA") == TRUE) {
        break
      }
    }
  }
  
  
  
  # 内含子小于500bp,大于300bp ------------------------------------------------------
  if (exists("gRNA") == FALSE) {
    KO_region3 <- KO_region_300(KO_region3)
    if (nrow(KO_region3) != 0) {
      for (a in 1:nrow(KO_region3)) {
        #外显子上游
        ko_start1 <- KO_region3[a, ]$start - KO_region3[a, ]$left + 1
        ko_end1 <- KO_region3[a, ]$start
        ko_seq1 <- substring(Gene, ko_start1, ko_end1)
        if (Gene_rev) {
          seq1 <- DNAString(ko_seq1)
          seq1_rev <- reverse(seq1)
          ko_seq1 <- as.character(seq1_rev)
        }
        #读取ko_seq的gRNA表格
        py$run(ko_seq1)
        gRNA.table1 <- read.csv("gRNA.csv", header = FALSE)
        #特异性得分80以下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        for (i in 1:length(gRNA.table1[, 1])) {
          if (gRNA.table1[i, 3] <= 80) {
            gRNA.del <- append(gRNA.del, i)
          }
          else if (grepl("Inefficient", gRNA.table1[i, 2])) {
            gRNA.del <- append(gRNA.del, i)
          }
        }
        gRNA.table1 <- gRNA.table1[-gRNA.del,]
        if (nrow(gRNA.table1) > 6) {
          gRNA.table1 <- gRNA.table1[1:6,]
        }
        #外显子下游
        ko_start2 <- KO_region3[a, ]$end
        ko_end2 <- KO_region3[a, ]$end + KO_region3[a, ]$right
        ko_seq2 <- substring(Gene, ko_start2, ko_end2)
        if (Gene_rev) {
          seq2 <- DNAString(ko_seq2)
          seq2_rev <- reverse(seq2)
          ko_seq2 <- as.character(seq2_rev)
        }
        #读取ko_seq的gRNA表格
        py$run(ko_seq2)
        gRNA.table2 <- read.csv("gRNA.csv", header = FALSE)
        #特异性得分80以下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        for (i in 1:length(gRNA.table2[, 1])) {
          if (gRNA.table2[i, 3] <= 80) {
            gRNA.del <- append(gRNA.del, i)
          }
          else if (grepl("Inefficient", gRNA.table2[i, 2])) {
            gRNA.del <- append(gRNA.del, i)
          }
        }
        gRNA.table2 <- gRNA.table2[-gRNA.del,]
        if (nrow(gRNA.table2) > 6) {
          gRNA.table2 <- gRNA.table2[1:6,]
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
                gRNA.table[i, 5] <- rev_start
                gRNA.table[i, 6] <- rev_end
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
        GC_avoid_region <- GC_analysis2(KO_region3[a,])
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
        source_python("C://Users//41518//Desktop//work//Ubigene//crispr_get_score(2)(1).py")
        py$reader_writer("C:\\Users\\41518\\Desktop\\work\\Ubigene\\CCTOP-predictor.csv")
        gRNA.table <- read.csv("CCTOP-predictor.csv", header = TRUE)
        
        #切割效率得分低于0.60的删除
        gRNA.table <-
          gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
        
        #上下游分开
        gRNA.table1 <-
          gRNA.table[which(gRNA.table$end <= KO_region3[a,]$start),]
        gRNA.table2 <-
          gRNA.table[which(gRNA.table$start >= KO_region3[a,]$end),]
        #table1是外显子上游的gRNA,table2是外显子下游的gRNA
        {
          if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
            next
          }
          else{
            gRNA <- Get_result1(gRNA.table,KO_region[a, ])        #相差0.05分以内优选
          }
        }
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA) == "NULL") {
          gRNA.table1 <-
            gRNA.table1[order(-gRNA.table1$crispr_score,-gRNA.table1$Score1),]
          gRNA.table2 <-
            gRNA.table2[order(-gRNA.table2$crispr_score,-gRNA.table2$Score1),]
          gRNA <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
        }
        
        #第二对gRNA
        if (class(gRNA) != "NULL" &
            nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
          #在gRNA表格中删除第一对gRNA
          n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
          gRNA.table_2 <- gRNA.table[-n, ]
          #相差0.05分以内优选
          gRNA2 <- Get_result1(gRNA.table_2, KO_region[a, ])
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA2) == "NULL") {
            if (gRNA.table1[1, ]$analysis_seq %in% gRNA$analysis_seq) {
              if (gRNA.table2[1, ]$analysis_seq %in% gRNA$analysis_seq) {
                gRNA2 <- rbind(gRNA.table1[2, ], gRNA.table2[2, ])
              }
              else{
                gRNA2 <- rbind(gRNA.table1[2, ], gRNA.table2[1, ])
              }
            }
            else{
              if (gRNA.table2[1, ]$analysis_seq %in% gRNA$analysis_seq) {
                gRNA2 <- rbind(gRNA.table1[1, ], gRNA.table2[2, ])
              }
              else{
                gRNA2 <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
              }
            }
          }
          
          if (class(gRNA2) == "NULL") {
            rm(gRNA2)
            #######换区域
          }
        }
        
        if (exists("gRNA") == TRUE) {
          if (KO_region3[a,]$times != transcript.count) {
            No.transcript <-
              ko.data[which(ko.data$start == KO_region3_1200[t, ]$start), ]$transcript
            which.transcript <-
              transcript.table$Name %in% No.transcript
            Not_KO_transcript <-
              transcript.table$Name[which(which.transcript == "FALSE")]
            sprintf("%s转录本可能无法被影响到", Not_KO_transcript)
          }
        }
        if (exists("gRNA") == TRUE) {
          judge3 <- "TRUE"
          break
        }
      }
    }
  }
  
  
  
  
  # 整个敲除 --------------------------------------------------------------------
  if (exists("gRNA") == FALSE) {
    KO_region2 <- KO_longregion(t_Exon_CDS)
    if (nrow(KO_region2) != 0) {
      #外显子上游
      ko_start1 <- KO_region2$start - 400
      ko_end1 <- KO_region2$start
      if(ko_start1<0){
        KO_region2$start<-min(t_Exon_CDS$start)
      }
      ko_start1 <- KO_region2$start - 400
      ko_end1 <- KO_region2$start
      ko_seq1 <- substring(Gene, ko_start1, ko_end1)
      if (Gene_rev) {
        seq1 <- DNAString(ko_seq1)
        seq1_rev <- reverse(seq1)
        ko_seq1 <- as.character(seq1_rev)
      }
      #读取ko_seq的gRNA表格
      py$run(ko_seq1)
      gRNA.table1 <- read.csv("gRNA.csv", header = FALSE)
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
      gRNA.table1 <- gRNA.table1[-gRNA.del,]
      if (nrow(gRNA.table1) > 7) {
        gRNA.table1 <- gRNA.table1[1:7, ]
      }
      
      #外显子下游
      ko_start2 <- KO_region2$end
      ko_end2 <- KO_region2$end + 400
      if(ko_end2>nchar(Gene)){
        KO_region2$end<-max(t_Exon_CDS$end)
      }
      ko_start2 <- KO_region2$end
      ko_end2 <- KO_region2$end + 400
      ko_seq2 <- substring(Gene, ko_start2, ko_end2)
      if (Gene_rev) {
        seq2 <- DNAString(ko_seq2)
        seq2_rev <- reverse(seq2)
        ko_seq2 <- as.character(seq2_rev)
      }
      #读取ko_seq的gRNA表格
      py$run(ko_seq2)
      gRNA.table2 <- read.csv("gRNA.csv", header = FALSE)
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
      gRNA.table2 <- gRNA.table2[-gRNA.del,]
      if (nrow(gRNA.table2) > 7) {
        gRNA.table2 <- gRNA.table2[1:7, ]
      }
      if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
        #筛选不到gRNA时要及时退出
        break
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
            if (gRNA.table[i, ]$strand == "fw") {
              gRNA_rev <- DNAString(gRNA.table[i, ]$gRNA_seq)
              gRNA_rev <- reverse(gRNA_rev)
              rev <-
                matchPattern(pattern = gRNA_rev, subject = gene)
              rev_start <- start(rev)
              rev_end <- end(rev)
              gRNA.table[i, 5] <- rev_start
              gRNA.table[i, 6] <- rev_end
            }
            else if (gRNA.table[i, ]$strand == "rev") {
              gRNA_fw <- DNAString(gRNA.table[i, ]$gRNA_seq)
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
            if (gRNA.table[i,]$strand == "rev") {
              gRNA_rev <- DNAString(gRNA.table[i,]$gRNA_seq)
              gRNA_rev <- reverseComplement(gRNA_rev)
              rev <-
                matchPattern(pattern = gRNA_rev, subject = gene)
              rev_start <- start(rev)
              rev_end <- end(rev)
              gRNA.table[i, 5] <- rev_start[1]
              gRNA.table[i, 6] <- rev_end[1]
            }
            else if (gRNA.table[i,]$strand == "fw") {
              fw <-
                matchPattern(pattern = gRNA.table[i,]$gRNA_seq,
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
      GC_avoid_region <- GC_analysis2(KO_region2)
      if (GC_avoid_region != FALSE) {
        GC_del <- numeric()
        for (i in 1:nrow(gRNA.table)) {
          for (j in 1:nrow(GC_avoid_region)) {
            if (gRNA.table[i, ]$start %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2]) |
                gRNA.table[i, ]$end %in% c(GC_avoid_region[j, 1]:GC_avoid_region[j, 2])) {
              GC_del <- append(GC_del, i)
            }
          }
        }
        if (length(GC_del) != 0) {
          gRNA.table <- gRNA.table[-GC_del,]
        }
      }
      
      #符合条件的gRNA进行切割效率预测
      write.csv(gRNA.table, file = "CCTOP-predictor.csv", row.names = FALSE)
      source_python("C://Users//41518//Desktop//work//Ubigene//crispr_get_score(2)(1).py")
      py$reader_writer("C:\\Users\\41518\\Desktop\\work\\Ubigene\\CCTOP-predictor.csv")
      gRNA.table <- read.csv("CCTOP-predictor.csv", header = TRUE)
      
      #切割效率得分低于0.60的删除
      gRNA.table <-
        gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]
      
      #上下游分开
      gRNA.table1 <-
        gRNA.table[which(gRNA.table$end <= KO_region2$start), ]
      gRNA.table2 <-
        gRNA.table[which(gRNA.table$start >= KO_region2$end), ]
      #table1是外显子上游的gRNA,table2是外显子下游的gRNA
      {
        if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
          next
        }
        else{
          gRNA <- Get_result1(gRNA.table, KO_region2)        #相差0.05分以内优选
        }
      }
      #如果没有相差0.05分以内的gRNA
      if (class(gRNA) == "NULL") {
        gRNA.table1 <-
          gRNA.table1[order(-gRNA.table1$crispr_score, -gRNA.table1$Score1), ]
        gRNA.table2 <-
          gRNA.table2[order(-gRNA.table2$crispr_score, -gRNA.table2$Score1), ]
        gRNA <- rbind(gRNA.table1[1,], gRNA.table2[1,])
      }
      
      #第二对gRNA
      if (class(gRNA) != "NULL" &
          nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
        #在gRNA表格中删除第一对gRNA
        n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
        gRNA.table_2 <- gRNA.table[-n, ]
        #相差0.05分以内优选
        gRNA2 <- Get_result1(gRNA.table_2, KO_region2)
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA2) == "NULL") {
          if (gRNA.table1[1, ]$analysis_seq %in% gRNA$analysis_seq) {
            if (gRNA.table2[1, ]$analysis_seq %in% gRNA$analysis_seq) {
              gRNA2 <- rbind(gRNA.table1[2, ], gRNA.table2[2, ])
            }
            else{
              gRNA2 <- rbind(gRNA.table1[2, ], gRNA.table2[1, ])
            }
          }
          else{
            if (gRNA.table2[1, ]$analysis_seq %in% gRNA$analysis_seq) {
              gRNA2 <- rbind(gRNA.table1[1, ], gRNA.table2[2, ])
            }
            else{
              gRNA2 <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
            }
          }
        }
        
        if (class(gRNA2) == "NULL") {
          rm(gRNA2)
          #######换区域
        }
      }
      if (exists("gRNA") == TRUE) {
        judge2 <- "TRUE"
      }
    }
  }
  
  
  
  # 画图 ----------------------------------------------------------------------
  {
    if (exists("gRNA") == TRUE) {
      y <- 0.5
      f <- data.frame(x = c(1:nchar(Gene)), y = y)
      p1 <-
        ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "grey", size = 4.5) + theme_bw() +
        theme(panel.grid = element_blank(), panel.border = element_blank()) +
        ylim(0, 1)
      
      for (i in 1:nrow(t_Exon_region)) {
        start <- as.numeric(t_Exon_region[i,]$Exon_start)
        end <- as.numeric(t_Exon_region[i,]$Exon_end)
        p1 <- p1 + annotate(
          "rect",
          xmin = start,
          xmax = end,
          ymin = y - 0.04,
          ymax = y + 0.04,
          colour = "orange2",
          alpha = .0
        ) + annotate(
          "text",
          label = i,
          x = mean(c(start, end)),
          y = y - 0.05,
          size = 3
        )
      }  #最长的转录本的展示
      for (i in 1:nrow(t_Exon_CDS)) {
        start <- as.numeric(t_Exon_CDS[i,]$start)
        end <- as.numeric(t_Exon_CDS[i,]$end)
        p1 <- p1 + annotate(
          "rect",
          xmin = start,
          xmax = end,
          ymin = y - 0.04,
          ymax = y + 0.04,
          fill = "orange2",
          alpha = .7
        )
      }
      
      ff <-
        data.frame(x = rep(min(gRNA[1, ]$start, gRNA[2, ]$start), 2), y = c(0.51, 0.57))
      fff <-
        data.frame(x = rep(max(gRNA[1, ]$end, gRNA[2, ]$end), 2), y = c(0.49, 0.43))
      p1 <-
        p1 + geom_line(
          data = ff,
          aes(x = x, y = y),
          arrow = arrow(
            length = unit(0.15, "cm"),
            ends = "first",
            type = "closed"
          ),
          color = "red"
        ) +
        geom_line(
          data = fff,
          aes(x = x, y = y),
          arrow = arrow(
            length = unit(0.15, "cm"),
            ends = "first",
            type = "closed"
          ),
          color = "red"
        )
      p1 <-
        p1 + labs(title = "Overview of the Targeting Strategy") + theme(plot.title = element_text(hjust = 0.5,                                                                                          size = 25))
      
      jpeg(file = "方案二.png")
      print(p1)
      dev.off()
      print(p1)
    }
    else{
      print("方案设计失败")
    }
  }
  
  
  # 输出 ----------------------------------------------------------------------
  #敲除大小
  if (exists("gRNA") == TRUE) {
    if (all(gRNA$strand == "rev") |
        all(gRNA$strand == "fw")) {
      KO_length <- abs(gRNA[1,]$end - gRNA[2,]$end)
    }
    else{
      if (gRNA[1, ]$strand == "rev" & gRNA[2, ]$strand == "fw") {
        pos1 <- gRNA[1, ]$start + 3
        pos2 <- gRNA[2, ]$end - 3
        KO_length <- abs(pos1 - pos2)
      }
      else{
        pos1 <- gRNA[1, ]$end - 3
        pos2 <- gRNA[2, ]$start + 3
        KO_length <- abs(pos1 - pos2)
      }
    }
    
    #敲除的CDS
    if (exists("judge")) {
      KO_length_CDS <- KO_region[t,]$Exon_length
    }
    else if (exists("judge2")) {
      KO_length_CDS <- KO_region2$Exon_length
    }
    else if (exists("judge3")) {
      KO_length_CDS <- KO_region3[t,]$Exon_length
    }
    else{
      KO_length_CDS <- KO_length
    }
  }
  
  if (exists("gRNA2") == TRUE) {
    if (all(gRNA2$strand == "rev") |
        all(gRNA2$strand == "fw")) {
      KO_length2 <- abs(gRNA2[1,]$end - gRNA2[2,]$end)
    }
    else{
      if (gRNA2[1, ]$strand == "rev" & gRNA2[2, ]$strand == "fw") {
        pos1 <- gRNA2[1, ]$start + 3
        pos2 <- gRNA2[2, ]$end - 3
        KO_length2 <- abs(pos1 - pos2)
      }
      else{
        pos1 <- gRNA2[1, ]$end - 3
        pos2 <- gRNA2[2, ]$start + 3
        KO_length2 <- abs(pos1 - pos2)
      }
    }
    
    #敲除的CDS
    if (exists("judge")) {
      KO_length_CDS2 <- KO_region[t,]$Exon_length
    }
    else if (exists("judge2")) {
      KO_length_CDS2 <- KO_region2$Exon_length
    }
    else if (exists("judge3")) {
      KO_length_CDS2 <- KO_region3[t,]$Exon_length
    }
    else{
      KO_length_CDS2 <- KO_length2
    }
  }
  
  if (exists("gRNA") == TRUE) {
    output <- read.csv("C://Users//41518//Desktop//测试//第四轮测试.csv")
    i <- nrow(output)
    output[i + 1, 1] <- term
    output[i + 1, 2] <- gRNA[1,]$analysis_seq
    output[i + 1, 3] <- gRNA[1,]$Score1
    output[i + 1, 4] <- gRNA[1,]$crispr_score
    output[i + 1, 5] <- gRNA[2,]$analysis_seq
    output[i + 1, 6] <- gRNA[2,]$Score1
    output[i + 1, 7] <- gRNA[2,]$crispr_score
    output[i + 1, 8] <- KO_length
    output[i + 1, 9] <- KO_length_CDS
    if(exists("gRNA2")==TRUE){
      output[i + 1, 10] <- gRNA2[1,]$analysis_seq
      output[i + 1, 11] <- gRNA2[1,]$Score1
      output[i + 1, 12] <- gRNA2[1,]$crispr_score
      output[i + 1, 13] <- gRNA2[2,]$analysis_seq
      output[i + 1, 14] <- gRNA2[2,]$Score1
      output[i + 1, 15] <- gRNA2[2,]$crispr_score
      output[i + 1, 16] <- KO_length2
      output[i + 1, 17] <- KO_length_CDS2
    }
    write.csv(output, file = "C://Users//41518//Desktop//测试//第四轮测试.csv", row.names = FALSE)
  }
  print(gRNA)
  print(gRNA$analysis_seq)
  print(KO_length)
  print(KO_length_CDS)
  print(Not_KO_transcript)
  
  # rm(gRNA)
  # rm(gRNA2)
  # rm(judge)
  # rm(judge2)
  # rm(judge3)
  # rm(KO_region)
  # rm(KO_region3)
  # rm(Not_KO_transcript)
  # rm(ko_seq)
  # rm(ko_seq1)
  # rm(ko_seq2)
}