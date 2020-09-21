#换区域找第二对gRNA
Get_gRNA2 <- function(KO_region) {
  if (nrow(KO_region) > t) {
    for (h in (t + 1):nrow(KO_region)) {
      # 设计在起始密码子后面 -----------------------------------------------------------------
      for (j in 1:nrow(Exon_CDS)) {
        if (all(KO_region[h, ] %in% Exon_CDS[j, ])) {
          ko_start <- KO_region[h, ]$start + 3
          ko_end <- KO_region[h, ]$end
          ko_seq <- substring(Gene, ko_start, ko_end)
          if (Gene_rev) {
            seq <- DNAString(ko_seq)
            seq_rev <- reverse(seq)
            ko_seq <- as.character(seq_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq,species,filepath,"1")
          gRNA.table <-read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          
          #特异性得分60以下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table<-gRNA.table[which(gRNA.table$V3!="No matches"),]
          for (i in 1:length(gRNA.table[, 1])) {
            if (as.numeric(gRNA.table[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table[i, 5]) < 40){
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table <- gRNA.table[-gRNA.del,]
          }
          if(nrow(gRNA.table)==0){
            next
          }
          #0-0-0(优化)
          count_0 <- numeric()
          for (p in 1:nrow(gRNA.table)) {
            if (gRNA.table[p,]$V9 == "0-0-0") {
              count_0 <- append(count_0, p)
            }
          }
          if (length(count_0) != 0) {
            gRNA.table_min <- gRNA.table[-count_0,]
            gRNA.table <-rbind(gRNA.table[count_0,], gRNA.table_min)
          }
          if (nrow(gRNA.table) < 2) {
            #筛选不到gRNA时要及时退出
            next
          }
          
          #获取每个gRNA在基因上的位置
          strand <- sub("[^a-zA-Z]+", "", gRNA.table[, 1])
          gRNA_seq <- substring(gRNA.table[, 2], 1, 20)
          Score1 <- gRNA.table[, 3]
          crispr_score<-as.numeric(gRNA.table[, 5])
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
          gRNA.table<-cbind(gRNA.table,crispr_score)
          #局部GC含量大于80%或小于25%，避免该区域
          GC_avoid_region <- GC_analysis2(KO_region[h, ])
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

          gRNA.table_85<-gRNA.table[which(gRNA.table$Score1>=85),]
          gRNA.table_70<-gRNA.table[which(gRNA.table$Score1>=70 & gRNA.table$Score1<85),]
          gRNA.table_60<-gRNA.table[which(gRNA.table$Score1<70),]
          
          gRNA.table<-rbind(gRNA.table_85[order(gRNA.table_85$crispr_score,decreasing = TRUE),],
                            gRNA.table_70[order(gRNA.table_70$crispr_score,decreasing = TRUE),],
                            gRNA.table_60[order(gRNA.table_60$crispr_score,decreasing = TRUE),])
        
          #70分
          gRNA2 <- Get_result2(gRNA.table)
          
          #60分
          if (class(gRNA2) == "NULL") {
            gRNA2 <- Get_result3(gRNA.table)
          }

          
          if (class(gRNA2) == "NULL") {
            rm(gRNA2)
          }
          
          if (exists("gRNA2") == TRUE) {
            break
          }
        }
      }
      if (exists("gRNA2") == TRUE) {
        break
      }
      
      # 设计在内含子上 -----------------------------------------------------------------
      for (j in 1:nrow(KO_region_500)) {
        if (all(KO_region[h, ] %in% KO_region_500[j, ])) {
          #外显子上游
          ko_start1 <- KO_region[h,]$start - 400
          ko_end1 <- KO_region[h, ]$start
          ko_seq1 <- substring(Gene, ko_start1, ko_end1)
          if (Gene_rev) {
            seq1 <- DNAString(ko_seq1)
            seq1_rev <- reverse(seq1)
            ko_seq1 <- as.character(seq1_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq1,species,filepath,"1")
          gRNA.table1 <-read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60下，Repeat和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table1<-gRNA.table1[which(gRNA.table1$V3!="No matches"),]
          for (i in 1:length(gRNA.table1[, 1])) {
            if (as.numeric(gRNA.table1[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table1[i, 5]) < 40){
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table1[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table1 <- gRNA.table1[-gRNA.del,]
          }
          
          if(nrow(gRNA.table1)==0){
            next
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

          
          #外显子下游
          ko_start2 <- KO_region[h, ]$end
          ko_end2 <- KO_region[h, ]$end + 400
          ko_seq2 <- substring(Gene, ko_start2, ko_end2)
          if (Gene_rev) {
            seq2 <- DNAString(ko_seq2)
            seq2_rev <- reverse(seq2)
            ko_seq2 <- as.character(seq2_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq2,species,filepath,"1")
          gRNA.table2 <-
            read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table2<-gRNA.table2[which(gRNA.table2$V3!="No matches"),]
          for (i in 1:length(gRNA.table2[, 1])) {
            if (as.numeric(gRNA.table2[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table2[i, 5]) < 40){
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table2[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table2 <- gRNA.table2[-gRNA.del,]
          }
          
          if(nrow(gRNA.table2)==0){
            next
          }
          #0-0-0(优化)
          count_0 <- numeric()
          for (p in 1:nrow(gRNA.table2)) {
            if (gRNA.table2[p,]$V9 == "0-0-0") {
              count_0 <- append(count_0, p)
            }
          }
          if (length(count_0) != 0) {
            gRNA.table_min <- gRNA.table2[-count_0,]
            gRNA.table2 <-rbind(gRNA.table2[count_0,], gRNA.table_min)
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
          crispr_score<-as.numeric(gRNA.table[, 5])
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
          gRNA.table<-cbind(gRNA.table,crispr_score)
          #局部GC含量大于80%或小于25%，避免该区域
          GC_avoid_region <- GC_analysis2(KO_region[h, ])
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
          
          gRNA.table_85<-gRNA.table[which(gRNA.table$Score1>=85),]
          gRNA.table_70<-gRNA.table[which(gRNA.table$Score1>=70 & gRNA.table$Score1<85),]
          gRNA.table_60<-gRNA.table[which(gRNA.table$Score1<70),]
          
          gRNA.table<-rbind(gRNA.table_85[order(gRNA.table_85$crispr_score,decreasing = TRUE),],
                            gRNA.table_70[order(gRNA.table_70$crispr_score,decreasing = TRUE),],
                            gRNA.table_60[order(gRNA.table_60$crispr_score,decreasing = TRUE),])
          
          
          #上下游分开
          gRNA.table1 <-
            gRNA.table[which(gRNA.table$end <= KO_region[h, ]$start), ]
          gRNA.table2 <-
            gRNA.table[which(gRNA.table$start >= KO_region[h, ]$end), ]
 
          
          #table1是外显子上游的gRNA,table2是外显子下游的gRNA
          {
            if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
              next
            }
            else{
              gRNA2 <- Get_result1(gRNA.table, KO_region[h,])        #相差0.05分以内优选
            }
          }
          
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA2) == "NULL") {
            gRNA2 <- rbind(gRNA.table1[1,], gRNA.table2[1,])
          }
          
          if (class(gRNA2) == "NULL") {
            rm(gRNA2)
          }
          
          if (exists("gRNA2") == TRUE) {
            break
          }
        }
      }
      if (exists("gRNA2") == TRUE) {
        judge_2 <- "TRUE"
        break
      }
      
      # 设计在外显子上 -----------------------------------------------------------------
      for (j in 1:nrow(Exon_300)) {
        if (all(KO_region[h, ] %in% Exon_300[j, ])) {
          ko_start <- KO_region[h, ]$start
          ko_end <- KO_region[h, ]$end
          ko_seq <- substring(Gene, ko_start, ko_end)
          if (Gene_rev) {
            seq <- DNAString(ko_seq)
            seq_rev <- reverse(seq)
            ko_seq <- as.character(seq_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq,species,filepath,"1")
          gRNA.table <-read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          
          #特异性得分60以下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table<-gRNA.table[which(gRNA.table$V3!="No matches"),]
          for (i in 1:length(gRNA.table[, 1])) {
            if (as.numeric(gRNA.table[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table[i, 5]) < 40){
              gRNA.del <- append(gRNA.del, i)
            }
            
            else if (grepl("Inefficient", gRNA.table[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table <- gRNA.table[-gRNA.del,]
          }
          
          if(nrow(gRNA.table)==0){
            next
          }
          #0-0-0(优化)
          count_0 <- numeric()
          for (p in 1:nrow(gRNA.table)) {
            if (gRNA.table[p,]$V9 == "0-0-0") {
              count_0 <- append(count_0, p)
            }
          }
          if (length(count_0) != 0) {
            gRNA.table_min <- gRNA.table[-count_0,]
            gRNA.table <-rbind(gRNA.table[count_0,], gRNA.table_min)
          }
          if (nrow(gRNA.table) < 2) {
            #筛选不到gRNA时要及时退出
            next
          }
          
          #获取每个gRNA在基因上的位置
          strand <- sub("[^a-zA-Z]+", "", gRNA.table[, 1])
          gRNA_seq <- substring(gRNA.table[, 2], 1, 20)
          Score1 <- gRNA.table[, 3]
          crispr_score<-as.numeric(gRNA.table[, 5])
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
          gRNA.table<-cbind(gRNA.table,crispr_score)
          #局部GC含量大于80%或小于25%，避免该区域
          GC_avoid_region <- GC_analysis2(KO_region[h, ])
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

          
          gRNA.table_85<-gRNA.table[which(gRNA.table$Score1>=85),]
          gRNA.table_70<-gRNA.table[which(gRNA.table$Score1>=70 & gRNA.table$Score1<85),]
          gRNA.table_60<-gRNA.table[which(gRNA.table$Score1<70),]
          
          gRNA.table<-rbind(gRNA.table_85[order(gRNA.table_85$crispr_score,decreasing = TRUE),],
                            gRNA.table_70[order(gRNA.table_70$crispr_score,decreasing = TRUE),],
                            gRNA.table_60[order(gRNA.table_60$crispr_score,decreasing = TRUE),])
          
          #70分
          gRNA2 <- Get_result2(gRNA.table)
          
          #60分
          if (class(gRNA2) == "NULL") {
            gRNA2 <- Get_result3(gRNA.table)
          }
          
          if (class(gRNA2) == "NULL") {
            rm(gRNA2)
          }
          if (exists("gRNA2") == TRUE) {
            break
          }
        }
      }
      if (exists("gRNA2") == TRUE) {
        break
      }
    }
  }
  
  # 内含子小于500bp,大于300bp ------------------------------------------------------
  if (exists("gRNA2") == FALSE) {
    if (nrow(KO_region_4) != 0) {
      KO_region3 <- KO_region_300(KO_region_4)
      if (exists("a") == FALSE) {
        a <- 0
      }
      if (nrow(KO_region3) > a) {
        for (h in (a + 1):nrow(KO_region3)) {
          #外显子上游
          ko_start1 <-
            KO_region3[h,]$start - KO_region3[h,]$left + 1
          ko_end1 <- KO_region3[h,]$start
          ko_seq1 <- substring(Gene, ko_start1, ko_end1)
          if (Gene_rev) {
            seq1 <- DNAString(ko_seq1)
            seq1_rev <- reverse(seq1)
            ko_seq1 <- as.character(seq1_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq1, species,filepath,"1")
          gRNA.table1 <-read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60以下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table1<-gRNA.table1[which(gRNA.table1$V3!="No matches"),]
          for (i in 1:length(gRNA.table1[, 1])) {
            if (as.numeric(gRNA.table1[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table1[i, 5]) < 40){
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
            next
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
            gRNA.table1 <-rbind(gRNA.table1[count_0, ], gRNA.table_min)
          }
          if (nrow(gRNA.table1) > 6) {
            gRNA.table1 <- gRNA.table1[1:6, ]
          }
          #外显子下游
          ko_start2 <- KO_region3[h,]$end
          ko_end2 <- KO_region3[h,]$end + KO_region3[h,]$right
          ko_seq2 <- substring(Gene, ko_start2, ko_end2)
          if (Gene_rev) {
            seq2 <- DNAString(ko_seq2)
            seq2_rev <- reverse(seq2)
            ko_seq2 <- as.character(seq2_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq2, species,filepath,"1")
          gRNA.table2 <-read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          #特异性得分60以下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table2<-gRNA.table2[which(gRNA.table2$V3!="No matches"),]
          for (i in 1:length(gRNA.table2[, 1])) {
            if (as.numeric(gRNA.table2[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table2[i, 5]) < 40){
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table2[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if(length(gRNA.del)!=0){
            gRNA.table2 <- gRNA.table2[-gRNA.del, ]
          }
          if (nrow(gRNA.table2) == 0) {
            next
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
            gRNA.table2 <-rbind(gRNA.table2[count_0, ], gRNA.table_min)
          }
          if (nrow(gRNA.table2) > 6) {
            gRNA.table2 <- gRNA.table2[1:6, ]
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
          crispr_score<-as.numeric(gRNA.table[, 5])
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
          gRNA.table<-cbind(gRNA.table,crispr_score)
          #局部GC含量大于80%或小于25%，避免该区域
          GC_avoid_region <- GC_analysis2(KO_region3[h, ])
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
          
          gRNA.table_85<-gRNA.table[which(gRNA.table$Score1>=85),]
          gRNA.table_70<-gRNA.table[which(gRNA.table$Score1>=70 & gRNA.table$Score1<85),]
          gRNA.table_60<-gRNA.table[which(gRNA.table$Score1<70),]
          
          gRNA.table<-rbind(gRNA.table_85[order(gRNA.table_85$crispr_score,decreasing = TRUE),],
                            gRNA.table_70[order(gRNA.table_70$crispr_score,decreasing = TRUE),],
                            gRNA.table_60[order(gRNA.table_60$crispr_score,decreasing = TRUE),])
          
          #上下游分开
          gRNA.table1 <-
            gRNA.table[which(gRNA.table$end <= KO_region3[h, ]$start), ]
          gRNA.table2 <-
            gRNA.table[which(gRNA.table$start >= KO_region3[h, ]$end), ]
  
          #table1是外显子上游的gRNA,table2是外显子下游的gRNA
          {
            if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
              next
            }
            else{
              gRNA2 <- Get_result1(gRNA.table, KO_region3[h, ])     #相差0.05分以内优选
            }
          }
          #如果没有相差0.05分以内的gRNA2
          if (class(gRNA2) == "NULL") {
            gRNA2 <- rbind(gRNA.table1[1,], gRNA.table2[1,])
          }
          
          if (class(gRNA2) == "NULL") {
            rm(gRNA2)
          }
          if (exists("gRNA2") == TRUE) {
            judge3_2 <- "TRUE"
            break
          }
        }
      }
    }
  }
  
  # 小基因全敲 -------------------------------------------------------------------
  if (exists("gRNA2") == FALSE) {
    if (exists("KO_region_all")) {
      if (exists("judge4") == FALSE) {
        if (all(KO_region_all %in% KO_region[nrow(KO_region),])) {
          print("all")
          KO_region_all$start <- KO_region_all$start + 500
          KO_region_all$end <- KO_region_all$end + 500
          #外显子上游
          ko_start1 <- KO_region_all$start - 400
          ko_end1 <- KO_region_all$start
          ko_seq1 <- substring(Gene2, ko_start1, ko_end1)
          if (Gene_rev) {
            seq1 <- DNAString(ko_seq1)
            seq1_rev <- reverse(seq1)
            ko_seq1 <- as.character(seq1_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq1, species, filepath,"1")
          gRNA.table1 <-
            read.csv(paste0(filepath, "//gRNA1.csv"),
                     header = FALSE,
                     encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table1 <-
            gRNA.table1[which(gRNA.table1$V3 != "No matches"), ]
          for (i in 1:length(gRNA.table1[, 1])) {
            if (as.numeric(gRNA.table1[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table1[i, 5]) < 40){
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
            next
          }
          #0-0-0(优化)
          count_0 <- numeric()
          for (p in 1:nrow(gRNA.table1)) {
            if (gRNA.table1[p, ]$V9 == "0-0-0") {
              count_0 <- append(count_0, p)
            }
          }
          if (length(count_0) != 0) {
            gRNA.table_min <- gRNA.table1[-count_0,]
            gRNA.table1 <-
              rbind(gRNA.table1[count_0,], gRNA.table_min)
          }

          
          #外显子下游
          ko_start2 <- KO_region_all$end
          ko_end2 <- KO_region_all$end + 400
          ko_seq2 <- substring(Gene2, ko_start2, ko_end2)
          if (Gene_rev) {
            seq2 <- DNAString(ko_seq2)
            seq2_rev <- reverse(seq2)
            ko_seq2 <- as.character(seq2_rev)
          }
          #读取ko_seq的gRNA表格
          source_python("crispor_table_download.py")
          py$run(ko_seq2, species, filepath,"1")
          gRNA.table2 <-read.csv(paste0(filepath, "//gRNA1.csv"),header = FALSE,encoding = "UTF-8")
          #特异性得分60下，和Inefficient的gRNA排除掉
          gRNA.del <- numeric()
          gRNA.table2 <-
            gRNA.table2[which(gRNA.table2$V3 != "No matches"), ]
          for (i in 1:length(gRNA.table2[, 1])) {
            if (as.numeric(gRNA.table2[i, 3]) < 60) {
              gRNA.del <- append(gRNA.del, i)
            }
            else if(as.numeric(gRNA.table2[i, 5]) < 40){
              gRNA.del <- append(gRNA.del, i)
            }
            else if (grepl("Inefficient", gRNA.table2[i, 2])) {
              gRNA.del <- append(gRNA.del, i)
            }
          }
          if (length(gRNA.del) != 0) {
            gRNA.table2 <- gRNA.table2[-gRNA.del,]
          }
          if (nrow(gRNA.table2) == 0) {
            next
          }
          #0-0-0(优化)
          count_0 <- numeric()
          for (p in 1:nrow(gRNA.table2)) {
            if (gRNA.table2[p, ]$V9 == "0-0-0") {
              count_0 <- append(count_0, p)
            }
          }
          if (length(count_0) != 0) {
            gRNA.table_min <- gRNA.table2[-count_0,]
            gRNA.table2 <-
              rbind(gRNA.table2[count_0,], gRNA.table_min)
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
          crispr_score<-as.numeric(gRNA.table[, 5])
          analysis_seq <-
            gsub(" ", "", substring(gRNA.table[, 2], 1, 24))
          gRNA.table <-
            cbind(strand, gRNA_seq, analysis_seq, Score1)
          gRNA.table <- as.data.frame(gRNA.table)
          gene <- DNAString(Gene2)
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
          gRNA.table<-cbind(gRNA.table,crispr_score)
          #局部GC含量大于80%或小于25%，避免该区域
          GC_avoid_region <- GC_analysis2(KO_region_all)
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
          
          
          gRNA.table_85<-gRNA.table[which(gRNA.table$Score1>=85),]
          gRNA.table_70<-gRNA.table[which(gRNA.table$Score1>=70 & gRNA.table$Score1<85),]
          gRNA.table_60<-gRNA.table[which(gRNA.table$Score1<70),]
          
          gRNA.table<-rbind(gRNA.table_85[order(gRNA.table_85$crispr_score,decreasing = TRUE),],
                            gRNA.table_70[order(gRNA.table_70$crispr_score,decreasing = TRUE),],
                            gRNA.table_60[order(gRNA.table_60$crispr_score,decreasing = TRUE),])
          
          
          #上下游分开
          gRNA.table1 <-
            gRNA.table[which(gRNA.table$end <= KO_region_all$start), ]
          gRNA.table2 <-
            gRNA.table[which(gRNA.table$start >= KO_region_all$end), ]
       
          #table1是外显子上游的gRNA,table2是外显子下游的gRNA
          {
            if (nrow(gRNA.table1) == 0 | nrow(gRNA.table2) == 0) {
              next
            }
            else{
              gRNA2 <- Get_result1(gRNA.table, KO_region_all)        #相差0.05分以内优选
            }
          }
          
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA2) == "NULL") {
            gRNA2 <- rbind(gRNA.table1[1,], gRNA.table2[1,])
          }
          if (class(gRNA2) == "NULL") {
            rm(gRNA2)
          }
          if (exists("gRNA2") == TRUE) {
            judge4_2 <- "TRUE"
          }
        }
      }
    }
  }
  
  
  # 输出 ----------------------------------------------------------------------
  #敲除大小
  if (exists("gRNA2") == TRUE) {
    if (Gene_rev) {
      if (all(gRNA2$strand == "rev") |
          all(gRNA2$strand == "fw")) {
        KO_length2 <- abs(gRNA2[1, ]$end - gRNA2[2, ]$end)
      }
      else{
        if (gRNA2[1,]$strand == "rev" & gRNA2[2,]$strand == "fw") {
          if (gRNA2[1,]$start > gRNA2[2,]$start) {
            pos1 <- gRNA2[1,]$end
            pos2 <- gRNA2[2,]$start
            KO_length2 <- abs(pos1 - pos2) + 1
          }
          else{
            pos1 <- gRNA2[1,]$end
            pos2 <- gRNA2[2,]$start
            KO_length2 <- abs(pos1 - pos2) - 1
          }
        }
        else{
          if (gRNA2[1,]$start > gRNA2[2,]$start) {
            pos1 <- gRNA2[1, ]$start
            pos2 <- gRNA2[2, ]$end
            KO_length2 <- abs(pos1 - pos2) - 1
          }
          else{
            pos1 <- gRNA2[1, ]$start
            pos2 <- gRNA2[2, ]$end
            KO_length2 <- abs(pos1 - pos2) + 1
          }
        }
      }
    }
    #正向
    else{
      if (all(gRNA2$strand == "rev") |
          all(gRNA2$strand == "fw")) {
        KO_length2 <- abs(gRNA2[1, ]$end - gRNA2[2, ]$end)
      }
      else{
        if (gRNA2[1,]$strand == "rev" & gRNA2[2,]$strand == "fw") {
          pos1 <- gRNA2[1,]$start
          pos2 <- gRNA2[2,]$end
          KO_length2 <- abs(pos1 - pos2- 1)
        }
        else{
          pos1 <- gRNA2[1,]$end
          pos2 <- gRNA2[2,]$start
          KO_length2 <- abs(pos1 - pos2+ 1) 
        }
      }
    }
    #敲除的CDS
    if (exists("judge_2")) {
      KO_length_CDS2 <- KO_region[h,]$Exon_length
    }
    else if (exists("judge3_2")) {
      KO_length_CDS2 <- KO_region3[h,]$Exon_length
    }
    else if (exists("judge4_2")) {
      KO_length_CDS2 <- KO_region_all$Exon_length
    }
    else{
      KO_length_CDS2 <- KO_length2
    }
    gRNA2[1, 8] <- KO_length2
    gRNA2[2, 8] <- KO_length_CDS2
    return(gRNA2)
  }
}
