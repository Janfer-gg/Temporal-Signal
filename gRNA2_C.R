Get_gRNA2_planC <- function(KO_region2) {
  if (nrow(KO_region2) > t) {
    for (h in (t + 1):nrow(KO_region2)) {
      # 长度不够 --------------------------------------------------------------------
      if(KO_region2[h,]$start<1200 | KO_region2[h,]$end+1200>nchar(Gene)){
        Gene3 <- Get_seq3(ID)
        if (Gene_rev) {
          gene3 <- DNAString(Gene3)
          rev3 <- reverse(gene3)
          Gene3 <- as.character(rev3)
        }
        KO_region2[h,]$start<-KO_region2[h,]$start+1200
        KO_region2[h,]$end<-KO_region2[h,]$end+1200
        ko_start1 <- KO_region2[h,]$start -400
        ko_end1 <- KO_region2[h,]$start
        ko_seq1 <- substring(Gene3, ko_start1, ko_end1)
        if (Gene_rev) {
          seq1 <- DNAString(ko_seq1)
          seq1_rev <- reverse(seq1)
          ko_seq1 <- as.character(seq1_rev)
        }
        source_python("crispor_table_download.py")
        py$run(ko_seq1, species,filepath)
        gRNA.table1 <-
          read.csv(paste0(filepath,"//gRNA.csv"), header = FALSE, encoding = "UTF-8")
        #特异性得分60下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        for (i in 1:length(gRNA.table1[, 1])) {
          if (gRNA.table1[i, 3] < 60) {
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
          rm(Gene3)
          next
        }
        
        #0-0-0(优化)
        count_0 <- numeric()
        for (p in 1:nrow(gRNA.table1)) {
          if (gRNA.table1[p,]$V9 == "0-0-0") {
            count_0 <- append(count_0, p)
          }
        }
        gRNA.table_min <- gRNA.table1[-count_0, ]
        gRNA.table1 <- rbind(gRNA.table1[count_0, ], gRNA.table_min)
        
        if (nrow(gRNA.table1) > 10) {
          gRNA.table1 <- gRNA.table1[1:10,]
        }
        
        #外显子下游
        ko_start2 <- KO_region2[h,]$end 
        ko_end2 <- KO_region2[h,]$end + 400
        ko_seq2 <- substring(Gene3, ko_start2, ko_end2)
        if (Gene_rev) {
          seq2 <- DNAString(ko_seq2)
          seq2_rev <- reverse(seq2)
          ko_seq2 <- as.character(seq2_rev)
        }
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq2, species,filepath)
        gRNA.table2 <-
          read.csv(paste0(filepath,"//gRNA.csv"), header = FALSE, encoding = "UTF-8")
        #特异性得分60下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        for (i in 1:length(gRNA.table2[, 1])) {
          if (gRNA.table2[i, 3] < 60) {
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
          rm(Gene3)
          next
        }
        
        #0-0-0(优化)
        count_0 <- numeric()
        for (p in 1:nrow(gRNA.table2)) {
          if (gRNA.table2[p,]$V9 == "0-0-0") {
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
        GC_avoid_region <- GC_analysis2(KO_region2[h,])
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
        write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor.csv"), row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer(paste0(filepath,"//CCTOP-predictor.csv"),species)
        gRNA.table <-
          read.csv(paste0(filepath,"//CCTOP-predictor.csv"), header = TRUE)
        print(gRNA.table)
        #切割效率得分低于0.60的删除
        gRNA.table <-
          gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
        
        #上下游分开
        gRNA.table1 <-
          gRNA.table[which(gRNA.table$end <= KO_region2[h,]$start),]
        gRNA.table2 <-
          gRNA.table[which(gRNA.table$start >= KO_region2[h,]$end),]
        
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
            gRNA2_planC <- Get_result1(gRNA.table, KO_region2[h,])        #相差0.05分以内优选
          }
        }
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA2_planC) == "NULL") {
          gRNA2_planC <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
        }
      }
      
      # 长度够 ---------------------------------------------------------------------
      else{
        ko_start1 <- KO_region2[h,]$start - 400
        ko_end1 <- KO_region2[h,]$start
        ko_seq1 <- substring(Gene, ko_start1, ko_end1)
        if (Gene_rev) {
          seq1 <- DNAString(ko_seq1)
          seq1_rev <- reverse(seq1)
          ko_seq1 <- as.character(seq1_rev)
        }
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq1, species,filepath)
        gRNA.table1 <-
          read.csv(paste0(filepath,"//gRNA.csv"), header = FALSE, encoding = "UTF-8")
        #特异性得分60下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        for (i in 1:length(gRNA.table1[, 1])) {
          if (gRNA.table1[i, 3] < 60) {
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
            py$run(ko_seq1, species,filepath)
            gRNA.table1 <-
              read.csv(paste0(filepath,"//gRNA.csv"), header = FALSE, encoding = "UTF-8")
            #特异性得分60下，和Inefficient的gRNA排除掉
            gRNA.del <- numeric()
            for (i in 1:length(gRNA.table1[, 1])) {
              if (gRNA.table1[i, 3] < 60) {
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
        gRNA.table_min <- gRNA.table1[-count_0, ]
        gRNA.table1 <- rbind(gRNA.table1[count_0, ], gRNA.table_min)
        
        if (nrow(gRNA.table1) > 10) {
          gRNA.table1 <- gRNA.table1[1:10,]
        }
        
        #外显子下游
        ko_start2 <- KO_region2[h,]$end
        ko_end2 <- KO_region2[h,]$end + 400
        ko_seq2 <- substring(Gene, ko_start2, ko_end2)
        if (Gene_rev) {
          seq2 <- DNAString(ko_seq2)
          seq2_rev <- reverse(seq2)
          ko_seq2 <- as.character(seq2_rev)
        }
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq2, species,filepath)
        gRNA.table2 <-
          read.csv(paste0(filepath,"//gRNA.csv"), header = FALSE, encoding = "UTF-8")
        #特异性得分60下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        for (i in 1:length(gRNA.table2[, 1])) {
          if (gRNA.table2[i, 3] < 60) {
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
            py$run(ko_seq2, species,filepath)
            gRNA.table2 <-
              read.csv(paste0(filepath,"//gRNA.csv"), header = FALSE, encoding = "UTF-8")
            #特异性得分60下，和Inefficient的gRNA排除掉
            gRNA.del <- numeric()
            for (i in 1:length(gRNA.table2[, 1])) {
              if (gRNA.table2[i, 3] < 60) {
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
        gRNA.table_min <- gRNA.table2[-count_0, ]
        gRNA.table2 <- rbind(gRNA.table2[count_0, ], gRNA.table_min)
        
        if (nrow(gRNA.table2) > 7) {
          gRNA.table2 <- gRNA.table2[1:7,]
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
        GC_avoid_region <- GC_analysis2(KO_region2[h,])
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
        write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor.csv"), row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer(paste0(filepath,"//CCTOP-predictor.csv"),species)
        gRNA.table <-
          read.csv(paste0(filepath,"//CCTOP-predictor.csv"), header = TRUE)
        print(gRNA.table)
        #切割效率得分低于0.60的删除
        gRNA.table <-
          gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
        
        #上下游分开
        gRNA.table1 <-
          gRNA.table[which(gRNA.table$end <= KO_region2[h,]$start),]
        gRNA.table2 <-
          gRNA.table[which(gRNA.table$start >= KO_region2[h,]$end),]
        
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
            gRNA2_planC <- Get_result1(gRNA.table, KO_region2[h,])        #相差0.05分以内优选
          }
        }
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA2_planC) == "NULL") {
          gRNA2_planC <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
        }
      }
      if (exists("gRNA2_planC") == TRUE) {
        break
      }
    }
  }
  # 输出 ----------------------------------------------------------------------
  #敲除大小
  if (exists("gRNA2_planC") == TRUE) {
    if (Gene_rev) {
      if (all(gRNA2_planC$strand == "rev") |
          all(gRNA2_planC$strand == "fw")) {
        KO_length <- abs(gRNA2_planC[1,]$end - gRNA2_planC[2,]$end) 
      }
      else{
        if (gRNA2_planC[1, ]$strand == "rev" & gRNA2_planC[2, ]$strand == "fw") {
          if(gRNA2_planC[1,]$start>gRNA2_planC[2,]$start){
            pos1 <- gRNA2_planC[1, ]$end
            pos2 <- gRNA2_planC[2, ]$start
            KO_length <- abs(pos1 - pos2) +1
          }
          else{
            pos1 <- gRNA2_planC[1, ]$end
            pos2 <- gRNA2_planC[2, ]$start
            KO_length <- abs(pos1 - pos2) -1
          }
        }
        else{
          if(gRNA2_planC[1,]$start>gRNA2_planC[2,]$start){
            pos1 <- gRNA2_planC[1, ]$start
            pos2 <- gRNA2_planC[2, ]$end
            KO_length <- abs(pos1 - pos2) - 1
          }
          else{
            pos1 <- gRNA2_planC[1, ]$start
            pos2 <- gRNA2_planC[2, ]$end
            KO_length <- abs(pos1 - pos2) + 1
          }
        }
      }
    }
    #正向
    else{
      if (all(gRNA2_planC$strand == "rev") |
          all(gRNA2_planC$strand == "fw")) {
        KO_length <- abs(gRNA2_planC[1,]$end - gRNA2_planC[2,]$end) 
      }
      else{
        if (gRNA2_planC[1, ]$strand == "rev" & gRNA2_planC[2, ]$strand == "fw") {
          pos1 <- gRNA2_planC[1, ]$start
          pos2 <- gRNA2_planC[2, ]$end
          KO_length <- abs(pos1 - pos2) -1
        }
        else{
          pos1 <- gRNA2_planC[1, ]$end
          pos2 <- gRNA2_planC[2, ]$start
          KO_length <- abs(pos1 - pos2) -1
        }
      }
    }
    #敲除的CDS
    KO_length_CDS <- KO_region2[h,]$Exon_length
    gRNA2_planC[1,8]<-KO_length
    gRNA2_planC[1,9]<-KO_length_CDS
    return(gRNA2_planC)
  }
}
