# 非移码方案的内侧靶位点 -------------------------------------------------------------------
# 整个敲除 --------------------------------------------------------------------
{
  if (exists("judge2")) {
    n <-which(t_Exon_CDS$start >= KO_region2[t, ]$start &
              t_Exon_CDS$end <= KO_region2[t, ]$end)
    PCR_seq1 <-
      substring(Gene, t_Exon_CDS[n[1],]$start, t_Exon_CDS[n[1],]$end)
    PCR_seq2 <-
      substring(Gene, t_Exon_CDS[n[length(n)],]$start, t_Exon_CDS[n[length(n)],]$end)
    if (Gene_rev) {
      seq1 <- DNAString(PCR_seq1)
      seq1_rev <- reverse(seq1)
      PCR_seq1 <- as.character(seq1_rev)
      seq2 <- DNAString(PCR_seq2)
      seq2_rev <- reverse(seq2)
      PCR_seq2 <- as.character(seq2_rev)
    }
    
    #上游的PCR引物
    py$run2(PCR_seq1,crispor.org)
    PCR.tab1 <- read.csv("PCR.csv", header = FALSE)
    
    #下游的PCR引物
    py$run2(PCR_seq2,crispor.org)
    PCR.tab2 <- read.csv("PCR.csv", header = FALSE)
    
    {
      # 反向基因 --------------------------------------------------------------------
      if (Gene_rev) {
        PCR.tab1 <-
          PCR.tab1[which(sub("[^a-zA-Z]+", "", PCR.tab1[, 1]) == "rev"), ]
        PCR.tab2 <-
          PCR.tab2[which(sub("[^a-zA-Z]+", "", PCR.tab2[, 1]) == "fw"), ]
        
        # 上游 ----------------------------------------------------------------------
        #删除连续4个碱基重复的序列
        del <- del_dup(PCR.tab1)
        if (length(del) != 0) {
          PCR.table1 <- PCR.tab1[-del,]
        }
        else{
          PCR.table1 <- PCR.tab1
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        pcr_seq1 <- DNAString(PCR_seq1)
        PCR.table_fw <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table1)) {
          PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
          PCR_fw <- reverseComplement(PCR_fw)
          fw <-
            matchPattern(pattern = PCR_fw, subject = pcr_seq1)
          if (length(fw) > 1) {
            next
          }
          fw_start <- start(fw)
          fw_end <- end(fw)
          PCR_pre <- substring(PCR_seq1, fw_start - 3, fw_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <-
                    100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_fw[m, 1] <- PCR_pre2
                    PCR.table_fw[m, 2] <- GC
                    PCR.table_fw[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        del <- numeric()
        del <- del_dup2(PCR.table_fw)
        if (length(del) != 0) {
          PCR.table_fw <- PCR.table_fw[-del,]
        }
        #去除重复
        PCR.table_fw <-
          PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
        #取20条序列
        if (nrow(PCR.table_fw) > 20) {
          PCR.table_fw <- PCR.table_fw[1:20,]
        }
        write.csv(PCR.table_fw,
                  file = "PCR1.csv",
                  row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        # 下游 ----------------------------------------------------------------------
        #去掉连续4个碱基相同的序列
        del <- del_dup(PCR.tab2)
        if (length(del) != 0) {
          PCR.table2 <- PCR.tab2[-del,]
        }
        else{
          PCR.table2 <- PCR.tab2
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        pcr_seq2 <- DNAString(PCR_seq2)
        PCR.table_rev <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table2)) {
          rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                              subject = pcr_seq2)
          if (length(rev) > 1) {
            next
          }
          rev_start <- start(rev)
          rev_end <- end(rev)
          PCR_pre <-
            substring(PCR_seq2, rev_start - 3, rev_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <-
                    100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_rev[m, 1] <- PCR_pre2
                    PCR.table_rev[m, 2] <- GC
                    PCR.table_rev[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        
        del <- numeric()
        del <- del_dup2(PCR.table_rev)
        if (length(del) != 0) {
          PCR.table_rev <- PCR.table_rev[-del,]
        }
        #去除重复
        PCR.table_rev <-
          PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
        #取20条
        if (nrow(PCR.table_rev) > 20) {
          PCR.table_rev <- PCR.table_rev[1:20,]
        }
        write.csv(PCR.table_rev,
                  file = "PCR2.csv",
                  row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
        #PCR1的反向互补序列
        for (i in 1:nrow(PCR2.table)) {
          pcr2 <- DNAString(PCR2.table[i,]$PCR)
          pcr2_revcom <- reverseComplement(pcr2)
          PCR2.table[i,]$PCR <- as.character(pcr2_revcom)
        }
        write.csv(PCR2.table,
                  file = "PCR3.csv",
                  row.names = FALSE)
        
        # Blast -------------------------------------------------------------------
        source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
        )
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
        )
        primer.table <-
          read.csv("primer_result1.csv", header = TRUE)
        print(primer.table)
        primer.table <-
          read.csv("primer_result2.csv", header = TRUE)
        print(primer.table)
      }
      
      # 正向基因 --------------------------------------------------------------------
      else{
        PCR.tab1 <-
          PCR.tab1[which(sub("[^a-zA-Z]+", "", PCR.tab1[, 1]) == "rev"), ]
        PCR.tab2 <-
          PCR.tab2[which(sub("[^a-zA-Z]+", "", PCR.tab2[, 1]) == "fw"), ]
        
        # 上游 ----------------------------------------------------------------------
        del <- del_dup(PCR.tab1)
        if (length(del) != 0) {
          PCR.table1 <- PCR.tab1[-del,]
        }
        else{
          PCR.table1 <- PCR.tab1
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        pcr_seq1 <- DNAString(PCR_seq1)
        PCR.table_fw <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table1)) {
          PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
          PCR_fw <- reverseComplement(PCR_fw)
          fw <-
            matchPattern(pattern = PCR_fw, subject = pcr_seq1)
          if (length(fw) > 1) {
            next
          }
          fw_start <- start(fw)
          fw_end <- end(fw)
          PCR_pre <- substring(PCR_seq1, fw_start - 3, fw_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <-
                    100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_fw[m, 1] <- PCR_pre2
                    PCR.table_fw[m, 2] <- GC
                    PCR.table_fw[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        del <- numeric()
        del <- del_dup2(PCR.table_fw)
        if (length(del) != 0) {
          PCR.table_fw <- PCR.table_fw[-del,]
        }
        #去除重复
        PCR.table_fw <-
          PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
        #取20条序列
        if (nrow(PCR.table_fw) > 20) {
          PCR.table_fw <- PCR.table_fw[1:20,]
        }
        write.csv(PCR.table_fw,
                  file = "PCR1.csv",
                  row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        # 下游 ----------------------------------------------------------------------
        #去掉连续4个碱基相同的序列
        del <- del_dup(PCR.tab2)
        if (length(del) != 0) {
          PCR.table2 <- PCR.tab2[-del,]
        }
        else{
          PCR.table2 <- PCR.tab2
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        pcr_seq2 <- DNAString(PCR_seq2)
        PCR.table_rev <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table2)) {
          rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                              subject = pcr_seq2)
          if (length(rev) > 1) {
            next
          }
          rev_start <- start(rev)
          rev_end <- end(rev)
          PCR_pre <-
            substring(PCR_seq2, rev_start - 3, rev_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <-
                    100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_rev[m, 1] <- PCR_pre2
                    PCR.table_rev[m, 2] <- GC
                    PCR.table_rev[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        
        del <- numeric()
        del <- del_dup2(PCR.table_rev)
        if (length(del) != 0) {
          PCR.table_rev <- PCR.table_rev[-del,]
        }
        #去除重复
        PCR.table_rev <-
          PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
        #取20条
        if (nrow(PCR.table_rev) > 20) {
          PCR.table_rev <- PCR.table_rev[1:20,]
        }
        write.csv(PCR.table_rev,
                  file = "PCR2.csv",
                  row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
        #PCR2的反向互补序列
        for (i in 1:nrow(PCR2.table)) {
          pcr2 <- DNAString(PCR2.table[i,]$PCR)
          pcr2_revcom <- reverseComplement(pcr2)
          PCR2.table[i,]$PCR <- as.character(pcr2_revcom)
        }
        write.csv(PCR2.table,
                  file = "PCR3.csv",
                  row.names = FALSE)
        
        
        # BLAST -------------------------------------------------------------------
        source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
        )
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
        )
        primer.table <-
          read.csv("primer_result1.csv", header = TRUE)
        print(primer.table)
        primer.table <-
          read.csv("primer_result2.csv", header = TRUE)
        print(primer.table)
      }
    }
    
  }
  
  # 内含子小于500 ----------------------------------------------------------------
  if (exists("judge3")) {
    PCR_seq <- substring(Gene, KO_region3[a, ]$start, KO_region3[a, ]$end)
    if (Gene_rev) {
      seq1 <- DNAString(PCR_seq)
      seq1_rev <- reverse(seq1)
      PCR_seq <- as.character(seq1_rev)
    }
    py$run2(PCR_seq,crispor.org)
    PCR.tab <- read.csv("PCR.csv", header = FALSE)
    {
      # 反向基因 --------------------------------------------------------------------
      if (Gene_rev) {
        PCR.tab1 <-
          PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "rev"),]
        PCR.tab2 <-
          PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "fw"),]
        
        # 上游 ----------------------------------------------------------------------
        #删除连续4个碱基重复的序列
        del <- del_dup(PCR.tab1)
        if (length(del) != 0) {
          PCR.table1 <- PCR.tab1[-del, ]
        }
        else{
          PCR.table1 <- PCR.tab1
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        pcr_seq <- DNAString(PCR_seq)
        PCR.table_fw <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table1)) {
          PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
          PCR_fw <- reverseComplement(PCR_fw)
          fw <- matchPattern(pattern = PCR_fw, subject = pcr_seq)
          if (length(fw) > 1) {
            next
          }
          fw_start <- start(fw)
          fw_end <- end(fw)
          PCR_pre <- substring(PCR_seq, fw_start - 3, fw_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_fw[m, 1] <- PCR_pre2
                    PCR.table_fw[m, 2] <- GC
                    PCR.table_fw[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        del <- numeric()
        del <- del_dup2(PCR.table_fw)
        if (length(del) != 0) {
          PCR.table_fw <- PCR.table_fw[-del, ]
        }
        #去除重复
        PCR.table_fw <- PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
        #取20条序列
        if (nrow(PCR.table_fw) > 20) {
          PCR.table_fw <- PCR.table_fw[1:20, ]
        }
        write.csv(PCR.table_fw, file = "PCR1.csv", row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        # 下游 ----------------------------------------------------------------------
        #去掉连续4个碱基相同的序列
        del <- del_dup(PCR.tab2)
        if (length(del) != 0) {
          PCR.table2 <- PCR.tab2[-del, ]
        }
        else{
          PCR.table2 <- PCR.tab2
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        PCR.table_rev <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table2)) {
          rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                              subject = pcr_seq)
          if (length(rev) > 1) {
            next
          }
          rev_start <- start(rev)
          rev_end <- end(rev)
          PCR_pre <- substring(PCR_seq, rev_start - 3, rev_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_rev[m, 1] <- PCR_pre2
                    PCR.table_rev[m, 2] <- GC
                    PCR.table_rev[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        
        del <- numeric()
        del <- del_dup2(PCR.table_rev)
        if (length(del) != 0) {
          PCR.table_rev <- PCR.table_rev[-del, ]
        }
        #去除重复
        PCR.table_rev <-
          PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
        #取20条
        if (nrow(PCR.table_rev) > 20) {
          PCR.table_rev <- PCR.table_rev[1:20, ]
        }
        write.csv(PCR.table_rev, file = "PCR2.csv", row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
        #PCR1的反向互补序列
        for (i in 1:nrow(PCR2.table)) {
          pcr2 <- DNAString(PCR2.table[i, ]$PCR)
          pcr2_revcom <- reverseComplement(pcr2)
          PCR2.table[i, ]$PCR <- as.character(pcr2_revcom)
        }
        write.csv(PCR2.table, file = "PCR3.csv", row.names = FALSE)
        
        # Blast -------------------------------------------------------------------
        source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
        )
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
        )
        
        primer.table <- read.csv("primer_result1.csv", header = TRUE)
        print(primer.table)
        primer.table <- read.csv("primer_result2.csv", header = TRUE)
        print(primer.table)
        
      }
      
      # 正向基因 --------------------------------------------------------------------
      else{
        PCR.tab1 <-
          PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "fw"),]
        PCR.tab2 <-
          PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "rev"),]
        
        # 上游 ----------------------------------------------------------------------
        del <- del_dup(PCR.tab1)
        if (length(del) != 0) {
          PCR.table1 <- PCR.tab1[-del, ]
        }
        else{
          PCR.table1 <- PCR.tab1
        }
        #匹配PCR的位置，然后前后扩大3个碱基
        pcr_seq <- DNAString(PCR_seq)
        PCR.table_fw <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table1)) {
          fw <- matchPattern(pattern = substring(PCR.table1[i, 2], 1, 20),
                             subject = pcr_seq)
          if (length(fw) > 1) {
            next
          }
          fw_start <- start(fw)
          fw_end <- end(fw)
          PCR_pre <- substring(PCR_seq, fw_start - 3, fw_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_fw[m, 1] <- PCR_pre2
                    PCR.table_fw[m, 2] <- GC
                    PCR.table_fw[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        del <- numeric()
        del <- del_dup2(PCR.table_fw)
        if (length(del) != 0) {
          PCR.table_fw <- PCR.table_fw[-del, ]
        }
        #去除重复
        PCR.table_fw <- PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
        #取20条
        if (nrow(PCR.table_fw) > 20) {
          PCR.table_fw <- PCR.table_fw[1:20, ]
        }
        write.csv(PCR.table_fw, file = "PCR1.csv", row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        
        
        
        # 下游 ----------------------------------------------------------------------
        #去掉连续4个碱基相同的序列
        del <- del_dup(PCR.tab2)
        if (length(del) != 0) {
          PCR.table2 <- PCR.tab2[-del, ]
        }
        else{
          PCR.table2 <- PCR.tab2
        }
        PCR.table2[, 1] <-
          as.numeric(sub("[^0-9]+", "", PCR.table2[, 1]))
        #匹配PCR的位置，然后前后扩大3个碱基
        PCR.table_rev <-
          data.frame(PCR = character(),
                     GC = numeric(),
                     TM = numeric())
        m <- 1
        for (i in 1:nrow(PCR.table2)) {
          PCR_rev <- DNAString(substring(PCR.table2[i, 2], 1, 20))
          PCR_rev <- reverseComplement(PCR_rev)
          rev <- matchPattern(pattern = PCR_rev, subject = pcr_seq)
          if (length(rev) > 1) {
            next
          }
          rev_start <- start(rev)
          rev_end <- end(rev)
          PCR_pre <- substring(PCR_seq, rev_start - 3, rev_end + 3)
          for (j in 20:25) {
            for (k in 1:(nchar(PCR_pre) - (j - 1))) {
              PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
              if (j > nchar(PCR_pre2)) {
                next
              }
              #如果以G或C结尾，则计算GC值和TM值
              if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                  strsplit(PCR_pre2, "")[[1]][j] == "C") {
                cal <- DNAString(PCR_pre2)
                GC_count <- sum(letterFrequency(cal, c("G", "C")))
                GC <- round(GC_count / j, 2)
                #如果GC大于0.4和小于0.6，计算TM值
                if (GC >= 0.4 & GC <= 0.6) {
                  TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                  #如果TM值也符合，则该PCR引物符合
                  if (TM >= 61 & TM <= 64) {
                    PCR.table_rev[m, 1] <- PCR_pre2
                    PCR.table_rev[m, 2] <- GC
                    PCR.table_rev[m, 3] <- TM
                    m <- m + 1
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
              else{
                next
              }
            }
          }
        }
        
        del <- numeric()
        del <- del_dup2(PCR.table_rev)
        if (length(del) != 0) {
          PCR.table_rev <- PCR.table_rev[-del, ]
        }
        #去除重复
        PCR.table_rev <-
          PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
        #取20条
        if (nrow(PCR.table_rev) > 20) {
          PCR.table_rev <- PCR.table_rev[1:20, ]
        }
        write.csv(PCR.table_rev, file = "PCR2.csv", row.names = FALSE)
        py$RunOligoCalc(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
        )
        PCR2.table <- read.csv(file = "PCR2.csv", header = TRUE)
        #PCR2的反向互补序列
        for (i in 1:nrow(PCR2.table)) {
          pcr2 <- DNAString(PCR2.table[i, ]$PCR)
          pcr2_revcom <- reverseComplement(pcr2)
          PCR2.table[i, ]$PCR <- as.character(pcr2_revcom)
        }
        write.csv(PCR2.table, file = "PCR3.csv", row.names = FALSE)
        
        
        # BLAST -------------------------------------------------------------------
        source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
        )
        py$run_primertool(
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
          "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
        )
        primer.table <- read.csv("primer_result1.csv", header = TRUE)
        print(primer.table)
        primer.table <- read.csv("primer_result2.csv", header = TRUE)
        print(primer.table)
      }
    }
    
  }
  
  
  # 内含子，外显子 ---------------------------------------------------------------------
  else{
    if (exists("judge")) {
      n <-which(t_Exon_CDS$start >= KO_region[t, ]$start &
                t_Exon_CDS$end <= KO_region[t, ]$end)
      # 单个外显子 -------------------------------------------------------------------
      if (length(n) == 1) {
        PCR_seq <- substring(Gene, t_Exon_CDS[n, ]$start, t_Exon_CDS[n, ]$end)
        if (Gene_rev) {
          seq1 <- DNAString(PCR_seq)
          seq1_rev <- reverse(seq1)
          PCR_seq <- as.character(seq1_rev)
        }
        py$run2(PCR_seq,crispor.org)
        PCR.tab <- read.csv("PCR.csv", header = FALSE)
        {
          # 反向基因 --------------------------------------------------------------------
          if (Gene_rev) {
            PCR.tab1 <-
              PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "rev"),]
            PCR.tab2 <-
              PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "fw"),]
            
            # 上游 ----------------------------------------------------------------------
            #删除连续4个碱基重复的序列
            del <- del_dup(PCR.tab1)
            if (length(del) != 0) {
              PCR.table1 <- PCR.tab1[-del, ]
            }
            else{
              PCR.table1 <- PCR.tab1
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            pcr_seq <- DNAString(PCR_seq)
            PCR.table_fw <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table1)) {
              PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
              PCR_fw <- reverseComplement(PCR_fw)
              fw <- matchPattern(pattern = PCR_fw, subject = pcr_seq)
              if (length(fw) > 1) {
                next
              }
              fw_start <- start(fw)
              fw_end <- end(fw)
              PCR_pre <- substring(PCR_seq, fw_start - 3, fw_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_fw[m, 1] <- PCR_pre2
                        PCR.table_fw[m, 2] <- GC
                        PCR.table_fw[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            del <- numeric()
            del <- del_dup2(PCR.table_fw)
            if (length(del) != 0) {
              PCR.table_fw <- PCR.table_fw[-del, ]
            }
            #去除重复
            PCR.table_fw <-
              PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
            #取20条序列
            if (nrow(PCR.table_fw) > 20) {
              PCR.table_fw <- PCR.table_fw[1:20, ]
            }
            write.csv(PCR.table_fw,
                      file = "PCR1.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            # 下游 ----------------------------------------------------------------------
            #去掉连续4个碱基相同的序列
            del <- del_dup(PCR.tab2)
            if (length(del) != 0) {
              PCR.table2 <- PCR.tab2[-del, ]
            }
            else{
              PCR.table2 <- PCR.tab2
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            PCR.table_rev <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table2)) {
              rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                                  subject = pcr_seq)
              if (length(rev) > 1) {
                next
              }
              rev_start <- start(rev)
              rev_end <- end(rev)
              PCR_pre <- substring(PCR_seq, rev_start - 3, rev_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_rev[m, 1] <- PCR_pre2
                        PCR.table_rev[m, 2] <- GC
                        PCR.table_rev[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            
            del <- numeric()
            del <- del_dup2(PCR.table_rev)
            if (length(del) != 0) {
              PCR.table_rev <- PCR.table_rev[-del, ]
            }
            #去除重复
            PCR.table_rev <-
              PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
            #取20条
            if (nrow(PCR.table_rev) > 20) {
              PCR.table_rev <- PCR.table_rev[1:20, ]
            }
            write.csv(PCR.table_rev,
                      file = "PCR2.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
            #PCR1的反向互补序列
            for (i in 1:nrow(PCR2.table)) {
              pcr2 <- DNAString(PCR2.table[i, ]$PCR)
              pcr2_revcom <- reverseComplement(pcr2)
              PCR2.table[i, ]$PCR <- as.character(pcr2_revcom)
            }
            write.csv(PCR2.table,
                      file = "PCR3.csv",
                      row.names = FALSE)
            
            # Blast -------------------------------------------------------------------
            source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
            )
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
            )
            
            primer.table <-
              read.csv("primer_result1.csv", header = TRUE)
            print(primer.table)
            primer.table <-
              read.csv("primer_result2.csv", header = TRUE)
            print(primer.table)
            
          }
          
          # 正向基因 --------------------------------------------------------------------
          else{
            PCR.tab1 <-
              PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "fw"),]
            PCR.tab2 <-
              PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "rev"),]
            
            # 上游 ----------------------------------------------------------------------
            del <- del_dup(PCR.tab1)
            if (length(del) != 0) {
              PCR.table1 <- PCR.tab1[-del, ]
            }
            else{
              PCR.table1 <- PCR.tab1
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            pcr_seq <- DNAString(PCR_seq)
            PCR.table_fw <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table1)) {
              fw <- matchPattern(pattern = substring(PCR.table1[i, 2], 1, 20),
                                 subject = pcr_seq)
              if (length(fw) > 1) {
                next
              }
              fw_start <- start(fw)
              fw_end <- end(fw)
              PCR_pre <- substring(PCR_seq, fw_start - 3, fw_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_fw[m, 1] <- PCR_pre2
                        PCR.table_fw[m, 2] <- GC
                        PCR.table_fw[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            del <- numeric()
            del <- del_dup2(PCR.table_fw)
            if (length(del) != 0) {
              PCR.table_fw <- PCR.table_fw[-del, ]
            }
            #去除重复
            PCR.table_fw <-
              PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
            #取20条
            if (nrow(PCR.table_fw) > 20) {
              PCR.table_fw <- PCR.table_fw[1:20, ]
            }
            write.csv(PCR.table_fw,
                      file = "PCR1.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            
            
            # 下游 ----------------------------------------------------------------------
            #去掉连续4个碱基相同的序列
            del <- del_dup(PCR.tab2)
            if (length(del) != 0) {
              PCR.table2 <- PCR.tab2[-del, ]
            }
            else{
              PCR.table2 <- PCR.tab2
            }
            PCR.table2[, 1] <-
              as.numeric(sub("[^0-9]+", "", PCR.table2[, 1]))
            #匹配PCR的位置，然后前后扩大3个碱基
            PCR.table_rev <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table2)) {
              PCR_rev <- DNAString(substring(PCR.table2[i, 2], 1, 20))
              PCR_rev <- reverseComplement(PCR_rev)
              rev <-
                matchPattern(pattern = PCR_rev, subject = pcr_seq)
              if (length(rev) > 1) {
                next
              }
              rev_start <- start(rev)
              rev_end <- end(rev)
              PCR_pre <- substring(PCR_seq, rev_start - 3, rev_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_rev[m, 1] <- PCR_pre2
                        PCR.table_rev[m, 2] <- GC
                        PCR.table_rev[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            
            del <- numeric()
            del <- del_dup2(PCR.table_rev)
            if (length(del) != 0) {
              PCR.table_rev <- PCR.table_rev[-del, ]
            }
            #去除重复
            PCR.table_rev <-
              PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
            #取20条
            if (nrow(PCR.table_rev) > 20) {
              PCR.table_rev <- PCR.table_rev[1:20, ]
            }
            write.csv(PCR.table_rev,
                      file = "PCR2.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            PCR2.table <- read.csv(file = "PCR2.csv", header = TRUE)
            #PCR2的反向互补序列
            for (i in 1:nrow(PCR2.table)) {
              pcr2 <- DNAString(PCR2.table[i, ]$PCR)
              pcr2_revcom <- reverseComplement(pcr2)
              PCR2.table[i, ]$PCR <- as.character(pcr2_revcom)
            }
            write.csv(PCR2.table,
                      file = "PCR3.csv",
                      row.names = FALSE)
            
            
            # BLAST -------------------------------------------------------------------
            source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
            )
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
            )
            primer.table <-
              read.csv("primer_result1.csv", header = TRUE)
            print(primer.table)
            primer.table <-
              read.csv("primer_result2.csv", header = TRUE)
            print(primer.table)
          }
        }
      }
      
      # 多个外显子 -------------------------------------------------------------------
      else{
        PCR_seq1 <-
          substring(Gene, t_Exon_CDS[n[1], ]$start, t_Exon_CDS[n[1], ]$end)
        PCR_seq2 <-
          PCR_seq <-
          substring(Gene, t_Exon_CDS[n[length(n)], ]$start, t_Exon_CDS[n[length(n)], ]$end)
        if (Gene_rev) {
          seq1 <- DNAString(PCR_seq1)
          seq1_rev <- reverse(seq1)
          PCR_seq1 <- as.character(seq1_rev)
          seq2 <- DNAString(PCR_seq2)
          seq2_rev <- reverse(seq2)
          PCR_seq2 <- as.character(seq2_rev)
        }
        
        #上游的PCR引物
        py$run2(PCR_seq1,crispor.org)
        PCR.tab1 <- read.csv("PCR.csv", header = FALSE)
        
        #下游的PCR引物
        py$run2(PCR_seq2,crispor.org)
        PCR.tab2 <- read.csv("PCR.csv", header = FALSE)
        
        {
          # 反向基因 --------------------------------------------------------------------
          if (Gene_rev) {
            PCR.tab1 <-
              PCR.tab1[which(sub("[^a-zA-Z]+", "", PCR.tab1[, 1]) == "rev"),]
            PCR.tab2 <-
              PCR.tab2[which(sub("[^a-zA-Z]+", "", PCR.tab2[, 1]) == "fw"),]
            
            # 上游 ----------------------------------------------------------------------
            #删除连续4个碱基重复的序列
            del <- del_dup(PCR.tab1)
            if (length(del) != 0) {
              PCR.table1 <- PCR.tab1[-del, ]
            }
            else{
              PCR.table1 <- PCR.tab1
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            pcr_seq1 <- DNAString(PCR_seq1)
            PCR.table_fw <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table1)) {
              PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
              PCR_fw <- reverseComplement(PCR_fw)
              fw <-
                matchPattern(pattern = PCR_fw, subject = pcr_seq1)
              if (length(fw) > 1) {
                next
              }
              fw_start <- start(fw)
              fw_end <- end(fw)
              PCR_pre <- substring(PCR_seq1, fw_start - 3, fw_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_fw[m, 1] <- PCR_pre2
                        PCR.table_fw[m, 2] <- GC
                        PCR.table_fw[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            del <- numeric()
            del <- del_dup2(PCR.table_fw)
            if (length(del) != 0) {
              PCR.table_fw <- PCR.table_fw[-del, ]
            }
            #去除重复
            PCR.table_fw <-
              PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
            #取20条序列
            if (nrow(PCR.table_fw) > 20) {
              PCR.table_fw <- PCR.table_fw[1:20, ]
            }
            write.csv(PCR.table_fw,
                      file = "PCR1.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            # 下游 ----------------------------------------------------------------------
            #去掉连续4个碱基相同的序列
            del <- del_dup(PCR.tab2)
            if (length(del) != 0) {
              PCR.table2 <- PCR.tab2[-del, ]
            }
            else{
              PCR.table2 <- PCR.tab2
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            pcr_seq2 <- DNAString(PCR_seq2)
            PCR.table_rev <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table2)) {
              rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                                  subject = pcr_seq2)
              if (length(rev) > 1) {
                next
              }
              rev_start <- start(rev)
              rev_end <- end(rev)
              PCR_pre <- substring(PCR_seq2, rev_start - 3, rev_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_rev[m, 1] <- PCR_pre2
                        PCR.table_rev[m, 2] <- GC
                        PCR.table_rev[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            
            del <- numeric()
            del <- del_dup2(PCR.table_rev)
            if (length(del) != 0) {
              PCR.table_rev <- PCR.table_rev[-del, ]
            }
            #去除重复
            PCR.table_rev <-
              PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
            #取20条
            if (nrow(PCR.table_rev) > 20) {
              PCR.table_rev <- PCR.table_rev[1:20, ]
            }
            write.csv(PCR.table_rev,
                      file = "PCR2.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
            #PCR1的反向互补序列
            for (i in 1:nrow(PCR2.table)) {
              pcr2 <- DNAString(PCR2.table[i, ]$PCR)
              pcr2_revcom <- reverseComplement(pcr2)
              PCR2.table[i, ]$PCR <- as.character(pcr2_revcom)
            }
            write.csv(PCR2.table,
                      file = "PCR3.csv",
                      row.names = FALSE)
            
            # Blast -------------------------------------------------------------------
            source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
            )
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
            )
            primer.table <-
              read.csv("primer_result1.csv", header = TRUE)
            print(primer.table)
            primer.table <-
              read.csv("primer_result2.csv", header = TRUE)
            print(primer.table)
          }
          
          # 正向基因 --------------------------------------------------------------------
          else{
            PCR.tab1 <-
              PCR.tab1[which(sub("[^a-zA-Z]+", "", PCR.tab1[, 1]) == "rev"),]
            PCR.tab2 <-
              PCR.tab2[which(sub("[^a-zA-Z]+", "", PCR.tab2[, 1]) == "fw"),]
            
            # 上游 ----------------------------------------------------------------------
            del <- del_dup(PCR.tab1)
            if (length(del) != 0) {
              PCR.table1 <- PCR.tab1[-del, ]
            }
            else{
              PCR.table1 <- PCR.tab1
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            pcr_seq1 <- DNAString(PCR_seq1)
            PCR.table_fw <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table1)) {
              PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
              PCR_fw <- reverseComplement(PCR_fw)
              fw <-
                matchPattern(pattern = PCR_fw, subject = pcr_seq1)
              if (length(fw) > 1) {
                next
              }
              fw_start <- start(fw)
              fw_end <- end(fw)
              PCR_pre <- substring(PCR_seq1, fw_start - 3, fw_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_fw[m, 1] <- PCR_pre2
                        PCR.table_fw[m, 2] <- GC
                        PCR.table_fw[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            del <- numeric()
            del <- del_dup2(PCR.table_fw)
            if (length(del) != 0) {
              PCR.table_fw <- PCR.table_fw[-del, ]
            }
            #去除重复
            PCR.table_fw <-
              PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
            #取20条序列
            if (nrow(PCR.table_fw) > 20) {
              PCR.table_fw <- PCR.table_fw[1:20, ]
            }
            write.csv(PCR.table_fw,
                      file = "PCR1.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            # 下游 ----------------------------------------------------------------------
            #去掉连续4个碱基相同的序列
            del <- del_dup(PCR.tab2)
            if (length(del) != 0) {
              PCR.table2 <- PCR.tab2[-del, ]
            }
            else{
              PCR.table2 <- PCR.tab2
            }
            #匹配PCR的位置，然后前后扩大3个碱基
            pcr_seq2 <- DNAString(PCR_seq2)
            PCR.table_rev <-
              data.frame(PCR = character(),
                         GC = numeric(),
                         TM = numeric())
            m <- 1
            for (i in 1:nrow(PCR.table2)) {
              rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                                  subject = pcr_seq2)
              if (length(rev) > 1) {
                next
              }
              rev_start <- start(rev)
              rev_end <- end(rev)
              PCR_pre <- substring(PCR_seq2, rev_start - 3, rev_end + 3)
              for (j in 20:25) {
                for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                  PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                  if (j > nchar(PCR_pre2)) {
                    next
                  }
                  #如果以G或C结尾，则计算GC值和TM值
                  if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                      strsplit(PCR_pre2, "")[[1]][j] == "C") {
                    cal <- DNAString(PCR_pre2)
                    GC_count <- sum(letterFrequency(cal, c("G", "C")))
                    GC <- round(GC_count / j, 2)
                    #如果GC大于0.4和小于0.6，计算TM值
                    if (GC >= 0.4 & GC <= 0.6) {
                      TM <- 100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                      #如果TM值也符合，则该PCR引物符合
                      if (TM >= 61 & TM <= 64) {
                        PCR.table_rev[m, 1] <- PCR_pre2
                        PCR.table_rev[m, 2] <- GC
                        PCR.table_rev[m, 3] <- TM
                        m <- m + 1
                      }
                      else{
                        next
                      }
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
              }
            }
            
            del <- numeric()
            del <- del_dup2(PCR.table_rev)
            if (length(del) != 0) {
              PCR.table_rev <- PCR.table_rev[-del, ]
            }
            #去除重复
            PCR.table_rev <-
              PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
            #取20条
            if (nrow(PCR.table_rev) > 20) {
              PCR.table_rev <- PCR.table_rev[1:20, ]
            }
            write.csv(PCR.table_rev,
                      file = "PCR2.csv",
                      row.names = FALSE)
            py$RunOligoCalc(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
            )
            
            PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
            #PCR2的反向互补序列
            for (i in 1:nrow(PCR2.table)) {
              pcr2 <- DNAString(PCR2.table[i, ]$PCR)
              pcr2_revcom <- reverseComplement(pcr2)
              PCR2.table[i, ]$PCR <- as.character(pcr2_revcom)
            }
            write.csv(PCR2.table,
                      file = "PCR3.csv",
                      row.names = FALSE)
            
            
            # BLAST -------------------------------------------------------------------
            source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
            )
            py$run_primertool(
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
              "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
            )
            # primer.table <-
            #   read.csv("primer_result1.csv", header = TRUE)
            # print(primer.table)
            # primer.table <-
            #   read.csv("primer_result2.csv", header = TRUE)
            # print(primer.table)
          }
        }
      }
    }
    
    # 外显子上 ---------------------------------------------------------------------
    else{
      PCR_seq <- substring(Gene,min(gRNA$end)+1,max(gRNA$start)-1)
      if (Gene_rev) {
        seq1 <- DNAString(PCR_seq)
        seq1_rev <- reverse(seq1)
        PCR_seq <- as.character(seq1_rev)
      }
      py$run2(PCR_seq,crispor.org)
      PCR.tab <- read.csv("PCR.csv", header = FALSE)
      {
        # 反向基因 --------------------------------------------------------------------
        if (Gene_rev) {
          PCR.tab1 <-
            PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "rev"), ]
          PCR.tab2 <-
            PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "fw"), ]
          
          # 上游 ----------------------------------------------------------------------
          #删除连续4个碱基重复的序列
          del <- del_dup(PCR.tab1)
          if (length(del) != 0) {
            PCR.table1 <- PCR.tab1[-del,]
          }
          else{
            PCR.table1 <- PCR.tab1
          }
          #匹配PCR的位置，然后前后扩大3个碱基
          pcr_seq <- DNAString(PCR_seq)
          PCR.table_fw <-
            data.frame(PCR = character(),
                       GC = numeric(),
                       TM = numeric())
          m <- 1
          for (i in 1:nrow(PCR.table1)) {
            PCR_fw <- DNAString(substring(PCR.table1[i, 2], 1, 20))
            PCR_fw <- reverseComplement(PCR_fw)
            fw <- matchPattern(pattern = PCR_fw, subject = pcr_seq)
            if (length(fw) > 1) {
              next
            }
            fw_start <- start(fw)
            fw_end <- end(fw)
            PCR_pre <- substring(PCR_seq, fw_start - 3, fw_end + 3)
            for (j in 20:25) {
              for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                if (j > nchar(PCR_pre2)) {
                  next
                }
                #如果以G或C结尾，则计算GC值和TM值
                if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                    strsplit(PCR_pre2, "")[[1]][j] == "C") {
                  cal <- DNAString(PCR_pre2)
                  GC_count <- sum(letterFrequency(cal, c("G", "C")))
                  GC <- round(GC_count / j, 2)
                  #如果GC大于0.4和小于0.6，计算TM值
                  if (GC >= 0.4 & GC <= 0.6) {
                    TM <-
                      100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                    #如果TM值也符合，则该PCR引物符合
                    if (TM >= 61 & TM <= 64) {
                      PCR.table_fw[m, 1] <- PCR_pre2
                      PCR.table_fw[m, 2] <- GC
                      PCR.table_fw[m, 3] <- TM
                      m <- m + 1
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
            }
          }
          del <- numeric()
          del <- del_dup2(PCR.table_fw)
          if (length(del) != 0) {
            PCR.table_fw <- PCR.table_fw[-del,]
          }
          #去除重复
          PCR.table_fw <-
            PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
          #取20条序列
          if (nrow(PCR.table_fw) > 20) {
            PCR.table_fw <- PCR.table_fw[1:20,]
          }
          write.csv(PCR.table_fw,
                    file = "PCR1.csv",
                    row.names = FALSE)
          py$RunOligoCalc(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
          )
          
          # 下游 ----------------------------------------------------------------------
          #去掉连续4个碱基相同的序列
          del <- del_dup(PCR.tab2)
          if (length(del) != 0) {
            PCR.table2 <- PCR.tab2[-del,]
          }
          else{
            PCR.table2 <- PCR.tab2
          }
          #匹配PCR的位置，然后前后扩大3个碱基
          PCR.table_rev <-
            data.frame(PCR = character(),
                       GC = numeric(),
                       TM = numeric())
          m <- 1
          for (i in 1:nrow(PCR.table2)) {
            rev <- matchPattern(pattern = substring(PCR.table2[i, 2], 1, 20),
                                subject = pcr_seq)
            if (length(rev) > 1) {
              next
            }
            rev_start <- start(rev)
            rev_end <- end(rev)
            PCR_pre <-
              substring(PCR_seq, rev_start - 3, rev_end + 3)
            for (j in 20:25) {
              for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                if (j > nchar(PCR_pre2)) {
                  next
                }
                #如果以G或C结尾，则计算GC值和TM值
                if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                    strsplit(PCR_pre2, "")[[1]][j] == "C") {
                  cal <- DNAString(PCR_pre2)
                  GC_count <- sum(letterFrequency(cal, c("G", "C")))
                  GC <- round(GC_count / j, 2)
                  #如果GC大于0.4和小于0.6，计算TM值
                  if (GC >= 0.4 & GC <= 0.6) {
                    TM <-
                      100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                    #如果TM值也符合，则该PCR引物符合
                    if (TM >= 61 & TM <= 64) {
                      PCR.table_rev[m, 1] <- PCR_pre2
                      PCR.table_rev[m, 2] <- GC
                      PCR.table_rev[m, 3] <- TM
                      m <- m + 1
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
            }
          }
          
          del <- numeric()
          del <- del_dup2(PCR.table_rev)
          if (length(del) != 0) {
            PCR.table_rev <- PCR.table_rev[-del,]
          }
          #去除重复
          PCR.table_rev <-
            PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
          #取20条
          if (nrow(PCR.table_rev) > 20) {
            PCR.table_rev <- PCR.table_rev[1:20,]
          }
          write.csv(PCR.table_rev,
                    file = "PCR2.csv",
                    row.names = FALSE)
          py$RunOligoCalc(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
          )
          
          PCR2.table <- read.csv(file = "PCR1.csv", header = TRUE)
          #PCR1的反向互补序列
          for (i in 1:nrow(PCR2.table)) {
            pcr2 <- DNAString(PCR2.table[i,]$PCR)
            pcr2_revcom <- reverseComplement(pcr2)
            PCR2.table[i,]$PCR <- as.character(pcr2_revcom)
          }
          write.csv(PCR2.table, file = "PCR3.csv", row.names = FALSE)
          
          # Blast -------------------------------------------------------------------
          source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
          
          py$run_primertool(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
          )
          py$run_primertool(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
          )
          
          primer.table <-
            read.csv("primer_result1.csv", header = TRUE)
          print(primer.table)
          primer.table <-
            read.csv("primer_result2.csv", header = TRUE)
          print(primer.table)
          
        }
        
        # 正向基因 --------------------------------------------------------------------
        else{
          PCR.tab1 <-
            PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "fw"), ]
          PCR.tab2 <-
            PCR.tab[which(sub("[^a-zA-Z]+", "", PCR.tab[, 1]) == "rev"), ]
          
          # 上游 ----------------------------------------------------------------------
          del <- del_dup(PCR.tab1)
          if (length(del) != 0) {
            PCR.table1 <- PCR.tab1[-del,]
          }
          else{
            PCR.table1 <- PCR.tab1
          }
          #匹配PCR的位置，然后前后扩大3个碱基
          pcr_seq <- DNAString(PCR_seq)
          PCR.table_fw <-
            data.frame(PCR = character(),
                       GC = numeric(),
                       TM = numeric())
          m <- 1
          for (i in 1:nrow(PCR.table1)) {
            fw <- matchPattern(pattern = substring(PCR.table1[i, 2], 1, 20),
                               subject = pcr_seq)
            if (length(fw) > 1) {
              next
            }
            fw_start <- start(fw)
            fw_end <- end(fw)
            PCR_pre <- substring(PCR_seq, fw_start - 3, fw_end + 3)
            for (j in 20:25) {
              for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                if (j > nchar(PCR_pre2)) {
                  next
                }
                #如果以G或C结尾，则计算GC值和TM值
                if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                    strsplit(PCR_pre2, "")[[1]][j] == "C") {
                  cal <- DNAString(PCR_pre2)
                  GC_count <- sum(letterFrequency(cal, c("G", "C")))
                  GC <- round(GC_count / j, 2)
                  #如果GC大于0.4和小于0.6，计算TM值
                  if (GC >= 0.4 & GC <= 0.6) {
                    TM <-
                      100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                    #如果TM值也符合，则该PCR引物符合
                    if (TM >= 61 & TM <= 64) {
                      PCR.table_fw[m, 1] <- PCR_pre2
                      PCR.table_fw[m, 2] <- GC
                      PCR.table_fw[m, 3] <- TM
                      m <- m + 1
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
            }
          }
          del <- numeric()
          del <- del_dup2(PCR.table_fw)
          if (length(del) != 0) {
            PCR.table_fw <- PCR.table_fw[-del,]
          }
          #去除重复
          PCR.table_fw <-
            PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
          #取20条
          if (nrow(PCR.table_fw) > 20) {
            PCR.table_fw <- PCR.table_fw[1:20,]
          }
          write.csv(PCR.table_fw,
                    file = "PCR1.csv",
                    row.names = FALSE)
          py$RunOligoCalc(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
          )
          
          
          
          # 下游 ----------------------------------------------------------------------
          #去掉连续4个碱基相同的序列
          del <- del_dup(PCR.tab2)
          if (length(del) != 0) {
            PCR.table2 <- PCR.tab2[-del,]
          }
          else{
            PCR.table2 <- PCR.tab2
          }
          PCR.table2[, 1] <-
            as.numeric(sub("[^0-9]+", "", PCR.table2[, 1]))
          #匹配PCR的位置，然后前后扩大3个碱基
          PCR.table_rev <-
            data.frame(PCR = character(),
                       GC = numeric(),
                       TM = numeric())
          m <- 1
          for (i in 1:nrow(PCR.table2)) {
            PCR_rev <- DNAString(substring(PCR.table2[i, 2], 1, 20))
            PCR_rev <- reverseComplement(PCR_rev)
            rev <-
              matchPattern(pattern = PCR_rev, subject = pcr_seq)
            if (length(rev) > 1) {
              next
            }
            rev_start <- start(rev)
            rev_end <- end(rev)
            PCR_pre <-
              substring(PCR_seq, rev_start - 3, rev_end + 3)
            for (j in 20:25) {
              for (k in 1:(nchar(PCR_pre) - (j - 1))) {
                PCR_pre2 <- substring(PCR_pre, k, k + (j - 1))
                if (j > nchar(PCR_pre2)) {
                  next
                }
                #如果以G或C结尾，则计算GC值和TM值
                if (strsplit(PCR_pre2, "")[[1]][j] == "G" |
                    strsplit(PCR_pre2, "")[[1]][j] == "C") {
                  cal <- DNAString(PCR_pre2)
                  GC_count <- sum(letterFrequency(cal, c("G", "C")))
                  GC <- round(GC_count / j, 2)
                  #如果GC大于0.4和小于0.6，计算TM值
                  if (GC >= 0.4 & GC <= 0.6) {
                    TM <-
                      100.5 + (41 * (GC_count) / j) - (820 / j) + 16.6 * log(0.05, 10)
                    #如果TM值也符合，则该PCR引物符合
                    if (TM >= 61 & TM <= 64) {
                      PCR.table_rev[m, 1] <- PCR_pre2
                      PCR.table_rev[m, 2] <- GC
                      PCR.table_rev[m, 3] <- TM
                      m <- m + 1
                    }
                    else{
                      next
                    }
                  }
                  else{
                    next
                  }
                }
                else{
                  next
                }
              }
            }
          }
          
          del <- numeric()
          del <- del_dup2(PCR.table_rev)
          if (length(del) != 0) {
            PCR.table_rev <- PCR.table_rev[-del,]
          }
          #去除重复
          PCR.table_rev <-
            PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
          #取20条
          if (nrow(PCR.table_rev) > 20) {
            PCR.table_rev <- PCR.table_rev[1:20,]
          }
          write.csv(PCR.table_rev,
                    file = "PCR2.csv",
                    row.names = FALSE)
          py$RunOligoCalc(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe"
          )
          PCR2.table <- read.csv(file = "PCR2.csv", header = TRUE)
          #PCR2的反向互补序列
          for (i in 1:nrow(PCR2.table)) {
            pcr2 <- DNAString(PCR2.table[i,]$PCR)
            pcr2_revcom <- reverseComplement(pcr2)
            PCR2.table[i,]$PCR <- as.character(pcr2_revcom)
          }
          write.csv(PCR2.table, file = "PCR3.csv", row.names = FALSE)
          
          
          # BLAST -------------------------------------------------------------------
          source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
          py$run_primertool(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR4.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result1.csv",species
          )
          py$run_primertool(
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR5.csv",
            "C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result2.csv",species
          )
          primer.table <-
            read.csv("primer_result1.csv", header = TRUE)
          print(primer.table)
          primer.table <-
            read.csv("primer_result2.csv", header = TRUE)
          print(primer.table)
        }
      }
    }
  }
}
