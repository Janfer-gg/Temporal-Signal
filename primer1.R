
source_python("C://Users//41518//Desktop//work//Ubigene//crispor_table_download1.py")

source_python("C://Users//41518//Desktop//work//Ubigene//OligoCalc.py")
source("del_dup.R")

# 外侧靶位点 -------------------------------------------------------------------
{
  if (min(gRNA$start) < 750 | max(gRNA$end) + 750 > nchar(Gene)) {
    if(exists("gRNA2")==TRUE){
      PCR_seq1 <- substring(Gene2, min(gRNA$start,gRNA2$start) -250, min(gRNA$start,gRNA2$start) + 250)
      PCR_seq2 <- substring(Gene2, max(gRNA$end,gRNA2$end) +750, max(gRNA$end,gRNA2$end) +1250)
    }
    else{
      PCR_seq1 <- substring(Gene2, min(gRNA$start) - 250, min(gRNA$start) + 250)
      PCR_seq2 <- substring(Gene2, max(gRNA$end) +750, max(gRNA$end) +1250)
    }
   }
  else{
    if(exists("gRNA2")==TRUE){
      PCR_seq1 <- substring(Gene, min(gRNA$start,gRNA2$start) - 750, min(gRNA$start,gRNA2$start) - 250)
      PCR_seq2 <- substring(Gene, max(gRNA$end,gRNA2$end) +250, max(gRNA$end,gRNA2$end) +750)
    }
    else{
      PCR_seq1 <- substring(Gene, min(gRNA$start) - 750, min(gRNA$start) - 250)
      PCR_seq2 <- substring(Gene, max(gRNA$end) +250, max(gRNA$end) +750)
    }
   }
}

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
PCR.tab1<-read.csv("PCR.csv", header = FALSE,encoding = "UTF-8")

#下游的PCR引物
py$run2(PCR_seq2,crispor.org)
PCR.tab2<-read.csv("PCR.csv", header = FALSE,encoding = "UTF-8")

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
    {
      if (length(del) != 0) {
        PCR.table1 <- PCR.tab1[-del, ]
      }
      else{
        PCR.table1 <- PCR.tab1
      }
    }
    #40分以上
    PCR.table1<-PCR.table1[which(PCR.table1$V3!="No matches"),]
    PCR.table1<-PCR.table1[which(as.numeric(PCR.table1$V3)>=40),]
    
    PCR.table1[,1]<-as.numeric(sub("[^0-9]+", "", PCR.table1[, 1]))
    PCR.table_300<-PCR.table1[which(PCR.table1[,1]<=300),]
    PCR.table_500<-PCR.table1[which(PCR.table1[,1]>300),]
    PCR.table1<-rbind(PCR.table_300,PCR.table_500)
    #匹配PCR的位置，然后前后扩大3个碱基
    pcr_seq1 <- DNAString(PCR_seq1)
    PCR.table_fw<-data.frame(PCR=character(),GC=numeric(),TM=numeric())
    m<-1
    for(i in 1:nrow(PCR.table1)){
      PCR_fw <- DNAString(substring(PCR.table1[i,2],1,20))
      PCR_fw <- reverseComplement(PCR_fw)
      fw <- matchPattern(pattern = PCR_fw,subject = pcr_seq1)
      if(length(fw)>1){
        next
      }
      fw_start <- start(fw)
      fw_end <- end(fw)
      PCR_pre<-substring(PCR_seq1,fw_start-3,fw_end+3)
      for(j in 20:25){
        for(k in 1:(nchar(PCR_pre)-(j-1))){
          PCR_pre2<-substring(PCR_pre,k,k+(j-1))
          if(j>nchar(PCR_pre2)){
            next
          }
          #如果以G或C结尾，则计算GC值和TM值
          if(strsplit(PCR_pre2,"")[[1]][j]=="G"|strsplit(PCR_pre2,"")[[1]][j]=="C"){
            cal<-DNAString(PCR_pre2)
            GC_count<-sum(letterFrequency(cal,c("G","C")))
            GC<-round(GC_count/j,2)
            #如果GC大于0.4和小于0.6，计算TM值
            if(GC >= 0.4 & GC <= 0.6){
              TM<-100.5 + (41*(GC_count)/j) - (820/j) + 16.6*log(0.05,10)
              #如果TM值也符合，则该PCR引物符合
              if(TM>=61&TM<=64){
                PCR.table_fw[m,1]<-PCR_pre2
                PCR.table_fw[m,2]<-GC
                PCR.table_fw[m,3]<-TM
                m<-m+1
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
    del<-numeric()
    del<-del_dup2(PCR.table_fw)
    if(length(del)!=0){
      PCR.table_fw<-PCR.table_fw[-del,]
    }
    #去除重复
    PCR.table_fw<-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
    #取20条序列
    if(nrow(PCR.table_fw)>15){
      PCR.table_fw<-PCR.table_fw[1:15,]
    }
    write.csv(PCR.table_fw,file = "PCR1.csv",row.names = FALSE)
    py$RunOligoCalc("C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe")
    
    

# 下游 ----------------------------------------------------------------------
    #去掉连续4个碱基相同的序列
    del<-del_dup(PCR.tab2)
    {
      if (length(del) != 0) {
        PCR.table2 <- PCR.tab2[-del, ]
      }
      else{
        PCR.table2 <- PCR.tab2
      }
    }
    #40分以上
    PCR.table2<-PCR.table2[which(PCR.table2$V3!="No matches"),]
    PCR.table2<-PCR.table2[which(as.numeric(PCR.table2$V3)>=40),]
    
    PCR.table2[,1]<-as.numeric(sub("[^0-9]+", "", PCR.table2[, 1]))
    PCR.table_300<-PCR.table2[which(PCR.table2[,1]<=300),]
    PCR.table_500<-PCR.table2[which(PCR.table2[,1]>300),]
    PCR.table2<-rbind(PCR.table_500,PCR.table_300)
    #匹配PCR的位置，然后前后扩大3个碱基
    pcr_seq2 <- DNAString(PCR_seq2)
    PCR.table_rev<-data.frame(PCR=character(),GC=numeric(),TM=numeric())
    m<-1
    for(i in 1:nrow(PCR.table2)){
      rev <- matchPattern(pattern = substring(PCR.table2[i,2],1,20),
                         subject = pcr_seq2)
      if(length(rev)>1){
        next
      }
      rev_start <- start(rev)
      rev_end <- end(rev)
      PCR_pre<-substring(PCR_seq2,rev_start-3,rev_end+3)
      for(j in 20:25){
        for(k in 1:(nchar(PCR_pre)-(j-1))){
          PCR_pre2<-substring(PCR_pre,k,k+(j-1))
          if(j>nchar(PCR_pre2)){
            next
          }
          #如果以G或C结尾，则计算GC值和TM值
          if(strsplit(PCR_pre2,"")[[1]][j]=="G"|strsplit(PCR_pre2,"")[[1]][j]=="C"){
            cal<-DNAString(PCR_pre2)
            GC_count<-sum(letterFrequency(cal,c("G","C")))
            GC<-round(GC_count/j,2)
            #如果GC大于0.4和小于0.6，计算TM值
            if(GC >= 0.4 & GC <= 0.6){
              TM<-100.5 + (41*(GC_count)/j) - (820/j) + 16.6*log(0.05,10)
              #如果TM值也符合，则该PCR引物符合
              if(TM>=61&TM<=64){
                PCR.table_rev[m,1]<-PCR_pre2
                PCR.table_rev[m,2]<-GC
                PCR.table_rev[m,3]<-TM
                m<-m+1
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
    
    del<-numeric()
    del<-del_dup2(PCR.table_rev)
    if(length(del)!=0){
      PCR.table_rev<-PCR.table_rev[-del,]
    }
    #去除重复
    PCR.table_rev<-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
    #取20条
    if(nrow(PCR.table_rev)>15){
      PCR.table_rev<-PCR.table_rev[1:15,]
    }
    write.csv(PCR.table_rev,file = "PCR2.csv",row.names = FALSE)
    py$RunOligoCalc("C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe")
    
    PCR2.table<-read.csv(file = "PCR1.csv",header = TRUE)
    #PCR2的反向互补序列
    for(i in 1: nrow(PCR2.table)){
      pcr2<-DNAString(PCR2.table[i,]$PCR)
      pcr2_revcom<-reverseComplement(pcr2)
      PCR2.table[i,]$PCR<-as.character(pcr2_revcom)
    }
    write.csv(PCR2.table,file = "PCR3.csv",row.names = FALSE)

# Blast -------------------------------------------------------------------
    source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
    py$run_primertool("C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result.csv",species)
    primer.table<-read.csv("primer_result.csv",header = TRUE)
    primer1<-data.frame(PCR=character())
    primer2<-data.frame(PCR=character())
    primer1[1,]<-primer.table[,2]
    primer2[1,]<-primer.table[,1]
    write.csv(primer1,file = "PCR4.csv",row.names = FALSE)
    write.csv(primer2,file = "PCR5.csv",row.names = FALSE)
    print(primer.table)
    {
      if (min(gRNA$start) < 750 | max(gRNA$end) + 750 > nchar(Gene)) {
        primer1.pos <-
          matchPattern(reverse(DNAString(primer1[1,])), subject = Gene2)
        primer2.pos <-
          matchPattern(complement(DNAString(primer2[1,])), subject = Gene2)
        WT.distance <-
          max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) +1
        
      }
      
      else{
        primer1.pos <-
          matchPattern(reverse(DNAString(primer1[1,])), subject = Gene)
        primer2.pos <-
          matchPattern(complement(DNAString(primer2[1,])), subject = Gene)
        WT.distance <-
          max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) +1
      }
    }
  }

# 正向基因 --------------------------------------------------------------------
  else{
    PCR.tab1 <-
      PCR.tab1[which(sub("[^a-zA-Z]+", "", PCR.tab1[, 1]) == "fw"), ]
    PCR.tab2 <-
      PCR.tab2[which(sub("[^a-zA-Z]+", "", PCR.tab2[, 1]) == "rev"), ]

# 上游 ----------------------------------------------------------------------
    del<-del_dup(PCR.tab1)
    {
      if (length(del) != 0) {
        PCR.table1 <- PCR.tab1[-del, ]
      }
      else{
        PCR.table1 <- PCR.tab1
      }
    }
    #40分以上
    PCR.table1<-PCR.table1[which(PCR.table1$V3!="No matches"),]
    PCR.table1<-PCR.table1[which(as.numeric(PCR.table1$V3)>=40),]
    
    PCR.table1[,1]<-as.numeric(sub("[^0-9]+", "", PCR.table1[, 1]))
    PCR.table_300<-PCR.table1[which(PCR.table1[,1]>=300),]
    PCR.table_500<-PCR.table1[which(PCR.table1[,1]<300),]
    PCR.table1<-rbind(PCR.table_300,PCR.table_500)
    #匹配PCR的位置，然后前后扩大3个碱基
    pcr_seq1 <- DNAString(PCR_seq1)
    PCR.table_fw<-data.frame(PCR=character(),GC=numeric(),TM=numeric())
    m<-1
    for(i in 1:nrow(PCR.table1)){
      fw <- matchPattern(pattern = substring(PCR.table1[i,2],1,20),
                         subject = pcr_seq1)
      if(length(fw)>1){
        next
      }
      fw_start <- start(fw)
      fw_end <- end(fw)
      PCR_pre<-substring(PCR_seq1,fw_start-3,fw_end+3)
      for(j in 20:25){
        for(k in 1:(nchar(PCR_pre)-(j-1))){
          PCR_pre2<-substring(PCR_pre,k,k+(j-1))
          if(j>nchar(PCR_pre2)){
            next
          }
          #如果以G或C结尾，则计算GC值和TM值
          if(strsplit(PCR_pre2,"")[[1]][j]=="G"|strsplit(PCR_pre2,"")[[1]][j]=="C"){
            cal<-DNAString(PCR_pre2)
            GC_count<-sum(letterFrequency(cal,c("G","C")))
            GC<-round(GC_count/j,2)
            #如果GC大于0.4和小于0.6，计算TM值
            if(GC >= 0.4 & GC <= 0.6){
              TM<-100.5 + (41*(GC_count)/j) - (820/j) + 16.6*log(0.05,10)
              #如果TM值也符合，则该PCR引物符合
              if(TM>=61&TM<=64){
                PCR.table_fw[m,1]<-PCR_pre2
                PCR.table_fw[m,2]<-GC
                PCR.table_fw[m,3]<-TM
                m<-m+1
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
    del<-numeric()
    del<-del_dup2(PCR.table_fw)
    if(length(del)!=0){
      PCR.table_fw<-PCR.table_fw[-del,]
    }
    #去除重复
    PCR.table_fw<-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
    #取20条
    if(nrow(PCR.table_fw)>15){
      PCR.table_fw<-PCR.table_fw[1:15,]
    }
    write.csv(PCR.table_fw,file = "PCR1.csv",row.names = FALSE)
    py$RunOligoCalc("C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe")
    
    

# 下游 ----------------------------------------------------------------------
    #去掉连续4个碱基相同的序列
    del<-del_dup(PCR.tab2)
    {
      if (length(del) != 0) {
        PCR.table2 <- PCR.tab2[-del, ]
      }
      else{
        PCR.table2 <- PCR.tab2
      }
    }
    #40分以上
    PCR.table2<-PCR.table2[which(PCR.table2$V3!="No matches"),]
    PCR.table2<-PCR.table2[which(as.numeric(PCR.table2$V3)>=40),]
    
    PCR.table2[,1]<-as.numeric(sub("[^0-9]+", "", PCR.table2[, 1]))
    PCR.table_300<-PCR.table2[which(PCR.table2[,1]>300),]
    PCR.table_500<-PCR.table2[which(PCR.table2[,1]<=300),]
    PCR.table2<-rbind(PCR.table_500,PCR.table_300)
    #匹配PCR的位置，然后前后扩大3个碱基
    pcr_seq2 <- DNAString(PCR_seq2)
    PCR.table_rev<-data.frame(PCR=character(),GC=numeric(),TM=numeric())
    m<-1
    for(i in 1:nrow(PCR.table2)){
      PCR_rev <- DNAString(substring(PCR.table2[i,2],1,20))
      PCR_rev <- reverseComplement(PCR_rev)
      rev <- matchPattern(pattern = PCR_rev,subject = pcr_seq2)
      if(length(rev)>1){
        next
      }
      rev_start <- start(rev)
      rev_end <- end(rev)
      PCR_pre<-substring(PCR_seq2,rev_start-3,rev_end+3)
      for(j in 20:25){
        for(k in 1:(nchar(PCR_pre)-(j-1))){
          PCR_pre2<-substring(PCR_pre,k,k+(j-1))
          if(j>nchar(PCR_pre2)){
            next
          }
          #如果以G或C结尾，则计算GC值和TM值
          if(strsplit(PCR_pre2,"")[[1]][j]=="G"|strsplit(PCR_pre2,"")[[1]][j]=="C"){
            cal<-DNAString(PCR_pre2)
            GC_count<-sum(letterFrequency(cal,c("G","C")))
            GC<-round(GC_count/j,2)
            #如果GC大于0.4和小于0.6，计算TM值
            if(GC >= 0.4 & GC <= 0.6){
              TM<-100.5 + (41*(GC_count)/j) - (820/j) + 16.6*log(0.05,10)
              #如果TM值也符合，则该PCR引物符合
              if(TM>=61&TM<=64){
                PCR.table_rev[m,1]<-PCR_pre2
                PCR.table_rev[m,2]<-GC
                PCR.table_rev[m,3]<-TM
                m<-m+1
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
    
    del<-numeric()
    del<-del_dup2(PCR.table_rev)
    if(length(del)!=0){
      PCR.table_rev<-PCR.table_rev[-del,]
    }
    #去除重复
    PCR.table_rev<-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
    #取20条
    if(nrow(PCR.table_rev)>15){
      PCR.table_rev<-PCR.table_rev[1:15,]
    }
    write.csv(PCR.table_rev,file = "PCR2.csv",row.names = FALSE)
    py$RunOligoCalc("C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR2.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\chromedriver.exe")
    PCR2.table<-read.csv(file = "PCR2.csv",header = TRUE)
    #PCR2的反向互补序列
    for(i in 1: nrow(PCR2.table)){
      pcr2<-DNAString(PCR2.table[i,]$PCR)
      pcr2_revcom<-reverseComplement(pcr2)
      PCR2.table[i,]$PCR<-as.character(pcr2_revcom)
    }
    write.csv(PCR2.table,file = "PCR3.csv",row.names = FALSE)
    

# BLAST -------------------------------------------------------------------
    source_python("C://Users//41518//Desktop//work//Ubigene//primertool.py") #人的
    #source_python("C://Users//41518//Desktop//work//Ubigene//primertool2.py") #鼠的
    py$run_primertool("C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR1.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\PCR3.csv","C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer_result.csv",species)
    primer.table<-read.csv("primer_result.csv",header = TRUE)
    primer1<-data.frame(PCR=character())
    primer2<-data.frame(PCR=character())
    primer1[1,]<-primer.table[,1]
    primer2[1,]<-primer.table[,2]
    write.csv(primer1,file = "PCR4.csv",row.names = FALSE)
    write.csv(primer2,file = "PCR5.csv",row.names = FALSE)
    print(primer.table)
    {
      if (min(gRNA$start) < 750 | max(gRNA$end) + 750 > nchar(Gene)) {
        primer1.pos <-
          matchPattern(DNAString(primer1[1,]), subject = Gene2)
        primer2.pos <-
          matchPattern(reverseComplement(DNAString(primer2[1,])), subject = Gene2)
        WT.distance <-
          max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) +1
        
      }
      
      else{
        primer1.pos <-
          matchPattern(DNAString(primer1[1,]), subject = Gene)
        primer2.pos <-
          matchPattern(reverseComplement(DNAString(primer2[1,])), subject = Gene)
        WT.distance <-
          max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) +1
      }
    }
  }
}



#matchPattern(reverse(DNAString(PCR.table_rev[21,]$PCR)),subject = gene)

# rm(PCR.tab1)
# rm(PCR.tab2)
# rm(PCR.table1)
# rm(PCR.table2)
# rm(PCR.table_fw)
# rm(PCR.table_rev)










