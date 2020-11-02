setwd("C://Users//41518//Desktop//work/ubigene")
#创建文件夹
filepath<-"C://Users//41518//Desktop//0925测试//FOS"
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
source("Get_dot_region.R")
source("Get_dot_region2.R")
source("Get_dot_region3.R")
source("GC_analysis.R")
source("Get_dot_plot.R")
source("Get_GC_image.R")
source("KO_longregion.R")
source("KO_region_300.R")
source("KO_region_image2.R")
source("Get_avoid_region.R")
source("Get_mark_region.R")
source("dot_analysis.R")


#进度条10%
write.table("1",paste0(filepath,"//","B10%.txt"),row.names = FALSE,col.names = FALSE)

term <- ("FOS")
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


print("ensembl success")
#进度条20%
write.table("1",paste0(filepath,"//","B20%.txt"),row.names = FALSE,col.names = FALSE)


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

if(transcript.count==0){
  write.table("1",paste0(filepath,"//","no protein.txt"),row.names = FALSE,col.names = FALSE)
}

transcript.name <-
  transcript.table[which(transcript.table$bp == max(transcript.table$bp)), ]$Name[1]    #最长的转录本
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
ATG_Exon<-t_Exon_CDS[1,]$Exon
stop_Exon<-t_Exon_CDS[nrow(t_Exon_CDS),]$Exon
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

# 参考转录本的KO区域 --------------------------------------------------------------
KO_region_20 <-ko.data[which(ko.data$transcript == transcript.name), ]


#删除不在前20%CDS的
l <- 0
CDS_length <- KO_region_20$end - KO_region_20$start + 1
CDS_length_20 <- floor(sum(CDS_length) * 0.2)
for (q in 1:length(CDS_length)) {
  l <- l + CDS_length[q]
  if (l > CDS_length_20) {
    break
  }
}
{
  if (Gene_rev) {
    {
      if (q == 1) {
        t_CDS_20 <- KO_region_20[q, ]$end -  CDS_length_20
      }
      else{
        t_CDS_20 <- KO_region_20[q, ]$end -  CDS_length_20 + sum(CDS_length[1:(q - 1)])
      }
    }
    
  }
  else{
    {
      if (q == 1) {
        t_CDS_20 <- KO_region_20[q, ]$start +  CDS_length_20
      }
      else{
        t_CDS_20 <- KO_region_20[q, ]$start +  CDS_length_20 - sum(CDS_length[1:(q - 1)])
      }
    }
  }
}

KO_region_20 <- KO_region_20[1:q,]

#在前20%的起始密码子
{
  if (Gene_rev) {
    n <- which(Exon_CDS$end > t_CDS_20)
    Exon_CDS <- Exon_CDS[n,]
  }
  else{
    n <- which(Exon_CDS$start < t_CDS_20)
    Exon_CDS <- Exon_CDS[n,]
  }
}


Exon_CDS <- Exon_CDS[!duplicated(Exon_CDS$start), ]
Exon_length <- Exon_CDS$end - Exon_CDS$start + 1
Exon_CDS <- cbind(Exon_CDS, Exon_length)

#判断该外显子是否为起始密码子所在的外显子
region_del <- numeric()
for (i in 1:nrow(KO_region_20)) {
  if (Gene_rev) {
    if (KO_region_20[i,]$end == t_CDS_start) {
      region_del <- append(region_del, i)
    }
  }
  else{
    if (KO_region_20[i,]$start == t_CDS_start) {
      region_del <- append(region_del, i)
    }
  }
}
if (length(region_del) != 0) {
  KO_region_20 <- KO_region_20[-region_del,]
}

Exon_length <- KO_region_20$end-KO_region_20$start+1
#KO区域排序
KO_region_20 <- cbind(KO_region_20, Exon_length)
{
  if (Gene_rev) {
    KO_region_20 <- KO_region_20[order(-KO_region_20$times,-KO_region_20$end),]
  }
  else{
    KO_region_20 <- KO_region_20[order(-KO_region_20$times, KO_region_20$end),]
  }
}

#times-1
if(nrow(Exon_CDS)!=0) {
  if (nrow(KO_region_20) != 0) {
    for (i in 1:nrow(KO_region_20)) {
      for (j in 1:nrow(Exon_CDS)) {
        if (Exon_CDS[j, ]$start >= KO_region_20[i, ]$start &
            Exon_CDS[j, ]$end <= KO_region_20[i, ]$end) {
          KO_region_20[i, ]$times <- KO_region_20[i, ]$times - 1
        }
      }
    }
  }
}

# if(nrow(Exon_CDS)!=0){
#   for (i in 1:nrow(Exon_CDS)) {
#     if(Exon_CDS[i,]$transcript!=transcript.name){
#       if(Exon_CDS[i,]$start %in% t_Exon_CDS$start){
#         Exon_CDS[i,]$end<-t_Exon_CDS[which(t_Exon_CDS$start==Exon_CDS[i,]$start),]$end
#       }
#       else if(Exon_CDS[i,]$end %in% t_Exon_CDS$end){
#         Exon_CDS[i,]$start<-t_Exon_CDS[which(t_Exon_CDS$end==Exon_CDS[i,]$end),]$start
#       }
#     }
#   }
# }

#其他转录本的ATG所在外显子与参考转录本不重叠
if(nrow(Exon_CDS)!=0){
  del<-numeric()
  for(i in 1:nrow(Exon_CDS)){
    times<-0
    for(j in 1:nrow(t_Exon_CDS)){
      if(Exon_CDS[i,]$start %in% c(t_Exon_CDS[j,]$start:t_Exon_CDS[j,]$end)){
        if(Exon_CDS[i,]$end>t_Exon_CDS[j,]$end){
          Exon_CDS[i,]$end<-t_Exon_CDS[j,]$end
        }
      }
      else if(Exon_CDS[i,]$end %in% c(t_Exon_CDS[j,]$start:t_Exon_CDS[j,]$end)){
        if(Exon_CDS[i,]$start<t_Exon_CDS[j,]$start){
          Exon_CDS[i,]$start<-t_Exon_CDS[j,]$start
        }
      }
      else{
        times<-times+1
      }
    }
    if(times==nrow(t_Exon_CDS)){
      del<-append(del,i)
    }
  }
  if(length(del)!=0){
    Exon_CDS<-Exon_CDS[-del,]
  }
}



Exon_CDS$Exon_length <- Exon_CDS$end - Exon_CDS$start + 1
KO_region_20 <- rbind(KO_region_20, Exon_CDS)
KO_region_20 <- KO_region_20[which(KO_region_20$Exon_length>=60),]

KO_region_20<-KO_region_20[!(duplicated(KO_region_20$Exon_length)&duplicated(KO_region_20$start)),]

#排序
{
  if (Gene_rev) {
    KO_region_20 <- KO_region_20[order(-KO_region_20$times,-KO_region_20$end),]
  }
  else{
    KO_region_20 <- KO_region_20[order(-KO_region_20$times, KO_region_20$end),]
  }
}
print(KO_region_20)

#敲除区域不在转录本内
if(nrow(KO_region_20)!=0){
  if(Gene_rev){
    KO_region_20<-KO_region_20[which(KO_region_20$end<=t_CDS_start),]
  }
  else{
    KO_region_20<-KO_region_20[which(KO_region_20$start>=t_CDS_start),]
  }
}

if(nrow(KO_region_20)==0){
  write.table("4",paste0(filepath,"//","result2.txt"),row.names = FALSE,col.names = FALSE)
}

if(nrow(KO_region_20)!=0){
  #如果与其他基因重叠
  if(nrow(avoid_region)!=0){
    avoid_ko_region<-data.frame()
    avoid_ko_del<-numeric()
    for (i in 1:nrow(KO_region_20)) {
      for(j in 1:nrow(avoid_region)){
        #如果重叠
        if(!(KO_region_20[i,]$end<avoid_region[j,]$start | KO_region_20[i,]$start>avoid_region[j,]$end)){
          avoid_ko_del<-append(avoid_ko_del,i)
          avoid_ko_region<-rbind(avoid_ko_region,KO_region_20[i,])
        }
      }
    }
    if(length(avoid_ko_del)!=0){
      KO_region_20 <- KO_region_20[-avoid_ko_del, ]
    }
    if(nrow(KO_region_20)==0){
      write.table("1",paste0(filepath,"//","result2.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}
if(nrow(KO_region_20)!=0){
  #GC含量分析:平均GC含量大于70%或小于30%，则删除该区域
  GC_del <- numeric()
  for (i in 1:nrow(KO_region_20)) {
    analysis_GC <- GC_analysis1(KO_region_20[i, ])
    if (analysis_GC == TRUE) {
      GC_del <- append(GC_del, i)
    }
  }
  if (length(GC_del) != 0) {
    KO_region_20 <- KO_region_20[-GC_del,]
  }
  if(nrow(KO_region_20)==0){
    write.table("2",paste0(filepath,"//","result2.txt"),row.names = FALSE,col.names = FALSE)
  }
}
if(nrow(KO_region_20)!=0){
  #点阵图分析
  dot_del <- numeric()
  for (i in 1:nrow(KO_region_20)) {
    analysis_dot1 <- Get_dot_region1(KO_region_20[i,])
    analysis_dot2 <- Get_dot_region2(KO_region_20[i,])
    analysis_dot4 <- Get_dot_region4(KO_region_20[i,])
    if (analysis_dot1 == TRUE | analysis_dot2 == TRUE |
        analysis_dot4 == TRUE ) {
      dot_del <- append(dot_del, i)
    }
  }
  if (length(dot_del) != 0) {
    KO_region_20 <- KO_region_20[-dot_del,]
  }
  if(nrow(KO_region_20)==0){
    write.table("2",paste0(filepath,"//","result2.txt"),row.names = FALSE,col.names = FALSE)
  }
}

#进度条60%
write.table("1",paste0(filepath,"//","B60%.txt"),row.names = FALSE,col.names = FALSE)
print(KO_region_20)


# 移码gRNA设计方案 --------------------------------------------------------------
gRNA<-data.frame(strand=character(),gRNA_seq=character(),analysis_seq=character(),Score1=numeric(),start=numeric(),end=numeric(),crispr_score=numeric())
if(nrow(KO_region_20)!=0){
  for(t in 1:nrow(KO_region_20)){
    ko_seq3<-substring(Gene,KO_region_20[t,]$start,KO_region_20[t,]$end)
    if (Gene_rev) {
      seq <- DNAString(ko_seq3)
      seq_rev <- reverse(seq)
      ko_seq3 <- as.character(seq_rev)
    }
    if(nchar(ko_seq3)>1500){
      ko_seq3<-substring(ko_seq3, 1, 1500)
    }
    source_python("crispor_table_download.py")
    py$run(ko_seq3,species,filepath,"2")
    gRNA.table <- read.csv(paste0(filepath,"//gRNA2.csv"), header = FALSE, encoding = "UTF-8")
    
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
    count_0<-numeric()
    for (p in 1:nrow(gRNA.table)){
      if (gRNA.table[p,]$V9 == "0-0-0") {
        count_0 <- append(count_0, p)
      }
    }
    if (length(count_0) != 0) {
      gRNA.table_min<-gRNA.table[-count_0,]
      gRNA.table<-rbind(gRNA.table[count_0,],gRNA.table_min)
    }
    #筛选不到gRNA时要及时退出
    if (nrow(gRNA.table) < 1) {
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
    GC_avoid_region <- GC_analysis2(KO_region_20[t, ])
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
    
    {
      if (Gene_rev) {
        gRNA.table <- gRNA.table[which(gRNA.table$start > t_CDS_20), ]
      }
      else{
        gRNA.table <- gRNA.table[which(gRNA.table$end < t_CDS_20), ]
      }
    }
    
    if(nrow(gRNA.table)==0){
      next
    }
    
    gRNA.table_85<-gRNA.table[which(gRNA.table$Score1>=85),]
    gRNA.table_70<-gRNA.table[which(gRNA.table$Score1>=70 & gRNA.table$Score1<85),]
    gRNA.table_60<-gRNA.table[which(gRNA.table$Score1<70),]
    
    gRNA.table<-rbind(gRNA.table_85[order(gRNA.table_85$crispr_score,decreasing = TRUE),],
                      gRNA.table_70[order(gRNA.table_70$crispr_score,decreasing = TRUE),],
                      gRNA.table_60[order(gRNA.table_60$crispr_score,decreasing = TRUE),])
    
    #删除重叠大于4bp的gRNA序列
    w<-1
    repeat{
      s<-which(abs(gRNA.table$end - gRNA.table[w,]$end) <= 20 & abs(gRNA.table$end - gRNA.table[w,]$end)!=0)
      if(length(s)!=0){
        gRNA.table<-gRNA.table[-s,]
      }
      w<-w+1
      if(w >= nrow(gRNA.table)){
        break
      }
    }
    
    #排除完全重叠的情况
    gRNA.table<-gRNA.table[!duplicated(gRNA.table$end),]
   
    {
      #如果有超过3条gRNA,直接取前3条，结束程序
      if (nrow(gRNA.table) >= 3) {
        gRNA <- gRNA.table[1:3, ]
        mark<-"FALSE"
        break
      }
      
      #如果不够3条gRNA,暂时存起来,换一个区域继续找
      else if (nrow(gRNA.table) < 3 & nrow(gRNA.table) > 0) {
        n <- nrow(gRNA)
        if (n == 3) {
          next
        }
        
        m <- nrow(gRNA.table)
        if (n + m > 3) {
          m <- 3 - n
        }
        {
          if (n == 0) {
            gRNA[1:(n + m), ] <- gRNA.table[1:m, ]
            #暂时只考虑前两条gRNA
            if(m==2){
              mark<-"FALSE"
            }
          }
          else{
            gRNA[(n+1):(n + m), ] <- gRNA.table[1:m, ]
          }
        }
      }
    }
  }
}

#进度条90%
write.table("1",paste0(filepath,"//","B90%.txt"),row.names = FALSE,col.names = FALSE)

if(exists("mark")==FALSE){
  mark<-"TRUE"
}

if(nrow(gRNA)!=0){
  if (KO_region_20[t, ]$times != transcript.count) {
    No.transcript <-ko.data[which(ko.data$start <= gRNA[1,]$start & 
                                    ko.data$end >= gRNA[1,]$end),]$transcript
    which.transcript <-
      transcript.table$Name %in% No.transcript
    Not_KO_transcript <-
      transcript.table$Name[which(which.transcript == "FALSE")]
  }
}


if(nrow(gRNA)!=0){
  # 敲了哪些外显子 -----------------------------------------------------------------
  t_Exon_region2 <-t_Exon_region[which(t_Exon_region$Exon_start <= gRNA[1,]$end &
                                         t_Exon_region$Exon_end >= gRNA[1,]$start),]
  t_Exon_CDS2 <- t_Exon_CDS[which(t_Exon_CDS$start <= gRNA[1,]$start &
                                    t_Exon_CDS$end >= gRNA[1,]$end),]
  which_ko <- t_Exon_CDS2$Exon
  
  if(nrow(gRNA)>=2){
    
    t_Exon_region3 <-t_Exon_region[which(t_Exon_region$Exon_start <= gRNA[2,]$end &
                                           t_Exon_region$Exon_end >= gRNA[2,]$start),]
    t_Exon_CDS3<- t_Exon_CDS[which(t_Exon_CDS$start <= gRNA[2,]$start &
                                     t_Exon_CDS$end >= gRNA[2,]$end),]
    
    which_ko2 <- t_Exon_CDS3$Exon
  }
  if(nrow(gRNA)==3){
    t_Exon_region4 <-t_Exon_region[which(t_Exon_region$Exon_start <= gRNA[3,]$end &
                                           t_Exon_region$Exon_end >= gRNA[3,]$start),]
    t_Exon_CDS4<- t_Exon_CDS[which(t_Exon_CDS$start <= gRNA[3,]$start &
                                     t_Exon_CDS$end >= gRNA[3,]$end),]
    
    which_ko3 <- t_Exon_CDS4$Exon
  }
  print(which_ko)
}


# 画图 ----------------------------------------------------------------------
{
  if (nrow(gRNA) != 0) {
    for (j in 1:length(allinfo$Transcript)) {
      if (allinfo$Transcript[[j]]$display_name == transcript.name) {
        start<-allinfo$start
        transcript.start <- allinfo$Transcript[[j]]$start
        transcript.end <-
          allinfo$Transcript[[j]]$end - start + 1
        
        y <- 0.5
        f <- data.frame(x = c(1:transcript.end), y = y)
        p1 <-
          ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#bdc4ca", size = 2) + theme_bw() +
          theme(panel.grid = element_blank(), panel.border = element_blank()) +
          scale_x_discrete(breaks = NULL) + scale_y_discrete(breaks = NULL) + xlab(NULL) +
          ylab(NULL)
        
        #最长的转录本的展示
        CDS_start <-
          allinfo$Transcript[[j]]$Translation$start - start + 1    #CDS起始位置
        CDS_end <-
          allinfo$Transcript[[j]]$Translation$end - start + 1   #CDS终止位置
        {
          if (Gene_rev) {
            label <- paste0(transcript.name, "<")
          }
          else{
            label <- paste0(transcript.name, ">")
          }
        }
        p1 <-
          p1 + annotate(
            "text",
            label = label,
            x = 1,
            y = y - 0.07,
            size = 4,
            hjust = 0,
            color = "orange2"
          )
        for (i in 1:length(allinfo$Transcript[[j]]$Exon)) {
          Exon_start <-
            allinfo$Transcript[[j]]$Exon[[i]]$start - start + 1   #外显子起始位置
          Exon_end <-
            allinfo$Transcript[[j]]$Exon[[i]]$end - start + 1       #外显子终止位置
          p1 <- p1 + annotate(
            "rect",
            xmin = Exon_start,
            xmax = Exon_end,
            ymin = y - 0.04,
            ymax = y + 0.04,
            colour = "orange2",
            alpha = .0
          )
          #完全在编码区
          if (length(CDS_start) != 0 & length(CDS_end) != 0) {
            if (Exon_start >= CDS_start & Exon_end <= CDS_end) {
              p1 <-
                p1 + annotate(
                  "rect",
                  xmin = Exon_start,
                  xmax = Exon_end,
                  ymin = y - 0.04,
                  ymax = y + 0.04,
                  fill = "orange2"
                )
            }
            #部分在编码区(上游)
            else if (Exon_start <= CDS_start &
                     Exon_end >= CDS_start & Exon_end <= CDS_end) {
              p1 <-
                p1 + annotate(
                  "rect",
                  xmin = CDS_start,
                  xmax = Exon_end,
                  ymin = y - 0.04,
                  ymax = y + 0.04,
                  fill = "orange2"
                )
            }
            #部分在编码区(下游)
            else if (Exon_start >= CDS_start &
                     Exon_start <= CDS_end & Exon_end >= CDS_end) {
              p1 <-
                p1 + annotate(
                  "rect",
                  xmin = Exon_start,
                  xmax = CDS_end,
                  ymin = y - 0.04,
                  ymax = y + 0.04,
                  fill = "orange2"
                )
            }
            #部分在编码区(中间)
            else if (Exon_start<CDS_start & Exon_end>CDS_start &Exon_end>CDS_end){
              p1 <-
                p1 + annotate(
                  "rect",
                  xmin = CDS_start,
                  xmax = CDS_end,
                  ymin = y - 0.04,
                  ymax = y + 0.04,
                  fill = "orange2"
                )
            }
          }
        }
      }
    }
    
    p1 <-
      p1 + annotate(
        "rect",
        xmin = min(gRNA$start),
        xmax = max(gRNA$end),
        ymin = y - 0.07,
        ymax = y + 0.07,
        alpha = .0,
        color = "#D01027"
      )
    
    
    
    # 敲除区域放大图 -----------------------------------------------------------------
    large_start <- min(gRNA$start)-300
    large_end <- max(gRNA$end)+300
    y <- 0.5
    f <- data.frame(x = c(large_start:large_end), y = y)
    p1_large <-
      ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "black", size = 2) + theme_bw() +
      theme(panel.grid = element_blank(), panel.border = element_blank()) +
      scale_x_discrete(breaks=NULL)+scale_y_discrete(breaks=NULL)+xlab(NULL)+ylab(NULL)
    {
      if(nrow(t_Exon_region2)==0) {
        start <- large_start
        end <- large_end
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
      else{
        for (i in 1:nrow(t_Exon_region2)) {
          start <- as.numeric(t_Exon_region2[i, ]$Exon_start)
          if (start < large_start) {
            start <- large_start
          }
          end <- as.numeric(t_Exon_region2[i, ]$Exon_end)
          if (end > large_end) {
            end <- large_end
          }
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
      }
      }
    for (i in 1:nrow(t_Exon_CDS2)) {
      start <- as.numeric(t_Exon_CDS2[i,]$start)
      if(start<large_start){
        start<-large_start
      }
      end <- as.numeric(t_Exon_CDS2[i,]$end)
      if(end>large_end){
        end<-large_end
      }
      p1_large <- p1_large + annotate(
        "rect",
        xmin = start,
        xmax = end,
        ymin = y - 0.04,
        ymax = y + 0.04,
        fill = "#D01027",
      )+ annotate(
        "text",
        label = sub("xon ", "", t_Exon_CDS2[i, ]$Exon),
        x = (start+end)/2,
        y = y - 0.07,
        size = 4,
        hjust = 0,
        color = "black"
      )
    }
    
    img<-'cut.png'
    img2<-"cut2.png"
    if(nrow(gRNA)==1){
      ff <-data.frame(x =gRNA[1, ]$start, y = 0.59)
      p1_large <-p1_large + geom_image(data = ff,aes(x = x, y = y),image=img)+annotate("text",label="g1",x = gRNA[1, ]$start,y = 0.67,size=6)
    }
    
    if(nrow(gRNA)==2){
      ff <-data.frame(x =gRNA[1, ]$start, y = 0.59)
      fff <-data.frame(x =gRNA[2, ]$start, y = 0.41)
      p1_large<-p1_large+geom_image(data = ff,aes(x = x, y = y),image=img)+annotate("text",label="g1",x = gRNA[1, ]$start,y = 0.67,size=6)
      +geom_image(data = fff,aes(x = x, y = y),image=img2)+annotate("text",label="g2",x = gRNA[2, ]$start,y = 0.33,size=6)
    }
    if(nrow(gRNA)==3){
      oo<-order(gRNA$start)
      ff <-data.frame(x =gRNA[oo[1], ]$start, y = 0.59)
      fff <-data.frame(x =gRNA[oo[2], ]$start, y = 0.41)
      ffff <-data.frame(x =gRNA[oo[3], ]$start, y = 0.59)
      p1_large<-p1_large+geom_image(data = ff,aes(x = x, y = y),image=img)+annotate("text",label=paste0("g",oo[1]),x = gRNA[oo[1], ]$start,y = 0.67,size=6)+
        geom_image(data = fff,aes(x = x, y = y),image=img2)+annotate("text",label=paste0("g",oo[2]),x = gRNA[oo[2], ]$start,y = 0.34,size=6)+
        geom_image(data = ffff,aes(x = x, y = y),image=img)+annotate("text",label=paste0("g",oo[3]),x = gRNA[oo[3], ]$start,y = 0.67,size=6)
    }
    png(file = paste0(filepath,"//","gRNA_position_large2.png"),width = 480*3,height = 480*2,res = 72*2)
    print(p1_large)
    dev.off()
    png(file = paste0(filepath,"//","gRNA_position2.png"),width = 480*3,height = 480*2,res = 72*2)
    print(p1)
    dev.off()
  }
  
  else{
    print("fail")
    if(nrow(KO_region_20)!=0){
      write.table("3",paste0(filepath,"//","result2.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}


# 输出 ----------------------------------------------------------------------
if (nrow(gRNA) != 0) {
  output <- data.frame(Gene=character(),gRNA1=character(),strand1=character(),
                       score1=numeric(),score2=numeric(),gRNA2=character(),
                       strand2=character(),score2_1=numeric(),score2_2=numeric(),
                       incomplete_transcript=character(),transcript=character(),
                       Exon_count=numeric(),start_condon=character(),stop_condon=character(),
                       ko_condon=character(),Not_KO_transcript=character(),overlap1=character(),
                       tip1=character(),tip2=character(),GC=numeric(),mark=character(),ko_condon2=character(),
                       gRNA3=character(),strand3=character(),score3_1=numeric(),score3_2=numeric(),
                       species=character())
  output[1, 1] <- term
  output[1, 2] <- gRNA[1,]$analysis_seq
  output[1, 3] <- gRNA[1,]$strand
  output[1, 4] <- gRNA[1,]$Score1
  output[1, 5] <- gRNA[1,]$crispr_score
  output[1, 27] <- species
  #有没有重叠的lncRNA,microRNA...
  if(nrow(mark_region)!=0){
    for(j in 1:nrow(mark_region)){
      if(min(gRNA$start)>mark_region[j,]$end | max(gRNA$end)<mark_region[j,]$start){
        overlap1<-"FALSE"
      }
      else{
        overlap1<-"TRUE"
      }
    }
    output[1, 17] <- overlap1
  }
  
  if(nrow(gRNA)>=2){
    output[1, 6] <- gRNA[2,]$analysis_seq
    output[1, 7] <- gRNA[2,]$strand
    output[1, 8] <- gRNA[2,]$Score1
    output[1, 9] <- gRNA[2,]$crispr_score
  }
  if(nrow(gRNA)==3){
    output[1, 23] <- gRNA[3,]$analysis_seq
    output[1, 24] <- gRNA[3,]$strand
    output[1, 25] <- gRNA[3,]$Score1
    output[1, 26] <- gRNA[3,]$crispr_score
  }
  output[1, 10] <- paste(incomplete.transcript,collapse = ",")
  output[1, 11] <- transcript.name
  output[1, 12] <- nrow(t_Exon_region)
  output[1, 13] <- ATG_Exon
  output[1, 14] <- stop_Exon
  output[1, 15] <- paste(which_ko,collapse = ",")
  if(exists("Not_KO_transcript")){
    output[1, 16] <- paste(Not_KO_transcript,collapse = ",")
  }
  
  #发卡结构
  tip1 <- dot_analysis1(KO_region_20[t,])
  #片段重复
  tip2 <- dot_analysis2(KO_region_20[t,])
  #GC含量
  GC <- GC_analysis3(KO_region_20[t,])
  output[1, 18] <- tip1
  output[1, 19] <- tip2
  output[1, 20] <- round(GC*100,2)
  output[1, 21] <- mark
  if(exists("which_ko2")){
    output[1, 22] <-which_ko2
  }
  print(gRNA)
  
  
  # 敲除区域的点阵图和GC含量图 ----------------------------------------------------------
  p3 <- Get_dot_plot(KO_region_20[t,])
  p4 <- Get_GC_image(KO_region_20[t,])
  
  
  png(file = paste0(filepath,"//","Lattice and diagram2.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p3)
  dev.off()
  
  png(file = paste0(filepath,"//","GC content2.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p4)
  dev.off()
  
  # 所有转录本的展示图 ---------------------------------------------------------------
  p2<-KO_region_image2(Gene,allinfo)
  image_write(p2,paste0(filepath,"//","transcript_region2.png"))
  
  write.csv(output, file = paste0(filepath,"//","planB.csv"), row.names = FALSE)
  print("success")
}








