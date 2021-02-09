setwd("C://Users//41518//Desktop//work/ubigene")
#创建文件夹
filepath<-"C://Users//41518//Desktop//IRF3"
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
source("Get_dot_region5.R")
source("Get_dot_region6.R")
source("Get_result1.R")
source("Get_result2.R")
source("Get_result3.R")
source("GC_analysis.R")
source("gRNA2.R")
source("Get_dot_plot.R")
source("Get_GC_image.R")
source("KO_longregion.R")
source("KO_region_300.R")
source("KO_region_image1.R")
source("KO_region_image2.R")
source("KO_region_image3.R")
source("Get_avoid_region.R")
source("Get_mark_region.R")
source_python("crispr_get_score.py")
source_python("crispor_table_download.py")

#进度条10%
write.table("1",paste0(filepath,"//","10%.txt"),row.names = FALSE,col.names = FALSE)

term <- ("IRF3")
species<-"Human"
ID <- Get_ID(term,species)

# 获取信息 --------------------------------------------------------------------
Gene <- Get_seq(ID)                    #基因序列
Gene2 <- Get_seq2(ID)                 #5'和3'端各增加500bp
allinfo <- Get_allinfo(ID)
start <- allinfo$start
source_python("ensembl_table_download.py")
py$download_csv(ID)
transcript.table <- read.csv("transcript.csv",header = TRUE)

#如果要敲除的基因与其他基因有重叠
othergene<-Get_othergene(species,allinfo$seq_region_name,allinfo$start,allinfo$end)
othergene.table<-othergene[which(othergene$gene_id!=ID),]
avoid_region <- data.frame(start = numeric(), end = numeric())
mark_region<-data.frame()
if(nrow(othergene.table)!=0){
  avoid_region<-Get_avoid_region(othergene.table)
  mark_region<-Get_mark_region(othergene.table)
}

print("ensembl数据调用成功")

#进度条20%
write.table("1",paste0(filepath,"//","20%.txt"),row.names = FALSE,col.names = FALSE)

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

t_Exon_CDS <-  ko.data[which(ko.data$transcript == transcript.name),]

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

#参考转录本的KO区域
KO_region <-
  ko.data[which(ko.data$transcript == transcript.name), ]

#删除不在前30%CDS的
l <- 0
CDS_length <- KO_region$end - KO_region$start + 1
CDS_length_30 <- floor(sum(CDS_length) * 0.33)
for (q in 1:length(CDS_length)) {
  l <- l + CDS_length[q]
  if (l > CDS_length_30) {
    break
  }
}
{
  if (Gene_rev) {
    {
      if (q == 1) {
        t_CDS_30 <- KO_region[q, ]$end -  CDS_length_30
      }
      else{
        t_CDS_30 <-
          KO_region[q, ]$end -  CDS_length_30 + sum(CDS_length[1:(q - 1)])
      }
    }
    
  }
  else{
    {
      if (q == 1) {
        t_CDS_30 <- KO_region[q, ]$start +  CDS_length_30
      }
      else{
        t_CDS_30 <-
          KO_region[q, ]$start +  CDS_length_30 - sum(CDS_length[1:(q - 1)])
      }
    }
  }
}
KO_region <- KO_region[1:q,]

#在前33%的起始密码子
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

KO_region_4 <- KO_region
#KO_region
Exon_150 <- data.frame()
Exon_length_500 <- numeric()
exon_del<-numeric()
k <- 1
Exon_length <- numeric()

if(nrow(KO_region)!=0){
  for (i in 1:nrow(KO_region)) {
    #判断该区域前后500bp是否有其他外显子，有则合并
    Exon.name <- KO_region[i, ]$Exon
    j <- which(t_Exon_region_sort$Exon_name == Exon.name)
    Exon_length[i] <- KO_region[i, ]$end - KO_region[i, ]$start + 1
    ko_start <- KO_region[i, ]$start - 500
    #Exon_150是长度大于150bp的外显子
    if (Exon_length[i] >= 150) {
      Exon_length_500[k] <- KO_region[i, ]$end - KO_region[i, ]$start + 1
      Exon_150 <- rbind(Exon_150, KO_region[i,])
      k <- k + 1
    }
    #往左合并
    repeat {
      if (j == 1) {
        break
      }
      if (ko_start <= t_Exon_region_sort[j - 1, ]$Exon_end ) {
        if(t_Exon_region_sort[j - 1, ]$Exon_name!=t_Exon_CDS[1,]$Exon){
          exon_del<-append(exon_del,i)
          break
        }
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
    #往右合并
    repeat {
      if (j == nrow(t_Exon_region_sort)) {
        break
      }
      if (ko_end >= t_Exon_region_sort[j + 1, ]$Exon_start) {
        if(t_Exon_region_sort[j + 1, ]$Exon_name==t_Exon_CDS[1,]$Exon){
          exon_del<-append(exon_del,i)
          break
        }
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
  if(length(exon_del)!=0){
    KO_region<-KO_region[-exon_del,]
    Exon_length<-Exon_length[-exon_del]
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
  
  if(nrow(Exon_150)!=0){
    Exon_150 <- cbind(Exon_150, Exon_length_500)          #单个大于300bp的外显子
    names(Exon_150)[6] <- "Exon_length"
  }
  
  #times-1
  if(nrow(KO_region_500)!=0){
    for (i in 1:nrow(KO_region_500)) {
      for (j in 1:nrow(Exon_CDS)) {
        if (Exon_CDS[j,]$start >= KO_region_500[i,]$start &
            Exon_CDS[j,]$end <= KO_region_500[i,]$end) {
          KO_region_500[i,]$times <- KO_region_500[i,]$times - 1
        }
      }
    }
  }
  
  if(nrow(Exon_150)!=0){
    for (i in 1:nrow(Exon_150)) {
      for (j in 1:nrow(Exon_CDS)) {
        if (Exon_CDS[j,]$start >= Exon_150[i,]$start &
            Exon_CDS[j,]$end <= Exon_150[i,]$end) {
          Exon_150[i,]$times <- Exon_150[i,]$times - 1
        }
      }
    }
  }
  
  #KO区域中包含了最后一个外显子时
  if(nrow(KO_region_500)!=0)
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
  
  if(nrow(KO_region_500)!=0){
    for(i in 1:nrow(KO_region_500)){
      #把ATG所在外显子合并了，去掉
      if(KO_region_500[i,]$start==min(t_Exon_region$Exon_start)){
        KO_region_500<-KO_region_500[-i,]
      }
    }
  }
  KO_region <- rbind(KO_region_500, Exon_150)
}

Exon_CDS <- Exon_CDS[which(Exon_CDS$Exon_length >= 150),]
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
#删除重复
KO_region<-KO_region[!duplicated(KO_region),]

#如果编码区域小于3KB,可设全敲
if (max(t_Exon_CDS$end) - min(t_Exon_CDS$start) <= 3000) {
  KO_region_all <- t_Exon_CDS[1, ]
  KO_region_all$start <- min(t_Exon_CDS$start)
  KO_region_all$end <- max(t_Exon_CDS$end)
  Exon_length <- sum(t_Exon_CDS$end-t_Exon_CDS$start+1)
  times<-1
  KO_region_all<-cbind(KO_region_all,times,Exon_length)
  KO_region<-rbind(KO_region,KO_region_all)
}

print(KO_region)


if(nrow(KO_region)!=0) {
  #如果与其他基因重叠
  if(nrow(avoid_region)!=0){
    avoid_ko_region<-data.frame()
    avoid_ko_del<-numeric()
    for (i in 1:nrow(KO_region)) {
      for(j in 1:nrow(avoid_region)){
        #如果重叠
        if(!(KO_region[i,]$end<avoid_region[j,]$start | KO_region[i,]$start>avoid_region[j,]$end)){
          avoid_ko_del<-append(avoid_ko_del,i)
          avoid_ko_region<-rbind(avoid_ko_region,KO_region[i,])
        }
      }
    }
    if(length(avoid_ko_del)!=0){
      KO_region <- KO_region[-avoid_ko_del, ]
    }
  }
  
  #GC含量分析:平均GC含量大于70%或小于30%，则删除该区域
  GC_del <- numeric()
  for (i in 1:nrow(KO_region)) {
    analysis_GC <- GC_analysis1(KO_region[i,])
    if (analysis_GC == TRUE) {
      GC_del <- append(GC_del, i)
    }
  }
  if (length(GC_del) != 0) {
    KO_region <- KO_region[-GC_del, ]
  }
  
  
  #点阵图分析
  dot_del <- numeric()
  for (i in 1:nrow(KO_region)) {
    analysis_dot1 <- Get_dot_region1(KO_region[i,])
    analysis_dot2 <- Get_dot_region2(KO_region[i,])
    analysis_dot3 <- Get_dot_region3(KO_region[i,])
    analysis_dot4 <- Get_dot_region4(KO_region[i,])
    if (analysis_dot1 == TRUE | analysis_dot2 == TRUE |
        analysis_dot3 == TRUE | analysis_dot4 == TRUE ) {
      dot_del <- append(dot_del, i)
    }
  }
  if (length(dot_del) != 0) {
    KO_region <- KO_region[-dot_del,]
    if(nrow(KO_region)==0){
      print("序列复杂，无法设计")
    }
  }
}

ko.seq<-data.frame(seq1=character(),seq2=character(),seq3=character(),seq4=character())
if (nrow(KO_region) != 0) {
  for (t in 1:nrow(KO_region)) {
    # 设计在起始密码子后面 -----------------------------------------------------------------
    for (j in 1:nrow(Exon_CDS)) {
      if (all(KO_region[t,] %in% Exon_CDS[j,])) {
        ko_start1 <- KO_region[t,]$start + 3
        ko_end1 <- KO_region[t,]$end
        ko_seq1 <- substring(Gene, ko_start1, ko_end1)
        if (Gene_rev) {
          seq1 <- DNAString(ko_seq1)
          seq1_rev <- reverse(seq1)
          ko_seq1 <- as.character(seq1_rev)
        }
        ko.seq[t,1]<-ko_seq1
      }
    }
    
    # 设计在内含子上 -----------------------------------------------------------------
    if (exists("KO_region_500") == TRUE) {
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
          ko.seq[t,1]<-ko_seq1
          #往前400bp
          if(t_Exon_region_sort[which(t_Exon_region_sort$Exon_start==KO_region[t, ]$start)-1,]$Exon_end+1000<=KO_region[t, ]$start){
            ko_start3 <- KO_region[t, ]$start - 800
            ko_end3 <- KO_region[t, ]$start-400
            ko_seq3 <- substring(Gene, ko_start3, ko_end3)
            if (Gene_rev) {
              seq3 <- DNAString(ko_seq3)
              seq3_rev <- reverse(seq3)
              ko_seq3 <- as.character(seq3_rev)
            }
            ko.seq[t,3]<-ko_seq3
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
          ko.seq[t,2]<-ko_seq2
          #往后400bp
          if(t_Exon_region_sort[which(t_Exon_region_sort$Exon_end==KO_region[t, ]$end)+1,]$Exon_start-1000>=KO_region[t, ]$end){
            print("往后400bp")
            ko_start4 <- KO_region[t, ]$end + 400
            ko_end4 <- KO_region[t, ]$end + 800
            ko_seq4 <- substring(Gene, ko_start4, ko_end4)
            if (Gene_rev) {
              seq4 <- DNAString(ko_seq4)
              seq4_rev <- reverse(seq4)
              ko_seq4 <- as.character(seq4_rev)
            }
            ko.seq[t,4]<-ko_seq4
          }
        }
      }
    }
    
    # 设计在外显子上 -----------------------------------------------------------------
    for (j in 1:nrow(Exon_150)) {
      if (all(KO_region[t, ] %in% Exon_150[j, ])) {
        ko_start1 <- KO_region[t, ]$start
        ko_end1 <- KO_region[t, ]$end
        ko_seq1 <- substring(Gene, ko_start1, ko_end1)
        if(nchar(ko_seq1)>1000){
          ko_seq1<-substring(Gene, ko_start1, ko_start1+1000)
        }
        if (Gene_rev) {
          seq1 <- DNAString(ko_seq1)
          seq1_rev <- reverse(seq1)
          ko_seq1 <- as.character(seq1_rev)
        }
        ko.seq[t,1]<-ko_seq1
      }
    }
  }
}

if (nrow(KO_region) != 0) {
  temp<-list.files(paste0(filepath,"//","gRNA"))
  for (t in 1:nrow(KO_region)) {
    # 设计在起始密码子后面 -----------------------------------------------------------------
    for (j in 1:nrow(Exon_CDS)) {
      if (all(KO_region[t,] %in% Exon_CDS[j,])) {
        gRNA.table <- read.csv(paste0(filepath,"//gRNA//",temp[t],"//gRNA1.csv"), header = FALSE,encoding = "UTF-8")
        
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
        gRNA.table <- gRNA.table[-gRNA.del, ]
        if(nrow(gRNA.table)==0){
          next
        }
        #0-0-0(优化)
        count_0 <- numeric()
        for (p in 1:nrow(gRNA.table)) {
          target <-
            str_extract_all(gRNA.table[p, ]$V9, "\\d+\\s\\-\\s\\d+\\s\\-\\s\\d+")[[1]][2]
          target <- gsub("\\s", "", target)
          if (target == "0-0-0") {
            count_0 <- append(count_0, p)
          }
        }
        gRNA.table_min <- gRNA.table[-count_0, ]
        gRNA.table <- rbind(gRNA.table[count_0, ], gRNA.table_min)
        
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
        GC_avoid_region <- GC_analysis2(KO_region[t,])
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
        
        if (nrow(gRNA.table) > 40) {
          gRNA.table <- gRNA.table[1:40, ]
        }
        
        #符合条件的gRNA进行切割效率预测
        write.csv(gRNA.table, file = "CCTOP-predictor.csv", row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer("CCTOP-predictor.csv",species)
        gRNA.table <-
          read.csv("CCTOP-predictor.csv", header = TRUE)
        print(gRNA.table)
        #切割效率得分低于0.60的删除
        gRNA.table <-
          gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
        
        #切割得分大于0.65的优先
        gRNA.table <-
          rbind(gRNA.table[which(gRNA.table$crispr_score >= 0.65), ], gRNA.table[which(gRNA.table$crispr_score <
                                                                                         0.65), ])
        
        #相差0.05分以内优选
        gRNA <- Get_result2(gRNA.table)
        
        
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA) == "NULL") {
          gRNA <- Get_result3(gRNA.table)
        }
        
        #第二对gRNA
        if (class(gRNA) != "NULL") {
          #在gRNA表格中删除第一对gRNA
          n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
          gRNA.table_2 <- gRNA.table[-n, ]
          #删除重叠的
          gRNA.table_3 <-
            gRNA.table_2[which(
              abs(gRNA.table_2$start - gRNA[1, ]$start) >= 20 &
                abs(gRNA.table_2$start - gRNA[2, ]$start) >= 20
            ), ]
          
          #相差0.05分以内优选
          gRNA2 <- Get_result2(gRNA.table_3)
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA2) == "NULL") {
            gRNA2 <- Get_result3(gRNA.table_3)
          }
          
          if (class(gRNA2) == "NULL") {
            print("换区域找第二对gRNA")
            rm(gRNA2)
          }
          
        }
        if (exists("gRNA") == TRUE) {
          if (KO_region[t, ]$times != transcript.count) {
            No.transcript <-
              ko.data[which(ko.data$end == KO_region[t,]$end),]$transcript
            which.transcript <-
              transcript.table$Name %in% No.transcript
            Not_KO_transcript <-
              transcript.table$Name[which(which.transcript == "FALSE")]
          }
        }
        if (class(gRNA) == "NULL") {
          rm(gRNA)
        }
      }
    }
    if (exists("gRNA") == TRUE) {
      break
    }
    
    # 设计在内含子上 -----------------------------------------------------------------
    if (exists("KO_region_500") == TRUE) {
      for (j in 1:nrow(KO_region_500)) {
        if (all(KO_region[t, ] %in% KO_region_500[j, ])) {

        }
      }
    }
    
    # 设计在外显子上 -----------------------------------------------------------------
    for (j in 1:nrow(Exon_150)) {
      if (all(KO_region[t, ] %in% Exon_150[j, ])) {

      }
    }
    
  }
}










