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
source("Get_dot_region.R")
source("Get_dot_region2.R")
source("Get_dot_region3.R")
source("Get_result1.R")
source("Get_result2.R")
source("Get_result3.R")
source("Get_result4.R")
source("GC_analysis.R")
source("gRNA2.R")
source("Get_dot_plot.R")
source("Get_GC_image.R")
source("KO_longregion.R")
source("KO_region_300.R")
source("KO_region_image1.R")
source("Get_avoid_region.R")
source("Get_mark_region.R")
source("dot_analysis.R")


#进度条10%
write.table("1",paste0(filepath,"//","A10%.txt"),row.names = FALSE,col.names = FALSE)

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


print("ensembl success")

#进度条20%
write.table("1",paste0(filepath,"//","A20%.txt"),row.names = FALSE,col.names = FALSE)

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
{
  if (Gene_rev) {
    n <- which(Exon_CDS$end > t_CDS_30)
    Exon_CDS <- Exon_CDS[n,]
  }
  else{
    n <- which(Exon_CDS$start < t_CDS_30)
    Exon_CDS <- Exon_CDS[n,]
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
print(KO_region)
KO_region_4 <- KO_region
#KO_region
Exon_300 <- data.frame()
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
    #Exon_300是长度大于150bp的外显子
    if (Exon_length[i] >= 150) {
      Exon_length_500[k] <- KO_region[i, ]$end - KO_region[i, ]$start + 1
      Exon_300 <- rbind(Exon_300, KO_region[i,])
      k <- k + 1
    }
    #往左合并
    repeat {
      if (j == 1) {
        break
      }
      if (ko_start <= t_Exon_region_sort[j - 1, ]$Exon_end ) {
        if(t_Exon_region_sort[j - 1, ]$Exon_name==t_Exon_CDS[1,]$Exon){
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
          Exon_length[i] + t_Exon_region_sort[j + 1, ]$Exon_end - t_Exon_region_sort[j + 1, ]$Exon_start +1
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
  
  KO_region <- KO_region[!(duplicated(KO_region$Exon_length)& duplicated(KO_region$start)), ]
  
  KO_region <-
    KO_region[which(KO_region$Exon_length %% 3 != 0), ]              #外显子非3的倍数
  KO_region_500 <-
    KO_region[which(KO_region$Exon_length >= 100),]         #外显子不小于100bp
  
  if(nrow(Exon_300)!=0){
    Exon_300 <- cbind(Exon_300, Exon_length_500)          #单个大于150bp的外显子
    names(Exon_300)[6] <- "Exon_length"
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
  
  if(nrow(Exon_300)!=0){
    for (i in 1:nrow(Exon_300)) {
      for (j in 1:nrow(Exon_CDS)) {
        if (Exon_CDS[j,]$start >= Exon_300[i,]$start &
            Exon_CDS[j,]$end <= Exon_300[i,]$end) {
          Exon_300[i,]$times <- Exon_300[i,]$times - 1
        }
      }
    }
  }
  
  #KO区域中包含了最后一个外显子时
  if(nrow(KO_region_500)!=0){
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
        KO_region_500[1,]$Exon_length<-sum((t_Exon_CDS[which(t_Exon_CDS$Exon==KO_region_500[1,]$Exon):nrow(t_Exon_CDS),]$end-t_Exon_CDS[which(t_Exon_CDS$Exon==KO_region_500[1,]$Exon):nrow(t_Exon_CDS),]$start+1))
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
  KO_region <- rbind(KO_region_500, Exon_300)
}

# 
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

#小片段小于3000bp
KO_region<-KO_region[which(KO_region$Exon_length<3000),]

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

#敲除区域不在转录本内
if(nrow(KO_region)!=0){
  if(Gene_rev){
    KO_region<-KO_region[which(KO_region$end<=t_CDS_start),]
  }
  else{
    KO_region<-KO_region[which(KO_region$start>=t_CDS_start),]
  }
}

if(nrow(KO_region)==0){
  write.table("4",paste0(filepath,"//","result1.txt"),row.names = FALSE,col.names = FALSE)
}

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
    if(nrow(KO_region)==0){
      write.table("1",paste0(filepath,"//","result1.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}
  

  #GC含量分析:平均GC含量大于70%或小于30%，则删除该区域
if(nrow(KO_region)!=0) {
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
  if(nrow(KO_region)==0){
    write.table("2",paste0(filepath,"//","result1.txt"),row.names = FALSE,col.names = FALSE)
  }
}
  
  
  #点阵图分析
if(nrow(KO_region)!=0) {
  dot_del <- numeric()
  for (i in 1:nrow(KO_region)) {
    analysis_dot1 <- Get_dot_region1(KO_region[i,])
    analysis_dot2 <- Get_dot_region2(KO_region[i,])
    analysis_dot4 <- Get_dot_region4(KO_region[i,])
    if (analysis_dot1 == TRUE | analysis_dot2 == TRUE |
        analysis_dot4 == TRUE ) {
      dot_del <- append(dot_del, i)
    }
  }
  if (length(dot_del) != 0) {
    KO_region <- KO_region[-dot_del,]
  }
  if(nrow(KO_region)==0){
    write.table("2",paste0(filepath,"//","result1.txt"),row.names = FALSE,col.names = FALSE)
  }
}

#进度条60%
write.table("1",paste0(filepath,"//","A60%.txt"),row.names = FALSE,col.names = FALSE)
print(KO_region)

# gRNA设计方案 ----------------------------------------------------------------
if (nrow(KO_region) != 0) {
  for (t in 1:nrow(KO_region)) {
    # 设计在起始密码子后面 -----------------------------------------------------------------
    for (j in 1:nrow(Exon_CDS)) {
      if (all(KO_region[t,] %in% Exon_CDS[j,])) {
        print("ATG")
        ko_start <- KO_region[t,]$start + 3
        ko_end <- KO_region[t,]$end
        ko_seq <- substring(Gene, ko_start, ko_end)
        if (Gene_rev) {
          seq <- DNAString(ko_seq)
          seq_rev <- reverse(seq)
          ko_seq <- as.character(seq_rev)
        }
        if(nchar(ko_seq)>1500){
          ko_seq<-substring(ko_seq, 1, 1500)
        }
        
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq, species,filepath,"1")
        gRNA.table <-
          read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
        
        #特异性得分60以下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        gRNA.table<-gRNA.table[which(gRNA.table$V3!="No matches"),]
        for (i in 1:length(gRNA.table[, 1])) {
          if (as.numeric(gRNA.table[i, 3]) < 60) {
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
          gRNA.table_min <- gRNA.table[-count_0, ]
          gRNA.table <- rbind(gRNA.table[count_0, ], gRNA.table_min)
        }
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
        write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor1.csv"), row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer(paste0(filepath,"//CCTOP-predictor1.csv"),species)
        gRNA.table <-
          read.csv(paste0(filepath,"//CCTOP-predictor1.csv"), header = TRUE)
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
        if (class(gRNA) == "NULL") {
          gRNA <- Get_result4(gRNA.table)
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
            gRNA2 <- Get_result4(gRNA.table_3)
          }
          if (class(gRNA2) == "NULL") {
            print("change region")
            rm(gRNA2)
          }
          
        }
        if (exists("gRNA") == TRUE) {
          if (KO_region[t, ]$times != transcript.count) {
            No.transcript <-ko.data[which(ko.data$start <= min(gRNA$start) & 
                                            ko.data$end >= max(gRNA$end)),]$transcript
            which.transcript <-
              transcript.table$Name %in% No.transcript
            Not_KO_transcript <-
              transcript.table$Name[which(which.transcript == "FALSE")]
          }
        }
        if (class(gRNA) == "NULL") {
          rm(gRNA)
        }
        
        if (exists("gRNA") == TRUE) {
          break
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
          print("intron")
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
          source_python("crispor_table_download.py")
          py$run(ko_seq1, species,filepath,"1")
          gRNA.table1 <-
            read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
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
            gRNA.table1 <- gRNA.table1[-gRNA.del,]
          }
          
          if (nrow(gRNA.table1) == 0) {
            #往前400bp
            if(t_Exon_region_sort[which(t_Exon_region_sort$Exon_start==KO_region[t, ]$start)-1,]$Exon_end+1000<=KO_region[t, ]$start){
              print("forward 400bp")
              ko_start1 <- KO_region[t, ]$start - 800
              ko_end1 <- KO_region[t, ]$start-400
              ko_seq1 <- substring(Gene, ko_start1, ko_end1)
              if (Gene_rev) {
                seq1 <- DNAString(ko_seq1)
                seq1_rev <- reverse(seq1)
                ko_seq1 <- as.character(seq1_rev)
              }
              #读取ko_seq的gRNA表格
              source_python("crispor_table_download.py")
              py$run(ko_seq1, species,filepath,"1")
              gRNA.table1 <-
                read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
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
                gRNA.table1 <- gRNA.table1[-gRNA.del,]
              }
              if(nrow(gRNA.table1) == 0){
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
            gRNA.table1 <-rbind(gRNA.table1[count_0, ], gRNA.table_min)
          }
          if (nrow(gRNA.table1) > 10) {
            gRNA.table1 <- gRNA.table1[1:10, ]
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
          source_python("crispor_table_download.py")
          py$run(ko_seq2, species,filepath,"1")
          gRNA.table2 <-
            read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
          
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
            gRNA.table2 <- gRNA.table2[-gRNA.del,]
          }
          
          if (nrow(gRNA.table2) == 0) {
            #往后400bp
            if(t_Exon_region_sort[which(t_Exon_region_sort$Exon_end==KO_region[t, ]$end)+1,]$Exon_start-1000>=KO_region[t, ]$end){
              print("back 400bp")
              ko_start2 <- KO_region[t, ]$end + 400
              ko_end2 <- KO_region[t, ]$end + 800
              ko_seq2 <- substring(Gene, ko_start2, ko_end2)
              if (Gene_rev) {
                seq2 <- DNAString(ko_seq2)
                seq2_rev <- reverse(seq2)
                ko_seq2 <- as.character(seq2_rev)
              }
              #读取ko_seq的gRNA表格
              source_python("crispor_table_download.py")
              py$run(ko_seq2, species,filepath,"1")
              gRNA.table2 <-
                read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
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
                gRNA.table2 <- gRNA.table2[-gRNA.del,]
              }
              
              if(nrow(gRNA.table2) == 0){
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
            gRNA.table_min <- gRNA.table2[-count_0,]
            gRNA.table2 <- rbind(gRNA.table2[count_0,], gRNA.table_min)
          }
          if (nrow(gRNA.table2) > 10) {
            gRNA.table2 <- gRNA.table2[1:10, ]
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
          write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor1.csv"), row.names = FALSE)
          source_python("crispr_get_score.py")
          py$reader_writer(paste0(filepath,"//CCTOP-predictor1.csv"), species)
          gRNA.table <-
            read.csv(paste0(filepath,"//CCTOP-predictor1.csv"), header = TRUE)
          print(gRNA.table)
          #切割效率得分低于0.60的删除
          gRNA.table <-
            gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]
          
          #上下游分开
          gRNA.table1 <-
            gRNA.table[which(gRNA.table$end <= KO_region[t, ]$start), ]
          gRNA.table2 <-
            gRNA.table[which(gRNA.table$start >= KO_region[t, ]$end), ]
          
          #切割得分大于0.65的优先
          gRNA.table1 <-
            rbind(gRNA.table1[which(gRNA.table1$crispr_score >= 0.65),], gRNA.table1[which(gRNA.table1$crispr_score <
                                                                                             0.65),])
          gRNA.table2 <-
            rbind(gRNA.table2[which(gRNA.table2$crispr_score >= 0.65),], gRNA.table2[which(gRNA.table2$crispr_score <
                                                                                             0.65),])
          
          #重新合并
          gRNA.table <- rbind(gRNA.table1, gRNA.table2)
          
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
            gRNA <- rbind(gRNA.table1[1,], gRNA.table2[1,])
          }
          {
            #第二对gRNA
            if (class(gRNA) != "NULL" &
                nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
              #在gRNA表格中删除第一对gRNA
              n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
              gRNA.table_2 <- gRNA.table[-n,]
              #删除重叠的
              gRNA.table_3 <-
                gRNA.table_2[which(
                  abs(gRNA.table_2$start - gRNA[1, ]$start) >= 20 &
                    abs(gRNA.table_2$start - gRNA[2, ]$start) >= 20
                ), ]
              #上下游分开
              gRNA.table1_2 <-
                gRNA.table_3[which(gRNA.table_3$end <= KO_region[t,]$start),]
              gRNA.table2_2 <-
                gRNA.table_3[which(gRNA.table_3$start >= KO_region[t,]$end),]
              {
                if (nrow(gRNA.table1_2) != 0 & nrow(gRNA.table2_2) != 0) {
                  #相差0.05分以内优选
                  gRNA2 <- Get_result1(gRNA.table_3, KO_region[t, ])
                  #如果没有相差0.05分以内的gRNA
                  if (class(gRNA2) == "NULL") {
                    gRNA2 <- rbind(gRNA.table1_2[1,], gRNA.table2_2[1,])
                  }
                }
                else{
                  print("change region")
                }
              }
            }
            
            else{
              print("change region")
            }
          }
          
          if (exists("gRNA") == TRUE) {
            if (KO_region[t,]$times != transcript.count) {
              No.transcript <-ko.data[which(ko.data$start >= min(gRNA$start) &
                                ko.data$end <= max(gRNA$end)), ]$transcript
              which.transcript <-
                transcript.table$Name %in% No.transcript
              Not_KO_transcript <-
                transcript.table$Name[which(which.transcript == "FALSE")]
            }
          }
          if (exists("gRNA") == TRUE) {
            judge <- "TRUE"
            break
          }
        }
      }
    }
    if (exists("gRNA") == TRUE) {
      break
    }
    
    # 设计在外显子上 -----------------------------------------------------------------
    for (j in 1:nrow(Exon_300)) {
      if (all(KO_region[t, ] %in% Exon_300[j, ])) {
        print("Exon")
        ko_start <- KO_region[t, ]$start
        ko_end <- KO_region[t, ]$end
        ko_seq <- substring(Gene, ko_start, ko_end)

        if (Gene_rev) {
          seq <- DNAString(ko_seq)
          seq_rev <- reverse(seq)
          ko_seq <- as.character(seq_rev)
        }
        if(nchar(ko_seq)>1500){
          ko_seq<-substring(ko_seq, 1, 1500)
        }
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq, species,filepath,"1")
        gRNA.table <-
          read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
        
        #特异性得分60以下，和Inefficient的gRNA排除掉
        gRNA.del <- numeric()
        gRNA.table<-gRNA.table[which(gRNA.table$V3!="No matches"),]
        for (i in 1:length(gRNA.table[, 1])) {
          if (as.numeric(gRNA.table[i, 3]) < 60) {
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
          gRNA.table <- rbind(gRNA.table[count_0,], gRNA.table_min)
        }
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
        if (nrow(gRNA.table) > 40) {
          gRNA.table <- gRNA.table[1:40,]
        }
        #符合条件的gRNA进行切割效率预测
        write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor1.csv"), row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer(paste0(filepath,"//CCTOP-predictor1.csv"),species)
        gRNA.table <-
          read.csv(paste0(filepath,"//CCTOP-predictor1.csv"), header = TRUE)
        print(gRNA.table)
        #切割效率得分低于0.60的删除
        gRNA.table <-
          gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]

        #切割得分大于0.65的优先
        gRNA.table <-
          rbind(gRNA.table[which(gRNA.table$crispr_score >= 0.65),], gRNA.table[which(gRNA.table$crispr_score <
                                                                                        0.65),])
        
        #相差0.05分以内优选(大于70分)
        gRNA <- Get_result2(gRNA.table)
        
        #如果没有相差0.05分以内的gRNA(大于70分)
        if (class(gRNA) == "NULL") {
          gRNA <- Get_result3(gRNA.table)
        }
        
        if (class(gRNA) == "NULL") {
          gRNA <- Get_result4(gRNA.table)
        }
        
        #第二对gRNA
        if (class(gRNA) != "NULL") {
          #在gRNA表格中删除第一对gRNA
          n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
          gRNA.table_2 <- gRNA.table[-n,]
          #删除重叠的
          gRNA.table_3 <-
            gRNA.table_2[which(
              abs(gRNA.table_2$start - gRNA[1,]$start) >= 20 &
                abs(gRNA.table_2$start - gRNA[2,]$start) >= 20
            ),]
          
          #相差0.05分以内优选
          gRNA2 <- Get_result2(gRNA.table_3)
          #如果没有相差0.05分以内的gRNA
          if (class(gRNA2) == "NULL") {
            gRNA2 <- Get_result3(gRNA.table_3)
          }
          if (class(gRNA2) == "NULL") {
            gRNA2 <- Get_result4(gRNA.table_3)
          }
          
          if (class(gRNA2) == "NULL") {
            print("change region")
            rm(gRNA2)
          }
        }
        
        if (class(gRNA) == "NULL") {
          rm(gRNA)
        }
        
        if (exists("gRNA") == TRUE) {
          if (KO_region[t, ]$times != transcript.count) {
            No.transcript <-ko.data[which(ko.data$start <= min(gRNA$start) & 
                              ko.data$end >= max(gRNA$end)),]$transcript
            which.transcript <-
              transcript.table$Name %in% No.transcript
            Not_KO_transcript <-
              transcript.table$Name[which(which.transcript == "FALSE")]
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
  if (nrow(KO_region_4) != 0) {
    KO_region3 <- KO_region_300(KO_region_4)
    if (nrow(KO_region3) != 0) {
      for (a in 1:nrow(KO_region3)) {
        #外显子上游
        ko_start1 <-
          KO_region3[a,]$start - KO_region3[a,]$left + 1
        ko_end1 <- KO_region3[a,]$start
        ko_seq1 <- substring(Gene, ko_start1, ko_end1)
        if (Gene_rev) {
          seq1 <- DNAString(ko_seq1)
          seq1_rev <- reverse(seq1)
          ko_seq1 <- as.character(seq1_rev)
        }
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq1, species,filepath,"1")
        gRNA.table1 <-
          read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
        #特异性得分60以下，和Inefficient的gRNA排除掉
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
          gRNA.table1 <- gRNA.table1[1:10, ]
        }
        
        #外显子下游
        ko_start2 <- KO_region3[a,]$end
        ko_end2 <- KO_region3[a,]$end + KO_region3[a,]$right
        ko_seq2 <- substring(Gene, ko_start2, ko_end2)
        if (Gene_rev) {
          seq2 <- DNAString(ko_seq2)
          seq2_rev <- reverse(seq2)
          ko_seq2 <- as.character(seq2_rev)
        }
        #读取ko_seq的gRNA表格
        source_python("crispor_table_download.py")
        py$run(ko_seq2, species,filepath,"1")
        gRNA.table2 <-
          read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
        #特异性得分60以下，和Inefficient的gRNA排除掉
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
          gRNA.table2 <- gRNA.table2[1:10, ]
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
        GC_avoid_region <- GC_analysis2(KO_region3[a, ])
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
        write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor1.csv"), row.names = FALSE)
        source_python("crispr_get_score.py")
        py$reader_writer(paste0(filepath,"//CCTOP-predictor1.csv"),species)
        gRNA.table <- read.csv(paste0(filepath,"//CCTOP-predictor1.csv"), header = TRUE)
        print(gRNA.table)
        #切割效率得分低于0.60的删除
        gRNA.table <-
          gRNA.table[which(gRNA.table$crispr_score >= 0.60), ]
        
        #上下游分开
        gRNA.table1 <-
          gRNA.table[which(gRNA.table$end <= KO_region3[a, ]$start), ]
        gRNA.table2 <-
          gRNA.table[which(gRNA.table$start >= KO_region3[a, ]$end), ]
        
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
            gRNA <- Get_result1(gRNA.table, KO_region3[a,])        #相差0.05分以内优选
          }
        }
        #如果没有相差0.05分以内的gRNA
        if (class(gRNA) == "NULL") {
          gRNA <- rbind(gRNA.table1[1,], gRNA.table2[1,])
        }
        
        #第二对gRNA
        if (class(gRNA) != "NULL" &
            nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
          #在gRNA表格中删除第一对gRNA
          n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
          gRNA.table_2 <- gRNA.table[-n,]
          #删除重叠的
          gRNA.table_3 <-
            gRNA.table_2[which(
              abs(gRNA.table_2$start - gRNA[1, ]$start) >= 20 &
                abs(gRNA.table_2$start - gRNA[2, ]$start) >= 20
            ), ]
          #上下游分开
          gRNA.table1_2 <-
            gRNA.table_3[which(gRNA.table_3$end <= KO_region3[a,]$start),]
          gRNA.table2_2 <-
            gRNA.table_3[which(gRNA.table_3$start >= KO_region3[a,]$end),]
          {
            if (nrow(gRNA.table1_2) != 0 & nrow(gRNA.table2_2) != 0) {
              #相差0.05分以内优选
              gRNA2 <- Get_result1(gRNA.table_3, KO_region3[a, ])
              #如果没有相差0.05分以内的gRNA
              if (class(gRNA2) == "NULL") {
                gRNA2 <- rbind(gRNA.table1_2[1,], gRNA.table2_2[1,])
              }
            }
            else{
              print("change region")
            }
          }
        }
        
        if (exists("gRNA") == TRUE) {
          if (KO_region3[a, ]$times != transcript.count) {
            No.transcript <-ko.data[which(ko.data$start >= min(gRNA$start) &
                                            ko.data$end <= max(gRNA$end)), ]$transcript
            which.transcript <-
              transcript.table$Name %in% No.transcript
            Not_KO_transcript <-
              transcript.table$Name[which(which.transcript == "FALSE")]
          }
        }
        if (exists("gRNA") == TRUE) {
          judge3 <- "TRUE"
          break
        }
      }
      
    }
    if(exists("gRNA") == FALSE){
      rm(KO_region3)
    }
  }
}

# 小基因全敲 -------------------------------------------------------------------
if (exists("gRNA") == FALSE) {
  if(exists("KO_region_all")){
    if (all(KO_region_all %in% KO_region[nrow(KO_region), ])) {
      print("all")
      t<-nrow(KO_region)
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
      py$run(ko_seq1, species,filepath,"1")
      gRNA.table1 <-
        read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
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
      if (nrow(gRNA.table1) != 0) {
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
        py$run(ko_seq2, species,filepath,"1")
        gRNA.table2 <-
          read.csv(paste0(filepath,"//gRNA1.csv"), header = FALSE, encoding = "UTF-8")
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
        
        if (nrow(gRNA.table2) != 0) {
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
          
          
          if (nrow(gRNA.table1) != 0 & nrow(gRNA.table2) != 0) {
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
            gene <- DNAString(Gene2)
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
            GC_avoid_region <- GC_analysis2(KO_region_all)
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
            
            if (nrow(gRNA.table) != 0) {
              #符合条件的gRNA进行切割效率预测
              write.csv(gRNA.table, file = paste0(filepath,"//CCTOP-predictor1.csv"), row.names = FALSE)
              source_python("crispr_get_score.py")
              py$reader_writer(paste0(filepath,"//CCTOP-predictor1.csv"), species)
              gRNA.table <-
                read.csv(paste0(filepath,"//CCTOP-predictor1.csv"), header = TRUE)
              print(gRNA.table)
              #切割效率得分低于0.60的删除
              gRNA.table <-
                gRNA.table[which(gRNA.table$crispr_score >= 0.60),]
              
              #上下游分开
              gRNA.table1 <-
                gRNA.table[which(gRNA.table$end <= KO_region_all$start),]
              gRNA.table2 <-
                gRNA.table[which(gRNA.table$start >= KO_region_all$end),]
              
              #切割得分大于0.65的优先
              gRNA.table1 <-
                rbind(gRNA.table1[which(gRNA.table1$crispr_score >= 0.65), ], gRNA.table1[which(gRNA.table1$crispr_score <
                                                                                                  0.65), ])
              gRNA.table2 <-
                rbind(gRNA.table2[which(gRNA.table2$crispr_score >= 0.65), ], gRNA.table2[which(gRNA.table2$crispr_score <
                                                                                                  0.65), ])
              #重新合并
              gRNA.table <- rbind(gRNA.table1, gRNA.table2)
              
              if (nrow(gRNA.table1) != 0 & nrow(gRNA.table2) != 0) {
                gRNA <- Get_result1(gRNA.table, KO_region_all)        #相差0.05分以内优选
                if (class(gRNA) == "NULL") {
                  gRNA <- rbind(gRNA.table1[1, ], gRNA.table2[1, ])
                }
                {
                  #第二对gRNA
                  if (class(gRNA) != "NULL" &
                      nrow(gRNA.table1) >= 2 & nrow(gRNA.table2) >= 2) {
                    #在gRNA表格中删除第一对gRNA
                    n = which(gRNA.table$analysis_seq %in% gRNA$analysis_seq)
                    gRNA.table_2 <- gRNA.table[-n, ]
                    #删除重叠的
                    gRNA.table_3 <-
                      gRNA.table_2[which(
                        abs(gRNA.table_2$start - gRNA[1,]$start) >= 20 &
                          abs(gRNA.table_2$start - gRNA[2,]$start) >= 20
                      ),]
                    #上下游分开
                    gRNA.table1_2 <-
                      gRNA.table_3[which(gRNA.table_3$end <= KO_region_all$start), ]
                    gRNA.table2_2 <-
                      gRNA.table_3[which(gRNA.table_3$start >= KO_region_all$end), ]
                    {
                      if (nrow(gRNA.table1_2) != 0 & nrow(gRNA.table2_2) != 0) {
                        #相差0.05分以内优选
                        gRNA2 <- Get_result1(gRNA.table_3, KO_region_all)
                        #如果没有相差0.05分以内的gRNA
                        if (class(gRNA2) == "NULL") {
                          gRNA2 <- rbind(gRNA.table1_2[1, ], gRNA.table2_2[1, ])
                        }
                      }
                      else{
                        print("change region")
                      }
                    }
                  }
                }
                if (exists("gRNA") == TRUE) {
                  judge4 <- "TRUE"
                }
              }
            }
          }
        }
      }
    }
  }
}

#进度条90%
write.table("1",paste0(filepath,"//","A90%.txt"),row.names = FALSE,col.names = FALSE)

if(exists("judge4")){
  gRNA$start<-gRNA$start-500
  gRNA$end<-gRNA$end-500
}

# 敲了哪些外显子 -----------------------------------------------------------------
if (exists("gRNA") == TRUE) {
  {
    if (exists("judge")) {
      t_Exon_region2 <-
        t_Exon_region[which(
          t_Exon_region$Exon_start >= min(gRNA$start) &
            t_Exon_region$Exon_end <= max(gRNA$end)
        ),]
      t_Exon_CDS2 <-
        t_Exon_CDS[which(t_Exon_CDS$start >= min(gRNA$start) &
                           t_Exon_CDS$end <= max(gRNA$end)),]
      which_ko <- t_Exon_CDS2$Exon
      
    }
    else if (exists("judge3")) {
      t_Exon_region2 <-
        t_Exon_region[which(
          t_Exon_region$Exon_start >= min(gRNA$start) &
            t_Exon_region$Exon_end <= max(gRNA$end)
        ),]
      t_Exon_CDS2 <-
        t_Exon_CDS[which(t_Exon_CDS$start >= min(gRNA$start) &
                           t_Exon_CDS$end <= max(gRNA$end)),]
      which_ko <- t_Exon_CDS2$Exon
      
    }
    else if (exists("judge4")){
      t_Exon_region2 <-
        t_Exon_region[which(
          t_Exon_region$Exon_start >= min(gRNA$start) &
            t_Exon_region$Exon_end <= max(gRNA$end)
        ),]
      t_Exon_CDS2 <-
        t_Exon_CDS[which(t_Exon_CDS$start >= min(gRNA$start) &
                           t_Exon_CDS$end <= max(gRNA$end)),]
      which_ko <- t_Exon_CDS2$Exon
    }
    else{
      t_Exon_region2 <-
        t_Exon_region[which(
          t_Exon_region$Exon_start <= max(gRNA$end) &
            t_Exon_region$Exon_end >= min(gRNA$start)
        ),]
      t_Exon_CDS2 <-
        t_Exon_CDS[which(t_Exon_CDS$start <= max(gRNA$end) &
                           t_Exon_CDS$end >= min(gRNA$start)),]
      which_ko <- t_Exon_CDS2$Exon
    }
  }
}


# 画图 ----------------------------------------------------------------------
{
  if (exists("gRNA") == TRUE) {
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
    ff <-data.frame(x = gRNA[1, ]$start, y = 0.59)
    fff <-data.frame(x = gRNA[2, ]$start, y = 0.59)
    large_start <- min(gRNA$start) - 300
    large_end <- max(gRNA$end) + 300
    y <- 0.5
    f <- data.frame(x = c(large_start:large_end), y = y)
    p1_large <-
      ggplot(data = f, aes(x = x, y = y)) + geom_path(color = "#333333", size = 2) + theme_bw() +
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
        fill = "#D01027"
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
    p1_large <-
      p1_large + geom_image(data = ff,aes(x = x, y = y),image=img,) +
      geom_image(data = fff,aes(x = x, y = y),image=img) + 
      annotate("text",label="g1",x = gRNA[1, ]$start,y = 0.67,size=6) +
      annotate("text",label="g2",x = gRNA[2, ]$start,y = 0.67,size=6)
    png(file = paste0(filepath,"//","gRNA_position_large.png"),width = 480*3,height = 480*2,res = 72*2)
    print(p1_large)
    dev.off()
    png(file = paste0(filepath,"//","gRNA_position.png"),width = 480*3,height = 480*2,res = 72*2)
    print(p1)
    dev.off()
  }
  
  else{
    print("fail")
    if(nrow(KO_region)!=0 | exists("gRNA.table")){
      write.table("3",paste0(filepath,"//","result1.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}


# 输出 ----------------------------------------------------------------------
#敲除大小
if (exists("gRNA") == TRUE) {
  if (Gene_rev) {
    if (all(gRNA$strand == "rev") |
        all(gRNA$strand == "fw")) {
      KO_length <- abs(gRNA[1,]$end - gRNA[2,]$end) 
    }
    else{
      if (gRNA[1, ]$strand == "rev" & gRNA[2, ]$strand == "fw") {
        if(gRNA[1,]$start>gRNA[2,]$start){
          pos1 <- gRNA[1, ]$end-3
          pos2 <- gRNA[2, ]$start+3
          KO_length <- abs(pos1 - pos2) +1
        }
        else{
          pos1 <- gRNA[1, ]$end-3
          pos2 <- gRNA[2, ]$start+3
          KO_length <- abs(pos1 - pos2) -1
        }
      }
      else{
        if(gRNA[1,]$start>gRNA[2,]$start){
          pos1 <- gRNA[1, ]$start+3
          pos2 <- gRNA[2, ]$end-3
          KO_length <- abs(pos1 - pos2) - 1
        }
        #正确
        else{
          pos1 <- gRNA[1, ]$start+3
          pos2 <- gRNA[2, ]$end-3
          KO_length <- abs(pos1 - pos2) + 1
        }
      }
    }
  }
  #正向
  else{
    if (all(gRNA$strand == "rev") |
        all(gRNA$strand == "fw")) {
      KO_length <- abs(gRNA[1,]$end - gRNA[2,]$end) 
    }
    else{
      if (gRNA[1, ]$strand == "rev" & gRNA[2, ]$strand == "fw") {
        pos1 <- gRNA[1, ]$start+3
        pos2 <- gRNA[2, ]$end-3
        KO_length <- abs(pos1 - pos2-1) 
      }
      else{
        pos1 <- gRNA[1, ]$end-3
        pos2 <- gRNA[2, ]$start+3
        KO_length <- abs(pos1 - pos2+1) 
      }
    }
  }
  
  #敲除的CDS
  if (exists("judge")) {
    KO_length_CDS <- KO_region[t,]$Exon_length
  }
  else if (exists("judge3")) {
    KO_length_CDS <- KO_region3[a,]$Exon_length
  }
  else if (exists("judge4")) {
    KO_length_CDS <- KO_region_all$Exon_length
  }
  else{
    KO_length_CDS <- KO_length
  }
  
  {
    if (exists("gRNA2")==TRUE) {
      if (Gene_rev) {
        if (all(gRNA2$strand == "rev") |
            all(gRNA2$strand == "fw")) {
          KO_length2 <- abs(gRNA2[1, ]$end - gRNA2[2, ]$end)
        }
        else{
          if (gRNA2[1,]$strand == "rev" & gRNA2[2,]$strand == "fw") {
            if(gRNA2[1,]$start>gRNA2[2,]$start){
              pos1 <- gRNA2[1,]$end-3
              pos2 <- gRNA2[2,]$start+3
              KO_length2 <- abs(pos1 - pos2) +1
            }
            else{
              pos1 <- gRNA2[1,]$end-3
              pos2 <- gRNA2[2,]$start+3
              KO_length2 <- abs(pos1 - pos2) -1
            }
          }
          else{
            if(gRNA2[1,]$start>gRNA2[2,]$start){
              pos1 <- gRNA2[1, ]$start+3
              pos2 <- gRNA2[2, ]$end-3
              KO_length2 <- abs(pos1 - pos2) - 1
            }
            else{
              pos1 <- gRNA2[1, ]$start+3
              pos2 <- gRNA2[2, ]$end-3
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
            pos1 <- gRNA2[1,]$start+3
            pos2 <- gRNA2[2,]$end-3
            KO_length2 <- abs(pos1 - pos2-1) 
          }
          else{
            pos1 <- gRNA2[1,]$end-3
            pos2 <- gRNA2[2,]$start+3
            KO_length2 <- abs(pos1 - pos2+1) 
          }
        }
      }
      
      #敲除的CDS
      if (exists("judge")) {
        KO_length_CDS2 <- KO_region[t,]$Exon_length
      }
      else if (exists("judge3")) {
        KO_length_CDS2 <- KO_region3[a,]$Exon_length
      }
      else if (exists("judge4")) {
        KO_length_CDS2 <- KO_region_all$Exon_length
      }
      else{
        KO_length_CDS2 <- KO_length2
      }
      mark<-"FALSE"
    }
    
    else{
      gRNA2<-Get_gRNA2(KO_region)
      mark<-"TRUE"        #判断两对是否在同一区域
    }
  }
}
if(exists("gRNA2")==TRUE){
  if(class(gRNA2)=="NULL"){
    rm(mark)
    rm(gRNA2)
  }
}

if (exists("gRNA") == TRUE) {
  output <- data.frame(Gene=character(),gRNA1=character(),strand1=character(),score1_1=numeric(),score1_2=numeric(),
                       gRNA2=character(),strand2=character(),score2_1=numeric(),score2_2=numeric(),region=numeric(),cds=numeric(),
                       gRNA3=character(),strand3=character(),score3_1=numeric(),score3_2=numeric(),
                       gRNA4=character(),strand4=character(),score4_1=numeric(),score4_2=numeric(),region2=numeric(),cds2=numeric(),
                       incomplete_transcript=character(),transcript=character(),Exon_count=numeric(),
                       start_condon=character(),stop_condon=character(),ko_condon=character(),Not_KO_transcript=character(),
                       overlap1=character(),overlap2=character(),tip1=character(),tip2=character(),GC=numeric(),
                       mark=character(),ko_condon2=character(),species=character())
  output[1, 1] <- term
  output[1, 2] <- gRNA[1,]$analysis_seq
  output[1, 3] <- gRNA[1,]$strand
  output[1, 4] <- gRNA[1,]$Score1
  output[1, 5] <- gRNA[1,]$crispr_score
  output[1, 6] <- gRNA[2,]$analysis_seq
  output[1, 7] <- gRNA[2,]$strand
  output[1, 8] <- gRNA[2,]$Score1
  output[1, 9] <- gRNA[2,]$crispr_score
  output[1, 10] <- KO_length
  output[1, 11] <- KO_length_CDS
  output[1, 36] <- species
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
    output[1, 29] <- overlap1
  }
  
  if(exists("gRNA2")=="TRUE"){
    output[1, 12] <- gRNA2[1,]$analysis_seq
    output[1, 13] <- gRNA2[1,]$strand
    output[1, 14] <- gRNA2[1,]$Score1
    output[1, 15] <- gRNA2[1,]$crispr_score
    output[1, 16] <- gRNA2[2,]$analysis_seq
    output[1, 17] <- gRNA2[2,]$strand
    output[1, 18] <- gRNA2[2,]$Score1
    output[1, 19] <- gRNA2[2,]$crispr_score
    if(mark=="TRUE"){
      output[1, 20] <- gRNA2[1,8]
      output[1, 21] <- gRNA2[2,8]
      output[1, 35] <- gRNA2[1,9]
    }
    else if(mark=="FALSE"){
      output[1, 20] <- KO_length2
      output[1, 21] <- KO_length_CDS2
    }
    output[1, 34] <- mark
    #有没有重叠的lncRNA,microRNA...
    if(nrow(mark_region)!=0){
      for(j in 1:nrow(mark_region)){
        if(min(gRNA2$start)>mark_region[j,]$end | max(gRNA2$end)<mark_region[j,]$start){
          overlap2<-"FALSE"
        }
        else{
          overlap2<-"TRUE"
        }
      }
      output[1, 30] <- overlap2
    }
  }
  output[1, 22] <- paste(incomplete.transcript,collapse = ",")
  output[1, 23] <- transcript.name
  output[1, 24] <- nrow(t_Exon_region)
  output[1, 25] <- ATG_Exon
  output[1, 26] <- stop_Exon
  output[1, 27] <- paste(which_ko,collapse = ",")
  if(exists("Not_KO_transcript")){
    output[1, 28] <-paste(Not_KO_transcript,collapse = ",")
  }
  
  print(gRNA)
  
  # 敲除区域的点阵图和GC含量图 ----------------------------------------------------------
  {
    if (exists("judge3")) {
      p3<-Get_dot_plot(KO_region3[a, ])
      p4<-Get_GC_image(KO_region3[a, ])
      #发卡结构
      tip1 <- dot_analysis1(KO_region3[a,])
      #片段重复
      tip2 <- dot_analysis2(KO_region3[a,])
      #GC含量
      GC <- GC_analysis3(KO_region3[a,])
    }
    else{
      p3<-Get_dot_plot(KO_region[t, ])
      p4<-Get_GC_image(KO_region[t, ])
      #发卡结构
      tip1 <- dot_analysis1(KO_region[t,])
      #片段重复
      tip2 <- dot_analysis2(KO_region[t,])
      #GC含量
      GC <- GC_analysis3(KO_region[t,])
    }
  }
  output[1, 31] <- tip1
  output[1, 32] <- tip2
  output[1, 33] <- round(GC*100,2)
  
  png(file = paste0(filepath,"//","Lattice and diagram.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p3)
  dev.off()
  
  png(file = paste0(filepath,"//","GC content.png"),width = 480*3,height = 480*3,res = 72*3)
  print(p4)
  dev.off()
  
  # 所有转录本的展示图 ---------------------------------------------------------------
  # 三种方案共用
  p2<-KO_region_image1(Gene,allinfo)
  image_write(p2,paste0(filepath,"//","transcript_region1.png"))
  
  write.csv(output, file = paste0(filepath,"//","planA.csv"), row.names = FALSE)
  print("success")
}

