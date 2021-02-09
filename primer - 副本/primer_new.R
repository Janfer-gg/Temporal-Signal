setwd("C://Users//41518//Desktop//work/ubigene/primer")
filepath<-"C://Users//41518//Desktop//靶位点测试3//KIF7"
filepath2<-"C://Users//41518//Desktop//靶位点测试3//KIF7//primer_A"
dir.create(filepath2)
library(httr)
library(jsonlite)
 library(xml2)
library(Biostrings)
library(reticulate)
library(stringr)
library(ggplot2)
# library(ggpubr)
source("Get_ID.R")
source("Get_allinfo.R")
source("Get_seq.R")
source("delete_noprotein.R")
source("Get_max_transcript.R")
source("del_dup.R")
source("Get_gRNA_pos.R")
source("primer_design.R")
source("primer_design_in.R")
source("Get_size.R")
source("primer_dis.R")
source("Get_primer.R")
source("primer_image.R")
source("analysis_GC.R")
source_python("crispor_table_download2.py")

PlanA.table<-read.csv(paste0(filepath,"//planA.csv"))
term<-PlanA.table$Gene

species<-"Human"

ID <- Get_ID(term,species)
Gene <- Get_seq(ID)                    #基因序列(左右两端各加1500bp)
allinfo <- Get_allinfo(ID)
start <- allinfo$start-1500

transcript.table <- read.csv(paste0(filepath,"//transcript.csv"),header = TRUE)

transcript.table <-
  delete_noprotein(transcript.table)           #删除非编码蛋白和不完整的转录本
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

t_Exon_CDS <- ko.data[which(ko.data$transcript == transcript.name),]

gene<-DNAString(Gene)

#如果mark是TRUE,两对gRNA不同区域；如果mark是NA,没有第二对gRNA；如果mark是FALSE,两对gRNA同区域

gRNA<-data.frame()
gRNA[1,1]<-PlanA.table$gRNA1
gRNA[2,1]<-PlanA.table$gRNA2
gRNA<-Get_gRNA_pos(gRNA)
size1<-Get_size(PlanA.table$region)


# 设计一组靶位点 -----------------------------------------------------------------
{
  #如果有两对gRNA
  if(!is.na(PlanA.table$mark)) {
    gRNA2 <- data.frame()
    gRNA2[1, 1] <- PlanA.table$gRNA3
    gRNA2[2, 1] <- PlanA.table$gRNA4
    gRNA2 <- Get_gRNA_pos(gRNA2)
    size2<-Get_size(PlanA.table$region2)
    
    #两对gRNA在同一个区域，一起设计靶位点
    if (PlanA.table$mark == "FALSE") {
      result1 <- primer_design2(Gene, gRNA, gRNA2 ,Gene_rev,size1)
      KO_length1 <- result1$WT - PlanA.table$region
      KO_length2 <- result1$WT - PlanA.table$region2
    }
  }
  
  
  #如果只有一对gRNA
  else{
    result1 <- primer_design1(Gene, gRNA, Gene_rev,size1)
    KO_length1 <- result1$WT - PlanA.table$region
  }
  
  #可删
  write.csv(result1,paste0(filepath2,"//primer_resultA.csv"),row.names = FALSE)
  if(class(result1)!="NULL"){
    #内侧
    ko_Exon1<-PlanA.table$ko_condon
    ko_Exon1<-strsplit(ko_Exon1,",")[[1]]
    #如果只敲一个外显子
    if(length(ko_Exon1)==1){
      Exon_region<-t_Exon_CDS[which(t_Exon_CDS$Exon %in% ko_Exon1),]
      PCR_seq <- substring(Gene, min(gRNA$end) + 1, max(gRNA$start) - 1)
      if(PlanA.table$region<30){
        result_in<-character()
      }
      #一个内侧靶位点
      if (PlanA.table$region >=30 & PlanA.table$region <=200) {
        result_in <- primer_design3_in(Gene, PCR_seq, Gene_rev, result1)
      }
      
      #两个内侧靶位点
      if (PlanA.table$region > 200) {
        result_in <- primer_design_in(Gene, PCR_seq, Gene_rev, result1)
      }
    }

    #敲多个外显子
    if(length(ko_Exon1)>1){
      Exon_region<-t_Exon_CDS[which(t_Exon_CDS$Exon %in% ko_Exon1),]
      #先找到外侧引物的位置
      {
        if (Gene_rev) {
          F1.pos <-matchPattern(reverse(DNAString(result1$primer1)), subject = Gene)
          R1.pos <-matchPattern(complement(DNAString(result1$primer2)), subject = Gene)
          {
            if (end(F1.pos) - 999 < Exon_region[1, ]$start) {
              PCR_seq1 <- substring(Gene, end(F1.pos) - 999, Exon_region[1, ]$end)
              Exon_pos<-data.frame(start=Exon_region[1, ]$end-Exon_region[which(Exon_region$end>=(end(F1.pos) - 999)),]$end+1,
                                   end=Exon_region[1, ]$end-Exon_region[which(Exon_region$end>=(end(F1.pos) - 999)),]$start+1)
            }
            else{
              PCR_seq1 <-substring(Gene, Exon_region[1, ]$start, Exon_region[1, ]$end)
              Exon_pos<-data.frame(start=1,end=2000)
            }
          }
          
          {
            if (start(R1.pos) + 999 > Exon_region[nrow(Exon_region), ]$end) {
              PCR_seq2 <-substring(Gene, Exon_region[nrow(Exon_region), ]$start, start(R1.pos) +999)
              Exon_pos2<-data.frame(start=(start(R1.pos) + 999)-Exon_region[which(Exon_region$start<=(start(R1.pos) + 999)),]$end+1,
                                   end=(start(R1.pos) + 999)-Exon_region[which(Exon_region$start<=(start(R1.pos) + 999)),]$start+1)
            }
            else{
              PCR_seq2 <-substring(Gene, Exon_region[nrow(Exon_region), ]$start, Exon_region[nrow(Exon_region), ]$end)
              Exon_pos2<-data.frame(start=1,end=2000)
            }
          }
        }
        else{
          F1.pos <-matchPattern(DNAString(result1$primer1), subject = Gene)
          R1.pos <-matchPattern(reverseComplement(DNAString(result1$primer2)), subject = Gene)
          {
            if(start(F1.pos)+999>Exon_region[1,]$end){
              PCR_seq1 <- substring(Gene,Exon_region[1,]$start,start(F1.pos)+999)
              Exon_pos<-data.frame(start=Exon_region[which(Exon_region$start<=(start(F1.pos)+999)),]$start-Exon_region[1,]$start+1,
                                   end=Exon_region[which(Exon_region$start<=(start(F1.pos)+999)),]$end-Exon_region[1,]$start+1)
            }
            else{
              PCR_seq1 <- substring(Gene,Exon_region[1,]$start,Exon_region[1,]$end)
              Exon_pos<-data.frame(start=1,end=2000)
            }
          }
          {
            if(end(R1.pos)-999<Exon_region[nrow(Exon_region),]$start){
              PCR_seq2 <- substring(Gene,end(R1.pos)-999,Exon_region[nrow(Exon_region),]$end)
              Exon_pos2<-data.frame(start=Exon_region[which(Exon_region$end>=(end(R1.pos)-999)),]$start-(end(R1.pos)-999)+1,
                                   end=Exon_region[which(Exon_region$end>=(end(R1.pos)-999)),]$end-(end(R1.pos)-999)+1)
            }
            else{
              PCR_seq2 <- substring(Gene,Exon_region[nrow(Exon_region),]$start,Exon_region[nrow(Exon_region),]$end)
              Exon_pos2<-data.frame(start=1,end=2000)
            }
          }
        }
      }
      result_in<-primer_design2_in(Gene,PCR_seq1,PCR_seq2,Gene_rev,result1,Exon_pos,Exon_pos2)
    }
  }
  
  
  # 输出 ----------------------------------------------------------------------
  {
    if(length(result1)!=0) {
      if(length(result_in)!=0){
        #两个内侧靶位点
        if(PlanA.table$region > 200){
          output <- WT.calculation(result_in)
          p1 <- image1(output,Exon_region)
        }
        #一个内侧靶位点
        else{
          output<-result_in
          p1 <- image3(output,Exon_region)
        }
      }
      else{
        output <- data.frame(F1=result1$primer1,R1=result1$primer2,R2=NA,F2=NA,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result1$GC,GC2=NA,GC3=NA)
        p1<-image2(result1,Exon_region)
      }
      write.csv(output,paste0(filepath2,"//resultA.csv"),row.names = FALSE)
      png(file = paste0(filepath2,"//","primer_image.png"),width = 480*3,height = 480*2,res = 72*2)
      print(p1)
      dev.off()
    }
    
    else{
      write.table("fail",paste0(filepath2,"//fail.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
}





# 设计两组靶位点 ---------------------------------------------------
if (!is.na(PlanA.table$mark) & PlanA.table$mark == "TRUE") {
  #第一对gRNA
  result1 <- primer_design1(Gene, gRNA, Gene_rev,size1)
  #可删
  write.csv(result1,paste0(filepath2,"//primer_resultA1.csv"),row.names = FALSE)
  KO_length1 <- result1$WT - PlanA.table$region
  #第二对gRNA
  size2<-Get_size(PlanA.table$region2)
  result2 <- primer_design1(Gene, gRNA2, Gene_rev,size2)
  write.csv(result1,paste0(filepath2,"//primer_resultA2.csv"),row.names = FALSE)
  KO_length2 <- result2$WT - PlanA.table$region2
  
  # 内侧 ----------------------------------------------------------------------
  # 第一对gRNA -----------------------------------------------------------------
  if(length(result1)!=0){
    ko_Exon1<-PlanA.table$ko_condon
    ko_Exon1<-strsplit(ko_Exon1,",")[[1]]
    #如果只敲一个外显子
    if(length(ko_Exon1)==1){
      Exon_region<-t_Exon_CDS[which(t_Exon_CDS$Exon %in% ko_Exon1),]
      PCR_seq <- substring(Gene, min(gRNA$end) + 1, max(gRNA$start) - 1)
      if(PlanA.table$region<30){
        result_in<-character()
      }
      #一个内侧靶位点
      if (PlanA.table$region >=30 & PlanA.table$region <=200) {
        result_in <- primer_design3_in(Gene, PCR_seq, Gene_rev, result1)
      }
      
      #两个内侧靶位点
      if (PlanA.table$region > 200) {
        result_in <- primer_design_in(Gene, PCR_seq, Gene_rev, result1)
      }
    }

    #敲多个外显子
    if(length(ko_Exon1)>1){
      Exon_region<-t_Exon_CDS[which(t_Exon_CDS$Exon %in% ko_Exon1),]
      #先找到外侧引物的位置
      {
        if (Gene_rev) {
          F1.pos <-matchPattern(reverse(DNAString(result1$primer1)), subject = Gene)
          R1.pos <-matchPattern(complement(DNAString(result1$primer2)), subject = Gene)
          {
            if (end(F1.pos) - 999 < Exon_region[1, ]$start) {
              PCR_seq1 <- substring(Gene, end(F1.pos) - 999, Exon_region[1, ]$end)
              Exon_pos<-data.frame(start=Exon_region[1, ]$end-Exon_region[which(Exon_region$end>=(end(F1.pos) - 999)),]$end+1,
                                   end=Exon_region[1, ]$end-Exon_region[which(Exon_region$end>=(end(F1.pos) - 999)),]$start+1)
            }
            else{
              PCR_seq1 <-substring(Gene, Exon_region[1, ]$start, Exon_region[1, ]$end)
              Exon_pos<-data.frame(start=1,end=2000)
            }
          }
          
          {
            if (start(R1.pos) + 999 > Exon_region[nrow(Exon_region), ]$end) {
              PCR_seq2 <-substring(Gene, Exon_region[nrow(Exon_region), ]$start, start(R1.pos) +999)
              Exon_pos2<-data.frame(start=(start(R1.pos) + 999)-Exon_region[which(Exon_region$start<=(start(R1.pos) + 999)),]$end+1,
                                    end=(start(R1.pos) + 999)-Exon_region[which(Exon_region$start<=(start(R1.pos) + 999)),]$start+1)
            }
            else{
              PCR_seq2 <-substring(Gene, Exon_region[nrow(Exon_region), ]$start, Exon_region[nrow(Exon_region), ]$end)
              Exon_pos2<-data.frame(start=1,end=2000)
            }
          }
        }
        else{
          F1.pos <-matchPattern(DNAString(result1$primer1), subject = Gene)
          R1.pos <-matchPattern(reverseComplement(DNAString(result1$primer2)), subject = Gene)
          {
            if(start(F1.pos)+999>Exon_region[1,]$end){
              PCR_seq1 <- substring(Gene,Exon_region[1,]$start,start(F1.pos)+999)
              Exon_pos<-data.frame(start=Exon_region[which(Exon_region$start<=(start(F1.pos)+999)),]$start-Exon_region[1,]$start+1,
                                   end=Exon_region[which(Exon_region$start<=(start(F1.pos)+999)),]$end-Exon_region[1,]$start+1)
            }
            else{
              PCR_seq1 <- substring(Gene,Exon_region[1,]$start,Exon_region[1,]$end)
              Exon_pos<-data.frame(start=1,end=2000)
            }
          }
          {
            if(end(R1.pos)-999<Exon_region[nrow(Exon_region),]$start){
              PCR_seq2 <- substring(Gene,end(R1.pos)-999,Exon_region[nrow(Exon_region),]$end)
              Exon_pos2<-data.frame(start=Exon_region[which(Exon_region$end>=(end(R1.pos)-999)),]$start-(end(R1.pos)-999)+1,
                                    end=Exon_region[which(Exon_region$end>=(end(R1.pos)-999)),]$end-(end(R1.pos)-999)+1)
            }
            else{
              PCR_seq2 <- substring(Gene,Exon_region[nrow(Exon_region),]$start,Exon_region[nrow(Exon_region),]$end)
              Exon_pos2<-data.frame(start=1,end=2000)
            }
          }
        }
      }
      result_in<-primer_design2_in(Gene,PCR_seq1,PCR_seq2,Gene_rev,result1,Exon_pos,Exon_pos2)
    }
  }
  
  
  # 第二对gRNA -----------------------------------------------------------------
  if(length(result2)!=0){
    ko_Exon2<-PlanA.table$ko_condon2
    ko_Exon2<-strsplit(ko_Exon2,",")[[1]]
    
    #如果只敲一个外显子
    if(length(ko_Exon2)==1){
      Exon_region2<-t_Exon_CDS[which(t_Exon_CDS$Exon %in% ko_Exon2),]
      PCR_seq <- substring(Gene,min(gRNA2$end)+1,max(gRNA2$start)-1)
      if(PlanA.table$region2<30){
        result2_in<-character()
      }
      #一个内侧靶位点
      if (PlanA.table$region2 >=30 & PlanA.table$region2 <=200) {
        result2_in <- primer_design3_in(Gene, PCR_seq, Gene_rev, result2)
      }
      
      #两个内侧靶位点
      if (PlanA.table$region2 > 200) {
        result2_in <- primer_design_in(Gene, PCR_seq, Gene_rev, result2)
      }
    }
    
    #敲多个外显子
    if(length(ko_Exon2)>1){
      Exon_region2<-t_Exon_CDS[which(t_Exon_CDS$Exon %in% ko_Exon2),]
      #先找到外侧引物的位置
      {
        if (Gene_rev) {
          F1.pos <-matchPattern(reverse(DNAString(result2$primer1)), subject = Gene)
          R1.pos <-matchPattern(complement(DNAString(result2$primer2)), subject = Gene)
          {
            if (end(F1.pos) - 999 < Exon_region2[1, ]$start) {
              PCR_seq1 <- substring(Gene, end(F1.pos) - 999, Exon_region2[1, ]$end)
              Exon_pos<-data.frame(start=Exon_region2[1, ]$end-Exon_region2[which(Exon_region2$end>=(end(F1.pos) - 999)),]$end+1,
                                   end=Exon_region2[1, ]$end-Exon_region2[which(Exon_region2$end>=(end(F1.pos) - 999)),]$start+1)
            }
            else{
              PCR_seq1 <-substring(Gene, Exon_region2[1, ]$start, Exon_region2[1, ]$end)
              Exon_pos<-data.frame(start=1,end=2000)
            }
          }
          
          {
            if (start(R1.pos) + 999 > Exon_region2[nrow(Exon_region2), ]$end) {
              PCR_seq2 <-substring(Gene, Exon_region2[nrow(Exon_region2), ]$start, start(R1.pos) +999)
              Exon_pos2<-data.frame(start=(start(R1.pos) + 999)-Exon_region2[which(Exon_region2$start<=(start(R1.pos) + 999)),]$end+1,
                                    end=(start(R1.pos) + 999)-Exon_region2[which(Exon_region2$start<=(start(R1.pos) + 999)),]$start+1)
            }
            else{
              PCR_seq2 <-substring(Gene, Exon_region2[nrow(Exon_region2), ]$start, Exon_region2[nrow(Exon_region2), ]$end)
              Exon_pos2<-data.frame(start=1,end=2000)
            }
          }
        }
        else{
          F1.pos <-matchPattern(DNAString(result2$primer1), subject = Gene)
          R1.pos <-matchPattern(reverseComplement(DNAString(result2$primer2)), subject = Gene)
          {
            if(start(F1.pos)+999>Exon_region2[1,]$end){
              PCR_seq1 <- substring(Gene,Exon_region2[1,]$start,start(F1.pos)+999)
              Exon_pos<-data.frame(start=Exon_region2[which(Exon_region2$start<=(start(F1.pos)+999)),]$start-Exon_region2[1,]$start+1,
                                   end=Exon_region2[which(Exon_region2$start<=(start(F1.pos)+999)),]$end-Exon_region2[1,]$start+1)
            }
            else{
              PCR_seq1 <- substring(Gene,Exon_region2[1,]$start,Exon_region2[1,]$end)
              Exon_pos<-data.frame(start=1,end=2000)
            }
          }
          {
            if(end(R1.pos)-999<Exon_region2[nrow(Exon_region2),]$start){
              PCR_seq2 <- substring(Gene,end(R1.pos)-999,Exon_region2[nrow(Exon_region2),]$end)
              Exon_pos2<-data.frame(start=Exon_region2[which(Exon_region2$end>=(end(R1.pos)-999)),]$start-(end(R1.pos)-999)+1,
                                    end=Exon_region2[which(Exon_region2$end>=(end(R1.pos)-999)),]$end-(end(R1.pos)-999)+1)
            }
            else{
              PCR_seq2 <- substring(Gene,Exon_region2[nrow(Exon_region2),]$start,Exon_region2[nrow(Exon_region2),]$end)
              Exon_pos2<-data.frame(start=1,end=2000)
            }
          }
        }
      }
      result2_in<-primer_design2_in(Gene,PCR_seq1,PCR_seq2,Gene_rev,result2,Exon_pos,Exon_pos2)
    }
  }
  
  # 输出 ----------------------------------------------------------------------
  output<-data.frame(F1=character(),R1=character(),R2=character(),F1=character(),WT_F1R1=numeric(),WT_F1R2=numeric(),WT_F2R1=numeric(),GC1=character(),GC2=character(),GC3=character())
  {
    #两组外侧靶位点存在
    if (length(result1) != 0 & length(result2) != 0) {
      #两组都有内侧靶位点
      if (length(result_in) != 0 & length(result2_in) != 0) {

        #两个内侧靶位点
        {
          if (PlanA.table$region > 200) {
            output <- WT.calculation(result_in)
            p1 <- image1(output, Exon_region)
          }
          #一个内侧靶位点
          else{
            output <- result_in
            p1 <- image3(output, Exon_region)
          }
        }

        
        {
          if (PlanA.table$region2 > 200) {
            output2 <- WT.calculation(result2_in)
            p2 <- image1(output2,Exon_region)
          }
          else{
            output2 <- result2_in
            p2 <- image3(output2,Exon_region2)
          }
        }
      }
      
      #第二组没有内侧靶位点
      else if (length(result_in) != 0 & length(result2_in) == 0) {
        #两个内侧靶位点
        {
          if (PlanA.table$region > 200) {
            output <- WT.calculation(result_in)
            p1 <- image1(output, Exon_region)
          }
          #一个内侧靶位点
          else{
            output <- result_in
            p1 <- image3(output, Exon_region)
          }
        }
        
        p2<-image2(result2,Exon_region2)
        output2 <- data.frame(F1=result2$primer1,R1=result2$primer2,R2=NA,F2=NA,WT_F1R1=result2$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result2$GC,GC2=NA,GC3=NA)
      }
      #第一组没有内侧靶位点
      else if (length(result_in) == 0 & length(result2_in) != 0) {
        p1<-image2(result1,Exon_region)
        output <- data.frame(F1=result1$primer1,R1=result1$primer2,R2=NA,F2=NA,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result1$GC,GC2=NA,GC3=NA)
        {
          if (PlanA.table$region2 > 200) {
            output2 <- WT.calculation(result2_in)
            p2 <- image1(output2,Exon_region2)
          }
          else{
            output2 <- result2_in
            p2 <- image3(output2,Exon_region2)
          }
        }
      }
      #都没有内侧靶位点
      else{
        p1<-image2(result1,Exon_region)
        p2<-image2(result2,Exon_region2)
        output <- data.frame(F1=result1$primer1,R1=result1$primer2,R2=NA,F2=NA,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result1$GC,GC2=NA,GC3=NA)
        output2 <- data.frame(F1=result2$primer1,R1=result2$primer2,R2=NA,F2=NA,WT_F1R1=result2$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result2$GC,GC2=NA,GC3=NA)
      }
      
      png(file = paste0(filepath2,"//","primer_image1.png"),width = 480*3,height = 480*2,res = 72*2)
      print(p1)
      dev.off()
      png(file = paste0(filepath2,"//","primer_image2.png"),width = 480*3,height = 480*2,res = 72*2)
      print(p2)
      dev.off()
      write.csv(rbind(output,output2), paste0(filepath2, "//resultA.csv"), row.names = FALSE)
    }
    
    
    #第二组没有外侧靶位点
    else if(length(result1) != 0 & length(result2) == 0){
      if(length(result_in) != 0){
        #两个内侧靶位点
        if(PlanA.table$region > 200){
          output[1,] <- WT.calculation(result_in)
          p1 <- image1(output[1,],Exon_region)
        }
        #一个内侧靶位点
        else{
          output[1,]<-result_in
          p1 <- image3(output[1,],Exon_region)
        }
      }
      else{
        output[1,]$F1<-result1$primer1
        output[1,]$R1<-result1$primer2
        output[1,]$WT_F1R1<-result1$WT
        output[1,]$GC1<-result1$GC
        p1<-image2(result1,Exon_region)
      }
      png(file = paste0(filepath2,"//","primer_image1.png"),width = 480*3,height = 480*2,res = 72*2)
      print(p1)
      dev.off()
      write.csv(output, paste0(filepath2, "//resultA.csv"), row.names = FALSE)
    }
    
    #第一组没有外侧靶位点
    else if(length(result1) == 0 & length(result2) != 0){
      if(length(result2_in) != 0){
        #两个内侧靶位点
        if(PlanA.table$region2 > 200){
          output[2,] <- WT.calculation(result2_in)
          p2 <- image1(output[2,],Exon_region2)
        }
        #一个内侧靶位点
        else{
          output[2,]<-result2_in
          p2 <- image3(output[2,],Exon_region2)
        }
      }
      else{
        output[2,]$F1<-result2$primer1
        output[2,]$R1<-result2$primer2
        output[2,]$WT_F1R1<-result2$WT
        output[2,]$GC1<-result2$GC
        p2<-image2(result2,Exon_region2)
      }
      png(file = paste0(filepath2,"//","primer_image2.png"),width = 480*3,height = 480*2,res = 72*2)
      print(p2)
      dev.off()
      write.csv(output, paste0(filepath2, "//resultA.csv"), row.names = FALSE)
    }
    
    else{
      write.table("fail", paste0(filepath2, "//fail.txt"))
    }
  }
}