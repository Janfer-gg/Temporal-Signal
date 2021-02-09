setwd("C://Users//41518//Desktop//work/ubigene/primer")
filepath<-"C://Users//41518//Desktop//0922测试//BMPR1A"
# library(ggplot2)
# library(ggpubr)
library(httr)
library(jsonlite)
library(xml2)
library(Biostrings)
library(reticulate)
library(stringr)
source("Get_ID.R")
source("Get_allinfo.R")
source("Get_seq.R")
source("delete_noprotein.R")
source("Get_max_transcript.R")

source_python("crispor_table_download2.py")

PlanA.table<-read.csv(paste0(filepath,"//planA.csv"))
term<-PlanA.table$Gene

species<-"Human"

ID <- Get_ID(term,species)
Gene <- Get_seq(ID)                    #基因序列(左右两端各加1500bp)
allinfo <- Get_allinfo(ID)
start <- allinfo$start

# source_python("ensembl_table_download.py")
# py$download_csv(filepath,ID)
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
 
gRNA<-data.frame()
# gRNA2<-data.frame()
# gRNA4<-data.frame()
gRNA[1,1]<-PlanA.table$gRNA1
gRNA[2,1]<-PlanA.table$gRNA2
# gRNA2[1,1]<-"GTTTAAAGTGTACTGTCTGCAGG"
# gRNA2[2,1]<-"CATAATCGTGCCCAGAATACTGG"
# gRNA4[1,1]<-"ACTCCCGGATGCTTGACCGTAGG"
# gRNA4[2,1]<-"CCGGCATTACCTTATGAAAGTGG"


# gRNA[1,1]<-tab1$g1
# gRNA[2,1]<-tab1$g2
# gRNA2[1,1]<-tab1$g3
# gRNA2[2,1]<-tab1$g4

{
  if (Gene_rev) {
    for (i in 1:2) {
      gRNA1 <- DNAString(gRNA[i, 1])
      pos <- matchPattern(pattern = gRNA1, subject = gene)
      if (length(pos) == 0) {
        gRNA3 <- complement(gRNA1)
        pos <- matchPattern(pattern = gRNA3, subject = gene)
      }
      if (length(pos) == 0) {
        gRNA3 <- reverse(gRNA1)
        pos <- matchPattern(pattern = gRNA3, subject = gene)
      }
      pos_start <- start(pos)
      pos_end <- end(pos)
      gRNA[i, 2] <- pos_start
      gRNA[i, 3] <- pos_end
    }
  }
  else{
    for (i in 1:2) {
      gRNA1 <- DNAString(gRNA[i, 1])
      pos <- matchPattern(pattern = gRNA1, subject = gene)
      if (length(pos) == 0) {
        gRNA3 <- reverseComplement(gRNA1)
        pos <- matchPattern(pattern = gRNA3, subject = gene)
      }
      pos_start <- start(pos)
      pos_end <- end(pos)
      gRNA[i, 2] <- pos_start
      gRNA[i, 3] <- pos_end
    }
  }
}
# 
# {
#   if (Gene_rev) {
#     for (i in 1:2) {
#       gRNA1 <- DNAString(gRNA2[i, 1])
#       pos <- matchPattern(pattern = gRNA1, subject = gene)
#       if (length(pos) == 0) {
#         gRNA3 <- complement(gRNA1)
#         pos <- matchPattern(pattern = gRNA3, subject = gene)
#       }
#       if (length(pos) == 0) {
#         gRNA3 <- reverse(gRNA1)
#         pos <- matchPattern(pattern = gRNA3, subject = gene)
#       }
#       pos_start <- start(pos)
#       pos_end <- end(pos)
#       gRNA2[i, 2] <- pos_start
#       gRNA2[i, 3] <- pos_end
#     }
#   }
#   else{
#     for (i in 1:2) {
#       gRNA1 <- DNAString(gRNA2[i, 1])
#       pos <- matchPattern(pattern = gRNA1, subject = gene)
#       if (length(pos) == 0) {
#         gRNA3 <- reverseComplement(gRNA1)
#         pos <- matchPattern(pattern = gRNA3, subject = gene)
#       }
#       pos_start <- start(pos)
#       pos_end <- end(pos)
#       gRNA2[i, 2] <- pos_start
#       gRNA2[i, 3] <- pos_end
#     }
#   }
# }

# {
#   if (Gene_rev) {
#     for (i in 1:2) {
#       gRNA1 <- DNAString(gRNA4[i, 1])
#       pos <- matchPattern(pattern = gRNA1, subject = gene)
#       if (length(pos) == 0) {
#         gRNA3 <- complement(gRNA1)
#         pos <- matchPattern(pattern = gRNA3, subject = gene)
#       }
#       if (length(pos) == 0) {
#         gRNA3 <- reverse(gRNA1)
#         pos <- matchPattern(pattern = gRNA3, subject = gene)
#       }
#       pos_start <- start(pos)
#       pos_end <- end(pos)
#       gRNA4[i, 2] <- pos_start
#       gRNA4[i, 3] <- pos_end
#     }
#   }
#   else{
#     for (i in 1:2) {
#       gRNA1 <- DNAString(gRNA4[i, 1])
#       pos <- matchPattern(pattern = gRNA1, subject = gene)
#       if (length(pos) == 0) {
#         gRNA3 <- reverseComplement(gRNA1)
#         pos <- matchPattern(pattern = gRNA3, subject = gene)
#       }
#       pos_start <- start(pos)
#       pos_end <- end(pos)
#       gRNA4[i, 2] <- pos_start
#       gRNA4[i, 3] <- pos_end
#     }
#   }
# }


names(gRNA)[2:3] <- c("start", "end")
# names(gRNA2)[2:3] <- c("start", "end")
# names(gRNA4)[2:3] <- c("start", "end")
print(gRNA)
# print(gRNA2)
# print(gRNA4)

#rm(gRNA)
# rm(gRNA2)
# rm(gRNA4)



