primer_tool<-function(term,species,gRNA1,gRNA2,filepath){
  setwd("C://Users//41518//Desktop//work/ubigene/primer")
  dir.create(filepath)
  library(httr)
  library(jsonlite)
  library(xml2)
  library(Biostrings)
  library(reticulate)
  library(stringr)
  library(ggplot2)
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
  source("inside_primer.R")
  source_python("primertool2.py") 
  source_python("crispor_table_download2.py")
  source_python("Oligo.py")
 
  gRNA<-data.frame()
  gRNA[1,1]<-gRNA1
  if(nchar(gRNA2)!=0){
    gRNA[2,1]<-gRNA2
  }
  term<-term
  species<-species

  ID <- Get_ID(term,species)
  Gene <- Get_seq(ID)
  allinfo <- Get_allinfo(ID)
  start <- allinfo$start-1500

  source_python("ensembl_table_download.py")
  py$download_csv(filepath,ID)
  transcript.table <- read.csv(paste0(filepath,"//transcript.csv"),header = TRUE)
  transcript.table <-delete_noprotein(transcript.table)
  transcript.name <-
    transcript.table[which(transcript.table$bp == max(transcript.table$bp)), ]$Name
  t_Exon_region <-
    Get_max_transcript(allinfo, transcript.table,start)
  t_Exon_region$Exon_start <- as.numeric(t_Exon_region$Exon_start)
  t_Exon_region$Exon_end <- as.numeric(t_Exon_region$Exon_end)
  t_Exon_region_sort <-
    t_Exon_region[order(t_Exon_region$Exon_start), ]



  Gene_rev <-
    all(t_Exon_region$Exon_start == sort(t_Exon_region$Exon_start, decreasing = TRUE))

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
  }



  ko.data <- data.frame()
  name <- paste(transcript.table$Name, collapse = "")
  for (j in 1:length(allinfo$Transcript)) {
    if (grepl(allinfo$Transcript[[j]]$display_name, name)) {
      Transcript <- allinfo$Transcript[[j]]
      if (length(Transcript$Translation) != 0) {

        CDS_start <-
          Transcript$Translation$start - start + 1
        CDS_end <-
          Transcript$Translation$end - start + 1
        Transcript_start <-
          Transcript$start - start + 1
        Transcript_end <-
          Transcript$end - start + 1
        Exon_info_CDS <- data.frame()
        k <- 1
        for (i in 1:length(Transcript$Exon)) {
          Exon_start <- Transcript$Exon[[i]]$start - start + 1
          Exon_end <- Transcript$Exon[[i]]$end - start + 1

          if (Exon_start >= CDS_start & Exon_end <= CDS_end) {
            Exon_info_CDS[k, 1] <- paste("Exon", i)
            Exon_info_CDS[k, 2] <- Exon_start
            Exon_info_CDS[k, 3] <- Exon_end
            k <- k + 1
          }

          else if (Exon_start <= CDS_start &
                   Exon_end >= CDS_start & Exon_end <= CDS_end) {
            Exon_info_CDS[k, 1] <- paste("Exon", i)
            Exon_info_CDS[k, 2] <- CDS_start
            Exon_info_CDS[k, 3] <- Exon_end
            k <- k + 1
          }

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
  gRNA<-Get_gRNA_pos(gRNA,gene,Gene_rev)
  size1<-Get_size(max(gRNA$start)-min(gRNA$end))
  result1 <- primer_design1(Gene, gRNA,Gene_rev,size1,filepath)
  write.csv(result1,paste0(filepath,"//resultA.csv"),row.names = FALSE)
  
  if(nchar(gRNA2)==0){
    if(class(result1)!="NULL"){
      names(result1)[1:3]<-c("F1","R1","WT_F1R1")
      write.csv(result1,paste0(filepath,"//result.csv"),row.names = FALSE)
    }
    else{
      write.table("fail",paste0(filepath, "//fail.txt"),row.names = FALSE,col.names = FALSE)
    }
  }
  else if(nchar(gRNA2)!=0){
    {
      if(class(result1)!="NULL") {
        
        Exon_region <-t_Exon_CDS[which(t_Exon_CDS$start <= min(gRNA$end) &t_Exon_CDS$end >= max(gRNA$start)), ]
        {
          
          if (nrow(Exon_region) != 0) {
            PCR_seq <- substring(Gene, min(gRNA$end) + 1, max(gRNA$start) - 1)
            {
              if (nchar(PCR_seq) < 30) {
                result_in <- NULL
              }
              else if (nchar(PCR_seq) >= 30 & nchar(PCR_seq) <= 200) {
                result_in <- primer_design3_in(Gene, PCR_seq, Gene_rev, result1,filepath)
              }
              
              else if (nchar(PCR_seq) > 200) {
                result_in <- primer_design_in(Gene, PCR_seq, Gene_rev, result1,filepath)
              }
            }
            
            {
              if(length(result_in)!=0){
                
                if(nchar(PCR_seq) > 200){
                  output <- WT.calculation(result_in,Gene,Gene_rev)
                  p1 <- image1(output,Exon_region,Gene,Gene_rev)
                }
                
                else{
                  output<-result_in
                  p1 <- image3(output,Exon_region,Gene,Gene_rev)
                }
              }
              
              else{
                output <- data.frame(F1=result1$primer1,R1=result1$primer2,R2=NA,F2=NA,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result1$GC,GC2=NA,GC3=NA)
                p1<-image2(result1,Exon_region,Gene,Gene_rev)
              }
            }
            write.csv(output,paste0(filepath,"//result.csv"),row.names = FALSE)
            png(file = paste0(filepath,"//","primer_image.png"),width = 480*3,height = 480*2,res = 72*2)
            print(p1)
            dev.off()
          }
          
          
          else{
            Exon_region <-t_Exon_CDS[which(t_Exon_CDS$start > min(gRNA$start) &t_Exon_CDS$end < max(gRNA$end)),]
            result_in <-inside_primer3(Gene, Gene_rev, result1, Exon_region,filepath)
            
            
            {
              if(length(result_in)!=0){
                output <- WT.calculation(result_in,Gene,Gene_rev)
                p1 <- image1(output,Exon_region,Gene,Gene_rev)
                
              }
              else{
                output <- data.frame(F1=result1$primer1,R1=result1$primer2,R2=NA,F2=NA,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=NA,GC1=result1$GC,GC2=NA,GC3=NA)
                p1<-image2(result1,Exon_region,Gene,Gene_rev)
              }
            }
            write.csv(output,paste0(filepath,"//result.csv"),row.names = FALSE)
            png(file = paste0(filepath,"//","primer_image.png"),width = 480*3,height = 480*2,res = 72*2)
            print(p1)
            dev.off()
          }
        }
      }
      else{
        write.table("fail",paste0(filepath, "//fail.txt"),row.names = FALSE,col.names = FALSE)
      }
    }
    
  }
}

# 
term<-"SBF1"
species<-"Human"
gRNA1<-"ACAGAGAGGCTGAACGACACAGG"
gRNA2<-"GGGGGGCCGAGAGCCTAAGTGGG"
filepath<-"C://Users//41518//Desktop//SBF1"
primer_tool(term,species,gRNA1,gRNA2,filepath)

