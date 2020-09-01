Get_max_transcript <- function(allinfo, transcript_table) {
  name <-
    transcript_table[which(transcript_table$bp == max(transcript_table$bp)),]$Name     #最长的转录本
  Exon_region <- data.frame()               #最长转录本的外显子位置
  for (j in 1:length(allinfo$Transcript)) {
    if (allinfo$Transcript[[j]]$display_name == name) {
      #起始密码子位点
      CDS_start <-
        allinfo$Transcript[[j]]$Translation$start - start + 1
      #前30%CDS的位点
      CDS_30 <-
        floor((allinfo[["Transcript"]][[j]][["Translation"]][["end"]] - allinfo[["Transcript"]][[j]][["Translation"]][["start"]]) *
                0.3) + CDS_start
      #终止位点
      CDS_end <- allinfo$Transcript[[j]]$Translation$end - start + 1
      #外显子位置
      for (i in 1:length(allinfo$Transcript[[j]]$Exon)) {
        Exon_start <- allinfo$Transcript[[j]]$Exon[[i]]$start - start + 1
        Exon_end <- allinfo$Transcript[[j]]$Exon[[i]]$end - start + 1
        Exon_name <- paste("Exon", i)
        Exon_region <-
          rbind(Exon_region, cbind(Exon_name, Exon_start, Exon_end))
      }
    }
  }
  return(Exon_region)
}
