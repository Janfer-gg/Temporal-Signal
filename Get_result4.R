Get_result4 <- function(gRNA,gRNA.table, KO_region) {
  for (i in 1:nrow(gRNA.table)) {
    j <- i + 1
    repeat {
      if (j > nrow(gRNA.table)) {
        break
      }
      #如果相差0.05以内
      if (gRNA.table[j, ]$crispr_score >= gRNA.table[i, ]$crispr_score - 0.05 & gRNA.table[j, ]$crispr_score<=gRNA.table[i, ]$crispr_score + 0.05) {
        #如果一个在上游，一个在下游，则符合
        if (gRNA.table[i, ]$end <= KO_region$start &
            gRNA.table[j, ]$start >= KO_region$end |
            gRNA.table[i, ]$start >= KO_region$end &
            gRNA.table[j, ]$end <= KO_region$start) {
          t1<-gRNA.table[i, ]$end %in% c((gRNA[1,]$start+3):(gRNA[1,]$end-3)) | gRNA.table[i, ]$start %in% c((gRNA[1,]$start+3):(gRNA[1,]$end-3))
          t2<-gRNA.table[i, ]$end %in% c((gRNA[2,]$start+3):(gRNA[2,]$end-3)) | gRNA.table[i, ]$start %in% c((gRNA[2,]$start+3):(gRNA[2,]$end-3))
          t3<-gRNA.table[j, ]$end %in% c((gRNA[1,]$start+3):(gRNA[1,]$end-3)) | gRNA.table[j, ]$start %in% c((gRNA[1,]$start+3):(gRNA[1,]$end-3))
          t4<-gRNA.table[j, ]$end %in% c((gRNA[2,]$start+3):(gRNA[2,]$end-3)) | gRNA.table[j, ]$start %in% c((gRNA[2,]$start+3):(gRNA[2,]$end-3))
          if(t1==FALSE & t2==FALSE & t3==FALSE & t4==FALSE){
            gRNA1 <- gRNA.table[i, ]
            gRNA2 <- gRNA.table[j, ]
            break
          }
        }
      }
      j <- j + 1
    }
    if (exists("gRNA1") == TRUE) {
      break
    }
  }
  if (exists("gRNA1") == TRUE) {
    return(rbind(gRNA1, gRNA2))
  }
}
