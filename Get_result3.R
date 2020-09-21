Get_result3 <- function(gRNA.table) {
  for (i in 1:nrow(gRNA.table)) {
    j <- i + 1
    repeat {
      if (j > nrow(gRNA.table)) {
        break
      }
      gRNA1 <- gRNA.table[i, ]
      gRNA2 <- gRNA.table[j, ]
      if (Gene_rev) {
        if (gRNA1$strand == "rev" & gRNA2$strand == "rev" |
            gRNA1$strand == "fw" & gRNA2$strand == "fw") {
          KO_length <- abs(gRNA1$end - gRNA2$end) 
        }
        else{
          if (gRNA1$strand == "rev" & gRNA2$strand == "fw") {
            if(gRNA1$start>gRNA2$start){
              KO_length <- abs(gRNA1$end - gRNA2$start) +1
            }
            else{
              KO_length <- abs(gRNA1$end - gRNA2$start) -1
            }
          }
          else{
            if(gRNA1$start>gRNA2$start){
              KO_length <- abs(gRNA1$start - gRNA2$end) -1
            }
            else{
              KO_length <- abs(gRNA1$start - gRNA2$end) +1
            }
          }
        }
      }
      #正向
      else{
        if (gRNA1$strand == "rev" & gRNA2$strand == "rev" |
            gRNA1$strand == "fw" & gRNA2$strand == "fw") {
          KO_length <- abs(gRNA1$end - gRNA2$end) 
        }
        else{
          if (gRNA1$strand == "rev" & gRNA2$strand == "fw") {
            pos1 <- gRNA1$start
            pos2 <- gRNA2$end
            KO_length <- abs(pos1 - pos2-1) 
          }
          else{
            pos1 <- gRNA1$end
            pos2 <- gRNA2$start
            KO_length <- abs(pos1 - pos2+1) 
          }
        }
      }
      #如果相距100bp以上且切割大小非3的倍数,且都大于70分,则符合
      if (KO_length>=100 & KO_length %% 3!=0 & gRNA1$Score1>=70 & gRNA2$Score1>=70) {
        {
          if (Gene_rev == "FALSE") {
            if (gRNA1$end < t_CDS_30 | gRNA2$end < t_CDS_30) {
              return(rbind(gRNA1, gRNA2))
              break
            }
          }
          else{
            if (gRNA1$start > t_CDS_30 | gRNA2$start > t_CDS_30) {
              return(rbind(gRNA1, gRNA2))
              break
            }
          }
        }
      }
      j <- j + 1
    }
  }
}
