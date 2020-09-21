Get_result1 <- function(gRNA.table, KO_region) {
  for (i in 1:nrow(gRNA.table)) {
    j <- i + 1
    repeat {
      if (j > nrow(gRNA.table)) {
        break
      }
      #如果一个在上游，一个在下游，则符合
      if (gRNA.table[i,]$end <= KO_region$start &
          gRNA.table[j,]$start >= KO_region$end |
          gRNA.table[i,]$start >= KO_region$end &
          gRNA.table[j,]$end <= KO_region$start) {
        if (gRNA.table[i, ]$Score1 >= 70 & gRNA.table[j, ]$Score1 >= 70)
        {
          gRNA1 <- gRNA.table[i,]
          gRNA2 <- gRNA.table[j,]
          break
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
