Get_mark_region <- function(othergene.table) {
  mark_region <- data.frame(start = numeric(), end = numeric())
  u <- 1
  for (i in 1:nrow(othergene.table)) {
    #如果是lncRNA，mircoRNA，circRNA基因,要加备注
    if (othergene.table[i,]$biotype=="lncRNA" | othergene.table[i,]$biotype=="mircoRNA" | othergene.table[i,]$biotype=="circRNA") {
      mark_region[u, 1] <- as.numeric(othergene.table[i, ]$start) - start + 1
      mark_region[u, 2] <-
        as.numeric(othergene.table[i, ]$end) - start + 1
      u <- u + 1
    }
  }
  return(mark_region)
}
