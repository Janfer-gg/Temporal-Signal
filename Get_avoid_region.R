Get_avoid_region <- function(othergene.table) {
  u <- 1
  for (i in 1:nrow(othergene.table)) {
    #如果是编码蛋白的基因，要避免重叠区域
    if (othergene.table[i, ]$biotype == "Protein coding") {
      avoid_region[u, 1] <- as.numeric(othergene.table[i, ]$start) - start + 1
      avoid_region[u, 2] <-
        as.numeric(othergene.table[i, ]$end) - start + 1
      u <- u + 1
    }
  }
  return(avoid_region)
}
