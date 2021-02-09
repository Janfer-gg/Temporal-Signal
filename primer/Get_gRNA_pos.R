Get_gRNA_pos<-function(gRNA,gene,Gene_rev){
  {
    if (Gene_rev) {
      for (i in 1:nrow(gRNA)) {
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
      for (i in 1:nrow(gRNA)) {
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
  names(gRNA)[2:3] <- c("start", "end")
  return(gRNA)
}
