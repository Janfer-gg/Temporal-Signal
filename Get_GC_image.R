Get_GC_image <- function(KO_region) {
  if (KO_region$start < 800 | KO_region$end + 800 > nchar(Gene)) {
    analysis_seq <-
      substring(Gene2, KO_region$start + 500 - 800, KO_region$end + 500 + 800)
  }
  else{
    analysis_seq <-
      substring(Gene, KO_region$start - 800, KO_region$end + 800)
  }
  analysis_seq <- DNAString(analysis_seq)
  GC <-
    rowSums(letterFrequencyInSlidingView(analysis_seq, 30, c("G", "C"))) / 30 
    GC.data <- data.frame(x = c(1:(nchar(analysis_seq) - 29)), y = GC)
  p4 <-
    ggplot(data = GC.data, aes(x = x, y = y)) + geom_line(color = "red") + theme_bw() +
    theme(panel.grid = element_blank())
  return(p4)
}
