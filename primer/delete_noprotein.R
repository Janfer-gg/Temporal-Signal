#排除掉不编码蛋白和不完整的转录本
delete_noprotein <- function(transcript.table) {
  rmi <- numeric()
  for (i in 1:nrow(transcript.table)) {
    if (!grepl("Protein coding", transcript.table[i, ]$Biotype)) {
      rmi <- append(rmi, i)
    }
    else if (grepl("incomplete", transcript.table[i, ]$Flags)) {
      rmi <- append(rmi, i)
    }
  }
  if (length(rmi) != 0) {
    transcript.table <- transcript.table[-rmi, ]
  }
  return(transcript.table)
}
