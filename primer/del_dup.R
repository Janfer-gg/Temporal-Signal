#去掉连续5个碱基重复
del_dup <- function(PCR.table) {
  del<-numeric()
  for(i in 1:nrow(PCR.table)){
    a<-substring(PCR.table[i,2],1,20)
    x<-strsplit(a,"")[[1]]
    times<-0
    for (j in 2:length(x)) {
      if (x[j] == x[j - 1]) {
        times <- times + 1
        if (times == 4) {
          del<-append(del,i)
          break
        }
      }
      else{
        times <- 0
      }
    }
  }
  {
    if (length(del) != 0) {
      del_result <- PCR.table[-del,]
    }
    else{
      del_result <- PCR.table
    }
  }
}


del_dup2 <- function(PCR.table) {
  del<-numeric()
  for(i in 1:nrow(PCR.table)){
    x<-strsplit(PCR.table[i,]$PCR,"")[[1]]
    times<-0
    for (j in 2:length(x)) {
      if (x[j] == x[j - 1]) {
        times <- times + 1
        if (times == 4) {
          del<-append(del,i)
          break
        }
      }
      else{
        times <- 0
      }
    }
  }
  if (length(del) != 0) {
    PCR.table <- PCR.table[-del, ]
  }
  return(PCR.table)
}
