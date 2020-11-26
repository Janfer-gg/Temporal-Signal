library(Biostrings)
setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")

gene.table<-read.csv("C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv")


empty<-numeric()
wrong<-numeric()

for(i in 1:nrow(gene.table)){
  if(is.na(gene.table[i,]$gRNA1)){
    empty<-append(empty,i)
  }
  else if(nchar(gene.table[i,]$gRNA1)!=0){
    if(gene.table[i,]$strand1=="fw"){
      pos=matchPattern(gene.table[i,]$gRNA1,gene.table[i,]$seq)
    }
    if(gene.table[i,]$strand1=="rev"){
      pos=matchPattern(reverseComplement(DNAString(gene.table[i,]$gRNA1)),gene.table[i,]$seq)
    }
    
    if(end(pos)!=gene.table[i,]$pos1){
      wrong<-append(wrong,i)
    }
  }
  print(i)
}

# write.csv(gene.table,"C://Users//41518//Desktop//微生物//micro.csv",row.names = FALSE)
write.csv(gene.table_emp,"C://Users//41518//Desktop//微生物//micro-无.csv",row.names=FALSE)
