table_select<-function(table){
  table<-table[which(!grepl("Inefficient", table$V2)),]
  
  table<-table[which(as.numeric(table$V3)==100 & as.numeric(table$V5)>=40),]
  table<-table[order(as.numeric(table$V8),decreasing = TRUE),]
  table<-rbind(table[which(table$V9=="0-0-0-0-0" & table$V10=="0-0-0-0-0"),],table[which((substring(table$V9,1,7)=="0-0-0-0" & table$V9!="0-0-0-0-0" )| (substring(table$V10,1,7)!="0-0-0-0" & table$V10!="0-0-0-0-0")),])
  if(nrow(table)==0){
    return(NULL)
  }
  table$V2 <-gsub(" ", "", substring(table$V2, 1, 24))
  strand <- sub("[^a-zA-Z]+", "", table$V1)
  table$V4<-strand
  pos <- as.numeric(sub("[^0-9]+", "", table$V1))
  for(i in 1:length(strand)){
    if(strand[i]=="fw"){
      pos[i]=pos[i]+2
    }
    else if(strand[i]=="rev"){
      pos[i]=pos[i]+22
    }
  }
  table$V1<-pos
  
  table<-rbind(table[grep("AAAA|TTTT",table$V2,invert = TRUE),],table[grep("AAAA|TTTT",table$V2),])
  
  return(table)
}

