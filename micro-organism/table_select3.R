table_select2<-function(table){
  table<-table[which(!grepl("Inefficient", table$V2)),]
  
  table<-table[which(as.numeric(table$V3)>70 & as.numeric(table$V8)>60),]
  if(nrow(table)==0){
    return(NULL)
  }
  
  table$V2 <-gsub(" ", "", substring(table$V2, 1, 24))
  strand <- sub("[^a-zA-Z]+", "", table$V1)
  table$V4<-strand
  pos <- as.numeric(sub("[^0-9]+", "", table$V1))
  for(j in 1:length(strand)){
    if(strand[j]=="fw"){
      pos[j]=pos[j]+2
    }
    else if(strand[j]=="rev"){
      pos[j]=pos[j]+22
    }
  }
  table$V1<-pos
  
  table1<-table[which(as.numeric(table$V1)>100 & as.numeric(table$V3)>=80 & as.numeric(table$V8)>=70 & table$V10=="0-0-0-0-0"),]
  if(nrow(table1)==0){
    table1<-table[which(as.numeric(table$V1)>=200 & as.numeric(table$V3)>=80 & as.numeric(table$V8)>=70 & table$V10=="0-0-0-0-0"),]
    if(nrow(table1)==0){
      table1<-table[which(as.numeric(table$V1)>100 & as.numeric(table$V3)>70 & as.numeric(table$V8)>60),]
      if(nrow(table1)==0){
        table1<-table[which(as.numeric(table$V1)>=200 & as.numeric(table$V3)>70 & as.numeric(table$V8)>60),]
      }
    }
  }
  
  return(table1)
}
