# a<-which(gene.table$up_dis<0 & gene.table$up_dis>(-200) & gene.table$pos1 < abs(gene.table$up_dis) )
# 
# b<-which(gene.table$down_dis<0 & gene.table$down_dis>(-200) & (nchar(gene.table$seq)-gene.table$pos2) < abs(gene.table$down_dis))
# 
# #重叠大于200直接删除方案
# a
# b
# gene.table[b,]$gRNA1<-NA
# gene.table[b,]$gRNA2<-NA
# gene.table[b,]$strand1<-NA
# gene.table[b,]$strand2<-NA
# gene.table[b,]$pos1<-NA
# gene.table[b,]$pos2<-NA
# gene.table[b,]$score1<-NA
# gene.table[b,]$score2<-NA
# 
# 
# b<-which(is.na(gene.table2$gRNA1))
# c<-b[which((is.na(gene.table2[b,]$up_dis)|gene.table2[b,]$up_dis>(-20)) & (is.na(gene.table2[b,]$down_dis)|gene.table2[b,]$down_dis>(-20)) ) ]
# d<-c[which(nchar(gene.table2[c,]$seq)>=400)]
gene.table<-read.csv("C://Users//41518//Desktop//微生物//Sac-single.csv")

a<-which(gene.table$up_dis<0 &gene.table$down_dis>=0)
b<-which(gene.table$down_dis<0 & gene.table$up_dis>=0)
c<-which(gene.table$up_dis<0 &gene.table$down_dis<0 )

a_wrong<-a[which(gene.table[a,]$pos1<abs(gene.table[a,]$up_dis))]

b_wrong<-b[which(gene.table[b,]$pos1>nchar(gene.table[b,]$seq)+gene.table[b,]$down_dis)]
b_wrong2
# nchar(gene.table[b_wrong,]$seq)+gene.table[b_wrong,]$down_dis

seq<-substring(gene.table[i,]$seq,329+1,2216)
# 
# for(i in a_wrong){
  
  # if((nchar(gene.table[i,]$seq)-3)<=abs(gene.table[i,]$up_dis)+1){
  #   print(paste(i,"delete"))
  #   next
  # }
  
  # if(nchar(gene.table[i,]$seq)<=2000){
  #   seq<-substring(gene.table[i,]$seq,abs(gene.table[i,]$up_dis)+1,nchar(gene.table[i,]$seq)-3)
  # }
  # if(nchar(gene.table[i,]$seq)>2000){
  #   seq<-substring(gene.table[i,]$seq,abs(gene.table[i,]$up_dis)+1,abs(gene.table[i,]$up_dis)+2000)
  # }
  
  py$run(seq,species,filepath)
  gRNA.table<-read.csv(paste0(filepath,"//gRNA.csv"),header = FALSE)
  
  gRNA.table<-gRNA.table[which(!grepl("Inefficient", gRNA.table$V2)),]
  
  gRNA.table<-gRNA.table[which(as.numeric(gRNA.table$V3)>=70 & as.numeric(gRNA.table$V8)>=60),]
  if(nrow(gRNA.table)==0){
    print(paste(i,"fail"))
    next
  }
  
  gRNA.table$V2 <-gsub(" ", "", substring(gRNA.table$V2, 1, 24))
  strand <- sub("[^a-zA-Z]+", "", gRNA.table$V1)
  gRNA.table$V4<-strand
  pos <- as.numeric(sub("[^0-9]+", "", gRNA.table$V1))
  for(j in 1:length(strand)){
    if(strand[j]=="fw"){
      pos[j]=pos[j]+2
    }
    else if(strand[j]=="rev"){
      pos[j]=pos[j]+22
    }
  }
  gRNA.table$V1<-pos
  
  gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=90 & as.numeric(gRNA.table$V8)>=80 & gRNA.table$V10=="0-0-0-0-0"),]
  if(nrow(gRNA.table1)==0){
    gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70 & gRNA.table$V10=="0-0-0-0-0"),]
    if(nrow(gRNA.table1)==0){
      gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=80 & as.numeric(gRNA.table$V8)>=70),]
      if(nrow(gRNA.table1)==0){
        gRNA.table1<-gRNA.table[which(as.numeric(gRNA.table$V3)>=70 & as.numeric(gRNA.table$V8)>=60),]
      }
    }
  }
  
  if(!is.null(gRNA.table1)){
    if(nrow(gRNA.table1)!=0){
      gRNA<-data.frame(strand=gRNA.table1[1,]$V4,gRNA=gRNA.table1[1,]$V2,pos=gRNA.table1[1,]$V1+329,score=gRNA.table1[1,]$V3)
    }
  }
  
  {
    if(exists("gRNA")){
      gene.table[i, ]$gRNA1 <- gRNA$gRNA
      gene.table[i,]$strand1 <- gRNA$strand
      gene.table[i,]$pos1 <- gRNA$pos
      gene.table[i,]$score1 <- gRNA$score
      
      print(paste(i,"success"))
      # write.csv(gene.table,"C://Users//41518//Desktop//微生物//Sac-single.csv",row.names = FALSE)
      rm(gRNA)
    }
    else{
      print(paste(i,"fail"))
    }
  }
  
  
  
# # }
gene.table[i, ]$gRNA1 <- NA
gene.table[i, ]$strand1 <- NA
gene.table[i, ]$pos1 <- NA
gene.table[i, ]$score1 <- NA
