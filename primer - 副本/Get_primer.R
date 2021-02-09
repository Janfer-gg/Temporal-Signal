Get_primer<-function(PCR_seq,PCR.table){
  #匹配PCR的位置，然后前后扩大3个碱基
  pcr_seq <- DNAString(PCR_seq)
  PCR.table_fw<-data.frame(PCR=character(),GC=numeric(),TM=numeric())
  m<-1
  for(i in 1:nrow(PCR.table)){
    if(sub("[^a-zA-Z]+", "", PCR.table[i, 1]) == "fw"){
      fw <- matchPattern(pattern = sub(" ","",substring(PCR.table[i,2],1,24)),subject = pcr_seq)
      if(length(fw)>1){
        next
      }
    }
    else if(sub("[^a-zA-Z]+", "", PCR.table[i, 1]) == "rev"){
      fw <- matchPattern(pattern =reverseComplement(DNAString(sub(" ","",substring(PCR.table[i,2],1,24)))),subject = pcr_seq)
      if(length(fw)>1){
        next
      }
    }
    fw_start <- start(fw)
    fw_end <- end(fw)
    PCR_pre<-substring(PCR_seq,fw_start-3,fw_end+3)
    for(j in 19:28){
      for(k in 1:(nchar(PCR_pre)-(j-1))){
        PCR_pre2<-substring(PCR_pre,k,k+(j-1))
        if(j>nchar(PCR_pre2)){
          next
        }
        #如果以G或C结尾，则计算GC值和TM值
        if(strsplit(PCR_pre2,"")[[1]][j]=="G"|strsplit(PCR_pre2,"")[[1]][j]=="C"){
          cal<-DNAString(PCR_pre2)
          GC_count<-sum(letterFrequency(cal,c("G","C")))
          GC<-round(GC_count/j,2)
          #如果GC大于0.35和小于0.65，计算TM值
          if(GC >= 0.35 & GC <= 0.65){
            TM<-100.5 + (41*(GC_count)/j) - (820/j) + 16.6*log(0.05,10)
            #如果TM值也符合，则该PCR引物符合
            if(TM>=61&TM<=64){
              PCR.table_fw[m,1]<-PCR_pre2
              PCR.table_fw[m,2]<-GC
              PCR.table_fw[m,3]<-TM
              m<-m+1
            }
            else{
              next
            }
          }
          else{
            next
          }
        }
        else{
          next
        }
      }
    }
  }
  return(PCR.table_fw)
}


Get_primer2<-function(PCR_seq,PCR.table){
  #匹配PCR的位置，然后前后扩大3个碱基
  pcr_seq <- DNAString(PCR_seq)
  PCR.table_fw<-data.frame(PCR=character(),GC=numeric(),TM=numeric())
  m<-1
  for(i in 1:nrow(PCR.table)){
    if(sub("[^a-zA-Z]+", "", PCR.table[i, 1]) == "fw"){
      fw <- matchPattern(pattern = sub(" ","",substring(PCR.table[i,2],1,24)),subject = pcr_seq)
      if(length(fw)>1){
        next
      }
    }
    else if(sub("[^a-zA-Z]+", "", PCR.table[i, 1]) == "rev"){
      fw <- matchPattern(pattern =reverseComplement(DNAString(sub(" ","",substring(PCR.table[i,2],1,24)))),subject = pcr_seq)
      if(length(fw)>1){
        next
      }
    }
    fw_start <- start(fw)
    fw_end <- end(fw)
    PCR_pre<-substring(PCR_seq,fw_start-3,fw_end+3)
    for(j in 19:28){
      for(k in 1:(nchar(PCR_pre)-(j-1))){
        PCR_pre2<-substring(PCR_pre,k,k+(j-1))
        if(j>nchar(PCR_pre2)){
          next
        }
        #如果以G或C结尾，则计算GC值和TM值
        if(strsplit(PCR_pre2,"")[[1]][1]=="G"|strsplit(PCR_pre2,"")[[1]][1]=="C"){
          cal<-DNAString(PCR_pre2)
          GC_count<-sum(letterFrequency(cal,c("G","C")))
          GC<-round(GC_count/j,2)
          #如果GC大于0.35和小于0.65，计算TM值
          if(GC >= 0.35 & GC <= 0.65){
            TM<-100.5 + (41*(GC_count)/j) - (820/j) + 16.6*log(0.05,10)
            #如果TM值也符合，则该PCR引物符合
            if(TM>=61&TM<=64){
              PCR.table_fw[m,1]<-PCR_pre2
              PCR.table_fw[m,2]<-GC
              PCR.table_fw[m,3]<-TM
              m<-m+1
            }
            else{
              next
            }
          }
          else{
            next
          }
        }
        else{
          next
        }
      }
    }
  }
  return(PCR.table_fw)
}
