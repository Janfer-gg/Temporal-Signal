#找两个内侧引物
primer_design_in <- function(Gene, PCR_seq, Gene_rev,result1) {
  if (Gene_rev) {
    seq1 <- DNAString(PCR_seq)
    seq1_rev <- reverse(seq1)
    PCR_seq <- as.character(seq1_rev)
  }
  py$run(PCR_seq, species, filepath2, "3")
  PCR.tab <- read.csv(paste0(filepath2,"//PCR3.csv"), header = FALSE,encoding = "UTF-8")
  #60分以上
  PCR.tab<-PCR.tab[which(PCR.tab$V3!="No matches"),]
  PCR.tab<-PCR.tab[which(as.numeric(PCR.tab$V3)>=60),]
  if(nrow(PCR.tab)==0){
    return(NULL)
  }
  
  {
    # 反向基因 --------------------------------------------------------------------
    if (Gene_rev) {
      #删除连续5个碱基重复的序列
      PCR.table1 <- del_dup(PCR.tab)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_fw<-Get_primer(PCR_seq,PCR.table1)
      PCR.table_rev<-Get_primer2(PCR_seq,PCR.table1)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续4个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }

      #去除重复
      PCR.table_fw <-PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
      #取30条序列
      if (nrow(PCR.table_fw) > 30) {
        PCR.table_fw <- PCR.table_fw[1:30, ]
      }
      write.csv(PCR.table_fw,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
      source_python("OligoCalc.py")
      py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
      
      PCR1.table<-read.csv(paste0(filepath2,"//PCR3.csv"),header = TRUE)
      primer_up<-PCR1.table$PCR
      
      
      
      #PCR1的反向互补序列
      for (i in 1:nrow(PCR1.table)) {
        pcr1 <- DNAString(PCR1.table[i, ]$PCR)
        pcr1_revcom <- reverseComplement(pcr1)
        PCR1.table[i, ]$PCR <- as.character(pcr1_revcom)
      }
      primer_down<-PCR1.table$PCR
      
      # Blast -------------------------------------------------------------------
      source_python("primertool2.py")       
      primer.table<-py$run_primertool(as.list(result1$primer1),primer_down,species,4000)
      
      primer.table2<-py$run_primertool(primer_up,as.list(result1$primer2),species,4000)
      
      primer.table <-rbind(primer.table,primer.table2)
      
    }
    
    # 正向基因 --------------------------------------------------------------------
    else{
      PCR.table1 <- del_dup(PCR.tab)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_fw<-Get_primer(PCR_seq,PCR.table1)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_fw <-PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
      #取30条
      if (nrow(PCR.table_fw) > 30) {
        PCR.table_fw <- PCR.table_fw[1:30, ]
      }
      write.csv(PCR.table_fw,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
      source_python("OligoCalc.py")
      py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
      
      primer_up<-read.csv(paste0(filepath2,"//PCR3.csv"), header = TRUE)$PCR
      primer_down<-character()
      for(i in 1:length(primer_up)){
        primer_down[i]<-as.character(reverseComplement(DNAString(primer_up[i])))
      }
      # BLAST -------------------------------------------------------------------
      source_python("primertool2.py")
      primer.table1<-py$run_primertool(as.list(result1$primer1),primer_down,species,4000)
      primer.table2<-py$run_primertool(primer_up,as.list(result1$primer2),species,4000)
      primer.table <-rbind(primer.table1,primer.table2)
    }
  }
  return(primer.table)
}



#多个外显子
primer_design2_in <- function(Gene, PCR_seq1, PCR_seq2, Gene_rev, result1,pos1,pos2) {

  if (Gene_rev) {
    seq1 <- DNAString(PCR_seq1)
    seq1_rev <- reverse(seq1)
    PCR_seq1 <- as.character(seq1_rev)
    seq2 <- DNAString(PCR_seq2)
    seq2_rev <- reverse(seq2)
    PCR_seq2 <- as.character(seq2_rev)
  }
  
  #上游的PCR引物
  py$run(PCR_seq1, species, filepath2, "3")
  PCR.tab1 <- read.csv(paste0(filepath2,"//PCR3.csv"), header = FALSE,encoding = "UTF-8")
  ####40分以上
  PCR.tab1<-PCR.tab1[which(PCR.tab1$V3!="No matches"),]
  PCR.tab1<-PCR.tab1[which(as.numeric(PCR.tab1$V3)>=40),]
  if(nrow(PCR.tab1)==0){
    return(NULL)
  }
  
  #下游的PCR引物
  py$run(PCR_seq2, species, filepath2, "4")
  PCR.tab2 <- read.csv(paste0(filepath2,"//PCR4.csv"), header = FALSE,encoding = "UTF-8")
  ####40分以上
  PCR.tab2<-PCR.tab2[which(PCR.tab2$V3!="No matches"),]
  PCR.tab2<-PCR.tab2[which(as.numeric(PCR.tab2$V3)>=40),]
  if(nrow(PCR.tab2)==0){
    return(NULL)
  }
  
  #找出在外显子上的序列
  PCR.exon1<-data.frame()
  PCR.exon2<-data.frame()
  
  {
    # 反向基因 --------------------------------------------------------------------
    if (Gene_rev) {
      # 上游 ----------------------------------------------------------------------
      #删除连续5个碱基重复的序列
      PCR.table1 <- del_dup(PCR.tab1)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #在外显子上的优先
      for(i in 1:nrow(pos1)){
        PCR.exon1<-rbind(PCR.exon1,PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1])) %in% c(pos1[i,]$start:pos1[i,]$end)),])
      }
      a<-rbind(PCR.exon1,PCR.table1)
      PCR.table1<-a[which(!duplicated.data.frame(a)),]
      
      PCR.table_fw<-Get_primer2(PCR_seq1,PCR.table1)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      
      #去除重复
      PCR.table_fw <-
        PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
      #取30条序列
      if (nrow(PCR.table_fw) > 30) {
        PCR.table_fw <- PCR.table_fw[1:30, ]
      }
      write.csv(PCR.table_fw,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
      source_python("OligoCalc.py")
      py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
      PCR1.table<-read.csv(paste0(filepath2,"//PCR3.csv"),header = TRUE)
      #PCR1的反向互补序列
      for (i in 1:nrow(PCR1.table)) {
        pcr1 <- DNAString(PCR1.table[i, ]$PCR)
        pcr1_revcom <- reverseComplement(pcr1)
        PCR1.table[i, ]$PCR <- as.character(pcr1_revcom)
      }
      primer_down<-PCR1.table$PCR
      
      # 下游 ----------------------------------------------------------------------
      #去掉连续5个碱基相同的序列
      PCR.table2 <- del_dup(PCR.tab2)
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      #在外显子上的优先
      for(i in nrow(pos2):1){
        PCR.exon2<-rbind(PCR.exon2,PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1])) %in% c(pos2[i,]$start:pos2[i,]$end)),])
      }
      b<-rbind(PCR.exon2,PCR.table2)
      PCR.table2<-b[which(!duplicated.data.frame(b)),]
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_rev<-Get_primer(PCR_seq2,PCR.table2)
      
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_rev <- del_dup2(PCR.table_rev)
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      
      #去除重复
      PCR.table_rev <-
        PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
      #取30条
      if (nrow(PCR.table_rev) > 30) {
        PCR.table_rev <- PCR.table_rev[1:30, ]
      }
      write.csv(PCR.table_rev,paste0(filepath2,"//PCR4.csv"),row.names = FALSE)
      source_python("OligoCalc.py")
      py$RunOligoCalc(paste0(filepath2,"//PCR4.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
      
      primer_up<-read.csv(paste0(filepath2,"//PCR4.csv"),header = TRUE)$PCR
      
      # Blast -------------------------------------------------------------------
      source_python("primertool2.py")       
      primer.table<-py$run_primertool(as.list(result1$primer1),primer_down,species,4000)
      primer.table2<-py$run_primertool(primer_up,as.list(result1$primer2),species,4000)
      primer.table <-rbind(primer.table,primer.table2)
    }
    
    # 正向基因 --------------------------------------------------------------------
    else{
      
      # 上游 ----------------------------------------------------------------------
      PCR.table1 <- del_dup(PCR.tab1)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #在外显子上的优先
      for(i in 1:nrow(pos1)){
        PCR.exon1<-rbind(PCR.exon1,PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1])) %in% c(pos1[i,]$start:pos1[i,]$end)),])
      }
      a<-rbind(PCR.exon1,PCR.table1)
      PCR.table1<-a[which(!duplicated.data.frame(a)),]
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_fw<-Get_primer2(PCR_seq1,PCR.table1)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      
      #去除重复
      PCR.table_fw <- PCR.table_fw[!duplicated(PCR.table_fw$PCR), ]
      #取30条序列
      if (nrow(PCR.table_fw) > 30) {
        PCR.table_fw <- PCR.table_fw[1:30, ]
      }
      write.csv(PCR.table_fw,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
      source_python("OligoCalc.py")
      py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
      PCR1.table<-read.csv(paste0(filepath2,"//PCR3.csv"),header = TRUE)
      #PCR1的反向互补序列
      for (i in 1:nrow(PCR1.table)) {
        pcr1 <- DNAString(PCR1.table[i, ]$PCR)
        pcr1_revcom <- reverseComplement(pcr1)
        PCR1.table[i, ]$PCR <- as.character(pcr1_revcom)
      }
      primer_down<-PCR1.table$PCR
      
      # 下游 ----------------------------------------------------------------------
      #去掉连续5个碱基相同的序列
      PCR.table2 <- del_dup(PCR.tab2)
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      #在外显子上的优先
      for(i in nrow(pos2):1){
        PCR.exon2<-rbind(PCR.exon2,PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1])) %in% c(pos2[i,]$start:pos2[i,]$end)),])
      }
      b<-rbind(PCR.exon2,PCR.table2)
      PCR.table2<-b[which(!duplicated.data.frame(b)),]
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_rev<-Get_primer(PCR_seq2,PCR.table2)
      
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #删除连续4个重复
      PCR.table_rev <- del_dup2(PCR.table_rev)
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      
      #去除重复
      PCR.table_rev <-PCR.table_rev[!duplicated(PCR.table_rev$PCR), ]
      #取30条
      if (nrow(PCR.table_rev) > 30) {
        PCR.table_rev <- PCR.table_rev[1:30, ]
      }
      write.csv(PCR.table_rev,paste0(filepath2,"//PCR4.csv"),row.names = FALSE)
      source_python("OligoCalc.py")
      py$RunOligoCalc(paste0(filepath2,"//PCR4.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
      
      primer_up<-read.csv(paste0(filepath2,"//PCR4.csv"),header = TRUE)$PCR
      
      
      # Blast -------------------------------------------------------------------
      source_python("primertool2.py")       
      primer.table<-py$run_primertool(as.list(result1$primer1),primer_down,species,4000)
      primer.table2<-py$run_primertool(primer_up,as.list(result1$primer2),species,4000)
      primer.table <-rbind(primer.table,primer.table2)
    }
  }
  return(primer.table)
}


#只找一个靶位点
primer_design3_in<-function(Gene,PCR_seq,Gene_rev,result1){
  
  if (Gene_rev) {
    seq <- DNAString(PCR_seq)
    seq_rev <- reverse(seq)
    PCR_seq <- as.character(seq_rev)
  }
  py$run(PCR_seq, species, filepath2, "3")
  PCR.tab <- read.csv(paste0(filepath2,"//PCR3.csv"), header = FALSE,encoding = "UTF-8")
  ####40分以上
  PCR.tab<-PCR.tab[which(PCR.tab$V3!="No matches"),]
  PCR.tab<-PCR.tab[which(as.numeric(PCR.tab$V3)>=40),]
  if(nrow(PCR.tab)==0){
    return(NULL)
  }
  
  {
    # 反向基因 --------------------------------------------------------------------
    if (Gene_rev) {
      #删除连续5个碱基重复的序列
      PCR.table1 <- del_dup(PCR.tab)
      if(nrow(PCR.table1)!=0){
        #匹配PCR的位置，然后前后扩大3个碱基
        PCR.table_fw<-Get_primer(PCR_seq,PCR.table1)
        PCR.table_rev<-Get_primer2(PCR_seq,PCR.table1)
        if(nrow(PCR.table_fw)!=0){
          #删除连续5个重复
          PCR.table_fw <- del_dup2(PCR.table_fw)
          if(nrow(PCR.table_fw)!=0){
            #去除重复
            PCR.table_fw <-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
            #取30条序列
            if (nrow(PCR.table_fw) > 30) {
              PCR.table_fw <- PCR.table_fw[1:30,]
            }
            write.csv(PCR.table_fw,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
            source_python("OligoCalc.py")
            py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
            primer_up<-read.csv(paste0(filepath2,"//PCR3.csv"),header = TRUE)$PCR
            if(length(primer_up)==0){
              rm(primer_up)
            }
          }
        }
        if(nrow(PCR.table_rev)!=0){
          #删除连续5个重复
          PCR.table_rev <- del_dup2(PCR.table_rev)
          if(nrow(PCR.table_rev)!=0){
            #去除重复
            PCR.table_rev <-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
            #取30条序列
            if (nrow(PCR.table_rev) > 30) {
              PCR.table_rev <- PCR.table_rev[1:30,]
            }
            
            for(i in 1: nrow(PCR.table_rev)){
              pcr1<-DNAString(PCR.table_rev[i,]$PCR)
              pcr1_revcom<-reverseComplement(pcr1)
              PCR.table_rev[i,]$PCR<-as.character(pcr1_revcom)
            }
            
            write.csv(PCR.table_rev,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
            source_python("OligoCalc.py")
            py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
            primer_down<-read.csv(paste0(filepath2,"//PCR3.csv"),header = TRUE)$PCR
            if(length(primer_down)==0){
              rm(primer_down)
            }
          }
        }
      }

      # Blast -------------------------------------------------------------------
      source_python("primertool2.py")      
      if(exists("primer_down")){
        primer.table1<-py$run_primertool(as.list(result1$primer1),primer_down,species,4000)
        if (length(primer.table1) != 0) {
          primer1<-primer.table1[1]
          primer4<-primer.table1[2]
          primer1.pos<-matchPattern(reverse(DNAString(primer1)), subject = Gene)
          primer4.pos<-matchPattern(complement(DNAString(primer4)),subject = Gene)
          analysis_seq<-substring(Gene,min(start(primer1.pos), start(primer4.pos)),max(end(primer1.pos), end(primer4.pos)))
          if(analysis_GC(analysis_seq)!="FALSE"){
            WT2.distance<-max(end(primer1.pos),end(primer4.pos))-min(start(primer1.pos),start(primer4.pos))+1
            if(analysis_GC(analysis_seq)=="Low"){
              primer.table<-data.frame(F1=primer1,R1=result1$primer2,R2=primer4,F2=NA,WT_F1R1=result1$WT,WT_F1R2=WT2.distance,WT_F2R1=NA,GC1=result1$GC,GC2="Low",GC3=NA)
            }
            if(analysis_GC(analysis_seq)=="High"){
              primer.table<-data.frame(F1=primer1,R1=result1$primer2,R2=primer4,F2=NA,WT_F1R1=result1$WT,WT_F1R2=WT2.distance,WT_F2R1=NA,GC1=result1$GC,GC2="High",GC3=NA)
            }
            if(analysis_GC(analysis_seq)=="Normal"){
              primer.table<-data.frame(F1=primer1,R1=result1$primer2,R2=primer4,F2=NA,WT_F1R1=result1$WT,WT_F1R2=WT2.distance,WT_F2R1=NA,GC1=result1$GC,GC2="Normal",GC3=NA)
            }
          }
          if(analysis_GC(analysis_seq)=="FALSE"){
            print("GC high or low")
          }
          return(primer.table)
        }
        
        if(exists("primer_up")){
          primer.table2<-py$run_primertool(primer_up,as.list(result1$primer2),species,4000)
          if(length(primer.table2)!=0){
            primer3<-primer.table2[1]
            primer2<-primer.table2[2]
            primer3.pos<-matchPattern(reverse(DNAString(primer3)), subject = Gene)
            primer2.pos<-matchPattern(complement(DNAString(primer2)),subject = Gene)
            analysis_seq<-substring(Gene,min(start(primer2.pos), start(primer3.pos)),max(end(primer2.pos), end(primer3.pos)))
            if(analysis_GC(analysis_seq)!="FALSE"){
              WT2.distance<-max(end(primer2.pos),end(primer3.pos))-min(start(primer2.pos),start(primer3.pos))+1
              if(analysis_GC(analysis_seq)=="Low"){
                primer.table<-data.frame(F1=result1$primer1,R1=primer2,R2=NA,F2=primer3,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=WT2.distance,GC1=result1$GC,GC2=NA,GC3="Low")
              }
              if(analysis_GC(analysis_seq)=="High"){
                primer.table<-data.frame(F1=result1$primer1,R1=primer2,R2=NA,F2=primer3,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=WT2.distance,GC1=result1$GC,GC2=NA,GC3="High")
              }
              if(analysis_GC(analysis_seq)=="Normal"){
                primer.table<-data.frame(F1=result1$primer1,R1=primer2,R2=NA,F2=primer3,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=WT2.distance,GC1=result1$GC,GC2=NA,GC3="Normal")
              }
            }
            if(analysis_GC(analysis_seq)=="FALSE"){
              print("GC high or low")
            }
            return(primer.table)
          }
        }
      }
    }
    
    
    # 正向基因 --------------------------------------------------------------------
    else{
      PCR.table1<-del_dup(PCR.tab)
      if(nrow(PCR.table1)!=0){
        #匹配PCR的位置，然后前后扩大3个碱基
        PCR.table_fw<-Get_primer(PCR_seq,PCR.table1)
        PCR.table_rev<-Get_primer2(PCR_seq,PCR.table1)
        if(nrow(PCR.table_fw)!=0){
          #删除连续5个重复
          PCR.table_fw <- del_dup2(PCR.table_fw)
          if(nrow(PCR.table_fw)!=0){
            #去除重复
            PCR.table_fw <-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
            #取30条
            if (nrow(PCR.table_fw) > 30) {
              PCR.table_fw <- PCR.table_fw[1:30,]
            }
            write.csv(PCR.table_fw,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
            source_python("OligoCalc.py")
            py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
            primer_up<-read.csv(paste0(filepath2,"//PCR3.csv"), header = TRUE)$PCR
            if(length(primer_up)==0){
              rm(primer_up)
            }
          }
        }
        if(nrow(PCR.table_rev)!=0){
          #删除连续5个重复
          PCR.table_rev <- del_dup2(PCR.table_rev)
          if(nrow(PCR.table_rev)!=0){
            #去除重复
            PCR.table_rev <-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
            #取30条序列
            if (nrow(PCR.table_rev) > 30) {
              PCR.table_rev <- PCR.table_rev[1:30,]
            }
            
            for(i in 1: nrow(PCR.table_rev)){
              pcr1<-DNAString(PCR.table_rev[i,]$PCR)
              pcr1_revcom<-reverseComplement(pcr1)
              PCR.table_rev[i,]$PCR<-as.character(pcr1_revcom)
            }
            
            write.csv(PCR.table_rev,paste0(filepath2,"//PCR3.csv"),row.names = FALSE)
            source_python("OligoCalc.py")
            py$RunOligoCalc(paste0(filepath2,"//PCR3.csv"),"C:\\Users\\41518\\Desktop\\work\\Ubigene\\primer\\chromedriver.exe")
            primer_down<-read.csv(paste0(filepath2,"//PCR3.csv"),header = TRUE)$PCR
            if(length(primer_down)==0){
              rm(primer_down)
            }
          }
        }
      }
      
      # BLAST -------------------------------------------------------------------
      source_python("primertool2.py")
      if(exists("primer_down")){
        primer.table1<-py$run_primertool(as.list(result1$primer1),primer_down,species,4000)
        if (length(primer.table1) != 0) {
          primer1<-primer.table1[1]
          primer4<-primer.table1[2]
          primer1.pos<-matchPattern(primer1, subject = Gene)
          primer4.pos<-matchPattern(reverseComplement(DNAString(primer4)),subject = Gene)
          
          analysis_seq<-substring(Gene,min(start(primer1.pos), start(primer4.pos)),max(end(primer1.pos), end(primer4.pos)))
          if(analysis_GC(analysis_seq)!="FALSE"){
            WT2.distance<-max(end(primer1.pos),end(primer4.pos))-min(start(primer1.pos),start(primer4.pos))+1
            if(analysis_GC(analysis_seq)=="Low"){
              primer.table<-data.frame(F1=primer1,R1=result1$primer2,R2=primer4,F2=NA,WT_F1R1=result1$WT,WT_F1R2=WT2.distance,WT_F2R1=NA,GC1=result1$GC,GC2="Low",GC3=NA)
            }
            if(analysis_GC(analysis_seq)=="High"){
              primer.table<-data.frame(F1=primer1,R1=result1$primer2,R2=primer4,F2=NA,WT_F1R1=result1$WT,WT_F1R2=WT2.distance,WT_F2R1=NA,GC1=result1$GC,GC2="High",GC3=NA)
            }
            if(analysis_GC(analysis_seq)=="Normal"){
              primer.table<-data.frame(F1=primer1,R1=result1$primer2,R2=primer4,F2=NA,WT_F1R1=result1$WT,WT_F1R2=WT2.distance,WT_F2R1=NA,GC1=result1$GC,GC2="Normal",GC3=NA)
            }
          }
          if(analysis_GC(analysis_seq)=="FALSE"){
            print("GC high or low")
          }
          return(primer.table)
        }
        
        if(exists("primer_up")){
          primer.table2<-py$run_primertool(primer_up,as.list(result1$primer2),species,4000)
          if(length(primer.table2)!=0){
            primer3<-primer.table2[1]
            primer2<-primer.table2[2]
            primer3.pos<-matchPattern(primer3, subject = Gene)
            primer2.pos<-matchPattern(reverseComplement(DNAString(primer2)),subject = Gene)
            
            analysis_seq<-substring(Gene,min(start(primer2.pos), start(primer3.pos)),max(end(primer2.pos), end(primer3.pos)))
            if(analysis_GC(analysis_seq)!="FALSE"){
              WT2.distance<-max(end(primer2.pos),end(primer3.pos))-min(start(primer2.pos),start(primer3.pos))+1
              if(analysis_GC(analysis_seq)=="Low"){
                primer.table<-data.frame(F1=result1$primer1,R1=primer2,R2=NA,F2=primer3,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=WT2.distance,GC1=result1$GC,GC2=NA,GC3="Low")
              }
              if(analysis_GC(analysis_seq)=="High"){
                primer.table<-data.frame(F1=result1$primer1,R1=primer2,R2=NA,F2=primer3,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=WT2.distance,GC1=result1$GC,GC2=NA,GC3="High")
              }
              if(analysis_GC(analysis_seq)=="Normal"){
                primer.table<-data.frame(F1=result1$primer1,R1=primer2,R2=NA,F2=primer3,WT_F1R1=result1$WT,WT_F1R2=NA,WT_F2R1=WT2.distance,GC1=result1$GC,GC2=NA,GC3="Normal")
              }
            }
            if(analysis_GC(analysis_seq)=="FALSE"){
              print("GC high or low")
            }
            
            return(primer.table)
          }
        }
      }
    }
  }
}
