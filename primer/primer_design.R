primer_design1<-function(Gene,gRNA,Gene_rev,size,filepath){
  PCR_seq1 <- substring(Gene, min(gRNA$start) - 800, min(gRNA$start) - 250)
  PCR_seq2 <- substring(Gene, max(gRNA$end) +250, max(gRNA$end) +800)
  
  if (Gene_rev) {
    seq1 <- DNAString(PCR_seq1)
    seq1_rev <- reverse(seq1)
    PCR_seq1 <- as.character(seq1_rev)
    seq2 <- DNAString(PCR_seq2)
    seq2_rev <- reverse(seq2)
    PCR_seq2 <- as.character(seq2_rev)
  }
  
  #上游的PCR引物
  source_python("crispor_table_download2.py")
  py$run(PCR_seq1,species,filepath,"1")
  PCR.tab1<-read.csv(paste0(filepath,"//PCR1.csv"), header = FALSE,encoding = "UTF-8")
  
  #下游的PCR引物
  py$run(PCR_seq2,species,filepath,"2")
  PCR.tab2<-read.csv(paste0(filepath,"//PCR2.csv"), header = FALSE,encoding = "UTF-8")
  
  {
    # 反向基因 --------------------------------------------------------------------
    if (Gene_rev) {
      # 上游 ----------------------------------------------------------------------
      #删除连续5个碱基重复的序列
      PCR.table1 <- del_dup(PCR.tab1)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #40分以上
      PCR.table1<-PCR.table1[which(PCR.table1$V3!="No matches"),]
      PCR.table1<-PCR.table1[which(as.numeric(PCR.table1$V3)>=40),]
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))<=250),]
      PCR.table_800<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))>250),]
      PCR.table1<-rbind(PCR.table_500,PCR.table_800)
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_rev<-Get_primer2(PCR_seq1,PCR.table1)
      
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_rev <- del_dup2(PCR.table_rev)
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_rev<-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
      #取40条序列
      if(nrow(PCR.table_rev)>40){
        PCR.table_rev<-PCR.table_rev[1:40,]
      }
      #PCR1的反向互补序列
      for(i in 1: nrow(PCR.table_rev)){
        pcr1<-DNAString(PCR.table_rev[i,]$PCR)
        pcr1_revcom<-reverseComplement(pcr1)
        PCR.table_rev[i,]$PCR<-as.character(pcr1_revcom)
      }
      
      primer_up<-py$Oligo(as.list(PCR.table_rev$PCR))
      write.csv(primer_up,paste0(filepath,"//PCR1.csv"),row.names = FALSE)

      # 下游 ----------------------------------------------------------------------
      #去掉连续5个碱基相同的序列
      PCR.table2 <- del_dup(PCR.tab2)
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      #40分以上
      PCR.table2<-PCR.table2[which(PCR.table2$V3!="No matches"),]
      PCR.table2<-PCR.table2[which(as.numeric(PCR.table2$V3)>=40),]
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      
      PCR.table_500<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))>=300),]
      PCR.table_800<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))<300),]
      PCR.table2<-rbind(PCR.table_500,PCR.table_800)
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_fw<-Get_primer(PCR_seq2,PCR.table2)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_fw<-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
      #取40条
      if(nrow(PCR.table_fw)>40){
        PCR.table_fw<-PCR.table_fw[1:40,]
      }
      
      primer_down<-py$Oligo(as.list(PCR.table_fw$PCR))
      write.csv(primer_down,paste0(filepath,"//PCR2.csv"),row.names = FALSE)
     
      if(length(primer_up)==0 | length(primer_down)==0){
        return(NULL)
      }
      # Blast -------------------------------------------------------------------
      
      primer.table<-py$run_primertool(as.list(primer_up),as.list(primer_down),species,size)
      
      #primer1是上游，primer2是下游
      primer1<-primer.table[2]
      primer2<-primer.table[1]
      
      primer1.pos <-matchPattern(reverse(DNAString(primer1)), subject = Gene)
      primer2.pos <-matchPattern(complement(DNAString(primer2)), subject = Gene)
      analysis_seq<-substring(Gene,min(start(primer1.pos), start(primer2.pos)),max(end(primer1.pos), end(primer2.pos)))
      if(analysis_GC(analysis_seq)!="FALSE"){
        WT.distance <-max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) + 1
        if(analysis_GC(analysis_seq)=="Low"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Low")
        }
        if(analysis_GC(analysis_seq)=="High"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="High")
        }
        if(analysis_GC(analysis_seq)=="Normal"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Normal")
        }
      }
      if(analysis_GC(analysis_seq)=="FALSE"){
        result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="FALSE")
      }
      
      #如果外侧引物超过500bp,在500bp内再设一条
      if(start(primer1.pos)-max(gRNA$end)-1>500){
        if(primer_down[1]!=primer1){
          primer1_500<-primer_down[1]
          primer1_500.pos <-matchPattern(reverse(DNAString(primer1_500)), subject = Gene)
          if(start(primer1_500.pos)-max(gRNA$end)-1<=500){
            result$primer1_500<-primer1_500
          }
        }
      }
      if(min(gRNA$start)-end(primer2.pos)-1>500){
        if(primer_up[1]!=primer2){
          primer2_500<-primer_up[1]
          primer2_500.pos<-matchPattern(complement(DNAString(primer2_500)), subject = Gene)
          if(min(gRNA$start)-end(primer2_500.pos)-1<=500){
            result$primer2_500<-primer2_500
          }
        }
      }
    }
    
    # 正向基因 --------------------------------------------------------------------
    else{
      
      # 上游 ----------------------------------------------------------------------
      PCR.table1<-del_dup(PCR.tab1)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #40分以上
      PCR.table1<-PCR.table1[which(PCR.table1$V3!="No matches"),]
      PCR.table1<-PCR.table1[which(as.numeric(PCR.table1$V3)>=40),]
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))>=250),]
      PCR.table_800<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))<250),]
      PCR.table1<-rbind(PCR.table_500,PCR.table_800)
      
      PCR.table_fw<-Get_primer(PCR_seq1,PCR.table1)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      
      #去除重复
      PCR.table_fw<-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
      #取40条
      if(nrow(PCR.table_fw)>40){
        PCR.table_fw<-PCR.table_fw[1:40,]
      }
      
      primer_up<-py$Oligo(as.list(PCR.table_fw$PCR))
      write.csv(primer_up,paste0(filepath,"//PCR1.csv"),row.names = FALSE)
  
      
      # 下游 ----------------------------------------------------------------------
      #去掉连续5个碱基相同的序列
      PCR.table2 <- del_dup(PCR.tab2)
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      #40分以上
      PCR.table2<-PCR.table2[which(PCR.table2$V3!="No matches"),]
      PCR.table2<-PCR.table2[which(as.numeric(PCR.table2$V3)>=40),]
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))<=250),]
      PCR.table_800<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))>250),]
      PCR.table2<-rbind(PCR.table_500,PCR.table_800)
      
      PCR.table_rev<-Get_primer2(PCR_seq2,PCR.table2)
      
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_rev <- del_dup2(PCR.table_rev)
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_rev<-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
      #取40条
      if(nrow(PCR.table_rev)>40){
        PCR.table_rev<-PCR.table_rev[1:40,]
      }
      for(i in 1:nrow(PCR.table_rev)){
        PCR.table_rev[i,]$PCR<-as.character(reverseComplement(DNAString(PCR.table_rev[i,]$PCR)))
      }
      
      primer_down<-py$Oligo(as.list(PCR.table_rev$PCR))
      write.csv(primer_down,paste0(filepath,"//PCR2.csv"),row.names = FALSE)

      if(length(primer_up)==0 | length(primer_down)==0){
        return(NULL)
      }
      # BLAST -------------------------------------------------------------------
      primer.table<-py$run_primertool(as.list(primer_up),as.list(primer_down),species,size)
      
      primer1<-primer.table[1]
      primer2<-primer.table[2]
      
      primer1.pos <-matchPattern(DNAString(primer1), subject = Gene)
      primer2.pos <-matchPattern(reverseComplement(DNAString(primer2)), subject = Gene)
      analysis_seq<-substring(Gene,min(start(primer1.pos), start(primer2.pos)),max(end(primer1.pos), end(primer2.pos)))
      if(analysis_GC(analysis_seq)!="FALSE"){
        WT.distance <-max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) + 1
        if(analysis_GC(analysis_seq)=="Low"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Low")
        }
        if(analysis_GC(analysis_seq)=="High"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="High")
        }
        if(analysis_GC(analysis_seq)=="Normal"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Normal")
        }
      }
      if(analysis_GC(analysis_seq)=="FALSE"){
        result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="FALSE")
      }
      
      #如果外侧引物超过500bp,在500bp内再设一条
      if(min(gRNA$start)-end(primer1.pos)-1>500){
        if(primer_up[1]!=primer1){
          primer1_500<-primer_up[1]
          primer1_500.pos <-matchPattern(DNAString(primer1_500), subject = Gene)
          if(min(gRNA$start)-end(primer1_500.pos)-1<=500){
            result$primer1_500<-primer1_500
          }
        }
      }
      if(start(primer2.pos)-max(gRNA$end)-1>500){
        if(primer_down[1]!=primer2){
          primer2_500<-primer_down[1]
          primer2_500.pos<-matchPattern(reverseComplement(DNAString(primer2_500)), subject = Gene)
          if(start(primer2_500.pos)-max(gRNA$end)-1<=500){
            result$primer2_500<-primer2_500
          }
        }
      }
    }
  }
  return(result)
}



primer_design2<-function(Gene,gRNA,gRNA2,Gene_rev,size,filepath){
  PCR_seq1 <- substring(Gene, min(gRNA$start,gRNA2$start) - 800, min(gRNA$start,gRNA2$start) - 250)
  PCR_seq2 <- substring(Gene, max(gRNA$end,gRNA2$end) + 250, max(gRNA$end,gRNA2$end) + 800)
  
  if (Gene_rev) {
    seq1 <- DNAString(PCR_seq1)
    seq1_rev <- reverse(seq1)
    PCR_seq1 <- as.character(seq1_rev)
    seq2 <- DNAString(PCR_seq2)
    seq2_rev <- reverse(seq2)
    PCR_seq2 <- as.character(seq2_rev)
  }
  
  #下游的PCR引物
  source_python("crispor_table_download2.py")
  py$run(PCR_seq1,species,filepath,"1")
  PCR.tab1<-read.csv(paste0(filepath,"//PCR1.csv"), header = FALSE,encoding = "UTF-8")
  
  #上游的PCR引物
  py$run(PCR_seq2,species,filepath,"2")
  PCR.tab2<-read.csv(paste0(filepath,"//PCR2.csv"), header = FALSE,encoding = "UTF-8")
  
  {
    # 反向基因 --------------------------------------------------------------------
    if (Gene_rev) {
      # 上游 ----------------------------------------------------------------------
      #删除连续5个碱基重复的序列
      PCR.table1 <- del_dup(PCR.tab1)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #40分以上
      PCR.table1<-PCR.table1[which(PCR.table1$V3!="No matches"),]
      PCR.table1<-PCR.table1[which(as.numeric(PCR.table1$V3)>=40),]
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))<=250),]
      PCR.table_800<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))>250),]
      PCR.table1<-rbind(PCR.table_500,PCR.table_800)
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_rev<-Get_primer2(PCR_seq1,PCR.table1)
      
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_rev <- del_dup2(PCR.table_rev)
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_rev<-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
      #取40条序列
      if(nrow(PCR.table_rev)>40){
        PCR.table_rev<-PCR.table_rev[1:40,]
      }
      #PCR1的反向互补序列
      for(i in 1: nrow(PCR.table_rev)){
        pcr1<-DNAString(PCR.table_rev[i,]$PCR)
        pcr1_revcom<-reverseComplement(pcr1)
        PCR.table_rev[i,]$PCR<-as.character(pcr1_revcom)
      }
      
      primer_up<-py$Oligo(as.list(PCR.table_rev$PCR))
      write.csv(primer_up,paste0(filepath,"//PCR1.csv"),row.names = FALSE)
      
      # 下游 ----------------------------------------------------------------------
      #去掉连续5个碱基相同的序列
      PCR.table2 <- del_dup(PCR.tab2)
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      #40分以上
      PCR.table2<-PCR.table2[which(PCR.table2$V3!="No matches"),]
      PCR.table2<-PCR.table2[which(as.numeric(PCR.table2$V3)>=40),]
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))>=300),]
      PCR.table_800<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))<300),]
      PCR.table2<-rbind(PCR.table_500,PCR.table_800)
      
      #匹配PCR的位置，然后前后扩大3个碱基
      PCR.table_fw<-Get_primer(PCR_seq2,PCR.table2)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_fw<-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
      #取40条
      if(nrow(PCR.table_fw)>40){
        PCR.table_fw<-PCR.table_fw[1:40,]
      }
      
      primer_down<-py$Oligo(as.list(PCR.table_fw$PCR))
      write.csv(primer_down,paste0(filepath,"//PCR2.csv"),row.names = FALSE)

      if(length(primer_up)==0 | length(primer_down)==0){
        return(NULL)
      }
      # Blast -------------------------------------------------------------------
      primer.table<-py$run_primertool(as.list(primer_up),as.list(primer_down),species,size)
      
      #primer1是上游，primer2是下游
      primer1<-primer.table[2]
      primer2<-primer.table[1]
      
      primer1.pos <-
        matchPattern(reverse(DNAString(primer1)), subject = Gene)
      primer2.pos <-
        matchPattern(complement(DNAString(primer2)), subject = Gene)
      analysis_seq<-substring(Gene,min(start(primer1.pos), start(primer2.pos)),max(end(primer1.pos), end(primer2.pos)))
      if(analysis_GC(analysis_seq)!="FALSE"){
        WT.distance <-max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) + 1
        if(analysis_GC(analysis_seq)=="Low"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Low")
        }
        if(analysis_GC(analysis_seq)=="High"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="High")
        }
        if(analysis_GC(analysis_seq)=="Normal"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Normal")
        }
      }
      if(analysis_GC(analysis_seq)=="FALSE"){
        result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="FALSE")
      }
      
      #如果外侧引物超过500bp,在500bp内再设一条
      if(start(primer1.pos)-max(gRNA$end,gRNA2$end)-1>500){
        if(primer_down[1]!=primer1){
          primer1_500<-primer_down[1]
          primer1_500.pos <-matchPattern(reverse(DNAString(primer1_500)), subject = Gene)
          if(start(primer1_500.pos)-max(gRNA$end,gRNA2$end)-1<=500){
            result$primer1_500<-primer1_500
          }
        }
      }
      if(min(gRNA$start,gRNA2$start)-end(primer2.pos)-1>500){
        if(primer_up[1]!=primer2){
          primer2_500<-primer_up[1]
          primer2_500.pos<-matchPattern(complement(DNAString(primer2_500)), subject = Gene)
          if(min(gRNA$start,gRNA2$start)-end(primer2_500.pos)-1<=500){
            result$primer2_500<-primer2_500
          }
        }
      }
    }
    
    # 正向基因 --------------------------------------------------------------------
    else{
      # 上游 ----------------------------------------------------------------------
      PCR.table1<-del_dup(PCR.tab1)
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      #40分以上
      PCR.table1<-PCR.table1[which(PCR.table1$V3!="No matches"),]
      PCR.table1<-PCR.table1[which(as.numeric(PCR.table1$V3)>=40),]
      if(nrow(PCR.table1)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))>=300),]
      PCR.table_800<-PCR.table1[which(as.numeric(sub("[^0-9]+", "", PCR.table1[,1]))<300),]
      PCR.table1<-rbind(PCR.table_500,PCR.table_800)
      
      PCR.table_fw<-Get_primer(PCR_seq1,PCR.table1)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_fw <- del_dup2(PCR.table_fw)
      
      if(nrow(PCR.table_fw)==0){
        return(NULL)
      }
      
      #去除重复
      PCR.table_fw<-PCR.table_fw[!duplicated(PCR.table_fw$PCR),]
      #取40条
      if(nrow(PCR.table_fw)>40){
        PCR.table_fw<-PCR.table_fw[1:40,]
      }
      
      primer_up<-py$Oligo(as.list(PCR.table_fw$PCR))
      write.csv(primer_up,paste0(filepath,"//PCR1.csv"),row.names = FALSE)
      

      # 下游 ----------------------------------------------------------------------
      #去掉连续5个碱基相同的序列
      PCR.table2 <- del_dup(PCR.tab2)
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      #40分以上
      PCR.table2<-PCR.table2[which(PCR.table2$V3!="No matches"),]
      PCR.table2<-PCR.table2[which(as.numeric(PCR.table2$V3)>=40),]
      if(nrow(PCR.table2)==0){
        return(NULL)
      }
      PCR.table_500<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))<=250),]
      PCR.table_800<-PCR.table2[which(as.numeric(sub("[^0-9]+", "", PCR.table2[,1]))>250),]
      PCR.table2<-rbind(PCR.table_500,PCR.table_800)
      
      PCR.table_rev<-Get_primer2(PCR_seq2,PCR.table2)
      
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #删除连续5个重复
      PCR.table_rev <- del_dup2(PCR.table_rev)
      if(nrow(PCR.table_rev)==0){
        return(NULL)
      }
      #去除重复
      PCR.table_rev<-PCR.table_rev[!duplicated(PCR.table_rev$PCR),]
      #取40条
      if(nrow(PCR.table_rev)>40){
        PCR.table_rev<-PCR.table_rev[1:40,]
      }
      for(i in 1:nrow(PCR.table_rev)){
        PCR.table_rev[i,]$PCR<-as.character(reverseComplement(DNAString(PCR.table_rev[i,]$PCR)))
      }
      
      primer_down<-py$Oligo(as.list(PCR.table_rev$PCR))
      write.csv(primer_down,paste0(filepath,"//PCR2.csv"),row.names = FALSE)

      if(length(primer_up)==0 | length(primer_down)==0){
        return(NULL)
      }
      # BLAST -------------------------------------------------------------------
      primer.table<-py$run_primertool(as.list(primer_up),as.list(primer_down),species,size)
      
      primer1<-primer.table[1]
      primer2<-primer.table[2]
      
      primer1.pos <-matchPattern(DNAString(primer1), subject = Gene)
      primer2.pos <-matchPattern(reverseComplement(DNAString(primer2)), subject = Gene)
      
      analysis_seq<-substring(Gene,min(start(primer1.pos), start(primer2.pos)),max(end(primer1.pos), end(primer2.pos)))
      if(analysis_GC(analysis_seq)!="FALSE"){
        WT.distance <-max(end(primer1.pos), end(primer2.pos)) - min(start(primer1.pos), start(primer2.pos)) + 1
        if(analysis_GC(analysis_seq)=="Low"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Low")
        }
        if(analysis_GC(analysis_seq)=="High"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="High")
        }
        if(analysis_GC(analysis_seq)=="Normal"){
          result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="Normal")
        }
      }
      if(analysis_GC(analysis_seq)=="FALSE"){
        result<-data.frame(primer1=primer1,primer2=primer2,WT=WT.distance,GC="FALSE")
      }
      
      #如果外侧引物超过500bp,在500bp内再设一条
      if(min(gRNA$start,gRNA2$start)-end(primer1.pos)-1>500){
        if(primer_up[1]!=primer1){
          primer1_500<-primer_up[1]
          primer1_500.pos <-matchPattern(DNAString(primer1_500), subject = Gene)
          if(min(gRNA$start,gRNA2$start)-end(primer1_500.pos)-1<=500){
            result$primer1_500<-primer1_500
          }
        }
      }
      if(start(primer2.pos)-max(gRNA$end,gRNA2$end)-1>500){
        if(primer_down[1]!=primer2){
          primer2_500<-primer_down[1]
          primer2_500.pos<-matchPattern(reverseComplement(DNAString(primer2_500)), subject = Gene)
          if(start(primer2_500.pos)-max(gRNA$end,gRNA2$end)-1<=500){
            result$primer2_500<-primer2_500
          }
        }
      }
    }
  }
  return(result)
}
