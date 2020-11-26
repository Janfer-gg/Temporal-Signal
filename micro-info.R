library(stringr)
library(Biostrings)
setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
aa<-readLines("C://Users//41518//Desktop//微生物//Bacillus//sequence6.gb")
l <- grep("\\s+gene\\s+ ", aa)


seq_all<-readLines("C://Users//41518//Desktop//微生物//Bacillus//sequence6.fasta")
seq_all<-seq_all[-1]
seq_all<-paste(seq_all,collapse = "")
gene.table6<-data.frame(species=character(),gene=character(),locus_tag=character(),start=numeric(),end=numeric(),seq=character(),
                        analysis_seq=character(),GC=numeric(),h_start=character(),h_end=character(),
                        gRNA1=character(),strand1=character(),pos1=numeric(),score1=numeric(),gRNA2=character(),
                        strand2=character(),pos2=numeric(),score2=numeric(),direction=character(),
                        destroy=character(),up=character(),up_dir=character(),up_dis=numeric(),up_locus_tag=character(),
                        down=character(),down_dir=character(),down_dis=numeric(),down_locus_tag=character(),locus=character()
                        )
j<-1
for(i in l){
  if(str_extract_all(aa[i+1],"\\w+")[[1]][1]=="gene"){
    gene.table6[j,]$gene<-str_extract_all(aa[i+1],"\\w+")[[1]][2]
    gene.table6[j,]$locus_tag<-str_extract_all(aa[i+2],"\\w+")[[1]][2]
    gene.table6[j,]$locus<-str_extract_all(aa[i+3],"\\w+")[[1]][2]
  }
  if(str_extract_all(aa[i+1],"\\w+")[[1]][1]=="locus_tag"){
    gene.table6[j,]$locus_tag<-str_extract_all(aa[i+1],"\\w+")[[1]][2]
    gene.table6[j,]$locus<-str_extract_all(aa[i+2],"\\w+")[[1]][2]
  }
  
  gene.table6[j,]$start<-str_extract_all(aa[i],"\\d+")[[1]][1]
  gene.table6[j,]$end<-str_extract_all(aa[i],"\\d+")[[1]][2]
  {
    if(grepl("complement",aa[i])){
      gene.table6[j,]$direction<-"rev"
    }
    else{
      gene.table6[j,]$direction<-"fw"
    }
  }
  
  j<-j+1
}


for(i in 1:(nrow(gene.table6))){
  gg<-character()
  if(i==1){
    if (gene.table6[i + 1, ]$direction == "fw") {
      gg <- paste(gg, gene.table6[i + 1, ]$gene, sep = "")
    }
    gene.table6[i,]$down<-gene.table6[i+1,]$gene
    gene.table6[i,]$down_dir<-gene.table6[i+1,]$direction
    gene.table6[i,]$down_dis<-as.numeric(gene.table6[i+1,]$start)-as.numeric(gene.table6[i,]$end)
    gene.table6[i,]$down_locus_tag<-gene.table6[i+1,]$locus_tag
  }
  
  else if(i==nrow(gene.table6)){
    if (gene.table6[i - 1, ]$direction == "rev") {
      gg <- paste(gg, gene.table6[i + 1, ]$gene, sep = "")
    }
    gene.table6[i,]$up<-gene.table6[i-1,]$gene
    gene.table6[i,]$up_dir<-gene.table6[i-1,]$direction
    gene.table6[i,]$up_dis<-as.numeric(gene.table6[i,]$start)-as.numeric(gene.table6[i-1,]$end)
    gene.table6[i,]$up_locus_tag<-gene.table6[i-1,]$locus_tag
  }
  else{
    if (gene.table6[i + 1, ]$direction == "fw") {
      gg <- paste(gg, gene.table6[i + 1, ]$gene, sep = "")
    }
    if (gene.table6[i - 1, ]$direction == "rev") {
      if(length(gg)==0){
        gg <- paste(gg, gene.table6[i - 1, ]$gene)
      }
      else{
        gg <- paste(gg, gene.table6[i - 1, ]$gene, sep = ",")
      }
    }
    gene.table6[i,]$up<-gene.table6[i-1,]$gene
    gene.table6[i,]$up_dir<-gene.table6[i-1,]$direction
    gene.table6[i,]$up_dis<-as.numeric(gene.table6[i,]$start)-as.numeric(gene.table6[i-1,]$end)
    gene.table6[i,]$up_locus_tag<-gene.table6[i-1,]$locus_tag
    gene.table6[i,]$down<-gene.table6[i+1,]$gene
    gene.table6[i,]$down_dir<-gene.table6[i+1,]$direction
    gene.table6[i,]$down_dis<-as.numeric(gene.table6[i+1,]$start)-as.numeric(gene.table6[i,]$end)
    gene.table6[i,]$down_locus_tag<-gene.table6[i+1,]$locus_tag
  }
  
  if(length(gg)!=0){
    gene.table6[i,]$destroy<-gg
  }
}

gene.table7<-rbind(gene.table1,gene.table2,gene.table3,gene.table4,gene.table5,gene.table6)
write.csv(gene.table7,"C://Users//41518//Desktop//微生物//Bacillus//Bacillus-info.csv",row.names = FALSE)
gene.table<-read.csv("C://Users//41518//Desktop//微生物//Bacillus//Bacillus.csv")




