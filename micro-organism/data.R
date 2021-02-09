library(stringr)
library(Biostrings)
setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
aa<-readLines("C://Users//41518//Desktop//微生物//BL21//sequence.gb")
l <- grep("\\s+gene\\s+ ", aa)


seq_all<-readLines("C://Users//41518//Desktop//微生物//BL21//sequence.fasta")
seq_all<-seq_all[-1]
seq_all<-paste(seq_all,collapse = "")
gene.table6<-data.frame(gene=character(),locus_tag=character(),start=numeric(),end=numeric(),seq=character(),
                        analysis_seq=character(),GC=numeric(),h_start=character(),h_end=character(),
                        gRNA1=character(),strand1=character(),pos1=numeric(),score1=numeric(),gRNA2=character(),
                        strand2=character(),pos2=numeric(),score2=numeric(),direction=character(),
                        destroy=character(),up=character(),up_dir=character(),up_dis=numeric(),up_locus_tag=character(),
                        down=character(),down_dir=character(),down_dis=numeric(),down_locus_tag=character(),
                        ID=numeric(),symbol=character(),description=character(),locus=character(),type=character())
j<-1
for(i in l){
  if(str_extract_all(aa[i+1],"\\w+")[[1]][1]=="gene"){
    gene.table6[j,]$gene<-str_extract_all(aa[i+1],"\\w+")[[1]][2]
    gene.table6[j,]$locus_tag<-str_extract_all(aa[i+2],"\\w+")[[1]][2]
    
  }
  if(str_extract_all(aa[i+1],"\\w+")[[1]][1]=="locus_tag"){
    gene.table6[j,]$locus_tag<-str_extract_all(aa[i+1],"\\w+")[[1]][2]
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


gene.table7<-gene.table6

gene.table1<-read.csv("C://Users//41518//Desktop//微生物//micro.csv")
gene.table2<-read.csv("C://Users//41518//Desktop//微生物//micro-补充.csv")


for(i in 1:nrow(gene.table2)){
  j<-which(gene.table6$locus_tag==gene.table2[i,]$locus_tag)
  if(length(j)!=0){
    print(i)
    gene.table6[j,]$seq<-gene.table2[i,]$seq
    gene.table6[j,]$analysis_seq<-gene.table2[i,]$analysis_seq
    gene.table6[j,]$start<-gene.table2[i,]$start
    gene.table6[j,]$end<-gene.table2[i,]$end
    gene.table6[j,]$GC<-gene.table2[i,]$GC
    gene.table6[j,]$h_start<-gene.table2[i,]$h_start
    gene.table6[j,]$h_end<-gene.table2[i,]$h_end
    gene.table6[j,]$gRNA1<-gene.table2[i,]$gRNA1
    gene.table6[j,]$strand1<-gene.table2[i,]$strand1
    gene.table6[j,]$pos1<-gene.table2[i,]$pos1
    gene.table6[j,]$score1<-gene.table2[i,]$score1
    gene.table6[j,]$gRNA2<-gene.table2[i,]$gRNA2
    gene.table6[j,]$strand2<-gene.table2[i,]$strand2
    gene.table6[j,]$pos2<-gene.table2[i,]$pos2
    gene.table6[j,]$score2<-gene.table2[i,]$score2
    gene.table6[j,]$ID<-gene.table2[i,]$ID
    gene.table6[j,]$symbol<-gene.table2[i,]$symbol
    gene.table6[j,]$description<-gene.table2[i,]$description
    gene.table6[j,]$locus<-gene.table2[i,]$locus
    gene.table6[j,]$type<-gene.table2[i,]$type
  }
}


for(i in 1:nrow(gene.table1)){
  j<-which(gene.table6$locus_tag==gene.table1[i,]$locus_tag)
  if(length(j)!=0){
    print(i)
    print(j)
    gene.table6[j,]$seq<-gene.table1[i,]$seq
    gene.table6[j,]$analysis_seq<-gene.table1[i,]$analysis_seq
    gene.table6[j,]$start<-gene.table1[i,]$start
    gene.table6[j,]$end<-gene.table1[i,]$end
    gene.table6[j,]$GC<-gene.table1[i,]$GC
    gene.table6[j,]$h_start<-gene.table1[i,]$h_start
    gene.table6[j,]$h_end<-gene.table1[i,]$h_end
    gene.table6[j,]$gRNA1<-gene.table1[i,]$gRNA1
    gene.table6[j,]$strand1<-gene.table1[i,]$strand1
    gene.table6[j,]$pos1<-gene.table1[i,]$pos1
    gene.table6[j,]$score1<-gene.table1[i,]$score1
    gene.table6[j,]$gRNA2<-gene.table1[i,]$gRNA2
    gene.table6[j,]$strand2<-gene.table1[i,]$strand2
    gene.table6[j,]$pos2<-gene.table1[i,]$pos2
    gene.table6[j,]$score2<-gene.table1[i,]$score2
    gene.table6[j,]$ID<-gene.table1[i,]$ID
    gene.table6[j,]$symbol<-gene.table1[i,]$symbol
    gene.table6[j,]$description<-gene.table1[i,]$description
    gene.table6[j,]$locus<-gene.table1[i,]$locus
    gene.table6[j,]$type<-gene.table1[i,]$type
  }
}

write.csv(gene.table,"C://Users//41518//Desktop//微生物//micro-全.csv",row.names = FALSE)
