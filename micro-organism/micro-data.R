library(stringr)
library(Biostrings)
# setwd("C://Users//41518//Desktop//work/ubigene/micro-organism")
# aa<-readLines("C://Users//41518//Desktop//微生物//Saccharomyces//sequence17.gb")
# l <- grep("\\s+gene\\s+ ", aa)
# 
# 
# seq_all<-readLines("C://Users//41518//Desktop//微生物//Saccharomyces//sequence1.fasta")
# seq_all<-seq_all[-1]
# seq_all<-paste(seq_all,collapse = "")
# gene.table67<-data.frame(gene=character(),locus_tag=character(),start=numeric(),end=numeric(),seq=character(),
#                        analysis_seq=character(),GC=numeric(),h_start=character(),h_end=character(),
#                        gRNA1=character(),strand1=character(),pos1=numeric(),score1=numeric(),gRNA2=character(),
#                        strand2=character(),pos2=numeric(),score2=numeric(),direction=character(),
#                        destroy=character(),up=character(),up_dir=character(),up_dis=numeric(),up_locus_tag=character(),
#                        down=character(),down_dir=character(),down_dis=numeric(),down_locus_tag=character(),ID=numeric())
# j<-1
# for(i in l){
#   if(str_extract_all(aa[i+1],"\\w+")[[1]][1]=="gene"){
#     gene.table67[j,]$gene<-str_extract_all(aa[i+1],"\\w+")[[1]][2]
#     gene.table67[j,]$locus_tag<-str_extract_all(aa[i+2],"\\w+")[[1]][2]
#   }
# 
#   if(str_extract_all(aa[i+1],"\\w+")[[1]][1]=="locus_tag"){
#     gene.table67[j,]$locus_tag<-str_extract_all(aa[i+1],"\\w+")[[1]][2]
# 
#   }
# 
#   for(k in i:(i+10)){
#     if(str_extract_all(aa[k],"\\w+")[[1]][1]=="db_xref"){
#       if(str_extract_all(aa[k],"\\w+")[[1]][2]=="GeneID"){
#         gene.table67[j,]$ID<-str_extract_all(aa[k],"\\d+")[[1]][1]
#         break
#       }
#     }
#   }
# 
#   gene.table67[j,]$start<-str_extract_all(aa[i],"\\d+")[[1]][1]
#   gene.table67[j,]$end<-str_extract_all(aa[i],"\\d+")[[1]][2]
#   {
#     if(grepl("complement",aa[i])){
#       gene.table67[j,]$direction<-"rev"
#     }
#     else{
#       gene.table67[j,]$direction<-"fw"
#     }
#   }
# 
#   j<-j+1
# }

gene.table1$seq<-substring(seq_all,as.numeric(gene.table1$start),as.numeric(gene.table6[NA.count1,]$end))
gene.table1$analysis_seq<-substring(seq_all,as.numeric(gene.table1$start)-500,as.numeric(gene.table1$end)+500)

for(i in 1:nrow(gene.table1)){
 
  GC<-sum(letterFrequency(DNAString(gene.table1[i,]$analysis_seq),c("G","C")))/nchar(gene.table1[i,]$analysis_seq)
  gene.table1[i,]$GC<-round(GC*100,2)
}

#发夹结构
for(m in 153:6400){
  seq<-DNAString(gene.table[m,]$analysis_seq)
  hairpin<-findPalindromes(seq,min.armlength = 4,max.looplength = 17)
  bb<-data.frame(start=start(hairpin),end=end(hairpin))
  bb<-bb[order(bb$start),]
  
  for(i in 1:nrow(bb)){
    for(j in 1:nrow(bb)){
      if(bb[j,]$start>=bb[i,]$start & bb[j,]$start<=bb[i,]$end){
        if(bb[j,]$end>bb[i,]$end){
          bb[i,]$end<-bb[j,]$end
        }
        else{
          bb[j,]$end<-bb[i,]$end
        }
      }
      else{
        next
      }
    }
  }
  bb<-bb[!duplicated(bb$end),]
  gene.table[m,]$h_start<-paste(bb$start,collapse = ";")
  gene.table[m,]$h_end<-paste(bb$end,collapse = ";")
  print(m)
}
write.csv(gene.table,"C://Users//41518//Desktop//微生物//Sac-single//Sac.csv",row.names = FALSE)

# for(i in 1:(nrow(gene.table67))){
#   gg<-character()
#   if(i==1){
#     if (gene.table67[i + 1, ]$direction == "fw") {
#       gg <- paste(gg, gene.table67[i + 1, ]$gene, sep = "")
#     }
#     gene.table67[i,]$down<-gene.table67[i+1,]$gene
#     gene.table67[i,]$down_dir<-gene.table67[i+1,]$direction
#     gene.table67[i,]$down_dis<-as.numeric(gene.table67[i+1,]$start)-as.numeric(gene.table67[i,]$end)
#     gene.table67[i,]$down_locus_tag<-gene.table67[i+1,]$locus_tag
#   }
# 
#   else if(i==nrow(gene.table67)){
#     if (gene.table67[i - 1, ]$direction == "rev") {
#       gg <- paste(gg, gene.table67[i + 1, ]$gene, sep = "")
#     }
#     gene.table67[i,]$up<-gene.table67[i-1,]$gene
#     gene.table67[i,]$up_dir<-gene.table67[i-1,]$direction
#     gene.table67[i,]$up_dis<-as.numeric(gene.table67[i,]$start)-as.numeric(gene.table67[i-1,]$end)
#     gene.table67[i,]$up_locus_tag<-gene.table67[i-1,]$locus_tag
#   }
#   else{
#     if (gene.table67[i + 1, ]$direction == "fw") {
#       gg <- paste(gg, gene.table67[i + 1, ]$gene, sep = "")
#     }
#     if (gene.table67[i - 1, ]$direction == "rev") {
#       if(length(gg)==0){
#         gg <- paste(gg, gene.table67[i - 1, ]$gene)
#       }
#       else{
#         gg <- paste(gg, gene.table67[i - 1, ]$gene, sep = ",")
#       }
#     }
#     gene.table67[i,]$up<-gene.table67[i-1,]$gene
#     gene.table67[i,]$up_dir<-gene.table67[i-1,]$direction
#     gene.table67[i,]$up_dis<-as.numeric(gene.table67[i,]$start)-as.numeric(gene.table67[i-1,]$end)
#     gene.table67[i,]$up_locus_tag<-gene.table67[i-1,]$locus_tag
#     gene.table67[i,]$down<-gene.table67[i+1,]$gene
#     gene.table67[i,]$down_dir<-gene.table67[i+1,]$direction
#     gene.table67[i,]$down_dis<-as.numeric(gene.table67[i+1,]$start)-as.numeric(gene.table67[i,]$end)
#     gene.table67[i,]$down_locus_tag<-gene.table67[i+1,]$locus_tag
#   }
# 
#   if(length(gg)!=0){
#     gene.table67[i,]$destroy<-gg
#   }
# }
# 
# gene.table6<-rbind(gene.table6,gene.table6,gene.table6,gene.table6,gene.table65,gene.table66,gene.table67,gene.table68,gene.table69,gene.table60,gene.table61,gene.table62,gene.table63,gene.table64,gene.table65,gene.table66,gene.table67)
# write.csv(gene.table6,"C://Users//41518//Desktop//微生物//Saccharomyces//Sac.csv",row.names = FALSE)

# gene.table6<-read.csv("C://Users//41518//Desktop//微生物//Bacillus//Bacillus.csv")


# 复杂性分析
gene.table<-read.csv("C://Users/41518/Desktop/micro.csv")

for(m in 22428:24529){
  analysis_seq<-gene.table[m,]$analysis_seq
  analysis_pos <- data.frame(start=numeric(),end=numeric())
  len <- nchar(analysis_seq) - 9
  for (i in 1:len) {
    #正向重复
    pattern <- substring(analysis_seq, i, i + 9)
    pos <- matchPattern(pattern, analysis_seq)
    pos_start <- start(pos)
    pos_end <- end(pos)
    if(length(pos)>1){
      for (j in 1:length(pos_start)){
        if(pos_start[j]!=i){
          analysis_pos<-rbind(analysis_pos,c(pos_start[j],pos_end[j]))
        }
      }
    }
  }

  #反向重复
  aa=findPalindromes(DNAString(analysis_seq),min.armlength = 10,max.looplength = len)
  
  if(length(aa)!=0){
    for(l in 1:length(aa)){
      analysis_pos<-rbind(analysis_pos,c(start(palindromeLeftArm(aa[l])),end(palindromeLeftArm(aa[l]))))
      analysis_pos<-rbind(analysis_pos,c(start(palindromeRightArm(aa[l])),end(palindromeRightArm(aa[l]))))
    }
  }
  
  
  #去重复,排序
  analysis_pos<-analysis_pos[!duplicated.data.frame(analysis_pos),]
  analysis_pos<-analysis_pos[order(analysis_pos[,1],analysis_pos[,2]),]
  
  if(nrow(analysis_pos)>1){
    #合并
    for(i in 1:nrow(analysis_pos)){
      for(j in 1:nrow(analysis_pos)){
        if(analysis_pos[j,1]>=analysis_pos[i,1] & analysis_pos[j,1]<=analysis_pos[i,2]){
          if(analysis_pos[j,2]>analysis_pos[i,2]){
            analysis_pos[i,2]<-analysis_pos[j,2]
          }
          else{
            analysis_pos[j,2]<-analysis_pos[i,2]
          }
        }
        else{
          next
        }
      }
    }
    
    analysis_pos<-analysis_pos[!duplicated(analysis_pos[,2]),]
    
    gene.table[m,]$r_start<-paste(analysis_pos[,1],collapse = ";")
    gene.table[m,]$r_end<-paste(analysis_pos[,2],collapse = ";")
  }
  print(m)
}

write.csv(gene.table,"C://Users/41518/Desktop/micro2.csv")

