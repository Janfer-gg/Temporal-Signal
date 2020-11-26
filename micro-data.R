# library(stringr)
# library(Biostrings)
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

# gene.table6[NA.count1,]$seq<-substring(seq_all,as.numeric(gene.table6[NA.count1,]$start),as.numeric(gene.table6[NA.count1,]$end))
# gene.table6[NA.count1,]$analysis_seq<-substring(seq_all,as.numeric(gene.table6[NA.count1,]$start)-1000,as.numeric(gene.table6[NA.count1,]$end)+1000)

seq_all<-readLines("C://Users//41518//Desktop//微生物//Saccharomyces//sequence6.fasta")
seq_all<-seq_all[-1]
seq_all<-paste(seq_all,collapse = "")
for(i in NA.count6){
  gene.table6[i,]$seq<-substring(seq_all,as.numeric(gene.table6[i,]$start),as.numeric(gene.table6[i,]$end))
  gene.table6[i,]$analysis_seq<-substring(seq_all,as.numeric(gene.table6[i,]$start)-1000,as.numeric(gene.table6[i,]$end)+1000)
  GC<-sum(letterFrequency(DNAString(gene.table6[i,]$analysis_seq),c("G","C")))/nchar(gene.table6[i,]$analysis_seq)
  gene.table6[i,]$GC<-round(GC*100,2)
}

#发夹结构
for(m in NA.count6){
  seq<-DNAString(gene.table6[m,]$analysis_seq)
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
  gene.table6[m,]$h_start<-paste(bb$start,collapse = ";")
  gene.table6[m,]$h_end<-paste(bb$end,collapse = ";")
  
  print(m)
}

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



# write.csv(gene.table6,"C://Users//41518//Desktop//微生物//Bacillus//Bacillus-发夹.csv",row.names = FALSE)





# #发夹结构
# library(Biostrings)
# aa<-read.csv("C://Users//41518//Desktop//微生物//micro.csv")
# # aa<-transform(aa,h_start=numeric(nrow(aa)),h_end=numeric(nrow(aa)))

# for(m in 1:nrow(aa)){
#   seq<-DNAString(aa[m,]$analysis_seq)
#   hairpin<-findPalindromes(seq,min.armlength = 4,max.looplength = 17)
#   bb<-data.frame(start=start(hairpin),end=end(hairpin))
#   bb<-bb[order(bb$start),]
#   
#   for(i in 1:nrow(bb)){
#     for(j in 1:nrow(bb)){
#       if(bb[j,]$start>=bb[i,]$start & bb[j,]$start<=bb[i,]$end){
#         if(bb[j,]$end>bb[i,]$end){
#           bb[i,]$end<-bb[j,]$end
#         }
#         else{
#           bb[j,]$end<-bb[i,]$end
#         }
#       }
#       else{
#         next
#       }
#     }
#   }
#   bb<-bb[!duplicated(bb$end),]
#   aa[m,]$h_start<-paste(bb$start,collapse = ";")
#   aa[m,]$h_end<-paste(bb$end,collapse = ";")
#   
#   print(m)
# }
# 
# 
# 
# 
# write.csv(aa,"C://Users//41518//Desktop//微生物//micro1.csv",row.names = FALSE)





# seq<-readLines("C://Users//41518//Desktop//微生物//sequence.fasta")
# seq<-seq[-1]
# seq<-paste(seq,collapse = "")
# 
# seq1<-readLines("C://Users//41518//Desktop//微生物//sequence(1).fasta")
# seq1<-seq1[-1]
# seq1<-paste(seq1,collapse = "")
# 
# seq2<-readLines("C://Users//41518//Desktop//微生物//sequence(2).fasta")
# seq2<-seq2[-1]
# seq2<-paste(seq2,collapse = "")
# 
# 
# seq3<-readLines("C://Users//41518//Desktop//微生物//sequence(3).fasta")
# seq3<-seq3[-1]
# seq3<-paste(seq3,collapse = "")
# gene.table67<-data.frame(gene=c("K-12 substr. MG1655","Escherichia coli BL21(DE3)","Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S","shigella flexneri 2a str. 301"),seq=c(seq,seq1,seq2,seq3))
# 



