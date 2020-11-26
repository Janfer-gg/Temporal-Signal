a<-which(gene.table$up_dis<0 & gene.table$up_dis>(-200) & gene.table$pos1 < abs(gene.table$up_dis) )

b<-which(gene.table$down_dis<0 & gene.table$down_dis>(-200) & (nchar(gene.table$seq)-gene.table$pos2) < abs(gene.table$down_dis))

#重叠大于200直接删除方案
a
b
gene.table[b,]$gRNA1<-NA
gene.table[b,]$gRNA2<-NA
gene.table[b,]$strand1<-NA
gene.table[b,]$strand2<-NA
gene.table[b,]$pos1<-NA
gene.table[b,]$pos2<-NA
gene.table[b,]$score1<-NA
gene.table[b,]$score2<-NA