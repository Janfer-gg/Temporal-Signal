setwd("C://Users/41518/Desktop/work/Ubigene/载体测序/")
library(sangerseqR)
library(Biostrings)

filepath1<-'C://Users/41518/Desktop/ICE/run_exprisement(1)/CK19-004-Npnt/PK200908-02-A05-PK593.ab1'
filepath2<-'C://Users/41518/Desktop/载体/载体测序结果/合格/YSH-LV001-hREST[shRNA1]/YSH-LV001-hREST[shRNA1].dna'
#读入ab1文件和dna文件
info1 = readsangerseq(filepath1)
info= read.abif(filepath1)
ref_dna = readLines(filepath2)


#从ab1文件获取测序序列
seq<-toString(info1@primarySeq)


{
  #全是N的无信号序列
  if(all(unlist(strsplit(seq,""))=="N")){
    print("无信号")
    result<-data.frame(seq=unlist(strsplit(seq,"")),fill="#2788ff",color="black")
    write.csv(result,"result.csv",row.names = FALSE)
    write.table("无信号","result.txt",row.names = FALSE,col.names = FALSE)
  }
  
  else{
    #从dna文件获取理论序列
    ref_seq<-ref_dna[grep("ORIGIN",ref_seq)+1]
    ref_seq = paste0(ref_seq,ref_seq)
    #从ab1文件获取质量分数
    score<-info@data[["PCON.1"]]
    
    #测序序列的峰值和位置
    peek_pos=info@data[["PLOC.1"]]+1
    
    #找双峰
    bc = makeBaseCalls(info1, ratio = 0.33)
    pa <- pairwiseAlignment(primarySeq(bc), secondarySeq(bc), type="global-local")
    #主峰和次峰错配的位置（双峰的位置）
    doublepeak<-pa@pattern@mismatch@unlistData
    
    peakPos<-as.data.frame(peakPosMatrix(bc))
    peakAMP<-as.data.frame(peakAmpMatrix(bc))
    colnames(peakPos)<-c("A","C","G","T")
    colnames(peakAMP)<-c("A","C","G","T")
    
    #找到原序列上的双峰的位置
    doublepeak_pos<-numeric()
    for(i in doublepeak){
      primary_nt<-toString(primarySeq(bc)[i])
      secondary_nt<-toString(secondarySeq(bc)[i])
      nt_pos1<-peakPos[i,primary_nt]
      nt_pos2<-peakPos[i,secondary_nt]
      primary_pos<-which.min(abs(peek_pos-nt_pos1))
      secondary_pos<-which.min(abs(peek_pos-nt_pos2))
      {
        if(primary_pos==secondary_pos){
          doublepeak_pos<-append(doublepeak_pos,primary_pos)
        }
        else{
          print("FALSE")
        }
      }
    }
    
    #测序序列与理论序列进行比对
    alignment<-pairwiseAlignment(DNAString(seq),DNAString(ref_seq),type="global-local",substitutionMatrix=nucleotideSubstitutionMatrix(3,2,TRUE),gapOpening=5, gapExtension=3)
    writePairwiseAlignments(alignment,block.width = 60)
    chromatogram(bc, width = 100, height = 2, showcalls = "both", filename = "chromatogram.pdf")
    
    # 完成比对后的序列
    aligned_seq1<-as.character(alignedPattern(alignment))
    aligned_seq2<-as.character(alignedSubject(alignment))
    write.table(aligned_seq2,"control.txt",quote = FALSE,row.names = FALSE,col.names = FALSE)
   
    
    #建立原序列位置映射到比对后序列位置的索引
    i<-1
    index<-data.frame(pos1=numeric(),pos2=numeric())
    for (j in (1:nchar(aligned_seq1))) {
      if(substring(aligned_seq1,j,j)!="-"){
        index[i,]<-c(i,j)
        i<-i+1
      }
      else{
        next
      }
    }
    
    #在序列上标出质量低，双峰，突变的位置
    mutable<-alignment@pattern@mismatch@unlistData       #突变
    indel<-gregexpr("-",aligned_seq1)[[1]]             #缺失
    if(length(indel)==1){
      if(indel==(-1)){
        indel<-numeric()
      }
    }
    
    insert<-gregexpr("-",aligned_seq2)[[1]]           #插入
    for(m in 1:length(mutable)){
      mutable[m]<-mutable[m]+length(which((mutable[m]>indel)==TRUE))
    }
    
    #通过原序列位置和索引找出在比对后序列上的位置
    peak_pos<-sapply(doublepeak_pos,function(x) index[which(index$pos1==x),]$pos2)    #双峰位置
    
    low_pos<-index[which(score<=20),]$pos2                #质量得分低于20
    middle_pos<-index[which(score<40&score>20),]$pos2     #质量得分低于40
    
    result<-data.frame(seq=character(),fill=character(),color=character())
    for (n in 1:nchar(aligned_seq1)) {
      
      result[n,1]<-substring(aligned_seq1,n,n)
      {

        if(n %in% low_pos){
          result[n,2]<-"#2788ff"
        }
        else if(n %in% middle_pos){
          result[n,2]<-"#16b71e"
        }
        else{
          result[n,2]<-"#cccccc"
        }
      }
      {
        if(n %in% peak_pos){
          result[n,3]<-"ff407a"
        }
        else if(n %in% mutable | n %in% indel | n %in% insert){
          result[n,3]<-"fed24c"
        }
        else{
          result[n,3]<-"black"
        }
        
      }
    }
    
    write.csv(result,"result.csv",row.names = FALSE)
    
    
    {
      if(grepl("S175",filepath1)){
        result2<-result[61:700,]
      }
      else{
        result2<-result[61:800,]
      }
    }

    {
      if(all(result2$color=="white")){
        print("合格")
        write.table("合格","result.txt",row.names = FALSE)
      }
      else{
        if(any(result2$color=="red")==TRUE & any(result2$color=="orange")==TRUE){
          print("双峰和错配都有")
        }
        else if(any(result2$color=="red")){
          print("双峰")
          # write.table("有双峰","result.txt",row.names = FALSE,col.names = FALSE)
        }
        write.table("不合格","result.txt",row.names = FALSE,col.names = FALSE)
      }
    }
  }
}





# data<-cbind(info@data[["DATA.10"]],info@data[["DATA.12"]],info@data[["DATA.9"]],info@data[["DATA.11"]])
# colnames(data)<-c("A","C","G","T")


