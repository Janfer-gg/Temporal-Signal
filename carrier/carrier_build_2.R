build_carrier1<-function(access,carname,number,gRNA1,scaffold,account,filepath1,filepath2){
  library(openxlsx)
  library(Biostrings)
  
  Mydata<-read.xlsx(filepath2,sheet=1,fillMergedCells=TRUE)
  Mydata2<-data.frame()
  for(i in 1:account){
    Mydata3<-Mydata[which(Mydata[,2]==scaffold[i]),]
    if(nrow(Mydata3)==0){
      return(FALSE)
    }
    Mydata2[i,1] <- access[i]
    Mydata2[i,2] <- carname[i]
    Mydata2[i,3] <- number[i]
    Mydata2[i,4] <- gRNA1[i]
    Mydata2[i,5] <- scaffold[i]
    Mydata2[i,6] <- paste0(number[i],"-F")
    Mydata2[i,7] <- paste0("CACCg", gRNA1[i])
    Mydata2[i,8] <- paste0(number[i],"-R")
    g <- as.character(reverseComplement(DNAString(paste0("g", gRNA1[i]))))
    Mydata2[i,9] <- paste0("AAAC", g)
    Mydata2[i,10] <-Mydata3[,7]
    Mydata2[i,11] <-Mydata3[,8]
    Mydata2[i,12] <-Mydata3[,9]
    names(Mydata2)<-c("协议号","载体名称","内部编号","gRNA序列","骨架信息","Oligo-F","Oligo-F","Oligo-R","Oligo-R","检菌引物","测序引物","酶切方案")
  }
  
  wb <- createWorkbook()
  addWorksheet(wb, 1)
  writeData(wb,1,Mydata2,borders = "all")
  
  mergeCells(wb, 1, rows =1,cols =6:7)
  mergeCells(wb, 1, rows =1,cols =8:9)
  
  setColWidths(wb,1,cols = 1, widths = "20")
  setColWidths(wb,1,cols = 2, widths = "32")
  setColWidths(wb,1,cols = 3, widths = "12")
  setColWidths(wb,1,cols = 4, widths = "28")
  setColWidths(wb,1,cols = 5, widths = "12")
  setColWidths(wb,1,cols = c(6,8), widths = "16")
  setColWidths(wb,1,cols = c(7,9), widths = "32")
  setColWidths(wb,1,cols = 10, widths = "30")
  setColWidths(wb,1,cols = 11, widths = "14")
  setColWidths(wb,1,cols = 12, widths = "38")
  setRowHeights(wb,1,rows = 1:(account+1),heights = "24")
  addStyle(wb,1,style = createStyle(fontSize = 11, halign = "center",valign = "center",border =c("top", "bottom", "left", "right")),rows=1,cols = 1:12)
  addStyle(wb,1,style = createStyle(fontSize = 11, valign = "center",border =c("top", "bottom", "left", "right")),rows=2:(account+1),cols = 1:12,gridExpand=TRUE)
  addStyle(wb,1,style = createStyle(fontSize = 11, halign = "center",valign = "center",border =c("top", "bottom", "left", "right")),rows=2:(account+1),cols = 1:5,gridExpand=TRUE)
  
  saveWorkbook(wb, paste0(filepath1,"//carrier.xlsx"), overwrite = TRUE)
  return(TRUE)
}

build_carrier2<-function(access,carname,number,gRNA1,gRNA2,scaffold,account,filepath1,filepath2){
  library(openxlsx)
  library(Biostrings)
  TM<-function(seq){
    GC_count<-length(gregexpr("[GCgc]",seq)[[1]])
    tm<-100.5 + (41*(GC_count)/nchar(seq)) - (820/nchar(seq)) + 16.6*log(0.05,10)
    return(tm)
  }
  Mydata<-read.xlsx(filepath2,sheet=2,fillMergedCells=TRUE)
  Mydata2<-data.frame()
  for(i in 1:account){
    Mydata3<-Mydata[which(Mydata[,2]==scaffold[i]),]
    if(nrow(Mydata3)==0){
      return(FALSE)
    }
    Mydata2[(i*2-1):(i*2),1] <- access[i]
    Mydata2[(i*2-1):(i*2),2] <- carname[i]
    Mydata2[(i*2-1):(i*2),3] <- number[i]
    Mydata2[(i*2-1):(i*2),4] <- gRNA1[i]
    Mydata2[(i*2-1):(i*2),5] <- gRNA2[i]
    Mydata2[(i*2-1):(i*2),6] <- scaffold[i]
    Mydata2[i*2-1,7] <- "一轮PCR"
    Mydata2[i*2,7] <- "二轮PCR"
    Mydata2[i*2-1,8] <- paste0(number[i],"-F1")
    Mydata2[i*2,8] <- paste0(number[i],"-F2")
    Mydata2[i*2-1,9] <- paste0("g",gRNA1[i],"GTTTTAGAGCTAGAAATAGCAAGTT")

    Mydata2[i*2-1,10] <- paste0(number[i],"-R1")
    Mydata2[i*2,10] <- paste0(number[i],"-R2")
    g <- as.character(reverseComplement(DNAString(gRNA2[[i]])))
    Mydata2[i*2-1,11] <- paste0(g,"c","GGTGTTTCGTCCTTTCCACAAG")

    Mydata2[i*2-1,12] <- "合成基因scaffold-U6"
    Mydata2[i*2,12] <- paste0(number[i],"-1")
    Mydata2[i*2-1,13] <- paste0(number[i],"-1 (381bp)")
    Mydata2[i*2,13] <- paste0(number[i],"-2 (420bp)")
    Mydata2[(i*2-1):(i*2),14] <- Mydata3[,10]
    Mydata2[(i*2-1):(i*2),15] <- Mydata3[,11]
    Mydata2[(i*2-1):(i*2),16] <- Mydata3[,12]

    seq1<-paste0("g",gRNA1[[i]],"gttttaga")
    seq2<-paste0(g,"cggtgtttc")
    m=n=21

    while(TM(substring(seq1,1,m))>70){
      m<-m-1
      if(TM(substring(seq1,1,m))<=70|m==18){
        break
      }
    }
    while(TM(substring(seq1,1,m))<60){
      m<-m+1
      if(TM(substring(seq1,1,m))>=60|m==29){
        break
      }
    }
    while(TM(substring(seq2,1,n))>70){
      n<-n-1
      if(TM(substring(seq2,1,n))<=70|n==18){
        break
      }
    }
    while(TM(substring(seq2,1,n))<60){
      n<-n+1
      if(TM(substring(seq2,1,n))>=60|n==29){
        break
      }
    }


    while(TM(substring(seq1,1,m))-TM(substring(seq2,1,n))>6){
      if(n<29){
        n<-n+1
        if(TM(substring(seq1,1,m))-TM(substring(seq2,1,n))<=6){
          break
        }
      }
      if(m>18){
        m<-m-1
        if(TM(substring(seq1,1,m))-TM(substring(seq2,1,n))<=6){
          break
        }
      }
      if(n==29&m==18){
        break
      }
    }

    while(TM(substring(seq2,1,n))-TM(substring(seq1,1,m))>6){
      if(n>18){
        n<-n-1
        if(TM(substring(seq2,1,n))-TM(substring(seq1,1,m))<=6){
          break
        }
      }
      if(m<29){
        m<-m+1
        if(TM(substring(seq2,1,n))-TM(substring(seq1,1,m))<=6){
          break
        }
      }
      if(m==29&n==18){
        break
      }
    }


    Mydata2[i*2,9] <- paste0("TGTGGAAAGGACGAAACACC",substring(seq1,1,m))
    Mydata2[i*2,11] <- paste0("GCTATTTCTAGCTCTAAAAC",substring(seq2,1,n))
    
    names(Mydata2)<-c("协议号","载体名称","内部编号","gRNA1序列","gRNA2序列","骨架信息","扩增轮次","Oligo-F","Oligo-F","Oligo-R","Oligo-R","PCR模板","片段名称(大小)","检菌引物","测序引物","酶切方案")
  }
  
  wb <- createWorkbook()
  addWorksheet(wb, 1)
  writeData(wb,1,Mydata2,borders = "all")
  mergeCells(wb, 1, rows=1, cols =8:9)
  mergeCells(wb, 1, rows=1, cols =10:11)
  for(i in 1:account){
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =1)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =2)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =3)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =4)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =5)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =6)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =14)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =15)
    mergeCells(wb, 1, rows=(i*2):(i*2+1), cols =16)
  }
  
  setColWidths(wb,1,cols = 1, widths = "20")
  setColWidths(wb,1,cols = 2, widths = "36")
  setColWidths(wb,1,cols = 3, widths = "12")
  setColWidths(wb,1,cols = 4:5, widths = "28")
  setColWidths(wb,1,cols = 6, widths = "12")
  setColWidths(wb,1,cols = 7, widths = "10")
  setColWidths(wb,1,cols = c(8,10), widths = "16")
  setColWidths(wb,1,cols = c(9,11), widths = "32")
  setColWidths(wb,1,cols = 12, widths = "20")
  setColWidths(wb,1,cols = 13, widths = "20")
  setColWidths(wb,1,cols = 14, widths = "30")
  setColWidths(wb,1,cols = 15, widths = "24")
  setColWidths(wb,1,cols = 16, widths = "38")
  setRowHeights(wb,1,rows = 1:(account*2+1),heights = "24")
  
  addStyle(wb,1,style = createStyle(fontSize = 11, halign = "center",valign = "center",border =c("top", "bottom", "left", "right")),rows=1,cols = 1:16)
  addStyle(wb,1,style = createStyle(fontSize = 11, halign = "center",valign = "center",border =c("top", "bottom", "left", "right")),rows=2:(account*2+1),cols = 1:7,gridExpand=TRUE)
  addStyle(wb,1,style = createStyle(fontSize = 11, valign = "center",border =c("top", "bottom", "left", "right")),rows=2:(account*2+1),cols = 8:16,gridExpand=TRUE)
  
  saveWorkbook(wb, paste0(filepath1,"//carrier.xlsx"), overwrite = TRUE)
  return(TRUE)
}


build_carrier3<-function(access,carname,number,gRNA1,scaffold,account,filepath1,filepath2){
  library(openxlsx)
  library(Biostrings)
  
  Mydata<-read.xlsx(filepath2,sheet=3,fillMergedCells=TRUE)
  Mydata2<-data.frame()
  for(i in 1:account){
    Mydata3<-Mydata[which(Mydata[,2]==scaffold[i]),]
    if(nrow(Mydata3)==0){
      return(FALSE)
    }
    Mydata2[i,1] <- access[i]
    Mydata2[i,2] <- carname[i]
    Mydata2[i,3] <- number[i]
    Mydata2[i,4] <- gRNA1[i]
    Mydata2[i,5] <- scaffold[i]
    Mydata2[i,6] <- paste0(number[i],"-F")
    g <- as.character(reverseComplement(DNAString(gRNA1[i])))
    Mydata2[i,7] <- paste0("CCGG",gRNA1[i],"CTCGAG",g,"TTTTTG")
    Mydata2[i,8] <- paste0(number[i],"-R")
    
    Mydata2[i,9] <- paste0("AATTCAAAAA",gRNA1[i],"CTCGAG",g)
    Mydata2[i,10] <-Mydata3[,7]
    Mydata2[i,11] <-Mydata3[,8]
    Mydata2[i,12] <-Mydata3[,9]
    names(Mydata2)<-c("协议号","载体名称","内部编号","shRNA序列","骨架信息","Oligo-F","Oligo-F","Oligo-R","Oligo-R","检菌引物","测序引物","酶切方案")
  }
  
  wb <- createWorkbook()
  addWorksheet(wb, 1)
  writeData(wb,1,Mydata2,borders = "all")
  
  mergeCells(wb, 1, rows =1,cols =6:7)
  mergeCells(wb, 1, rows =1,cols =8:9)
  
  setColWidths(wb,1,cols = 1, widths = "20")
  setColWidths(wb,1,cols = 2, widths = "32")
  setColWidths(wb,1,cols = 3, widths = "12")
  setColWidths(wb,1,cols = 4, widths = "28")
  setColWidths(wb,1,cols = 5, widths = "12")
  setColWidths(wb,1,cols = c(6,8), widths = "16")
  setColWidths(wb,1,cols = c(7,9), widths = "32")
  setColWidths(wb,1,cols = 10, widths = "30")
  setColWidths(wb,1,cols = 11, widths = "20")
  setColWidths(wb,1,cols = 12, widths = "38")
  setRowHeights(wb,1,rows = 1:(account+1),heights = "24")
  addStyle(wb,1,style = createStyle(fontSize = 11, halign = "center",valign = "center",border =c("top", "bottom", "left", "right")),rows=1,cols = 1:12)
  addStyle(wb,1,style = createStyle(fontSize = 11, valign = "center",border =c("top", "bottom", "left", "right")),rows=2:(account+1),cols = 1:12,gridExpand=TRUE)
  addStyle(wb,1,style = createStyle(fontSize = 11, halign = "center",valign = "center",border =c("top", "bottom", "left", "right")),rows=2:(account+1),cols = 1:5,gridExpand=TRUE)
  
  saveWorkbook(wb, paste0(filepath1,"//carrier.xlsx"), overwrite = TRUE)
  return(TRUE)
}

