# library(RSQLite)
# con <- dbConnect(SQLite(),"test.db")
# dbWriteTable(con, "test", Mydata)
# Mydata2<-dbReadTable(con,"test")
library(openxlsx)

build_carrier<-function(type,scaffoldname,genename,gRNA1,gRNA2=NULL,number=1){
  if(type=="双gRNA"){
    if(is.null(gRNA2)){
      stop("please input gRNA2")
    }
    Mydata<-read.xlsx("敲除载体构建方案.xlsx",sheet=2)
    Mydata$载体类型<-type
    for (i in 1:nrow(Mydata)) {
      if(is.na(Mydata[i,]$骨架信息)){
        Mydata[i,]$骨架信息<-Mydata[i-1,]$骨架信息
        Mydata[i,]$检菌引物<-Mydata[i-1,]$检菌引物
        Mydata[i,]$测序引物<-Mydata[i-1,]$测序引物
      }
    }
    Mydata2<-Mydata[which(Mydata$骨架信息==scaffoldname),]
    Mydata2$F引物 <- c(paste0(genename,"[gRNA",number,"]-F1"),paste0(genename,"[gRNA1",number,"]-F2"))
    Mydata2$R引物 <- c(paste0(genename,"[gRNA2",number,"]-R1"),paste0(genename,"[gRNA2",number,"]-R2"))
    Mydata2[1,]$X5 <- paste0("g",gRNA1,"GTTTTAGAGCTAGAAATAGCAAGTT")
    Mydata2[2,]$X5 <- paste0("TGTGGAAAGGACGAAACACC","g",gRNA1)
    g <- as.character(reverseComplement(DNAString(gRNA2)))
    Mydata2[1,]$X7 <- paste0(g,"c","GGTGTTTCGTCCTTTCCACAAG")
    Mydata2[2,]$X7 <- paste0("ACTTGCTATTTCTAGCTCTAAAAC",g)
    #return(Mydata2)
    #carrier.table<-read.xlsx("C:\\Users\\41518\\Desktop\\测试\\敲除载体构建双gRNA.xlsx")
    #n<-nrow(carrier.table)
    #carrier.table[n+1,] <- Mydata2
    table <- Mydata2
    write.xlsx(table,"C:\\Users\\41518\\Desktop\\测试\\敲除载体构建双gRNA.xlsx")
  }
  
  if (type == "单gRNA") {
    Mydata <- read.xlsx("敲除载体构建方案.xlsx", sheet = 1)
    
    Mydata$载体类型 <- type
    for (i in 1:nrow(Mydata)) {
      if (is.na(Mydata[i, ]$骨架信息)) {
        Mydata[i, ]$骨架信息 <- Mydata[i - 1, ]$骨架信息
        Mydata[i, ]$检菌引物 <- Mydata[i - 1, ]$检菌引物
        Mydata[i, ]$测序引物 <- Mydata[i - 1, ]$测序引物
      }
    }
    Mydata2 <- Mydata[which(Mydata$骨架信息 == scaffoldname), ]
    Mydata2$`Oligo-F` <- c(paste0(genename, "[gRNA", number, "]-F"))
    Mydata2$`Oligo-R` <- c(paste0(genename, "[gRNA", number, "]-R"))
    Mydata2$X4 <- paste0("CACCg", gRNA1)
    names(Mydata2)[4]<-"Oligo-F sequence"
    g <- as.character(reverseComplement(DNAString(paste0("g", gRNA1))))
    Mydata2$X6 <- paste0("AAAC", g)
    names(Mydata2)[6]<-"Oligo-R sequence"
    carrier.table<-read.xlsx("C:\\Users\\41518\\Desktop\\测试\\敲除载体构建单gRNA.xlsx")
    n<-nrow(carrier.table)
    carrier.table[n+1,] <- Mydata2
    #table <- Mydata2
    write.xlsx(carrier.table, "C:\\Users\\41518\\Desktop\\测试\\敲除载体构建单gRNA.xlsx")
  }
}

build_carrier(type = "双gRNA",scaffoldname = "YKO-RP005",genename = "Cd44 ",gRNA1 ="AATGTAACCTGCCGCTACGC",gRNA2 ="CAACTTCATTTGGTCCATGG" ,number = 1)





