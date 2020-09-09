KO_longregion <- function(t_Exon_CDS) {
  region_del <- numeric()
  KO_region2 <- t_Exon_CDS
  
  #标记在前30%的外显子
  l <- 0
  CDS_length <- KO_region2$end - KO_region2$start + 1
  CDS_length_30 <- floor(sum(CDS_length) * 0.33)
  for (q in 1:length(CDS_length)) {
    l <- l + CDS_length[q]
    if (l > CDS_length_30) {
      break
    }
  }
  {
    if (Gene_rev) {
      {
        if (q == 1) {
          t_CDS_30 <- KO_region2[q,]$end -  CDS_length_30
        }
        else{
          t_CDS_30 <-
            KO_region2[q,]$end -  CDS_length_30 + sum(CDS_length[1:(q - 1)])
        }
      }
      
    }
    else{
      {
        if (q == 1) {
          t_CDS_30 <- KO_region2[q,]$start +  CDS_length_30
        }
        else{
          t_CDS_30 <-
            KO_region2[q,]$start +  CDS_length_30 - sum(CDS_length[1:(q - 1)])
        }
      }
    }
  }
  for (i in 1:nrow(KO_region2)) {
    #删除起始密码子所在的外显子
    if (Gene_rev) {
      if (KO_region2[i,]$end == t_CDS_start) {
        region_del <- append(region_del, i)
        Exon_start_t <- KO_region2[i,]$start
        for (j in (i + 1):nrow(KO_region2)) {
          if (KO_region2[j,]$end + 500 > Exon_start_t) {
            region_del <- append(region_del, j)
            Exon_start_t <- KO_region2[j,]$start
          }
        }
      }
    }
    
    else{
      if (KO_region2[i,]$start == t_CDS_start) {
        region_del <- append(region_del, i)
        Exon_end_t <- KO_region2[i,]$end
        for (j in (i + 1):nrow(KO_region2)) {
          if (KO_region2[j,]$start - 500 < Exon_end_t) {
            region_del <- append(region_del, j)
            Exon_end_t <- KO_region2[j,]$end
          }
        }
      }
    }
  }
  if (length(region_del) != 0) {
    KO_region2 <- KO_region2[-region_del,]
  }
  
  k <- 1
  Exon_length <- numeric()
  for (i in 1:nrow(KO_region2)) {
    #判断该区域前后500bp是否有其他外显子，有则合并
    Exon.name <- KO_region2[i, ]$Exon
    j <- which(t_Exon_CDS_sort$Exon == Exon.name)
    Exon_length[i] <- KO_region2[i, ]$end - KO_region2[i, ]$start + 1
    #往左边合并
    ko_start <- t_Exon_CDS_sort[j, ]$start - 500
    repeat {
      if (j == 1) {
        break
      }
      if (ko_start <= t_Exon_CDS_sort[j - 1, ]$end) {
        ko_start <- t_Exon_CDS_sort[j - 1, ]$start - 500
        # if()
        Exon_length[i] <-
          Exon_length[i] + t_Exon_CDS_sort[j - 1, ]$end - t_Exon_CDS_sort[j - 1, ]$start +
          1
        j <- j - 1
        if (j == 1) {
          break
        }
      }
      else{
        break
      }
    }
    
    j <- which(t_Exon_CDS_sort$Exon == Exon.name)
    if (j == nrow(t_Exon_CDS_sort)) {
      break
    }
    
    #往右边合并
    ko_end <- t_Exon_CDS_sort[j, ]$end + 500
    repeat {
      if (ko_end > t_Exon_CDS_sort[j + 1, ]$start) {
        ko_end <- t_Exon_CDS_sort[j + 1, ]$end + 500
        Exon_length[i] <-
          Exon_length[i] + t_Exon_CDS_sort[j + 1, ]$end - t_Exon_CDS_sort[j + 1, ]$start +
          1
        j <- j + 1
        if (j == nrow(t_Exon_CDS_sort)) {
          break
        }
      }
      else{
        break
      }
    }
    
    ko_start<-ko_start + 500
    ko_end<-ko_end - 500
    
    if(ko_start<min(t_Exon_CDS$start)){
      ko_start <- min(t_Exon_CDS$start)
    }
    if(ko_start>max(t_Exon_CDS$end)){
      ko_end<-max(t_Exon_CDS$end)
    }
    
    KO_region2[i, ]$start <- ko_start 
    KO_region2[i, ]$end <- ko_end

  }
  
  KO_region2 <- cbind(KO_region2, Exon_length)
  KO_region2 <- KO_region2[!duplicated(KO_region2$Exon_length), ]
  KO_region4 <- data.frame()
  KO_region5 <- data.frame()
  
  {
    if (Gene_rev) {
      #选起点（这里有多个起点选择）
      for (i in 1:nrow(KO_region2)) {
        if (KO_region2[i, ]$end >= t_CDS_30) {
          KO_region4 <- KO_region2[i, ]
          length <- KO_region4$Exon_length
          if (nrow(KO_region4) == 0) {
            return(KO_region4)
          }
          
          #如果只有一个区域
          if (nrow(KO_region2) == 1) {
            if (length %% 3 != 0) {
              KO_region4$start <- KO_region2[1, ]$start
              KO_region4$Exon_length <- length
              # if (KO_region4$start == t_Exon_CDS_sort[1, ]$Exon_start) {
              #   KO_region4$start <- t_Exon_CDS[nrow(t_Exon_CDS), ]$start
              # }
              return(KO_region4)
            }
          }
          
          #加起来非3的倍数且大于3kb，小于10kb
          else{
            length<-0
            for (j in i:nrow(KO_region2)) {
              length <- length + KO_region2[j,]$Exon_length
              if (length %% 3 != 0) {
                KO_region4$start <- KO_region2[j, ]$start
                if (KO_region4$end - KO_region4$start >= 3000 &
                    KO_region4$end - KO_region4$start <= 20000) {
                  KO_region4$Exon_length <-
                    sum(KO_region2[i:j, ]$Exon_length)
                  # if (KO_region4$start == t_Exon_CDS_sort[1, ]$Exon_start) {
                  #   KO_region4$start <- t_Exon_CDS[nrow(t_Exon_CDS), ]$start
                  # }
                  KO_region5 <- rbind(KO_region5, KO_region4)
                }
              }
            }
          }
        }
      }
    }
    
    else{
      #选起点（这里有多个起点选择）
      for (i in 1:nrow(KO_region2)) {
        if (KO_region2[i, ]$start <= t_CDS_30) {
          KO_region4 <- KO_region2[i, ]
          length <- KO_region4$Exon_length
          if (nrow(KO_region4) == 0) {
            return(KO_region4)
          }
          
          #如果只有一个区域
          if (nrow(KO_region2) == 1) {
            length <- KO_region4$Exon_length
            if (length %% 3 != 0) {
              KO_region4$end <- KO_region2[1, ]$end
              KO_region4$Exon_length <- length
              # if (KO_region4$end == t_Exon_region[nrow(t_Exon_region), ]$end) {
              #   KO_region4$end <- t_Exon_CDS[nrow(t_Exon_CDS), ]$end
              # }
              return(KO_region4)
            }
          }
          
          #加起来非3的倍数且大于3kb，小于10kb
          else{
            length<-0
            for (j in i:nrow(KO_region2)) {
              length <- length + KO_region2[j, ]$Exon_length
              if (length %% 3 != 0) {
                KO_region4$end <- KO_region2[j, ]$end
                if (KO_region4$end - KO_region4$start >= 3000 &
                    KO_region4$end - KO_region4$start <= 20000) {
                  KO_region4$Exon_length <-
                    sum(KO_region2[i:j, ]$Exon_length)
                  #末尾
                  # if (KO_region4$end == t_Exon_region[nrow(t_Exon_region), ]$Exon_end) {
                  #   KO_region4$end <- t_Exon_CDS[nrow(t_Exon_CDS), ]$end
                  # }
                  KO_region5 <- rbind(KO_region5, KO_region4)
                }
              }
            }
          }
        }
      }
    }
  }
  if (max(t_Exon_CDS$end) - min(t_Exon_CDS$start) <= 20000) {
    KO_region6 <- t_Exon_CDS[1, ]
    KO_region6$start <- min(t_Exon_CDS$start)
    KO_region6$end <- max(t_Exon_CDS$end)
    Exon_length <- sum(t_Exon_CDS$end-t_Exon_CDS$start+1)
    KO_region6<-cbind(KO_region6,Exon_length)
    KO_region5<-rbind(KO_region5,KO_region6)
  }
  return(KO_region5)
}
