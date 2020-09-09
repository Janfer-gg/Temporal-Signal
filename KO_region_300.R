# 内含子小于500bp,大于300bp
KO_region_300 <- function(KO_region3) {
  region_del <- numeric()
  left <- numeric()
  right <- numeric()
  for (i in 1:nrow(KO_region3)) {
    Exon.name <- KO_region3[i,]$Exon
    j <- which(t_Exon_region_sort$Exon_name == Exon.name)
    Exon_length[i] <- KO_region3[i,]$end - KO_region3[i,]$start + 1
    {
      if(j==1){
        left[i] <- 400
        #右侧也小于500,大于300
        if (t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end < 500&
            t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end >= 300) {
          right[i] <-
            t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end - 200
        }
        #右侧大于500
        else if (t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end >= 500) {
          right[i] <- 400
        }
        else{
          region_del <- append(region_del, i)
        }
      }
      
      #左侧小于500
      else if (KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end < 500 &
          KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end >= 300) {
        left[i] <-
          KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end - 200
        if(j==nrow(t_Exon_region_sort)){
          right[i] <- 400
        }
        #右侧也小于500,大于300
        else if (t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end < 500&
            t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end >= 300) {
          right[i] <-
            t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end - 200
        }
        #右侧大于500
        else if (t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end >= 500) {
          right[i] <- 400
        }
        else{
          region_del <- append(region_del, i)
        }
      }
      
      else if(j==nrow(t_Exon_region_sort)){
        right[i] <- 400
        #左侧也小于500,大于300
        if (KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end < 500 &
            KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end >= 300) {
          left[i] <-
            KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end - 200
        }
        #左侧大于500
        else if (KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end >= 500) {
          left[i] <- 400
        }
        else{
          region_del <- append(region_del, i)
        }
      }
      
      #右侧小于500
      else if (t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end < 500 &
               t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end >= 300) {
        right[i] <-
          t_Exon_region_sort[j + 1, ]$Exon_start - KO_region3[i, ]$end - 200
        if(j==1){
          left[i]<-400
        }
        #左侧也小于500,大于300
        else if (KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end < 500 &
            KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end >= 300) {
          left[i] <-
            KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end - 200
        }
        #左侧大于500
        else if (KO_region3[i,]$start - t_Exon_region_sort[j - 1,]$Exon_end >= 500) {
          left[i] <- 400
        }
        else{
          region_del <- append(region_del, i)
        }
      }
      else{
        region_del <- append(region_del, i)
      }
    }
  }
  if (length(region_del) != 0) {
    KO_region3 <- KO_region3[-region_del,]
    left <- left[-region_del]
    right <- right[-region_del]
    Exon_length<- Exon_length[-region_del]
  }
  if(nrow(KO_region3)==0){
    return(KO_region3)
  }
  
  #KO区域排序
  KO_region3 <- cbind(KO_region3, Exon_length, left, right)
  {
    if (Gene_rev) {
      KO_region3 <- KO_region3[order(-KO_region3$times,-KO_region3$end),]
    }
    else{
      KO_region3 <- KO_region3[order(-KO_region3$times, KO_region3$end),]
    }
  }
  
  KO_region3 <- KO_region3[!duplicated(KO_region3$Exon_length), ]
  
  KO_region3 <-
    KO_region3[which(KO_region3$Exon_length %% 3 != 0), ]              #外显子非3的倍数
  KO_region3 <-
    KO_region3[which(KO_region3$Exon_length >= 100),]         #外显子不小于100bp
  
  if(nrow(Exon_CDS)!=0) {
    for (i in 1:nrow(KO_region3)) {
      for (j in 1:nrow(Exon_CDS)) {
        if (Exon_CDS[j, ]$start >= KO_region3[i, ]$start &
            Exon_CDS[j, ]$end <= KO_region3[i, ]$end) {
          KO_region3[i, ]$times <- KO_region3[i, ]$times - 1
        }
      }
    }
  }
  if(nrow(Exon_CDS)==0){
    return(KO_region3)
  }
  if (nrow(KO_region3) != 0) {
    #GC含量分析:平均GC含量大于70%或小于30%，则退出该区域
    GC_del <- numeric()
    for (i in 1:nrow(KO_region3)) {
      analysis_GC <- GC_analysis1(KO_region3[i, ])
      if (analysis_GC == TRUE) {
        GC_del <- append(GC_del, i)
      }
    }
    if (length(GC_del) != 0) {
      KO_region3 <- KO_region3[-GC_del,]
    }
    dot_del <- numeric()
    for (i in 1:nrow(KO_region3)) {
      analysis_dot1 <- Get_dot_region1(KO_region3[i, ])
      analysis_dot2 <- Get_dot_region2(KO_region3[i, ])
      analysis_dot3 <- Get_dot_region3(KO_region3[i, ])
      analysis_dot4 <- Get_dot_region4(KO_region3[i, ])
      if (analysis_dot1 == TRUE |
          analysis_dot2 == TRUE |
          analysis_dot3 == TRUE | analysis_dot4 == TRUE) {
        dot_del <- append(dot_del, i)
      }
    }
    if (length(dot_del) != 0) {
      KO_region3 <- KO_region3[-dot_del, ]
    }
  }
  return(KO_region3)
}
