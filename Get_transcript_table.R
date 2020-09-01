#下载转录本表格
Get_transcript_table<-function(ID,species){
  IMUFE_url <-
    paste0("https://asia.ensembl.org/",species,"/Gene/Summary?db=core;g=",
           ID)
  allSourceCode <- IMUFE_url %>%
    read_html(encoding = "UTF-8")
  allSourceCode <- allSourceCode %>%
    html_node("div[class=transcripts_table]")
  t <- allSourceCode %>%
    html_node("table") %>%
    html_table()
  return(t)
}
