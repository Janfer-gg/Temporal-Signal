#下载基因的所有信息
Get_allinfo<-function(ID){
  server <- "https://rest.ensembl.org"
  ext <- paste0("/lookup/id/", ID, "?expand=1")
  r <-
    GET(paste(server, ext, sep = ""),
        content_type("application/json"))
  stop_for_status(r)
  allinfo <- content(r)
  return(allinfo)
} 
