Get_ID<-function(term,species){
  server <- "https://rest.ensembl.org"
  ext <-paste0("/lookup/symbol/",species,"/",term,"?")
  # r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  r<-RETRY("GET", paste(server, ext, sep = ""))
  stop_for_status(r)
  ID=content(r)$id
  return(ID)
}
