#下载基因组序列
Get_seq <- function(ID) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", ID, "?")
  # r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  r<-RETRY("GET", paste(server, ext, sep = ""))
  stop_for_status(r)
  Gene <- content(r)$seq
  return(Gene)
}

Get_seq2 <- function(ID) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", ID, "?","expand_5prime=500;","expand_3prime=500")
  # r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  r<-RETRY("GET", paste(server, ext, sep = ""))
  stop_for_status(r)
  Gene <- content(r)$seq
  return(Gene)
}

Get_seq3 <- function(ID) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", ID, "?","expand_5prime=1200;","expand_3prime=1200")
  # r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  r<-RETRY("GET", paste(server, ext, sep = ""))
  stop_for_status(r)
  Gene <- content(r)$seq
  return(Gene)
}
