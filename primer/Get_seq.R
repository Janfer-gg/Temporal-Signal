Get_seq <- function(ID) {
  server <- "https://rest.ensembl.org"
  ext <- paste0("/sequence/id/", ID, "?","expand_5prime=1500;","expand_3prime=1500")
  # r <- GET(paste(server, ext, sep = ""), content_type("text/plain"))
  r<-RETRY("GET", paste(server, ext, sep = ""))
  stop_for_status(r)
  Gene <- content(r)$seq
  return(Gene)
}