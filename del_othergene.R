#某个区域有无其他基因重叠
Get_othergene<-function(species,chr,start,end){
  chr<-as.numeric(chr)
  server <- "https://rest.ensembl.org"
  #ext <- "/overlap/region/human/6:11093834-11138733?feature=gene"
  ext <- paste0("/overlap/region/",species,"/",chr,":",start,"-",end,"?feature=gene")
  r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  a<-head(fromJSON(toJSON(content(r))))
  return(a)
}
