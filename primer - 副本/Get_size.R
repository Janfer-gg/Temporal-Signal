Get_size<-function(length){
  {
    if (length <= 2000) {
      size<-4000
    }
    else if (length <= 8000){
      size<-10000
    }
    else{
      size<-20000
    }
  }
  return(size)
}
