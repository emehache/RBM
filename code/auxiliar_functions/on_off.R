on_off <- function(traye, de1, de2, n_datos){
  
  traye <- copy(traye)
  T <- nrow(traye)
  p <- floor((T+de1-1)/(de1+de2))
  
  # if (is.null(n_datos)) {
  #   n_datos <- p*de1
  # }
  
  if (missing(n_datos)) {
    n_datos <- p*de1
  }
  
  periodo <- rep(c(1,0), c(de1,de2)) 
  i <- rep(periodo, length.out = T)
  
  indices <- which(i==1)[seq_len(n_datos)]
  # indices <- indices[indices<=nrow(traye)]
  
  traye[,on:=0]
  traye[indices,on:=1]
  traye[,p:=c(1,diff(on))]
  traye[p==-1,p:=0]
  traye[,p:=cumsum(p)][]
  
  return(traye)
}
