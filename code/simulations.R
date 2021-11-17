#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  args <- c(50, 1e5)
} else if (length(args)==1) {
  args[2] <- 1e5
}

args <- as.numeric(args)


## ----setup,echo=F,warning=F,message=F-----------------------------------------------------------------------------------------
paquetes <- c("data.table","magrittr", 'Rcpp')
sapply(paquetes,require,character.only=T)

# setwd("~/Dropbox/IngeMat/RBM/manuel")
sourceCpp('rbmd.cpp')

# parametros --------------------------------------------------------------

a<-1.5 # eje mayor elipse
b<-1 # eje menor elipse
r<-.5 # radio del circulo
xc<-0.8 # abscisa del centro del círculo
yc<-0 # ordenadad del centro del círculo
ini <- c(0,0) # inicio del proceso

th<-seq(0,2*pi,len=1000)
mg<-0.0000001 # margin for the plot
elipse <- cbind(a * cos(th),b * sin(th))
circulo <- cbind((r-mg)*cos(th)+xc,(r-mg)*sin(th)+yc)

# importante, el h tiene que tener relacion con el radio del circulo
h <- c(0.001,0.002,0.003)
# N <- 1e5
N <- args[2]

param <- expand.grid(N=N,h=h)

# n_rep <- 50
n_rep <- args[1]

cat('Numero de replicas: ', n_rep, '\n')
cat('Largo de cadena: ', N, '\n')

set.seed(1234)
t0 <- Sys.time()
trayectorias <- 
  apply(param,1,function(par){
    cat(par[2], '\n')
    replicate(n_rep, {
      rbmd_cpp(N=par[1],
               h=par[2],
               a=a, b=b, r=r,
               xc=xc, yc=yc,
               sigma=1,
               ini =ini,
               eps=10^-6) %>% data.table
    }, simplify = F)
  })

names(trayectorias) <- apply(param,1,paste,collapse='_')

invisible(sapply(names(trayectorias), function(nm){
  readr::write_rds(trayectorias[[nm]], sprintf('objetos_r/trayectorias/trayectoria_%s.rds',nm))
}))

t1 <- Sys.time()
t1 - t0
