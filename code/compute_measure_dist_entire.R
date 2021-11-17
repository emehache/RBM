## ----setup,echo=F,warning=F,message=F-----------------------------------------------------------------------------------------
paquetes <- c("data.table","magrittr")
sapply(paquetes,require,character.only=TRUE)

# setwd("~/Dropbox/Documentos/Tesis Ingemat/RBM Manuel/manuel")
setwd("~/Dropbox/RBM/manuel")

# trayectorias <- readRDS('objetos_r/sim.rds')
# trayectorias <- readRDS('objetos_r/sim_replicas.rds')
# 
# sapply(names(trayectorias), function(nm){
#   readr::write_rds(trayectorias[[nm]], sprintf('objetos_r/trayectorias/trayectoria_%s.rds',nm))
#   })
# 
# rm(trayectorias)

source('code/auxiliar_functions/measure_dist.R')
source('code/auxiliar_functions/on_off.R')


# -------------------------------------------------------------------------

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

h <- c(0.001,0.002,0.003)
N <- 1e5

param <- setDT(expand.grid(N=N,h=h))
param[,param:=apply(param,1,paste,collapse=', ')]

# distancia en medida --------------------------------------------------

set.seed(1234)
m <- 1e5
cuadrado <- cbind(runif(m,-1.5,1.5),runif(m,-1,1))
index <- (cuadrado[,1]-xc)^2+(cuadrado[,2]-yc)^2 >= r^2
index2 <- (cuadrado[,1])^2/a^2+(cuadrado[,2])^2/b^2 <= 1
conjunto <- cuadrado[as.logical(index*index2),]

ralpha <- .4

archivos <- list.files('data/trajectories', full.names = TRUE)

t0 <- Sys.time()
sapply(archivos, function(ar){
  
  param_on_off <- setDT(expand.grid(de1=c(100,250,500),de2=c(100,250,500)))
  
  h <- tstrsplit(gsub('^\\D*|.rds','',ar),', ')[[2]] %>% as.numeric
  
  t <- readRDS(ar)
  
  res <- lapply(t, function(tray){
    
    TT <- nrow(tray)
    apply(param_on_off[,.(de1,de2)],1,function(par) {
      
      de1 <- par[1]; de2 <- par[2]
      p <- floor((TT+de1-1)/(de1+de2))
      length <- p*de1
      
      cat(N, h, de1, p, '\n')
      data.table(N, h, TT, de1, de2, p, length = p*de1, dist = dist_medida(tray[seq_len(length)], ralpha))
      
    }) %>% rbindlist
  }) %>% rbindlist(idcol = 'rep')
  
  
  readr::write_rds(res, sprintf('r_objects/measure_distances/entire_model_%s.rds',h))
  c
})
t1 <- Sys.time()
t1-t0

