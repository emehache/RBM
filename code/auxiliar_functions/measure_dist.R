
dist_medida <- function(trayectoria, ralpha) {
  cr_tray_sim <- RcppAlphahull::ahull(trayectoria, alpha = ralpha)
  
  Var1 <- runif(1e5, -1.5, 1.5)
  Var2 <- runif(length(Var1), -1, 1)
  xy <- data.table(Var1, Var2)
  
  inS <- as.logical( xy[, (Var1/a)^2 + (Var2/b)^2 < 1] * xy[, (Var1 - xc)^2 + (Var2 - yc)^2 > r^2] )
  xy <- xy[inS]
  
  inhull <- RcppAlphahull::inahull(cr_tray_sim, x = as.matrix(xy)[, 1], y = as.matrix(xy)[, 2], alpha = ralpa)
  
  areaS <- 1.5 * pi - .5^2*pi
  (1-mean(inhull)) * areaS
}
