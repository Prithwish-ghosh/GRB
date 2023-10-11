library(Directional)
library(circular)
library(CircStats)
library(CircSpaceTime)
library(readr)
library(ggplot2)
library(SPADAR)

grb = read_csv('GRB.csv')
head(grb)
grb = data.frame(grb)
watson.test(grb$X_RA, alpha = 0.01, dist = "vonmises")
watson.test(grb$X_DE, alpha = 0.01, dist = "vonmises")

library(movMF)

d = cbind(grb$X_DE, grb$X_RA)
Evmf <- function(K){
  movMF(d, k= K, control = list(nruns = 20))
}

set.seed(122)
Esd = lapply(1:10, Evmf)
Esd
sapply(Esd, BIC)



library(ggplot2)
library(cooltools)

x= grb$X_RA
y <- rvonmises(n=1000, mu=mean(grb$X_RA), kappa=est.kappa(grb$X_RA))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(grb$X_RA , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Gamma Ray Burst RA data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

x= grb$X_DE
y <- rvonmises(n=1000, mu=mean(grb$X_DE), kappa=est.kappa(grb$X_DE))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(grb$X_DE , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Gamma Ray Burst DE data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)
dim(grb)


xray = read_csv('Xray.csv')
xray = data.frame(xray)
watson.test(xray$X_RA,  alpha = 0.01, dist = "vonmises")
watson.test(xray$X_DE,  alpha = 0.01, dist = "vonmises")

d = cbind(xray$X_DE, xray$X_RA)
Evmf <- function(K){
  movMF(d, k= K, control = list(nruns = 20))
}


set.seed(1227)
Esd = lapply(1:10, Evmf)
Esd
sapply(Esd, BIC)

x= xray$X_RA
y <- rvonmises(n=1000, mu=mean(xray$X_RA), kappa=est.kappa(xray$X_RA))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(xray$X_RA , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for X-Ray Burst RA data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

x= xray$X_DE
y <- rvonmises(n=1000, mu=mean(xray$X_DE), kappa=est.kappa(xray$X_DE))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(xray$X_DE , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for X-Ray Burst DE data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)


radio = read_csv('Radio.csv')
radio = data.frame(radio)
watson.test(radio$RAJ2000, alpha = 0.01, dist = "vonmises")
watson.test(radio$DEJ2000, alpha = 0.01, dist = "vonmises")

d = cbind(radio$DEJ2000, radio$RAJ2000)
Evmf <- function(K){
  movMF(d, k= K, control = list(nruns = 20))
}


set.seed(1227)
Esd = lapply(1:10, Evmf)
Esd
sapply(Esd, BIC)


x= radio$RAJ2000
y <- rvonmises(n=1000, mu=mean(radio$RAJ2000), kappa=est.kappa(radio$RAJ2000))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(radio$RAJ2000 , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Radio Burst RA data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

x= radio$DEJ2000
y <- rvonmises(n=1000, mu=mean(radio$DEJ2000), kappa=est.kappa(radio$DEJ2000))
resx <- density(x, bw=25)
resy <- density(y, bw=25)
pp=plot(density.circular(xray$X_DE , bw = 30) ,points.plot=F, pty = 10, lty=1, lwd=2, col="blue",xlim=c(-1.2,1), ylim=c(-1.1, 1.2), main="Comparison of  estimated sample density\n with fitted Von Mises\n for Radio Burst DE data", cex.main=0.25)
lines(resy, points.plot=F, col="red", points.col=2,lwd=2, lty=4, plot.info=pp)

legend("topleft", legend=c("estimated  \n kernel density", "von mises density"),
       col=c("blue", "red"), lty=c(1,4),lwd=c(2,2), cex=0.59,
       box.lty=0)

summary(grb)
summary(xray)
summary(radio)

library(SPADAR)
SPADAR::createAllSkyGridChart(longitude = grb$X_RA, latitude = grb$X_DE,mainGrid = "galactic",
                              galCol = "black")
library(sphereplot)

?SPADAR::createAllSkyScatterPlotChart()

createAllSkyScatterPlotChart(grb$X_RA, grb$X_DE, mainGrid = "galactic", eqDraw = F,
                             eclDraw = F, galCol = "black", main = "GRB plot", pointcol = "red")


createAllSkyScatterPlotChart(xray$X_RA, xray$X_DE, mainGrid = "galactic", eqDraw = F,
                             eclDraw = F, galCol = "black", main = "X-Ray plot", pointcol = "red")

createAllSkyScatterPlotChart(radio$RAJ2000, radio$DEJ2000, mainGrid = "galactic", eqDraw = F,
                             eclDraw = F, galCol = "black", main = "Radio Burst plot", pointcol = "red")


library(Directional)
vmf_density_grid = function(u, ngrid = 100){
  u[,1] <- u[,1] + 90
  u[,2] <- u[,2] 
  res <- vmf.kerncontour(u, thumb = "none", den.ret = T, full = T,
                         ngrid = ngrid)
  ret <- expand.grid(DEC = res$lat - 90 , RA = res$long )
  ret$Density <- c(res$den)
  ret
}


vmf.kerncontour

u = cbind(grb$X_DE, grb$X_RA)
u2 = cbind(xray$X_DE, xray$X_RA)
u3 = cbind(radio$DEJ2000, radio$RAJ2000)

grb.dens = vmf_density_grid(u , 300)
xray.dens = vmf_density_grid(u2, 300)
radio.dens = vmf_density_grid(u3, 300)
summary(grb.dens)

ggplot(grb) + 
  geom_point(aes(X_RA,X_DE)) +
  coord_map(projection="aitoff",orientation=c(90,180,0)) + 
  scale_y_continuous(breaks=(-2:2)*30,limits=c(-90,90)) + 
  scale_x_continuous(breaks=(0:8)*45,limits=c(0,360), labels=c("","","","","","","","",""))+
  geom_density_2d(data = grb,
                  aes(x = X_RA, y = X_DE),
                  color = "green", alpha = 2) +
  geom_contour(data = grb.dens, aes(x = RA , y = DEC, z = Density),
               color = "blue") +
  labs(x="R.A.(°)", y="Decl. (°)",title="Map of the GRB")+
  theme_void()



ggplot(xray) + 
  geom_point(aes(X_RA,X_DE)) +
  coord_map(projection="aitoff",orientation=c(90,180,0)) + 
  scale_y_continuous(breaks=(-2:2)*30,limits=c(-90,90)) + 
  scale_x_continuous(breaks=(0:8)*45,limits=c(0,360), labels=c("","","","","","","","",""))+
  geom_density_2d(data = grb,
                  aes(x = X_RA, y = X_DE),
                  color = "green", alpha = 2) +
  geom_contour(data = xray.dens, aes(x = RA , y = DEC, z = Density),
               color = "blue") +
  labs(x="R.A.(°)", y="Decl. (°)",title="Map of the X-Ray")+
  theme_void()

ggplot(radio) + 
  geom_point(aes(RAJ2000,DEJ2000)) +
  coord_map(projection="aitoff",orientation=c(90,180,0)) + 
  scale_y_continuous(breaks=(-2:2)*30,limits=c(-90,90)) + 
  scale_x_continuous(breaks=(0:8)*45,limits=c(0,360), labels=c("","","","","","","","",""))+
  geom_density_2d(data = radio,
                  aes(x = RAJ2000, y = DEJ2000),
                  color = "green", alpha = 2) +
  geom_contour(data = radio.dens, aes(x = RA , y = DEC, z = Density),
               color = "blue") +
  labs(x="R.A.(°)", y="Decl. (°)",title="Map of the Radio Burst")+
  theme_void()
