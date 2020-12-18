# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------
#
# Script para el análisis descriptivo de series de tiempo
#
# ---------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------

# CARGA DE LIBRERIAS

# A

library(ade4)
library(astsa)
library(ape)

# C
library(car) 
library(changepoint)
library(cluster)
library(coda)

#D
library(data.table)
library(dlm)
library(dplyr)
library(dygraphs)

# E
library(e1071)
library(expsmooth)

# F
library(fma)
library(FinTS)
library(forecast)
library(foreign)
library(fUnitRoots)
library(fftwtools)
library(fpp)

# G
library(ggcorrplot)
library(ggplot2)

# H
library(Hmisc)

# I
library(ipred)
library(influence.ME)
library(influence.SEM)

#L
library(lattice)
library(lme4)
library(lmtest)
library(lubridate)
library(longitudinal)

# M
library(MASS)
library(Matrix)
library(mlbench)

# O
library(openxlsx)

#P
library(partsm)  
library(pbkrtest)

# R
library(randomForest)
library(readr)
library(readxl)
library(relimp)
library(rpart)
library(RLRsim)
library(reshape)

# S
library(signal)
library(scales)
library(st)
library(StatMeasures)

# T
library(tseries)
library(tsoutliers)
library(TSrepr)
library(TSA)
library(tidyverse)
library(tuneR)

#U
library(UsingR)

# X
library(XML)
library(xts)

# Z
library(zoo)
  
# ----------------------------------------------------------------------------------------------
#
# DATOS EJEMPLOS
#
# ----------------------------------------------------------------------------------------------


# MEDICINA
data(ecg)
is.ts(ecg)

## 1. Crear la serie de tiempo
#tsdata <- ts(ecg, start = c(1), end = (1000))
#is.ts(tsdata)
# 
ts.plot(ecg)

# ECOLOGÍA

#OZONO

# Importacion de los datos
ozono <- read.table("c:/ozono.txt",header=FALSE)

# Crear la serie de tiempo

serie_ozono <- ts(ozono,start=c(1955),end=(1970),frequency=12)
is.ts(serie_ozono)
serie_ozono

# Representación gráfica la Serie Ozono

ts.plot(serie_ozono)
plot(serie_ozono,xlab= "Muestra",ylab="Ozono (pphm)")

# CLIMA
# VARIABLES: Temperatura minima (oC)
#            Temperatura maxima (oC)
#            Presipitacion (mm)

# Importacion de los datos

Crecimiento <- read_excel("c:/Users/Usuario/Desktop/series de tiempo/serie2.xlsx")
View(Crecimiento)

Clima <- read_excel("c:/Users/Usuario/Desktop/series de tiempo/serie.xlsx")
View(Clima)

# Crear la serie de tiempo
Temperatura_minima <- ts(Clima$tmmn,start=c(1995,5),end=c(1996,6),frequency=12)
is.ts(Temperatura_minima)

# Representación gráfica la Serie 
plot(Temperatura_minima,xlab= "Time",ylab="Temperatura minima (oC)")

Temperatura_maxima <- ts(Clima$tmmx,start=c(1995,5),end=c(1996,6),frequency=12)
is.ts(Temperatura_maxima)

# Representación gráfica la Serie 
plot(Temperatura_maxima,xlab= "Time",ylab="Temperatura maxima (oC)")

Precipitacion <- ts(Clima$Rain,start=c(1995,5),end=c(1996,6),frequency=12)
is.ts(Precipitacion)

# Representación gráfica la Serie 
plot(Precipitacion,xlab= "Time",ylab="Temperatura maxima (oC)")

# ----------------------------------------------------------------------------------------------
#
# DATOS FALTANTES
#
# ----------------------------------------------------------------------------------------------

Sp1 <- na.interp(Temperatura_minima)

# ----------------------------------------------------------------------------------------------
#
# VERIFICACION DE ESTACIONALIDAD
#
# ----------------------------------------------------------------------------------------------

# OZONO

adf <- adf.test(Temperatura_minima, alternative="stationary", k=0)
adf 

if(adf$p.value < 0.05) {
  Dickey.Fuller.Test <- "Estacionaria"
} else {
  Dickey.Fuller.Test <- "No Estacionaria"
}

Dickey.Fuller.Test

Serie_Diferencia <- diff(Temperatura_minima)
adf.test(Serie_Diferencia, alternative="stationary", k=0)
 
ppt <- pp.test(serie_ozono, alternative="stationary")

if(ppt$p.value < 0.05) {
  Phillips.Perron.Unit.Root.Test <- "Estacionaria"
} else {
  Phillips.Perron.Unit.Root.Test  <- "No Estacionaria"
}
Phillips.Perron.Unit.Root.Test 


# ecg

adf <- adf.test(ecg, alternative="stationary", k=0)

if(adf$p.value < 0.05) {
  Dickey.Fuller.Test <- "Estacionaria"
} else {
  Dickey.Fuller.Test <- "No Estacionaria"
}
Dickey.Fuller.Test
ppt <- pp.test(ecg, alternative="stationary")

if(ppt$p.value < 0.05) {
  Phillips.Perron.Unit.Root.Test <- "Estacionaria"
} else {
  Phillips.Perron.Unit.Root.Test  <- "No Estacionaria"
}
Phillips.Perron.Unit.Root.Test 


# ----------------------------------------------------------------------------------------------
#
# VERIFICACION DE DATOS ATÍPICOS
#
# ----------------------------------------------------------------------------------------------


# Obtención de los posibles datos outliers


outlier.serie_ozono <- tsoutliers::tso(serie_ozono, types = c("AO","LS","TC","IO","SLS"),maxit.iloop=10)
typo <- outlier.serie_ozono$outliers
typo <- typo$type

Total.outliers
NTC  <- length(typo[typo == "TC"])
NTC
NAO  <- length(typo[typo == "AO"])
NAO
NIO  <- length(typo[typo == "IO"])
NIO
NLS  <- length(typo[typo == "LS"])
NLS
NSLS  <- length(typo[typo == "SLS"])
NSLS
Total.outliers <- (NTC + NAO + NIO + NLS + NSLS)
Total.outliers

# Corrección del efecto de los datos atípicos

outliers <- tso(serie_ozono) 
newserie.ozono <- outliers$yadj 
ts.plot(serie_ozono)
plot(outliers)

# ----------------------------------------------------------------------------------------------
#
# DESCOMPOSICION DE LA SERIE DE TIEMPO
#
# ----------------------------------------------------------------------------------------------
ozono.fit.aditiva<-decompose(serie_ozono,type="additive")
plot(ozono.fit.aditiva)

ozono.fit.multi<-decompose(serie_ozono,type="multiplicative")
plot(ozono.fit.multi)

newozono.fit.aditiva<-decompose(newserie.ozono,type="additive")
plot(newozono.fit.aditiva)

ozono.diferencia <- diff(serie_ozono)
fit.ozono.diferencia<-decompose(ozono.diferencia,type="additive")
plot(fit.ozono.diferencia)
# ----------------------------------------------------------------------------------------------
#
# GENERACION DE MODELOS
#
# ----------------------------------------------------------------------------------------------

## Estudio en la parte regular de la serie.

par(mfrow=c(1,2))
acf(serie_ozono, 45, ylim=c(-1,1))
pacf(serie_ozono,45, ylim=c(-1,1))
lags <- length(serie_ozono)
lags

Bt2 <- Box.test(coredata(serie_ozono),type = 'Ljung-Box', lag =12)
if(Bt2$p.value < 0.05) {
  Ljung.Box.Test <- "Existe relación lineal, se rechaza la hipótesis de independencia"
} else {
  Ljung.Box.Test <- "No existe relación lineal, se acepta la hipótesis de independencia"
}
Ljung.Box.Test

## Estudio en la parte estacional de la serie. 

seasonplot(serie_ozono,12,ylab = "Ozono(pphm)",col=c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6))

At <- ArchTest(serie_ozono, lag=5)
if(At$p.value < 0.05) {
  Arch.Test <- "Existe varianza heterocedástica"
} else {
  Arch.Test<- "No existe varianza heterocedástica"
}
Arch.Test

# Modelos de ajuste 

plot(ozono.diferencia,xlab= "muestra",ylab="Ozono (pphm)")

y <- as.vector(ozono.diferencia)
x <- as.vector(time(ozono.diferencia))
r <- lm( y ~ poly(x,1) + cos(2*pi*x) + sin(2*pi*x) )
plot(y~x, type='l', xlab="time", ylab="Ozono(pphm")
lines(predict(r)~x, lty=3, col='red', lwd=3)

# Estimación de los coeficietes del Modelo ARIMA(p,d,q)(P,D,Q)12

ARIMA_A <- arima(ozono.diferencia, order=c(0, 0, 1))
ARIMA_A
summary(ARIMA_A)

ARIMA_B <- arima(ozono.diferencia, order=c(1, 0, 1))
ARIMA_B
summary(ARIMA_B)

ARIMA_C <- arima(ozono.diferencia, order=c(1, 1, 1))
ARIMA_C
summary(ARIMA_C)


Modelos.ARIMA <- c("ARIMA_A","ARIMA_B","ARIMA_C")
AIC.ARIMA <- c(ARIMA_A$aic, ARIMA_B$aic, ARIMA_C$aic)
Compara.ARIMA <- data.frame(Modelos.ARIMA,AIC.ARIMA)
View(Compara.ARIMA)


SARIMA_A <-arima(ozono.diferencia,c(0,0,1),list(order=c(0,1,1)))
SARIMA_A
summary(SARIMA_A)

SARIMA_B<-arima(ozono.diferencia,c(1,0,1),list(order=c(0,1,1)))
SARIMA_B
summary(SARIMA_B)

## Contraste y validez del modelo

tsdiag(SARIMA_A)

tsdiag(SARIMA_B)

## Análisis de los residuos
par(mfrow=c(3,1))
residuos_arima <-residuals(SARIMA_B)
acf(residuos_arima,24, ylim=c(-1,1))
pacf(residuos_arima,24, ylim=c(-1,1))
qqnorm(residuos_arima)
shapiro.test(residuos_arima)
shapiro.test

# Predicciones del modelo

par(mfrow=c(1,1))
SARIMA_A.f<-predict(SARIMA_A,n.ahead=150)
plot(ozono.diferencia,xlim=c(1955,1971),ylim=range(c(ozono.diferencia,SARIMA_A.f$pred)))
lines(SARIMA_A.f$pred,col='red')
lines(SARIMA_A.f$pred+qnorm(.025)*SARIMA_A.f$se,col='blue',lty=3)
lines(SARIMA_A.f$pred+qnorm(.975)*SARIMA_A.f$se,col='blue',lty=3)

par(mfrow=c(1,1))
SARIMA_B.f<-predict(SARIMA_B,n.ahead=150)
plot(ozono.diferencia,xlim=c(1955,1971),ylim=range(c(ozono.diferencia,SARIMA_B.f$pred)))
lines(SARIMA_B.f$pred,col='red')
lines(SARIMA_B.f$pred+qnorm(.025)*SARIMA_B.f$se,col='blue',lty=3)
lines(SARIMA_B.f$pred+qnorm(.975)*SARIMA_B.f$se,col='blue',lty=3)

## PREDICCIONES CON FORECAST

par(mfrow=c(1,1))
moda_1<- forecast(SARIMA_B,h=20)
moda_1 
plot(moda_1)

#TRANSFORMADA DE FOURIER

spectrum(serie_ozono)
abline(v=1:10, lty=3)
ozono.fourier <- fft(ozono$V1)
plot(abs(ozono.fourier), type="h")
plot(Arg(ozono.fourier), type="h")
obj = stft(ozono$V1)
plot(obj, cex = 5)

dev.new(); plot(obj, mode = "pval", reassign = FALSE, topthresh = Inf, log = "y") 

# MODELO DE EFECTOS MIXTOS

ggplot(data = Crecimiento, aes(x = Semana, y = Altura, color = Nplant)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~ Nplant) + 
  theme(legend.position = "none")

# MODELO DE INTERCEPTO ALEATORIO CON MEDIA FIJA

M0 <- lmer(NHojas ~ 1 + Semana + Trat + TT + Altura + (1|Nplant), Crecimiento, REML = FALSE) 
summary(M0)

# MODELO CON PENDIENTE E INTERCEPTO ALEATORIO CONDICIONADO 
M1 <- lmer(NHojas ~ 1 + Semana + Trat + TT + Altura + (Semana|Nplant), Crecimiento, REML = FALSE) 
summary(M1)

# MODELO CON PENDIENTE E INTERCEPTO ALEATORIO NO CONDICIONADO
M2 <- lmer(NHojas ~ 1 + Semana + Trat + TT + Altura + (1|Nplant) + (0 + Semana|Nplant), Crecimiento, REML = FALSE)
summary(M2)

# Validación (Analisis de residuos)
plot(fitted(M0),resid(M0,type="pearson"),col="blue")
densityplot(~ resid(M0)| factor(Semana), Crecimiento, groups = Trat,
            plot.points = TRUE, auto.key = TRUE,  xlim=c(-8,8), ylim=c(0.0,0.8))

par(mfrow = c(1, 2))
xyplot( NHojas ~Altura|Semana, data = Crecimiento, group = Trat, color='red',  grid = TRUE, xlim=c(0,8), ylim=c(0,20))
xyplot( NHojas ~predict(M0)|Semana, data = Crecimiento,   grid = TRUE,xlim=c(0,8), ylim=c(0,20))



