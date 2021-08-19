
install.packages("aTSA")
install.packages("TSA")
install.packages("forecast")
library("forecast")
library("aTSA")
library("TSA")
install.packages("tseries")

donnees = read.csv("FluDeaths.csv", header = TRUE, sep = ';')
weeks = as.data.frame(donnees$Weeks)
deaths = donnees$Deaths

checkupRes <- function(Res){
  Res = Res[!is.na(Res)]
  layout(matrix(c(1:9), nrow=3, ncol=3, byrow=TRUE))
  plot(Res,type = 'l',xlab = "time", main = "Evolution des résidus" )
  acf(Res, plot = TRUE,na.action = na.pass, lag.max = 10)
  pacf(Res, plot = TRUE,na.action = na.pass, lag.max = 10)
  x = c(0,Res[-length(Res)])
  plot(Res,x,xlab="Res[i-1]",ylab = 'Res[i]',type = 'p', col="green",
       main = "Résidus en fonction de leurs retards directs")
  hist(Res, col = "red",freq = FALSE)
  qqplot(1:length(Res),Res,plot.it = TRUE,xlab = "time", main = 
           "QQ plot par rapport aux quantiles gaussiens")
  ResNorm = (Res-mean(Res))/sd(Res)
  plot(ResNorm, type = 'p', main = "Résidus normalisés")
  abline(h = c(-1.96,1.96), col = 'red',lty = 2)
  plot()
  print(shapiro.test(ResNorm))
}


dev.off()
layout(matrix(c(1:6),3,2, byrow=TRUE))
plot(deaths, type = 'l',main="série originale",col = 'blue')
plot(log(deaths), type = 'l',main="série log")
acf(deaths, plot = TRUE,ci.col = "magenta",na.action = na.pass,lag.max =100,main="acf de la série originale")
acf(log(deaths), plot = TRUE,na.action = na.pass,main="acf de la série log")
pacf(deaths, plot = TRUE,na.action = na.pass ,main="pacf de la série originale")
pacf(log(deaths), plot = TRUE,na.action = na.pass,main="pacf de la série log")

# La périodicité est bien visible sur la série originale, avec une intensité constante.
# La série suit un modèle additif.C'est d'ailleurs la raison pour laquelle il n' y a pas
# de différence remarquable entre la série originale et son logarithme.
# On ne prend donc pas en considération la série log
# La pacf de la série originale est quasi nulle après le rang 3. On envisagera donc tester dans la  
# suite une modélisation AR(3)

# Tests de stationnarité

kpss.test(deaths) # stationnaire
adf.test(deaths) # non stationnaire jusqu'au seuil de 4%
per = periodogram(deaths,plot = TRUE, ylab = "Périodogramme", xlab = frequency)

# La fonction a utilisé n=300 pour éstimer les fréquences

Box.test(deaths,type = "Ljung") # série corrélée


# Tests de stationnarité sur la série différenciée

Ideaths = diff(deaths)
dev.off()
plot(Ideaths, type = 'l',main="série originale différenciée",col = 'blue')
acf(Ideaths, plot = TRUE,na.action = na.pass,main="acf de la série différenciée")
pacf(Ideaths, plot = TRUE,na.action = na.pass,main="pacf de la série différenciée")

# Au regard des acf qui sont presque toutes nulles après le rang 1 , la série différenciée
# conviendrait d'être modélisée par un MA(1)
# De même on envisagera de tester un AR(2) sur la série différenciée( compte tenu les pacf)

kpss.test(Ideaths) # stationnaire
adf.test(Ideaths) # aussi stationnaire

#On se limite aux premiers incréments pour nous faciliter l'interprétation, d'autant 
# plus que la série obtenue est parfaitement stationnaire.


# Voyons aussi les différenciations saisonnières (avec comme péride=52,4,8, et 12)

# période = 52
DSdeaths52 = diff(deaths, lag = 52)
plot(DSdeaths52, type = 'l',main="série saisonnièrement différenciée",col = 'blue')
acf(DSdeaths52, plot = TRUE,na.action = na.pass,main="acf de la série différenciée")
pacf(DSdeaths52, plot = TRUE,na.action = na.pass,main="pacf de la série différenciée")
# On pourra tester un AR(3) avec une diff saisonnière de période 52

kpss.test(DSdeaths52) # stationnaire
adf.test(DSdeaths52) # stationnaire

# période =12
DSdeaths12 = diff(deaths, lag = 12)
plot(DSdeaths12, type = 'l',main="série saisonnièrement différenciée",col = 'blue')
acf(DSdeaths12, plot = TRUE,na.action = na.pass,lag.max =100,main="acf de la série différenciée")
pacf(DSdeaths12, plot = TRUE,na.action = na.pass,lag.max=100,main="pacf de la série différenciée")
# Les acf se comportent très mal ,les pacf ont des pics saisonniers en 0,12,24,36 et 60
# on pourra tester un AR(3) ou un AR(5)

# période =13
DSdeaths13 = diff(deaths, lag = 13)
plot(DSdeaths13, type = 'l',main="série saisonnièrement différenciée",col = 'blue')
acf(DSdeaths13, plot = TRUE,na.action = na.pass,lag.max =100,main="acf de la série différenciée")
pacf(DSdeaths13, plot = TRUE,na.action = na.pass,lag.max =100,main="pacf de la série différenciée")
# On voit aussi des pics saisonniers à 13 et 26. On testera le modèle AR(3)

kpss.test(DSdeaths13,lag.short = FALSE) # stationnaire
adf.test(DSdeaths13, nlag = 1) # stationnaire


# période =8
DSdeaths8 = diff(deaths, lag = 8)
plot(DSdeaths8, type = 'l',main="série saisonnièrement différenciée",col = 'blue')
acf(DSdeaths8, plot = TRUE,na.action = na.pass,lag.max =100,main="acf de la série différenciée")
pacf(DSdeaths8, plot = TRUE,na.action = na.pass,lag.max =100,main="pacf de la série différenciée")
#De même ici , les acf et pacf se comportent très mal. les pics saisonniers jusqu'à 40.
# On testera un AR(3)

kpss.test(DSdeaths8,lag.short = FALSE)# stationnaire
adf.test(DSdeaths8, nlag = 1)# stationnaire

# période =4
DSdeaths4 = diff(deaths, lag = 4)
plot(DSdeaths4, type = 'l',main="série saisonnièrement différenciée",col = 'blue')
acf(DSdeaths4, plot = TRUE,na.action = na.pass,lag.max =100,main="acf de la série différenciée")
pacf(DSdeaths4, plot = TRUE,na.action = na.pass,lag.max =100,main="pacf de la série différenciée")
# idem
kpss.test(DSdeaths4,lag.short = FALSE) # stationnaire
adf.test(DSdeaths4, nlag = 1)# stationnaire

#  MODELISATIONS

# premier modèle : AR(3) sur la série originale

mod1 = Arima(deaths, order = c(3,0,0),include.drift =TRUE)
summary(mod1)
checkupRes(mod1$residuals)

# Seule l'estimation de la pente de la tendance linéaire est significative,bien que le modèle 
# semble tenir la route au regard de ses  résidus .
# effet , les acf et pacf se comportent comme un parfait bruit blanc stationnaire,le test 
# de shapiro ne rejette pas l'hypothèse de normalité et l'histogramme
# colle bien avec une densité gaussienne.BIC=3043.7

# deuxième modèle : MA(1) sur la série différenciée cad ARIMA(0,1,1)

mod2 = Arima(deaths, order = c(0,1,1),include.drift =TRUE)
summary(mod2)
checkupRes(mod2$residuals)
# Les conclusions sont similaires, seul le coéfficient angulaire de la tendance est significatif
# et les résidus se comportent très bien. le test de normalité de Shapiro est correct
#BIC=3034.2

# Troisième modèle : AR(2) sur les données différenciées

mod3 = Arima(deaths, order = c(2,1,0),include.drift =TRUE)
summary(mod3)
checkupRes(mod3$residuals)
# idem, BIC=3034.21

# Quatrième modèle : Un AR(3) sur les données avec une diff saisonnière de période 52

mod4 = Arima(deaths, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 52))
summary(mod4)
checkupRes(mod4$residuals)
# Modèle inadapté, coéfficients non significatifs et résidus non gaussiens d'aprèsle
# le test de Shapiro. Néanmoins, les acf et pacf des résidus ressemblent bien à un 
# bruit blanc Il a par un BIC plus petit que tous les modèles que nous avons déjà testés.
#On verra sa qualité prédictive dans le choix du modèle final.BIC=2670.74

# Cinquième modèle : Un AR(3) sur les données avec une diff saisonnière de période 12

mod5 = Arima(deaths, order = c(5,0,0), seasonal = list(order = c(0,1,0), period = 12))
summary(mod5)
checkupRes(mod5$residuals)
#Deux des 5 coéfficients du polynôme autoregressif sont significatifs, BIC=3084.55

# Sixième modèle : Un AR(3) sur les données avec une diff saisonnière de période 8

mod6 = Arima(deaths, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 8))
summary(mod6)
checkupRes(mod6$residuals)
# Pas de coéfficient significatif, BIC=3137.36

# Septième modèle : Un AR(3) sur les données avec une diff saisonnière de période 4

mod7 = Arima(deaths, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 4))
summary(mod7)
checkupRes(mod7$residuals)
# Pas de coéfficient significatif, BIC=3137.36 et les résidus non gaussiens( tant sur le 
# test de Shapiro que sur les acf et pacf qui montrent des corrélations )

# Huitième modèle : Un AR(3) sur les données avec une diff saisonnière de période 13

mod8 = Arima(deaths, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 13))
summary(mod8)
checkupRes(mod8$residuals)
# Pas de coéf significatif, BIC=3075.93 et les résidus parfaitement gaussiens et bruit blanc

# Neuvième modèle : Recherche automatique du meilleur modèle sur critère du BIC
mod10 <- auto.arima(deaths,seasonal= TRUE, allowdrift = TRUE ,allowmean = TRUE, ic = "bic",xreg = time(deaths))
summary(mod10)
dev.off()
hist(mod10$residuals)
# L'auto arima confirme que le meilleur SARIMA est le modèle mod2
# Pour choisir un modèle parmi les 5, on omettra la dernière saison pour la reprédire.
bic = c()
for (p in 0:2){
  for (q in 0:2){
    modele =  Arima(deaths, order = c(p,1,q),include.drift =TRUE)
    bic = c(bic,modele$bic)
  }
}
mod9 = Arima(deaths, order = c(0,1,1),include.drift =TRUE)
summary(mod9)
checkupRes(mod9$residuals) # c'est le modèle 2

Temps = c(1:247)
deathsTtronque = deaths[1:247]

mod1T = Arima(deathsTtronque, order = c(3,0,0),include.drift =TRUE)
mod2T = Arima(deathsTtronque, order = c(0,1,1),include.drift =TRUE,xreg = Temps)
mod3T = Arima(deathsTtronque, order = c(2,1,0),include.drift =TRUE)
mod4T = Arima(deathsTtronque, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 52))
mod5T = Arima(deathsTtronque, order = c(5,0,0), seasonal = list(order = c(0,1,0), period = 12))
mod6T = Arima(deathsTtronque, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 8))
mod7T = Arima(deathsTtronque, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 4))
mod8T = Arima(deathsTtronque, order = c(3,0,0), seasonal = list(order = c(0,1,0), period = 13))
lhw = HoltWinters(ts(deathsTtronque, frequency = 52),seasonal = "additive")
mod10T = Arima(deathsTtronque, order = c(1,0,1) , include.mean = TRUE ,xreg = Temps)

Pred1T <- forecast::forecast(mod1T, h=52, level=95)
Pred2T <- forecast::forecast(mod2T, h=52, level=95)
Pred3T <- forecast::forecast(mod3T, h=52, level=95)
Pred4T <- forecast::forecast(mod4T, h=52, level=95)
Pred5T <- forecast::forecast(mod5T, h=52, level=95)
Pred6T <- forecast::forecast(mod6T, h=52, level=95)
Pred7T <- forecast::forecast(mod7T, h=52, level=95)
Pred8T <- forecast::forecast(mod8T, h=52, level=95)
PredHT <- forecast::forecast(lhw, h=52, level=95)
Pred10T <- forecast::forecast(mod10T, h=52, level=95,xreg = tempsTronque)

tempsTronque = c(248:299)

dev.off()
plot(deaths, type = 'l', xlim = c(1,299), ylim = c(650,1000),
     main = "Superposition avec les prédictions de l'Arima(3,0,0)",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred1T$upper[,1],rev(Pred1T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred1T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred1T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred1T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,299), ylim = c(500,1050),
     main = "Superposition avec les prédictions de l'Arima(0,1,1)",col.main = "magenta")
#polygon(c(tempsTronque,rev(tempsTronque)), c(Pred2T$upper[,1],rev(Pred2T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred2T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred2T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred2T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,299), ylim = c(500,1100),
     main = "Superposition avec les prédictions de l'Arima(2,1,0)",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred3T$upper[,1],rev(Pred3T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred3T$mean, col = 'blue',lwd = 1)
lines(tempsTronque, Pred3T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred3T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,299),ylim = c(600,1100),
     main = "Superposition avec les prédictions du Sarima(3,0,0)x(0,1,0) de période 52",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred4T$upper[,1],rev(Pred4T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred4T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred4T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred4T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,299),ylim = c(350,1150),
     main = "Superposition avec les prédictions du Sarima(5,0,0)x(0,1,0) 12",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred5T$upper[,1],rev(Pred5T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred5T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred5T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred5T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,299),ylim = c(350,1150),
     main = "Superposition avec les prédictions du Sarima(3,0,0)x(0,1,0) 8",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred6T$upper[,1],rev(Pred6T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred6T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred6T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred6T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,299),ylim = c(400,1300),
     main = "Superposition avec les prédictions du Sarima(3,0,0)x(0,1,0) 4",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred7T$upper[,1],rev(Pred7T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred7T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred7T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred7T$lower[,1], col = 'blue',lty = 1)

dev.off()
plot(deaths, type = 'l', xlim = c(1,299),ylim = c(350,1100),
     main = "Superposition avec les prédictions du Sarima(3,0,0)x(0,1,0) 13",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(Pred8T$upper[,1],rev(Pred8T$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, Pred8T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred8T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred8T$lower[,1], col = 'blue',lty = 1)


dev.off()
plot(deaths, type = 'l', xlim = c(1,300),ylim = c(600,1100),
     main = "Superposition avec les prédictions du lissage de Holt Winters",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)), c(PredHT$upper[,1],rev(PredHT$lower[,1])), col = 'gray', border=FALSE )
lines(tempsTronque, PredHT$mean, col = 'blue',lwd = 2)
lines(tempsTronque, PredHT$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, PredHT$lower[,1], col = 'blue',lty = 1)

dev.off()
plot(deaths, type = 'l', xlim = c(1,299), ylim = c(600,1100),
     main = "Superposition avec les prédictions de l'ARIMA(1,0,1)",col.main = "magenta")
polygon(c(tempsTronque,rev(tempsTronque)),
        c(Pred10T$upper[,1],rev(Pred10T$lower[,1])),
        col = 'gray', border=FALSE )
lines(tempsTronque, Pred10T$mean, col = 'blue',lwd = 2)
lines(tempsTronque, Pred10T$upper[,1], col = 'blue',lty = 1)
lines(tempsTronque, Pred10T$lower[,1], col = 'blue',lty = 1)

n = length(deaths)

MSE1 =sum((Pred1T$mean - deaths[(n-51):n])**2)/52
MSE2 =sum((Pred2T$mean - deaths[(n-51):n])**2)/52
MSE3 =sum((Pred3T$mean - deaths[(n-51):n])**2)/52
MSE4 =sum((Pred4T$mean - deaths[(n-51):n])**2)/52
MSE5 =sum((Pred5T$mean - deaths[(n-51):n])**2)/52
MSE6 =sum((Pred6T$mean - deaths[(n-51):n])**2)/52
MSE7 =sum((Pred7T$mean - deaths[(n-51):n])**2)/52
MSE8 =sum((Pred8T$mean - deaths[(n-51):n])**2)/52
MSE9 =sum((PredHT$mean - deaths[(n-51):n])**2)/52
MSE10 =sum((Pred10T$mean - deaths[(n-51):n])**2)/52
mse = c(MSE1,MSE2,MSE3,MSE4,MSE5,MSE6,MSE7,MSE8,MSE9,MSE10)

# Meilleurs modèles choisis sur un compromis de pouvoir prédictif
# et du bic( le modèle 4 et le 9 )

Pred1 <- forecast::forecast(mod4, h=78, level=95)
LHW = HoltWinters(ts(deaths, frequency = 52),seasonal = "additive")
Pred2 <- forecast::forecast(LHW, h=78, level=95)
temps = c(300:377)

dev.off()
plot(deaths, type = 'l', xlim = c(1,378),ylim = c(550,1250),
     main = "Superposition avec les prédictions du Sarima(3,0,0)x(0,1,0) de période 52",col.main = "magenta")
legend(0,1250,xjust = 0,yjust = 1, legend = c("prédictions","données présentes"),lwd = 1,
       col = c('blue','black'),text.col=c('blue','black'),lty=rep(1,2),cex=0.7 )
polygon(c(temps,rev(temps)), c(Pred1$upper[,1],rev(Pred1$lower[,1])), col = 'gray', border=FALSE )
lines(temps, Pred1$mean, col = 'blue',lwd = 2)
lines(temps, Pred1$upper[,1], col = 'gray',lty = 1)
lines(temps, Pred1$lower[,1], col = 'gray',lty = 1)

dev.off()
plot(deaths, type = 'l', xlim = c(1,378),ylim = c(600,1250),
     main = "Superposition avec les prédictions du lissage exponentiel de Holtwinters",col.main = "magenta")
polygon(c(temps,rev(temps)), c(Pred2$upper[,1],rev(Pred2$lower[,1])), col = 'gray', border=FALSE )
legend(0,1250,xjust = 0,yjust = 1, legend = c("prédictions","données présentes"),lwd = 1,
       col = c('blue','black'),text.col=c('blue','black'),lty=rep(1,2),cex=0.7 )
lines(temps, Pred2$mean, col = 'blue',lwd = 2)
lines(temps, Pred2$upper[,1], col = 'gray',lty = 1)
lines(temps, Pred2$lower[,1], col = 'gray',lty = 1)
