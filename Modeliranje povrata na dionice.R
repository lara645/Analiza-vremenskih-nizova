library("fGarch")
library("TSA")
library("tseries")
library("aTSA")
library("astsa")
library("forecast")
library("FinTS")
library("fBasics")

acf <- stats::acf
adf.test <- tseries::adf.test

setwd("C:/Users/LARA/Desktop/AVN - R")
baza <- read.csv2('Ericsson.csv',header=TRUE,sep = ';')
baza1 <- baza[765:1253,] #12.2.2019 - 12.2.2021
baza1$Date <- as.Date(baza1$Date, "%d.%m.%Y")
detach(baza1)
attach(baza1)
str(baza1)

najvisa <- ts(High.Price)
summary(najvisa)
baza1[High.Price == 992,]
baza1[High.Price == 1530,]
sd(najvisa)
par(mfrow=c(1,1))
plot(High.Price ~ Date, xaxt = "n", type = "l", xlab = "Vrijeme", ylab = "Najviša dnevna cijena dionice Ericsson Nikola Tesla d.d.", col = "gray55", lwd = 1)           
axis(1, Date, format(Date, "%b %Y"), cex.axis = .7)
points(High.Price~ Date, cex = 1.5, pch = 20, col = adjustcolor(4, .4))

par(mfrow=c(1,1))
acf(najvisa,col = "gray57", lwd = 2.5, main = "")
adf.test(najvisa) #p-value = 0.4136
log_povrati <- (diff(log(najvisa))*100)
adf.test(log_povrati) #p-value = 0.01


### Log-povrati

par(mfrow=c(1,1))
plot(log_povrati, xlab = "Vrijeme", ylab = "Log-povrati", type = "l", col = "darkgrey", lwd = 2.5)
points(log_povrati, cex = 1, pch = 19, col = adjustcolor(4, .4))
abline(h=0,col="deepskyblue3",lwd = 2.5)

par(mfrow=c(1,2))
hist(log_povrati,prob = T, col = c("azure3"),main="",xlab="Log-povrati",ylab="Relativna frekvencija",ylim = c(0,0.5))
lines(density(log_povrati),col="blue",lwd = 1.8)
curve(dnorm(x,mean(log_povrati),sd(log_povrati)),col="red",add=T)
legend("topleft",c("Uzoraèka funkcija gustoæe", 
                   "teorijska funkcija gustoæe normalne distribucije"), col=c("deepskyblue3",
                    "red"), lwd=2, cex=0.6)
qqnorm(log_povrati,main="",xlab="Teorijski kvantili",ylab="Uzoraèki kvantili", cex = 1, pch = 19,
       col = adjustcolor(4, .4)) 
qqline(log_povrati,col="darkgrey", lwd=3)

kurtosis(log_povrati)

par(mfrow=c(1,3))
acf(log_povrati,col = "gray57", lwd = 2.5, main = "") 
acf(log_povrati^2,col = "gray57", lwd = 2.5, main = "")
acf(abs(log_povrati),col = "gray57", lwd = 2.5, main = "") 

par(mfrow=c(1,1))
McLeod.Li.test(y =log_povrati)
ArchTest(log_povrati)


### Identifikacija modela

najboljiAIC_norm <- function(pmax, qmax, log_povrati) {
  pred <- c()
  qred <- c()
  AICvrijednost <- c()
  for (i in 1:pmax) {
    for (j in 0:qmax) {
      model <- garchFit(substitute(~ garch(p,q), list(p=i, q=j)), data = log_povrati, trace = FALSE)
      fitaic <- tryCatch(model@fit$ics[1], error = function(e) NaN) 
      pred <- c(pred,i)
      qred <- c(qred,j)
      AICvrijednost <- c(AICvrijednost, fitaic)
    }
  }
  rez <- data.frame(p = pred, q = qred, AICvrijednosti = AICvrijednost)
  rez <- rez[order(rez$AICvrijednosti), ]
  return(rez[1:min(5,length(pred)), ])
}

najboljiAIC_st <- function(pmax, qmax, log_povrati) {
  pred <- c()
  qred <- c()
  AICvrijednost <- c()
  for (i in 1:pmax) {
    for (j in 0:qmax) {
      model <- garchFit(substitute(~ garch(p,q), list(p=i, q=j)), data = log_povrati, trace = FALSE, cond.dist = "sstd")
      fitaic <- tryCatch(model@fit$ics[1], error = function(e) NaN) 
      pred <- c(pred,i)
      qred <- c(qred,j)
      AICvrijednost <- c(AICvrijednost, fitaic)
    }
  }
  rez <- data.frame(p = pred, q = qred, AICvrijednosti = AICvrijednost)
  rez <- rez[order(rez$AICvrijednosti), ]
  return(rez[1:min(5,length(pred)), ])
}

najboljiBIC_norm <- function(pmax, qmax, log_povrati) {
  pred <- c()
  qred <- c()
  BICvrijednost <- c()
  for (i in 1:pmax) {
    for (j in 0:qmax) {
      model <- garchFit(substitute(~ garch(p,q), list(p=i, q=j)), data = log_povrati, trace = FALSE)
      fitbic <- tryCatch(model@fit$ics[2], error = function(e) NaN)
      pred <- c(pred,i)
      qred <- c(qred,j)
      BICvrijednost <- c(BICvrijednost, fitbic)
    }
  }
  rez <- data.frame(p = pred, q = qred, BICvrijednosti = BICvrijednost)
  rez <- rez[order(rez$BICvrijednosti), ]
  return(rez[1:min(5,length(pred)), ])
}

najboljiBIC_st <- function(pmax, qmax, log_povrati) {
  pred <- c()
  qred <- c()
  BICvrijednost <- c()
  for (i in 1:pmax) {
    for (j in 0:qmax) {
      model <- garchFit(substitute(~ garch(p,q), list(p=i, q=j)), data = log_povrati, trace = FALSE, cond.dist = "sstd")
      fitbic <- tryCatch(model@fit$ics[2], error = function(e) NaN)
      pred <- c(pred,i)
      qred <- c(qred,j)
      BICvrijednost <- c(BICvrijednost, fitbic)
    }
  }
  rez <- data.frame(p = pred, q = qred, BICvrijednosti = BICvrijednost)
  rez <- rez[order(rez$BICvrijednosti), ]
  return(rez[1:min(5,length(pred)), ])
}

najboljiAIC_norm(4,4,log_povrati)
# p q AICvrijednosti
# 2  1 1       3.254942
# 7  2 1       3.256454
# 3  1 2       3.259169
# 8  2 2       3.260011
# 12 3 1       3.260560
najboljiBIC_norm(4,4,log_povrati)
# p q BICvrijednosti
# 2  1 1       3.289289
# 7  2 1       3.299388
# 3  1 2       3.302102
# 8  2 2       3.311532
# 12 3 1       3.312080
najboljiAIC_st(4,4,log_povrati)
# p q AICvrijednosti
# 2  1 1       3.160323
# 7  2 1       3.164053
# 3  1 2       3.164501
# 12 3 1       3.167785
# 8  2 2       3.168151
najboljiBIC_st(4,4,log_povrati)
# p q BICvrijednosti
# 2  1 1       3.211843
# 7  2 1       3.224160
# 3  1 2       3.224608
# 12 3 1       3.236479
# 8  2 2       3.236845


### Modeliranje

# GARCH(1,1)
model1 <- garchFit(~garch(1,1),data=log_povrati,trace=FALSE)
model1
summary(model1)
par(mfrow=c(2,2))
plot(model1)

# GARCH(2,1)
model2 <- garchFit(~garch(2,1),data=log_povrati,trace=FALSE)
model2
summary(model2)
plot(model2)

# GARCH(1,2)
model3 <- garchFit(~garch(1,2),data=log_povrati,trace=FALSE)
model3 
summary(model2)
plot(model2)

# GARCH(2,2)
model4 <- garchFit(~garch(2,2),data=log_povrati,trace=FALSE)
model4

# GARCH(1,1) sa Studentovom distr.
model4 <- garchFit(~garch(1,1),data=log_povrati,trace=FALSE,cond.dist = "sstd")
model4
summary(model4) 
plot(model4)

model5 <- garchFit(~garch(1,2),data=log_povrati,trace=FALSE,cond.dist = "sstd")
model5
model6 <- garchFit(~garch(2,1),data=log_povrati,trace=FALSE,cond.dist = "sstd")
model6
model7 <- garchFit(~garch(2,2),data=log_povrati,trace=FALSE,cond.dist = "sstd")
model7

### Predikcija
par(mfrow=c(1,2))
#[765:1253,]

sve <- (diff(log(ts(baza$High.Price[765:1380])))*100)
#prediction with confidence interval
predict(model1, n.ahead = 127, trace = FALSE, mse = c("cond","uncond"),
        plot=TRUE, nx=488,main = "")
lines(sve,col="red")
legend("topleft", c("Stvarni log-povrati", "Predviðeni log-povrati"), cex=0.8, col = c("red", "black"), pch = c(16, 16))
predict(model4, n.ahead = 127, trace = FALSE, mse = c("cond","uncond"),
        plot=TRUE, nx=488,main = "")
lines(sve,col="red")
legend("topleft", c("stvarni log-povrati", "predviðeni log-povrati"), cex=0.8, col = c("red", "black"), pch = c(16, 16))
