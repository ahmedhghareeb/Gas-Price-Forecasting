
#setwd("~/Winter 2015/STAT 443/Project")
#gas <- read.table('Absent_data2.txt');

# Function for residual diagnostics
resdiags <- function(res) # you give this function a vector containing residuals from a model, as well as vector of fitteds
{
  par(mfcol=c(2,2)) # splits the view to show 4 plots
  ts.plot(res) # time series plot of residuals
  points(res) # points to make counting runs easier
  abline(h=mean(res)) # mean line
  qqnorm(res) #qq plot
  qqline(res)
  acf(res) #acf
  acf(res, type="partial") #pacf
}

# Data broken up into different training and testing sets
# For all data, we take the last year as testing set
gas.all <- ts(gas, start = c(2005, 1), frequency = 12)
gas.train.all <- ts(gas[1:108,1], start = c(2005, 1), frequency = 12)
gas.test.all <- ts(gas[109:120,1], start = c(2014, 1), frequency = 12)
# We also consider the data from 2009 onwards due to recession
gas.train.new <- ts(gas[49:108,1], start = c(2009, 1), frequency = 12)
gas.test.new <- ts(gas[109:120,1], start = c(2014, 1), frequency = 12)

# Complete plot of all the data
plot(as.ts(gas.all))



# Try to get rid of 2008 spike for easier linear regression


# Fitting pure seasonal model to training data
time <- time(gas.train.all)
months <- as.factor(cycle(time))
seas1 <- lm(gas.train.all~months)
points(time, seas1$fitted, col="red", type = "l")
resdiags(seas1$res)

# Fit a model with linear time component, quadratic time component
tim.train.all <- time(gas.train.all) # For linear time term
tim2.train.all <- tim.train.all^2 # For quadratic time term
seas2 <- lm(gas.train.all~tim.train.all+tim2.train.all+months)
points(time, seas2$fitted, col="red", type = "l")

par(mfcol=c(1,1))
plot(as.ts(log(gas.all)))

# Fitting pure seasonal model to log(training data)
seas3 <- lm(log(gas.train.all)~months)
q<-points(time, seas3$fitted, col="red", type = "l")
resdiags(seas3$res)

# Fit a log model with linear time component, quadratic time component
seas4 <- lm(log(gas.train.all)~tim.train.all+tim2.train.all+months)
points(time, seas4$fitted, col="red", type = "l")
resdiags(seas4$res)

#tadj <- seq(1:108)
#tadj2 <- tadj^2
tadj <- model.matrix(seas4)[,2] -2004
tadj2 <- tadj^2
seas4.t.adj <- lm(gas.train.all~tadj+tadj2+months)

# We prepare the testing set for prediction intervals for seas4, our quadratic time, log, seasonality
tim.test.all <- time(gas.test.all)
tim2.test.all <- tim.test.all^2
month.test <- as.factor(cycle(gas.test.all))
pred.test <- predict.lm(seas4, newdata=data.frame(tim.train.all = tim.test.all, tim2.train.all = tim2.test.all, months = month.test), interval = "prediction")
pred.test
log(gas.test.all)

#####################################################################33

# For seas4, our residuals look like AR(1), so we fit this to residuals
seas4.adjx <- model.matrix(seas4.t.adj)
seas4.adjx2 <- seas4.adjx[,-1]
seas2.ar1 <- arima(gas.train.all, order=c(1,0,0), xreg=seas4.adjx2, method="ML")
resdiags(seas2.ar1$res)

#####################MA Filter
#General MA filter function
MAsmooth <- function(series, q)
{
  c = 1/(2*q+1) # constant to multiply by
  series.MA <- series # starting point Yt
  for(i in 1:q){ series.MA <- series.MA + lag(series, k=-i) + lag(series, k=i) } # adding Yt-i and Yt+i to the averaged series for each i from 1 to q
  series.MA <- series.MA*c # multiplying by the constant
  return(series.MA)
}

#training set
plot(gas.train.all) # Data plot
lines(MAsmooth(gas.train.all,10),col="red") #q=10
lines(MAsmooth(gas.train.all,5),col="green") #q=5
lines(MAsmooth(gas.train.all,3),col="orange") #q=3
lines(MAsmooth(gas.train.all,1),col="blue") #q=1
#testing set
plot(gas.test.all) # Data plot
lines(MAsmooth(gas.test.all,3),col="orange") #q=3
lines(MAsmooth(gas.test.all,2),col="green") #q=2
lines(MAsmooth(gas.test.all,1),col="blue") #q=1

###################ARIMA
plot(gas.train.all) # clear upward trend and not stationary
acf(gas.train.all) # clearly correlated and not stationary  

#Diff=1
plot(diff(gas.train.all)) # trend seems to be removed
acf(diff(gas.train.all)) 
acf(diff(gas.train.all),type="partial")
#MA1, d=1
gas.train.all.diff1ma1<-arima(gas.train.all,order=c(0,1,1),method="ML")
gas.train.all.diff1ma1
#AR1, d=1
gas.train.all.diff1ar1<-arima(gas.train.all,order=c(1,1,0),method="ML")
gas.train.all.diff1ar1
#ARMA11, d=1
gas.train.all.diff1arma11<-arima(gas.train.all,order=c(1,1,1),method="ML")
gas.train.all.diff1arma11

##Diff=2
plot(diff(gas.train.all,difference=2))
acf(diff(gas.train.all,difference=2))
acf(diff(gas.train.all,difference=2),type="partial")
#MA1, d=2
gas.train.all.diff2ma1<-arima(gas.train.all,order=c(0,2,1),method="ML")
gas.train.all.diff2ma1
#AR1, d=2
gas.train.all.diff2ar1<-arima(gas.train.all,order=c(1,2,0),method="ML")
gas.train.all.diff2ar1
#ARMA11, d=2
gas.train.all.diff2arma11<-arima(gas.train.all,order=c(1,2,1),method="ML")
gas.train.all.diff2arma11

#According to AIC, we look at MA1, d=1
resdiags(gas.train.all.diff1ma1$res)
par(mfcol=c(1,1))
#################SARIMA
plot(gas.train.all)
plot(diff(gas.train.all,lag=12))
#D=1
plot(diff(gas.train.all,lag=12)) # seasonal diff 1
acf(diff(gas.train.all,lag=12),lag.max=60) 
acf(diff(gas.train.all,lag=12),lag.max=60,type="partial")
# From the acf and pacf, we see that the series is not stationary hence we move on to the following.

##diff 1 normal 1 seasonal

plot(diff(diff(gas.train.all,lag=12))) # plot normal diff1 sesonaldiff 1
acf(diff(diff(gas.train.all,lag=12)),lag.max=60)# acf
acf(diff(diff(gas.train.all,lag=12)),lag.max=60,type="partial") #pacf
# Acf and pacf shows stationarity.Hence, we can move on to fit the data into SARIMA models.

#SARIMA0021, d=1 D=1
gas.train.all.diff1sarma0021<-arima(diff(diff(gas.train.all,lag=12)),order=c(0,0,0),seasonal=list(order=c(2,0,1), period=12),method="ML")
gas.train.all.diff1sarma0021
#SARIMA0011, d=1 D=1
gas.train.all.diff1sarma0011<-arima(diff(diff(gas.train.all,lag=12)),order=c(0,0,0),seasonal=list(order=c(1,0,1), period=12),method="ML")
gas.train.all.diff1sarma0011
#SARIMA P=2 Q=1 is preferred with a lower AIC value
#Now we will determine p and q by trial and error.
#p=1 q=0
gas.train.all.diff1sarma1021<-arima(diff(diff(gas.train.all,lag=12)),order=c(1,0,0),seasonal=list(order=c(2,0,1), period=12),method="ML")
gas.train.all.diff1sarma1021
#p=0 q=1
gas.train.all.diff1sarma0121<-arima(diff(diff(gas.train.all,lag=12)),order=c(0,0,1),seasonal=list(order=c(2,0,1), period=12),method="ML")
gas.train.all.diff1sarma0121
#p=2, q=0
gas.train.all.diff1sarma2021<-arima(diff(diff(gas.train.all,lag=12)),order=c(2,0,0),seasonal=list(order=c(2,0,1), period=12),method="ML")
gas.train.all.diff1sarma2021
#p=0, q=2
gas.train.all.diff1sarma0221<-arima(diff(diff(gas.train.all,lag=12)),order=c(0,0,2),seasonal=list(order=c(2,0,1), period=12),method="ML")
gas.train.all.diff1sarma0221
#p=0, q=3
gas.train.all.diff1sarma0321<-arima(diff(diff(gas.train.all,lag=12)),order=c(0,0,3),seasonal=list(order=c(2,0,1), period=12),method="ML")
gas.train.all.diff1sarma0321
# By observation and trial/error, we conclude that SARIMA(0,1,2)x(2,1,1)_12 is the best model out of all SARIMA models

install.packages("tseries")
library("tseries")

## d=1
acf(diff(gas.train.all))
acf(diff(gas.train.all)^2)
acf(abs(diff(gas.train.all)))
#q=0,p=1
garch01i<-garch(diff(gas.train.all),order=c(0,1),trace=F)
summary(garch01i)
resdiags(garch01i$res[3:93])

## d=1 D=1
acf(diff(diff(gas.train.all,lag=12)))
acf(diff(diff(gas.train.all,lag=12))^2)
acf(abs(diff(diff(gas.train.all,lag=12))))
#q=0,p=1
garch01ii<-garch(diff(diff(gas.train.all,lag=12)),order=c(0,1),trace=F)
summary(garch01ii)
resdiags(garch01ii$res[3:93])
par(mfcol=c(1,1))
AIC(garch01ii)
#q=1,p=1
garch11ii<-garch(diff(diff(gas.train.all,lag=12)),order=c(1,1),trace=F)
summary(garch11ii)
resdiags(garch11ii$res[3:93])
AIC(garch11ii)
#q=1,p=0
garch10ii<-garch(diff(diff(gas.train.all,lag=12)),order=c(1,0),trace=F)
summary(garch10ii)
resdiags(garch10ii$res[3:93])
AIC(garch10ii)
#q=1,p=2
garch12ii<-garch(diff(diff(gas.train.all,lag=12)),order=c(1,2),trace=F)
summary(garch12ii)
resdiags(garch12ii$res[3:93])
AIC(garch12ii)
#q=2,p=2
garch22ii<-garch(diff(diff(gas.train.all,lag=12)),order=c(2,2),trace=F)
summary(garch22ii)
resdiags(garch22ii$res[3:93])
AIC(garch22ii)

#Although garch p=1 q=1 seems to have a better residual plot,
#garch p=1 q=0 seems to be the best model out of all garch models,
#by comparing the AIC values and it has a better QQ-plot.

