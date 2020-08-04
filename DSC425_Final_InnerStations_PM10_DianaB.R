
library(lubridate)
library(moments)
library(ggplot2)
library(astsa)
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)
library(fpp)
library(fpp2)
library(dplyr) # needs to be run every time you start R and want to use %>%
library(ggfortify)
library(stats)
library(fUnitRoots)
library(urca)
library(forecast)
library(lmtest)
library(astsa)
library(leaps)
library(zoo)

#================================================================
# Monthly data set: Explore the trend
#================================================================

pm10_monthly = read.csv('PM10_monthly.csv')
pollutants = read.csv('Madrid_pollutants.csv')


head(pm10_monthly)
tail(pm10_monthly)

summary(pm10_monthly)
date = pm10_monthly[,2]
pm10 = pm10_monthly[,3]
str(pm10_monthly)

#Parse_date_time(coil$date, orders = "dby")

pm10_monthly$date = parse_date_time(pm10_monthly$dateParsed, orders = "ymd")
head(pm10_monthly)
str(pm10_monthly)


qplot(pm10_monthly$date, pm10_monthly$PM10, geom="line",xlab = "Date", ylab = "PM10", 
      main = "Monthly averaged PM10 across all inner stations") + 
  scale_x_datetime(date_breaks = "1 year", date_labels = "%Y")+theme(axis.text.x = element_text(angle = 90))


#### pollutannts

head(pollutants)
summary(pollutants)

date1 = pollutants[,2]
PM10 = pollutants[,3]
NO2 = pollutants[,4]
O3 =pollutants[,5]
pollutants$date = parse_date_time(pollutants$dateParsed, orders = "ymd")

autoplot(pollutants[,c("PM10","NO2","O3")]) +
  ylab("concentration") + xlab("Year")


pollutants %>%
  as.data.frame() %>%
  ggplot(aes(x=PM10, y=NO2)) +
  ylab("N02 (monthly concentration)") +
  xlab("PM10 (monthly concentration)") +
  ggtitle("A scatter plot of PM10 changes against NO2 changes")+
  geom_point() +
  geom_smooth(method="lm", se=FALSE)

pollutants %>%
  as.data.frame() %>%
  ggplot(aes(x=PM10, y=O3)) +
  ylab("O3 (monthly concentration)") +
  xlab("PM10 (monthly concentration)") +
  ggtitle("A scatter plot of PM10 changes against O3 changes")+
  geom_point() +
  geom_smooth(method="lm", se=FALSE)

pollutants %>%
  as.data.frame() %>%
  ggplot(aes(x=NO2, y=O3)) +
  ylab("N02 (monthly concentration)") +
  xlab("PM10 (monthly concentration)") +
  ggtitle("A scatter plot of PM10 changes against NO2 changes")+
  geom_point() +
  geom_smooth(method="lm", se=FALSE)

#================================================================
# 2. Monthly data set: Histogram, QQPlot, Normality
#================================================================

qplot(pm10, main = "Histogram of average PM10 across all inner stations: monthly agg")
qqnorm(pm10)
qqline(pm10)

mean(pm10)
sd(pm10)
skewness(pm10)  
kurtosis(pm10)


# # Test the null hypothesis of perfect symmetry 
# for the distribution of rate at 5% significance level

# number of samples:
ss =skewness(pm10)/sqrt(6/length(pm10)) # Compute t-statistic
pss = 2*(1-pnorm(ss))
print(ss)
print(pss)

#Test the null hypothesis of excess kurtosis 
# equal to zero (normal tails) at 5% significance level

kk = (kurtosis(pm10)-3)/sqrt(24/length(pm10)) # Compute t-statistic
pkk = 2*(1-pnorm(kk))
print(kk)
print(pkk)

#================================================================
#3. COMPUTE ACF and PACF AND PLOT CORRELOGRAM: whole data set
#================================================================


AutoCorrelation <- acf(pm10, lag.max=60)
plot(AutoCorrelation, main = "Autocorrelation at lags 1-60 of the PM10 ")
Box.test(pm10, lag=15, type="Ljung-Box")

PAutoCorrelation <- pacf(pm10, lag.max=60)
plot(PAutoCorrelation, main = "Partial autocorrelation at lags 1-60 of the PM10 ")

acf2(pm10, max.lag = 60)
# acf decreases slowly with oscillations, we should think that this data has to be differenced


#================================================================
# 4. Train-Test Split_by month
#================================================================


# first 13 years is train data set
# last 5 years is test data set

# 18*12 = 216
# (216*0.3)/12 = 5

pm10_ts = ts(pm10_monthly$PM10, start=c(2001,1), frequency = 12)
test = window(pm10_ts, start = c(2014,01))
train = window(pm10_ts, start= c(2001,01),end=c(2013,12))

no2_ts = ts(pollutants$NO2, start=c(2001,1), frequency = 12)
no2_t = window(no2_ts, start = c(2014,01))
no2 = window(no2_ts, start= c(2001,01),end=c(2013,12))

AutoCorrelation <- acf(train, lag.max=60)
plot(AutoCorrelation, main = "ACF at lags 1-60 of the PM10 train_data ")
Box.test(train, lag=15, type="Ljung-Box")

PAutoCorrelation <- pacf(train, lag.max=60)
plot(PAutoCorrelation, main = "PACF at lags 1-60 of the PM10 train_data ")
test


#================================================================
#5 Checking stationarity : Original data
#================================================================

# The null hypothesis for the KPSS test is that the data are stationary. 
# For this test, we do NOT want to reject the null hypothesis. 
# In other words, we want the p-value to be greater than 0.05 not less than 0.05.

kpss_train <- kpss.test(train)
kpss_train # p = 0.01 # data is not stationary p<0.05

#================================================================
#5.1 Data transormation: First order of Diference by month
#================================================================

#  HOW TO MAKE DATA STATIONARY

# First Oreder of DIFF
train_diff <- diff(train)

AutoCorrelation <- acf(train_diff, lag.max=60)
plot(AutoCorrelation, main = "ACF : PM10 train_data first order of difference")
Box.test(train_diff, lag=15, type="Ljung-Box")

PAutoCorrelation <- pacf(train_diff, lag.max=60)
plot(PAutoCorrelation, main = "PACF: PM10 train_data first order of difference")

plot(diff(train), geom="line",xlab = "Date", ylab = "PM10 diff", 
      main = "Time plot-first order of difference of PM10") 


m = ar(train, methode = 'mle')
m
adf_train <- urca::ur.df(train, type="trend", lags=11)
summary(adf_train)
#================================================================
#5.2 Checking stationarity First Order difference
#================================================================

# The null hypothesis for the KPSS test is that the data are stationary. 
# For this test, we do NOT want to reject the null hypothesis. 
# In other words, we want the p-value to be greater than 0.05 not less than 0.05.

kpss_train1 <- kpss.test(diff(train))
kpss_train1 # p = 0.1 # data is stationary

train %>% diff()%>% acf()
train %>% diff()%>% pacf()


m = ar(diff(train), methode = 'mle')
m
adf_train <- urca::ur.df(diff(train), type="trend", lags=10)
summary(adf_train)

#================================================================
#5.2 Checking stationarity Second Order difference
#================================================================

# The null hypothesis for the KPSS test is that the data are stationary. 
# For this test, we do NOT want to reject the null hypothesis. 
# In other words, we want the p-value to be greater than 0.05 not less than 0.05.

kpss_train2 <- kpss.test(diff(train_diff))
kpss_train2 # p = 0.1 # data is stationary

train %>% diff()%>%diff()%>% acf()
train %>% diff()%>%diff()%>% pacf()


m = ar(diff(train_diff), methode = 'mle')
m
adf_train <- urca::ur.df(diff(train_diff), type="trend", lags=12)
summary(adf_train)


AutoCorrelation <- acf(diff(train_diff), lag.max=60)
plot(AutoCorrelation, main = "ACF : PM10 train_data second order of difference")
Box.test(diff(train_diff), lag=15, type="Ljung-Box")

PAutoCorrelation <- pacf(diff(train_diff), lag.max=60)
plot(PAutoCorrelation, main = "PACF: PM10 train_data second order of difference")

plot(diff(train_diff), geom="line",xlab = "Date", ylab = "PM10 diff", 
     main = "Time plot-second order of difference of PM10") 

#================================================================
#Model building: 
#================================================================

# ARIMA

m = auto.arima(train, ic=c("bic"), seasonal=F)
coeftest(m)

library(TSA)
M0 = Arima(train, order = c(5,1,1))
coeftest(M0)
M0_1 = Arima(train, order = c(5,1,1), fixed = c(NA, NA, 0, 0, 0, NA))
coeftest(M0_1)
M0_1
#AIC = 1016 BIC = 1028.77


M1 = Arima(train, order = c(3,1,2))
coeftest(M1)
M1_1 = Arima(train, order = c(3,1,2), fixed = c(NA, 0, NA, NA, NA))
coeftest(M1_1)
M1_1  # AIC = 1005.44, BIC = 1020.65


M2 = Arima(train, order = c(5,2,1))
coeftest(M2)
M2_1 = Arima(train, order = c(5,2,1), fixed = c(NA, 0, NA, NA, NA, NA))
M2_1  # AIC = 1028, BIC = 1046.97


plot(M2_1$resid, type='l', main = "ARIMA(5,2,1) Residuals")
plot(M1_1$resid, type='l', main = "ARIMA(3,1,2) Residuals")

# LINEAR REGRESSION

M3 <- tslm(train~season+trend)
M3
summary(M3)

M4 <- tslm(train ~ season + trend + no2)
M4
summary(M4)

M5 <- tslm(train~season+trend, lambda = 0)
M5
summary(M5)

M6 <- tslm(train~season+trend+no2, lambda = 0)
M6
summary(M6)


autoplot(BoxCox(train~season+trend, lambda = 2))


#================================================================
# Residuals: Ljung-Box
#================================================================

AutoCorrelation <- acf(M0_1$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ARIMA(5,1,1)")
Box.test(M0_1$resi, lag=15, type="Ljung-Box")
Box.test(M0_1$resi, lag=12, type="Ljung-Box")
Box.test(M0_1$resi, lag=7, type="Ljung-Box")
Box.test(M0_1$resi, lag=5, type="Ljung-Box")

AutoCorrelation <- acf(M1_1$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ARIMA(3,1,2)")
Box.test(M1_1$resi, lag=15, type="Ljung-Box")
Box.test(M1_1$resi, lag=12, type="Ljung-Box")
Box.test(M1_1$resi, lag=7, type="Ljung-Box")
Box.test(M1_1$resi, lag=5, type="Ljung-Box")

AutoCorrelation <- acf(M2_1$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ARIMA(5,2,1)")
Box.test(M2_1$resi, lag=15, type="Ljung-Box")
Box.test(M2_1$resi, lag=12, type="Ljung-Box")
Box.test(M2_1$resi, lag=7, type="Ljung-Box")
Box.test(M2_1$resi, lag=5, type="Ljung-Box")

AutoCorrelation <- acf(M3$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ~season+trend")

AutoCorrelation <- acf(M4$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ~season+trend+no2")

AutoCorrelation <- acf(M5$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ~season+trend, lambda = 0")
Box.test(M5$resi, lag=5, type="Ljung-Box")

AutoCorrelation <- acf(M6$resid, lag.max=15)
plot(AutoCorrelation, main = "Autocorrelation of residuals ~season+trend+no2, lambda = 0")



#================================================================
# Residuals Normality: Shapiro-Wilk's test for normality
#================================================================

#Reject the null hypothesis H0 that data come from normal distribution if p-value 
# is less than 0.05. 
# Since p value is higher than 0.05 for all three models we fail to reject H0. 
# Residuals do come from normal distribution.

normalTest(M0_1$residuals, method=c("sw")) 
normalTest(M1_1$residuals, method=c("sw"))
normalTest(M2_1$residuals, method=c("sw")) 

normalTest(M3$residuals, method=c("sw")) 
normalTest(M4$residuals, method=c("sw"))
normalTest(M5$residuals, method=c("sw"))
normalTest(M6$residuals, method=c("sw"))

#================================================================
# Residuals Normality: Perform Jarque-Bera normality test
#================================================================

normalTest(M0_1$residuals, method=c("jb")) 
normalTest(M1_1$residuals, method=c("jb")) 
normalTest(M2_1$residuals, method=c("jb")) 


normalTest(M3$residuals, method=c("jb")) 
normalTest(M4$residuals, method=c("jb")) 
normalTest(M5$residuals, method=c("jb"))
normalTest(M6$residuals, method=c("jb"))

#================================================================
# Residuals Normality: QQ plot Histogram
#================================================================

qplot(M0_1$residuals, main = "Histogram of ARIMA(5,1,1) residuals")
qqnorm(M0_1$residuals)
qqline(M0_1$residuals, col = 2) 

gghistogram(M1_1$residuals) + ggtitle("Histogram of ARIMA(3,1,2) residuals")
gghistogram(M2_1$residuals) + ggtitle("Histogram of ARIMA(5,2,1) residuals")

qplot(M1_1$residuals, main = "Histogram of ARIMA(3,1,2) residuals")
qqnorm(M1_1$residuals, main = "Q-Q plot of ARIMA(3,1,2) residuals")
qqline(M1_1$residuals, col = 2)

qplot(M2_1$residuals, main = "Histogram of ARIMA(5,2,1) residuals")
qqnorm(M2_1$residuals, main = "Q-Q plot of ARIMA(5,2,1) residuals")
qqline(M2_1$residuals, col = 2) 

qplot(M3$residuals, main = "Histogram of ~season+trend")
qqnorm(M3$residuals, main = "Q-Q plot of ~season+trend")
qqline(M3$residuals, col = 2) 

qplot(M4$residuals, main = "Histogram of ~season+trend+no2")
qqnorm(M4$residuals, main = "Q-Q plot of ~season+trend+no2")
qqline(M4$residuals, col = 2) 

qplot(M5$residuals, main = "Histogram of ~season+trend, lambda =0")
qqnorm(M5$residuals, main = "Q-Q plot of ~season+trend, lambda = 0 ")
qqline(M5$residuals, col = 2) 

qplot(M6$residuals, main = "Histogram of ~season+trend+no2, lambda =0")
qqnorm(M6$residuals, main = "Q-Q plot of ~season+trend+no2, lambda = 0 ")
qqline(M6$residuals, col = 2) 

#================================================================
#Plot Residuals
#================================================================

plot(M2_1$resid, type='p', main = "ARIMA(5,2,1) Residuals", xlab = "Date", ylab = "Residuals")
abline(h=0, col="grey")
plot(M1_1$resid, type='p', main = "ARIMA(3,1,2) Residuals", xlab = "Date", ylab = "Residuals")
abline(h=0, col="grey")

plot(M3$resid, type='l', main = "~season+trend", xlab = "Date", ylab = "Residuals")
abline(h=0, col="grey")
plot(M4$resid, type='l', main = "~season+trend+no2", xlab = "Date", ylab = "Residuals")
abline(h=0, col="grey")
plot(M5$resid, type='l', main = "~season+trend, lambda =0", xlab = "Date", ylab = "Residuals")
abline(h=0, col="grey")
plot(M6$resid, type='l', main = "~season+trend+no2, lambda =0", xlab = "Date", ylab = "Residuals")
abline(h=0, col="grey")

#================================================================
#Checking Stationarity of Residuals
#================================================================


adf_train <- urca::ur.df(M1_1$residuals, type="trend", lags=0)
summary(adf_train)

# The intercept and tt estimates indicate where 
# there is a non-zero level (intercept) or linear trend (tt).
# We can see that trend and intercept p values are higher than 0.05
# therfore they are not significant. Only z.lag1 is significant
# and it's value is gamma = -1.053 where we want  -1< gamma+1 <1
# that is true therfore data is stationary

adf_train <- urca::ur.df(M3$residuals, type="trend", lags=5)
summary(adf_train)

adf_train <- urca::ur.df(M4$residuals, type="trend", lags=5)
summary(adf_train)

adf_train <- urca::ur.df(M5$residuals, type="trend", lags=5)
summary(adf_train)

adf_train <- urca::ur.df(M6$residuals, type="trend", lags=5)
summary(adf_train)

checkresiduals(M3$residuals)
checkresiduals(M4$residuals)
checkresiduals(M5$residuals)
checkresiduals(M6$residuals)



#================================================================
#Forecast
#================================================================


M1_test <- Arima(test, model=M1_1) # RMSE 4.72, MAPE: 24.6
accuracy(M1_test)
plot(M1_test$residuals, main = 'Residuals on test data: ARIMA (3,1,2)', xlab = 'Time', ylab = 'Errors')

M2_test <- Arima(test, model=M2_1) # RMSE 4.94, MAE: 23.9
accuracy(M2_test)

f1 <- forecast(M1_1,h=53)
forecast(M1_1, h = 53) 
lines(test)
accuracy(f1$mean,test) # RMSE_train = 5.9, test = 5.6 MAPE_train = 17.6, MAPE_test = 29.37

f2 <- forecast(M2_1,h=53)
forecast(M2_1, h = 53) 
lines(test)
accuracy(f2$mean,test)  # RMSE_train = 6.37, test = 5.2 MAPE_train = 17.89, MAPE_test = 28.35


### FORECAST LINEAR REGRESSION

f3 <- forecast(M3,h=53)
forecast(M3, h = 53) 
lines(test)

f4 <- forecast(M4,h=53,newdata = data.frame(no2=no2_t))
forecast(M4, h = 53,newdata = data.frame(no2=no2_t)) 
lines(test)

f5 <- forecast(M5,h=53)
forecast(M5, h = 53) 
lines(test)

f6 <- forecast(M6,h=53,newdata = data.frame(no2=no2_t))
forecast(M6, h = 53,newdata = data.frame(no2=no2_t)) 
lines(test)




plot(f6, main = 'PM10 forecast(Linear Regression ~season+trend_no2, lambda = 0)', xlab="Year", ylab="PM10 values")
lines(test)


plot(f3, main = 'PM10 forecast (Linear Regression ~season+trend)', xlab="Year", ylab="PM10 values")
lines(test)


plot(f5, main ='PM10 forecast(Linear Regression ~season+trend, log)', xlab="Year", ylab="PM10 values")
lines(test)

plot(f4, main ='PM10 forecast(Linear Regression ~season +trend+no2)', xlab="Year", ylab="PM10 values")
lines(test)


autoplot(f1) + labs(title = 'PM10 forecast(ARIMA(3,2,1))', x='Year', y='PM10 values')
plot(f1, main = 'PM10 forecast(ARIMA(3,1,2))', xlab="Year", ylab="PM10 values")
lines(test)

plot(forecast(auto.arima(train),h=53))

plot(f2, main = 'PM10 forecast(ARIMA(5,2,1))', xlab="Year", ylab="PM10 values")
lines(test)

accuracy(f1, test)
accuracy(f2, test)

my_fc$residuals
test.RMSE<-(sqrt(mean((my_fc$residuals)^2))) # 5.9
test.RMSE
test.MAE<-(mean(abs(my_fc$residuals))) # 4.62
test.MAE


my_fc$residuals
test.RMSE<-(sqrt(mean((my_fc$residuals)^2))) # 5.9
test.RMSE
test.MAE<-(mean(abs(my_fc$residuals))) # 4.62
test.MAE

# residuals in original scale
errortestM3 <- test - f3$mean
plot(errortestM3, type="o", bty="l", main="Residuals in test part: ~trend, season", xlab = "Time", ylab = "Errors")

errortestM5 <- test - f5$mean
plot(errortestM5, type="o", bty="l", main="Residuals in test part: ~trend, season, lambda = 0", xlab = "Time", ylab = "Errors")

errortestM6 <- test - f6$mean
plot(errortestM6, type="o", bty="l", main="Residuals in test part: ~trend, season,No2 lambda = 0", xlab = "Time", ylab = "Errors")


error_test <- test - my_fc$mean
plot(error_test, type="o", bty="l", main="Residuals in test part: ~trend, season", xlab = "Time", ylab = "Errors")
accuracy(my_fc$mean, test)


accuracy(f3, test) # RMSE_train = 5.9, RMSE_test = 4.98, MAPE_train = 16.85, MAPE_test = 26.5
accuracy(f4, test) # RMSE_train = 5.3, RMSE_test = 5.69, MAPE_train = 15.11, MAPE_test = 34.57
accuracy(f5, test) # RMSE_train = 6.11, RMSE_test = 4.94, MAPE_train = 16.84, MAPE_test = 26.99
accuracy(f6, test) # RMSE_train = 5.52, RMSE_test = 5.3, MAPE_train = 14.89, MAPE_test = 32.06









