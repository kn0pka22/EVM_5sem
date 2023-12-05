library(forecast)
library(tseries)
library(openxlsx)
library(ggplot2)
library(quantmod)
library(rusquant)
library(HARModel)
library(midasr)
library(lubridate)

library(vars)


#library(urca)



x <- read.csv("C:\\Users\\maria\\Downloads\\NFLX (5).csv")
str(x)
x <- ts(x$Close)
View(x)

y <- read.csv("C:\\Users\\maria\\Downloads\\ZM (4).csv")
#View(y)
str(y)
y <- ts(y$Close)

m <- lm(y ~ x)
summary(m)
adf.test(m$residuals)

summary(ca.jo(data.frame(x, y))) # Johansen test

# #model=VECM(data.frame(Close),lag=2,r=1)


# x <- arima.sim(model = list(order = c(0,1,0)), 100)
# y <- arima.sim(model = list(order = c(0,1,0)), 100)

# plot(ts(x))
# lines(ts(y), col = "red")

# summary(lm(y ~ x))


# Real data


x <- read.csv("C:\\Users\\maria\\Downloads\\ZM (4).csv")
str(x)
#x <- ts(x$Close, start = c(2019,12), end = c(2022,12),frequency = 758)
#View(x)

 y <- read.csv("C:\\Users\\maria\\Downloads\\NFLX (5).csv")
 str(y)

# z <- read.csv("C:\\Users\\maria\\Downloads\\GOOGL.csv")
# str(z)

# u <- read.csv("C:\\Users\\maria\\Downloads\\UBER (1).csv")
# str(u)

a <- read.csv("C:\\Users\\maria\\Downloads\\AMZN (3).csv")
str(a)


#j <- read.csv("C:\\Users\\maria\\Downloads\\JET.L (4).csv")
#str(j)
#View(u)
#y <- ts(y$Close, start = c(2019,12), end = c(2022,12),frequency = 758)
#View(y)
#dd <- read.xlsx('data_ts.xlsx')
   
     #    netflix = exp(diff(log(y$Close))) - 1,
                #    google = exp(diff(log(z$Close))) - 1,
                #    uber = exp(diff(log(u$Close))) - 1,
                #    eat= exp(diff(log(j$Close))) - 1, 

# zoom = exp(diff(log(x$Close))) - 1, 

dUse <- data.frame(
                    amazon = exp(diff(log(a$Close))) - 1,
                    zoom = exp(diff(log(x$Close))) - 1,
                    netflix = exp(diff(log(y$Close))) - 1 )

m <- VAR(dUse, p = 1)
summary(m)

plot(irf(m, runs = 1000))  # Plots orthogonal by default!
plot(irf(m, cumulative = T))  
plot(irf(m, ortho = F))  # not orthogonal
plot(irf(m, impulse = "zoom"))  

fevd(m)






