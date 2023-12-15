library(forecast)
library(tseries)
library(vars)
library(openxlsx)
library(urca)
library(ggplot2)
library(midasr)
library(strucchange)
library(quantmod)
library(rusquant)
library(dplyr)
library(tsbox)
library(MSwM)


macro_data <- read.csv("gdp_russia.csv")
macro_ts <- xts(macro_data$DOLLARS, order.by = as.yearmon(macro_data$YEAR))
macro_ts <- ts_ts(macro_ts)


stock_prices <- read.csv("JNJ_1.csv")
stock_prices <- xts(stock_prices$Close, order.by = as.yearmon(stock_prices$Date))
stock_prices <- ts_ts(stock_prices)

CUSUM_result <- efp(macro_ts ~ stock_prices)
MOSUM_result <- efp(macro_ts ~ stock_prices, type = "Rec-MOSUM")

plot(CUSUM_result)
plot(MOSUM_result)