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


#Загрузка данных высокочастотных акций
#stock_prices <- read.csv("JNJ.xlsx")
stock_prices <- read.csv("C:\\Users\\maria\\Downloads\\JNJ (1).csv")
str(stock_prices)
stock_prices <- ts(stock_prices$Close)
View(stock_prices)

# Загрузка данных низкочастотных макропараметров (например, ВВП)
macro_data <- read.csv("C:\\Users\\maria\\Downloads\\gdp_russia.csv")
str(macro_data)
#macro_data <- ts(macro_data$DOLLARS)
# Преобразование данных акций во временной ряд
#stock_ts <- ts(stock_prices$Close, frequency = 252) # 252 рабочих дня в году

# Преобразование данных ВВП во временной ряд
macro_ts <- ts(macro_data$DOLLARS, frequency = 12) # 12 месяцев в году
View(macro_data)














