
#..................................................
# 0) Set-up ----
rm(list = ls())

library(dplyr) # for data manipulation
library(lmtest) # for robust inference
library(sandwich)
library(dynlm) # for estimation
library(vars)

#..................................................
# 1) Load data and estimated factors ----

macro.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/us_macro_data.txt",
                        header = TRUE,
                        sep = ",",
                        colClasses = c("character","numeric","numeric","numeric"))

macro.ts <- ts(macro.dat[,-c(1)], frequency = 4, start = c(1955, 2), end = c(2017, 4))

#..................................................
# 2) Merge time series ----

# use some date sequence as basis for ts
date <- ts(seq.Date(from = as.Date("1950-01-01"), to = as.Date("2025-01-01"), by = "quarter"), frequency = 4, start = c(1950, 1), end = c(2025,1))
TT <- length(date)

# merge
data.all.ts <- cbind(date, macro.ts)

# extract time period of interest
data.all.ts <- window(data.all.ts, start = c(1959, 3), end = c(2017, 4))

# transform to data frame
data.all.df <- data.frame(date = as.Date(data.all.ts[,1]),
                          year = lubridate::year(as.Date(data.all.ts[,1])),
                          quarter = lubridate::quarter(as.Date(data.all.ts[,1])),
                          GDP     = as.numeric(data.all.ts[,2]),
                          GDPGR   = as.numeric(data.all.ts[,3]),
                          TSpread = as.numeric(data.all.ts[,4]))

#..................................................
# 3) Transform time series ----

# data manipulation using tidyverse (transform to tibble)
data.all.tb <- data.all.df %>%
  mutate(GDPGR_h01 = (400/1) * (log(GDP / lag(GDP, 1))),
         GDPGR_h04 = (400/4) * (log(GDP / lag(GDP, 4))),
         GDPGR_h08 = (400/8) * (log(GDP / lag(GDP, 8))))

# transform back to time series
data.all.ts <- ts(data.all.tb, frequency = 4, start = c(1959, 3), end = c(2017, 4))

#..................................................
# 4) VAR Model: 1981-Q1 - 2017-Q3 ----

y <- window(data.all.ts, start = c(1980,3), end = c(2017,3))[,c(5,6)]
var.res <- VAR(y = y, p = 2, type = "const")

print("--------------------------------------------------")
print("Estimates VAR model (see: S&W, 2020, p. 653)")
var.res$varresult$GDPGR
var.res$varresult$TSpread

#..................................................
# 5) Forecasts ----

# 5.1) Multiple-period forecast (iterated) (VAR Model) ----

y <- window(data.all.ts, start = c(1980,3), end = c(2017,3))[,c(5,6)]
var.res <- VAR(y = y, p = 2, type = "const")

print("--------------------------------------------------")
print("Estimates VAR model (see: S&W, 2020, p. 570)")
var.res$varresult$GDPGR
var.res$varresult$TSpread

print("--------------------------------------------------")
print("Multi-period forecast based on iterated forecast (see: S&W, 2020, p. 655)")
predict(var.res, n.ahead = 10)

# 5.2) Multiple-period forecast (direct) (ADL Model) ----

# Auxiliary regression over complete sample to get complete model data frame
tmp <-  dynlm(GDPGR ~ L(GDPGR,1) + L(GDPGR,2) + L(TSpread,1) + L(TSpread,2),
              data = data.all.ts,
              start = c(1959, 3), end = c(2017, 4))

data.mod.01 <- ts(cbind(tmp$model[,1], 1 , tmp$model[,-1]),
                  frequency = 4, start = c(1960, 1), end = c(2017, 4))

# Regression for direct forecast
adl.dynlm.tmp <- dynlm(GDPGR ~ L(GDPGR,2) + L(GDPGR,3) + L(TSpread,2) + L(TSpread,3),
                       data = data.all.ts,
                       start = c(1981,1), end = c(2017, 3))

print("--------------------------------------------------")
print("Estimates auxiliary regression for direct forecast (see: S&W, 2020, p. 657)")
adl.dynlm.tmp

# Coefficients
bet.tmp <- matrix(coefficients(adl.dynlm.tmp), nrow = 1)
# Value predictors
dat.tmp <- matrix(window(data.mod.01[,-1], start = c(2017, 4), end = c(2017, 4)), ncol = 1)
# Forecast
y.hat <- bet.tmp %*% dat.tmp

print("--------------------------------------------------")
print("Multi-period forecast based on direct forecast (see: S&W, 2020, p. 657)")
y.hat

# Actual
y.act <- matrix(window(data.mod.01[,1], start = c(2017, 4), end = c(2017, 4)), ncol = 1)
# Prediction error
u.tmp <- y.act -  y.hat

print("--------------------------------------------------")
print("Multi-period forecast error based on direct forecast")
u.tmp
