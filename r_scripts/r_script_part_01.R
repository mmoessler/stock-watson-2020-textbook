
rm(list = ls())

library(dynlm)

library(lmtest) # for robust inference
library(sandwich)

library(vars)

library(dplyr)



#..................................................
# 1) Load data and estimated factors ----

# Observable macro variables
# macro.dat <- read.table("../data/us_macro_data.txt",
#                         header = TRUE,
#                         sep = ",",
#                         colClasses = c("character","numeric","numeric","numeric"))

macro.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/us_macro_data.txt",
                        header = TRUE,
                        sep = ",",
                        colClasses = c("character","numeric","numeric","numeric"))

# head(macro.dat)

macro.ts <- ts(macro.dat[,-c(1)], frequency = 4, start = c(1955, 2), end = c(2017, 4))
# dim(macro.ts)



#..................................................
# 2) Merge time series ----

# use some date sequence as basis for ts
date <- ts(seq.Date(from = as.Date("1950-01-01"), to = as.Date("2025-01-01"), by = "quarter"), frequency = 4, start = c(1950, 1), end = c(2025,1))
TT <- length(date)

# merge
data.all.ts <- cbind(date, macro.ts)
# head(data.all.ts)
# tail(data.all.ts)

# extract time period of interest
data.all.ts <- window(data.all.ts, start = c(1959, 3), end = c(2017, 4))
# head(data.all.ts)

# transform to data frame
data.all.df <- data.frame(date = as.Date(data.all.ts[,1]),
                          year = lubridate::year(as.Date(data.all.ts[,1])),
                          quarter = lubridate::quarter(as.Date(data.all.ts[,1])),
                          GDP     = as.numeric(data.all.ts[,2]),
                          GDPGR   = as.numeric(data.all.ts[,3]),
                          TSpread = as.numeric(data.all.ts[,4]))
# head(data.all.df)



#..................................................
# 3) Transform time series ----

# data manipulation using tidyverse (transform to tibble)
data.all.tb <- data.all.df %>%
  mutate(GDPGR_h01 = (400/1) * (log(GDP / lag(GDP, 1))),
         GDPGR_h04 = (400/4) * (log(GDP / lag(GDP, 4))),
         GDPGR_h08 = (400/8) * (log(GDP / lag(GDP, 8))))
# head(data.all.tb,10)

# transform back to time series
data.all.ts <- ts(data.all.tb, frequency = 4, start = c(1959, 3), end = c(2017, 4))
# head(data.all.ts)



#..................................................
# 4) Estimation (Part 1) ----

# 4.1) AR(2) Model: 1962-Q1 - 2017-Q3 ----
ar02.dynlm <- dynlm(GDPGR ~ L(GDPGR,1) + L(GDPGR,2),
                    data = data.all.ts,
                    start = c(1962, 1), end = c(2017, 3))
ct.ar02.dynlm <- coeftest(ar02.dynlm, vcov=vcovHC(ar02.dynlm, type="HC0"))

print("--------------------------------------------------")
print("Results AR(2) Model (see: S&W, 2020, p. 567)")
ct.ar02.dynlm
print("--------------------------------------------------")
# -> see: S&W, 2020, p. 567

# 4.2) ADL Model: 1962-Q1 - 2017-Q3 ----
adl.dynlm <- dynlm(GDPGR ~ L(GDPGR,1) + L(GDPGR,2) + L(TSpread,1) + L(TSpread,2),
                   data = data.all.ts,
                   start = c(1962, 1), end = c(2017, 3))
ct.adl.dynlm <- coeftest(adl.dynlm, vcov=vcovHC(adl.dynlm, type="HC1"))

print("--------------------------------------------------")
print("Results ADL Model (see: S&W, 2020, p. 570)")
ct.adl.dynlm
print("--------------------------------------------------")
# -> see: S&W, 2020, p. 570

# Compare with results from VAR for same periods
y <- window(data.all.ts, start = c(1961,3), end = c(2017,3))[,c(5,6)] # you loose 1961-Q3 and 1961-Q4 due to two lags!
var.res <- VAR(y = y, p = 2, type = "const")

print("--------------------------------------------------")
print("Results VAR Model (compare with ADL model above)")
var.res$varresult$GDPGR
print("--------------------------------------------------")
# -> see: adl.dynlm 

# head(model.frame(var.res$varresult$GDPGR))
# tail(model.frame(var.res$varresult$GDPGR))
# nrow(model.frame(var.res$varresult$GDPGR))
# head(model.frame(adl.dynlm))
# tail(model.frame(adl.dynlm))
# nrow(model.frame(adl.dynlm))

# 4.3) VAR Model: 1981-Q1 - 2017-Q3 ----
y <- window(data.all.ts, start = c(1980,3), end = c(2017,3))[,c(5,6)]
var.res <- VAR(y = y, p = 2, type = "const")

print("--------------------------------------------------")
print("Results ADL Model (see: S&W, 2020, p. 653)")
var.res$varresult$GDPGR
print("--------------------------------------------------")
# -> see: S&W, 2020, p. 653 (type for TSpread t-2 coefficient?)



#..................................................
# 5) Predictions (Part 1) ----

# 5.1) One-step ahead forecast (ADL Model) ----

# Select GDP growth as ts variable
gdp.grw.ts <- macro.ts[,2]
# Select term spread  as ts variable
ter.spr.ts <- macro.ts[,3]

# Use 2017 Q3 and 2017 Q2
X1 <- matrix(rev(window(gdp.grw.ts,start=c(2017,2),end=c(2017,3))),ncol=1) # I reverse the elements s.t. they fit to the order of the coefficients!
X2 <- matrix(rev(window(ter.spr.ts,start=c(2017,2),end=c(2017,3))),ncol=1)
XX <- rbind(1,X1,X2)
# Use coefficients of ar01.dynlm
bet <- matrix(adl.dynlm$coefficients,nrow=1)
# # Check
# XX
# bet
# Prediction of 2017 Q4
for.17.q4.adl <- bet %*% XX
# for.17.q4.adl
# -> see: S&W, 2020, p. 570

# 5.2) Multiple-step ahead forecast (iterated) (VAR Model) ----

y <- window(data.all.ts, start = c(1980,3), end = c(2017,3))[,c(5,6)]
var.res <- VAR(y = y, p = 2, type = "const")
# var.res$varresult$GDPGR
# var.res$varresult$TSpread
# -> see: S&W, 2020, p. 570
# predict(var.res, n.ahead = 10)
# -> see: S&W, 2020, p. 655

# 5.3) Multiple-step ahead forecast (direct) (ADL Model) ----

# Auxiliary regression over complete sample to get complete model data frame
# tmp <-  dynlm(GDPGR ~ L(GDPGR,2) + L(GDPGR,3) + L(TSpread,2) + L(TSpread,3),
#               data = data.all.ts,
#               start = c(1959, 3), end = c(2017, 4))
tmp <-  dynlm(GDPGR ~ L(GDPGR,1) + L(GDPGR,2) + L(TSpread,1) + L(TSpread,2),
              data = data.all.ts,
              start = c(1959, 3), end = c(2017, 4))

data.mod.01 <- ts(cbind(tmp$model[,1], 1 , tmp$model[,-1]),
                  frequency = 4, start = c(1960, 1), end = c(2017, 4))

# head(data.mod.01)
# class(data.mod.01)
# tsp(data.mod.01)

# Regression for direct forecast
adl.dynlm.tmp <- dynlm(GDPGR ~ L(GDPGR,2) + L(GDPGR,3) + L(TSpread,2) + L(TSpread,3),
                       data = data.all.ts,
                       start = c(1981,1), end = c(2017, 3))
adl.dynlm.tmp
# -> see: S&W, 2020, p. 657

# Coefficients
bet.tmp <- matrix(coefficients(adl.dynlm.tmp), nrow = 1)
bet.tmp

# Value predictors
dat.tmp <- matrix(window(data.mod.01[,-1], start = c(2017, 4), end = c(2017, 4)), ncol = 1)
# dat.tmp

# Prediction
y.hat <- bet.tmp %*% dat.tmp
# y.hat
# -> see: S&W, 2020 p. 657

# Actual
y.act <- matrix(window(data.mod.01[,1], start = c(2017, 4), end = c(2017, 4)), ncol = 1)
# y.act

# Prediction error
u.tmp <- y.act -  y.hat
# u.tmp
