
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
# 4) Estimation (Part 1) ----

# 4.1) AR(2) Model: 1962-Q1 - 2017-Q3 ----
ar02.dynlm <- dynlm(GDPGR ~ L(GDPGR,1) + L(GDPGR,2),
                    data = data.all.ts,
                    start = c(1962, 1), end = c(2017, 3))
ct.ar02.dynlm <- coeftest(ar02.dynlm, vcov=vcovHC(ar02.dynlm, type="HC0"))

print("--------------------------------------------------")
print("Estimates AR(2) model (see: S&W, 2020, p. 567)")
ct.ar02.dynlm

# 4.2) ADL Model: 1962-Q1 - 2017-Q3 ----
adl.dynlm <- dynlm(GDPGR ~ L(GDPGR,1) + L(GDPGR,2) + L(TSpread,1) + L(TSpread,2),
                   data = data.all.ts,
                   start = c(1962, 1), end = c(2017, 3))
ct.adl.dynlm <- coeftest(adl.dynlm, vcov=vcovHC(adl.dynlm, type="HC1"))

print("--------------------------------------------------")
print("Estimates ADL model (see: S&W, 2020, p. 570)")
ct.adl.dynlm

# Compare with results from VAR for same periods
y <- window(data.all.ts, start = c(1961,3), end = c(2017,3))[,c(5,6)] # you loose 1961-Q3 and 1961-Q4 due to two lags!
var.res <- VAR(y = y, p = 2, type = "const")

print("--------------------------------------------------")
print("Estimates VAR model (compare with ADL model above)")
var.res$varresult$GDPGR

#..................................................
# 5) One-step ahead forecast (ADL Model) ----

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
# Prediction of 2017 Q4
for.17.q4.adl <- bet %*% XX

print("--------------------------------------------------")
print("Results forecast based on ADL Model (see: S&W, 2020, p. 570)")
for.17.q4.adl
