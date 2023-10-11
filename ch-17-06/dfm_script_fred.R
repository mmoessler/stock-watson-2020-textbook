
#..................................................
# Set-up ----
rm(list = ls())

#..................................................
# Load functions ----
source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/ch-17-06/dfm_functions.R")

#..................................................
# Load data and estimated factors ----

# Load data
dat.00 <- read.csv(file = "https://files.stlouisfed.org/files/htdocs/fred-md/quarterly/current.csv")

dat.01 <- dat.00[-c(1, 2),]

# transform to time series
dat.01.ts <- ts(dat.01[,-c(1)], frequency = 4, start = c(1959, 1), end = c(2023, 3))
# corresponding date
date <- seq.Date(from = as.Date("1959-01-01"), to = as.Date("2023-07-01"), by = "quarter")

#..................................................
# A 0) Simple AR(2) model (for checks/comparisons) ----

# A 0.1) h=1 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(dat.01.ts[,1])
h <- 1

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

AR2.POOS.h1 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 1, n.fac = 4, fre = "quarter",
                                  all.per.sta = as.Date("1960-01-01"), all.per.end = as.Date("2022-10-01"),
                                  est.per.sta = as.Date("1960-01-01"), est.per.end = as.Date("2007-10-01"),
                                  pre.per.sta = as.Date("2007-10-01"), pre.per.end = as.Date("2022-10-01"),
                                  mis.val = TRUE, print = TRUE)
AR2.POOS.h1$RMSFE.POOS

# A 0.2) h=4 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(dat.01.ts[,1])
h <- 4

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

AR2.POOS.h4 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 4, n.fac = 4, fre = "quarter",
                                  all.per.sta = as.Date("1960-01-01"), all.per.end = as.Date("2022-10-01"),
                                  est.per.sta = as.Date("1960-01-01"), est.per.end = as.Date("2007-10-01"),
                                  pre.per.sta = as.Date("2007-10-01"), pre.per.end = as.Date("2022-10-01"),
                                  mis.val = TRUE, print = TRUE)
AR2.POOS.h4$RMSFE.POOS

# A 0.3) h=8 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(dat.01.ts[,1])
h <- 8

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

AR2.POOS.h8 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 8, n.fac = 4, fre = "quarter",
                                  all.per.sta = as.Date("1960-01-01"), all.per.end = as.Date("2022-10-01"),
                                  est.per.sta = as.Date("1960-01-01"), est.per.end = as.Date("2007-10-01"),
                                  pre.per.sta = as.Date("2007-10-01"), pre.per.end = as.Date("2022-10-01"),
                                  mis.val = TRUE, print = TRUE)
AR2.POOS.h8$RMSFE.POOS


AR2.POOS.h1$RMSFE.POOS
AR2.POOS.h4$RMSFE.POOS
AR2.POOS.h8$RMSFE.POOS
