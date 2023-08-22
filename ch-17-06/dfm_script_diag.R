
#..................................................
# Set-up ----
rm(list = ls())

# to save the results into .json
library(jsonlite)

#..................................................
# Load functions ----
source("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/ch-17-06/dfm_functions.R")

#..................................................
# Load data and estimated factors ----

# Load observable macro variables
macro.01.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/us_macro_data.txt",
                           header = TRUE,
                           sep = ",",
                           colClasses = c("character", "numeric", "numeric", "numeric"))

macro.01.ts <- ts(macro.01.dat[,-c(1)], frequency = 4, start = c(1955, 2), end = c(2017, 4))

# Load all observable macro variables
macro.02.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/us_all_macro_data.txt",
                           header = TRUE,
                           sep = ",")

macro.02.ts <- ts(macro.02.dat[,-c(1)], frequency = 4, start = c(1955, 2), end = c(2017, 4))
colnames(macro.02.ts) <- paste0("V", seq(1,ncol(macro.02.ts)))

#..................................................
# Estimate unobserved factors based on all macro variables ----

# determine the start and end of period for factor estimation
dat.seq <- seq.Date(from = as.Date("1955-04-01"), to = as.Date("2017-10-01"), by = "quarter")
dat.seq[3]   # start used
dat.seq[236] # end used

dat.tmp <- window(macro.02.ts, start = c(1955, 4), end = c(2014, 1))
xdata <- matrix(as.vector(dat.tmp), nrow = nrow(dat.tmp), ncol = ncol(dat.tmp)) # transform to "matrix"/"array"

# factor estimation using variables without missing values based on PCA

# fa based on ordinary pca
X <- as.matrix(scale(packr(xdata), center = TRUE, scale = FALSE))
n.fac <- 4

# # 1) based on princomp (i.e., eigenvalue decomposition)
# XX <- X %*% t(X)
# Sig <- dim(X)[2]^(-1) * XX
# pca.res <- princomp(Sig)
# F.est <- pca.res$scores[,1:n.fac]

# # 2) based on prcomp (i.e., svd)
# pca.res <- prcomp(X, retx = TRUE, rank. = n.fac)
# F.est <- pca.res$x[,1:n.fac]

# 3) based on svd directly
sv <- svd(X)
scores <- X %*% sv$v
F.est <- scores[,1:n.fac]

factor.est.pca <- F.est
factor.est.pca.ts <- ts(factor.est.pca, frequency = 4, start = c(1959, 3), end = c(2017, 4))

# factor estimation using all variables based on iterative LS
factor.est.ite <- factor_estimation_ls(xdata = xdata, n.fac = 4, nt.min = 20, init = "svd")
factor.est.ite.ts <- ts(factor.est.ite$fac, frequency = 4, start = c(1959, 3), end = c(2017, 4))

#..................................................
# Merge time series ----

# use some date sequence as basis for ts
date <- ts(seq.Date(from = as.Date("1950-01-01"), to = as.Date("2025-01-01"), by = "quarter"), frequency = 4, start = c(1950, 1), end = c(2025,1))
TT <- length(date)

# merge all variables/factors
data.all.ts <- cbind(date, macro.01.ts, macro.02.ts, factor.est.pca.ts, factor.est.ite.ts) # based on estimated factors

# extract time period of interest
data.all.ts <- window(data.all.ts, start = c(1959, 3), end = c(2017, 4))
# corresponding date
date <- seq.Date(from = as.Date("1959-07-01"), to = as.Date("2017-10-01"), by = "quarter")

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
GDP <- as.numeric(data.all.ts[,2])
h <- 1

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/AR2_POOS_h1_diag.txt")
AR2.POOS.h1 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 1, n.fac = 4, fre = "quarter",
                                  all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                  est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                  pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                  mis.val = TRUE, print = TRUE)
AR2.POOS.h1$RMSFE.POOS
sink()

jsonlite::write_json(AR2.POOS.h1, "./ch-17-06/results/AR2_POOS_h1.json")

# A 0.2) h=4 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 4

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/AR2_POOS_h4_diag.txt")
AR2.POOS.h4 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 4, n.fac = 4, fre = "quarter",
                                  all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                  est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                  pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                  mis.val = TRUE, print = TRUE)
AR2.POOS.h4$RMSFE.POOS
sink()

jsonlite::write_json(AR2.POOS.h4, "./ch-17-06/results/AR2_POOS_h4.json")

# A 0.3) h=8 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 8

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/AR2_POOS_h8_diag.txt")
AR2.POOS.h8 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 8, n.fac = 4, fre = "quarter",
                                  all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                  est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                  pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                  mis.val = TRUE, print = TRUE)
AR2.POOS.h8$RMSFE.POOS
sink()

jsonlite::write_json(AR2.POOS.h8, "./ch-17-06/results/AR2_POOS_h8.json")

#..................................................
# A 1) Unfair analysis on variables without missing values ----

# A 1.1) h=1 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 1

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
F1 <- as.numeric(data.all.ts[,136])
F2 <- as.numeric(data.all.ts[,137])
F3 <- as.numeric(data.all.ts[,138])
F4 <- as.numeric(data.all.ts[,139])

data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      F1_h_l1 = lag_fun(F1, h),
                      F2_h_l1 = lag_fun(F2, h),
                      F3_h_l1 = lag_fun(F3, h),
                      F4_h_l1 = lag_fun(F4, h),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h1_01_diag.txt")
FVAR.POOS.r4.h1.01 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 1, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h1.01$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h1.01, "./ch-17-06/results/FVAR_POOS_r4_h1_01.json")

# A 1.2) h=4 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 4

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
F1 <- as.numeric(data.all.ts[,136])
F2 <- as.numeric(data.all.ts[,137])
F3 <- as.numeric(data.all.ts[,138])
F4 <- as.numeric(data.all.ts[,139])

data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      F1_h_l1 = lag_fun(F1, h),
                      F2_h_l1 = lag_fun(F2, h),
                      F3_h_l1 = lag_fun(F3, h),
                      F4_h_l1 = lag_fun(F4, h),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h4_01_diag.txt")
FVAR.POOS.r4.h4.01 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 4, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h4.01$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h4.01, "./ch-17-06/results/FVAR_POOS_r4_h4_01.json")

# A 1.3) h=8 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 8

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
F1 <- as.numeric(data.all.ts[,136])
F2 <- as.numeric(data.all.ts[,137])
F3 <- as.numeric(data.all.ts[,138])
F4 <- as.numeric(data.all.ts[,139])

data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      F1_h_l1 = lag_fun(F1, h),
                      F2_h_l1 = lag_fun(F2, h),
                      F3_h_l1 = lag_fun(F3, h),
                      F4_h_l1 = lag_fun(F4, h),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h8_01_diag.txt")
FVAR.POOS.r4.h8.01 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 8, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h8.01$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h8.01, "./ch-17-06/results/FVAR_POOS_r4_h8_01.json")

#..................................................
# A 2) Unfair analysis based on all observations ----

# A 2.1) h=1 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 1

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
F1 <- as.numeric(data.all.ts[,140])
F2 <- as.numeric(data.all.ts[,141])
F3 <- as.numeric(data.all.ts[,142])
F4 <- as.numeric(data.all.ts[,143])

data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      F1_h_l1 = lag_fun(F1, h),
                      F2_h_l1 = lag_fun(F2, h),
                      F3_h_l1 = lag_fun(F3, h),
                      F4_h_l1 = lag_fun(F4, h),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h1_02_diag.txt")
FVAR.POOS.r4.h1.02 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 1, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h1.02$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h1.02, "./ch-17-06/results/FVAR_POOS_r4_h1_02.json")

# A 2.2) h=4 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 4

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
F1 <- as.numeric(data.all.ts[,140])
F2 <- as.numeric(data.all.ts[,141])
F3 <- as.numeric(data.all.ts[,142])
F4 <- as.numeric(data.all.ts[,143])

data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      F1_h_l1 = lag_fun(F1, h),
                      F2_h_l1 = lag_fun(F2, h),
                      F3_h_l1 = lag_fun(F3, h),
                      F4_h_l1 = lag_fun(F4, h),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h4_02_diag.txt")
FVAR.POOS.r4.h4.02 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 4, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h4.02$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h4.02, "./ch-17-06/results/FVAR_POOS_r4_h4_02.json")

# A 2.3) h=8 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 8

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
F1 <- as.numeric(data.all.ts[,140])
F2 <- as.numeric(data.all.ts[,141])
F3 <- as.numeric(data.all.ts[,142])
F4 <- as.numeric(data.all.ts[,143])

data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      F1_h_l1 = lag_fun(F1, h),
                      F2_h_l1 = lag_fun(F2, h),
                      F3_h_l1 = lag_fun(F3, h),
                      F4_h_l1 = lag_fun(F4, h),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- NULL

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h8_02_diag.txt")
FVAR.POOS.r4.h8.02 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 8, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h8.02$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h8.02, "./ch-17-06/results/FVAR_POOS_r4_h8_02.json")

#..................................................
# A 3) Fair analysis based on variables without missing values ----

# A 3.1) h=1 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 1

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- data.all.ts[,which(colnames(data.all.ts) %in% paste0("macro.02.ts.V", seq(1,ncol(macro.02.ts))))]
colnames(data.03) <- paste0("V", seq(1,ncol(data.03)))

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h1_03_diag.txt")
FVAR.POOS.r4.h1.03 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 1, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = FALSE, print = TRUE)
FVAR.POOS.r4.h1.03$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h1.03, "./ch-17-06/results/FVAR_POOS_r4_h1_03.json")

# A 3.2) h=4 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 4

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- data.all.ts[,which(colnames(data.all.ts) %in% paste0("macro.02.ts.V", seq(1,ncol(macro.02.ts))))]
colnames(data.03) <- paste0("V", seq(1,ncol(data.03)))

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h4_03_diag.txt")
FVAR.POOS.r4.h4.03 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 4, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = FALSE, print = TRUE)
FVAR.POOS.r4.h4.03$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h4.03, "./ch-17-06/results/FVAR_POOS_r4_h4_03.json")

# A 3.3) h=8 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 8

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- data.all.ts[,which(colnames(data.all.ts) %in% paste0("macro.02.ts.V", seq(1,ncol(macro.02.ts))))]
colnames(data.03) <- paste0("V", seq(1,ncol(data.03)))

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h8_03_diag.txt")
FVAR.POOS.r4.h8.03 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 8, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = FALSE, print = TRUE)
FVAR.POOS.r4.h8.03$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h8.03, "./ch-17-06/results/FVAR_POOS_r4_h8_03.json")

#..................................................
# A 4) Fair analysis based on all variables ----

# A 4.1) h=1 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 1

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- data.all.ts[,which(colnames(data.all.ts) %in% paste0("macro.02.ts.V", seq(1,ncol(macro.02.ts))))]
colnames(data.03) <- paste0("V", seq(1,ncol(data.03)))

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h1_04_diag.txt")
FVAR.POOS.r4.h1.04 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 1, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h1.04$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h1.04, "./ch-17-06/results/FVAR_POOS_r4_h1_04.json")

# A 4.2) h=4 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 4

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- data.all.ts[,which(colnames(data.all.ts) %in% paste0("macro.02.ts.V", seq(1,ncol(macro.02.ts))))]
colnames(data.03) <- paste0("V", seq(1,ncol(data.03)))

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h4_04_diag.txt")
FVAR.POOS.r4.h4.04 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 4, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h4.04$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h4.04, "./ch-17-06/results/FVAR_POOS_r4_h4_04.json")

# A 4.3) h=8 ----

# dates data frame (date column is important!)
date.df <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      month = format(date, "%m"),
                      day = format(date, "%d"))

# target variable (data frame)
GDP <- as.numeric(data.all.ts[,2])
h <- 8

data.01 <- data.frame(GDPGR_h = (400/h) * (log(GDP / lag_fun(GDP, h))))

# observed predictors (variables) (data frame)
data.02 <- data.frame(GDPGR_h_l1 = lag_fun(log(GDP / lag_fun(GDP, 1)), h),
                      GDPGR_h_l2 = lag_fun(log(GDP / lag_fun(GDP, 1)), h + 1),
                      CONST = 1)

# observed predictors for (fair) construction of unobserved factors (matrix/array)
data.03 <- data.all.ts[,which(colnames(data.all.ts) %in% paste0("macro.02.ts.V", seq(1,ncol(macro.02.ts))))]
colnames(data.03) <- paste0("V", seq(1,ncol(data.03)))

sink(file = "./ch-17-06/diagnostics/FVAR_POOS_r4_h8_04_diag.txt")
FVAR.POOS.r4.h8.04 <- FVAR_POOS_function(date.df = date.df, data.01 = data.01, data.02 = data.02, data.03 = data.03, h = 8, n.fac = 4, fre = "quarter",
                                         all.per.sta = as.Date("1981-01-01"), all.per.end = as.Date("2017-10-01"),
                                         est.per.sta = as.Date("1981-01-01"), est.per.end = as.Date("2002-10-01"),
                                         pre.per.sta = as.Date("2002-10-01"), pre.per.end = as.Date("2017-10-01"),
                                         mis.val = TRUE, print = TRUE)
FVAR.POOS.r4.h8.04$RMSFE.POOS
sink()

jsonlite::write_json(FVAR.POOS.r4.h8.04, "./ch-17-06/results/FVAR_POOS_r4_h8_04.json")

#..................................................
# Comparison of results ----

sink(file = "./ch-17-06/results/results_dfm_all.txt")

  print("..................................................")
  print("AR(2) model (compare with S&W, 2020, table 17.3, row 1)")
  
  print("AR2.POOS.h1$RMSFE.POOS")
  AR2.POOS.h1$RMSFE.POOS
  print("AR2.POOS.h4$RMSFE.POOS")
  AR2.POOS.h4$RMSFE.POOS
  print("AR2.POOS.h8$RMSFE.POOS")
  AR2.POOS.h8$RMSFE.POOS

  print("..................................................")
  print("Non-Recursive Estimation of the Factors based on variables without missing values")
  
  print("FVAR.POOS.r4.h1.01$RMSFE.POOS")
  FVAR.POOS.r4.h1.01$RMSFE.POOS
  print("FVAR.POOS.r4.h4.01$RMSFE.POOS")
  FVAR.POOS.r4.h4.01$RMSFE.POOS
  print("FVAR.POOS.r4.h8.01$RMSFE.POOS")
  FVAR.POOS.r4.h8.01$RMSFE.POOS

  print("..................................................")
  print("Non-Recursive Estimation of the Factors based on all variables (compare with S&W, 2020, table 17.3, row 3)")

  print("FVAR.POOS.r4.h1.02$RMSFE.POOS")
  FVAR.POOS.r4.h1.02$RMSFE.POOS
  print("FVAR.POOS.r4.h4.02$RMSFE.POOS")
  FVAR.POOS.r4.h4.02$RMSFE.POOS
  print("FVAR.POOS.r4.h8.02$RMSFE.POOS")
  FVAR.POOS.r4.h8.02$RMSFE.POOS

  print("..................................................")
  print("Recursive Estimation of the Factors based on variables without missing values")
  
  print("FVAR.POOS.r4.h1.03$RMSFE.POOS")
  FVAR.POOS.r4.h1.03$RMSFE.POOS
  print("FVAR.POOS.r4.h4.03$RMSFE.POOS")
  FVAR.POOS.r4.h4.03$RMSFE.POOS
  print("FVAR.POOS.r4.h8.03$RMSFE.POOS")
  FVAR.POOS.r4.h8.03$RMSFE.POOS
  
  print("..................................................")
  print("Recursive Estimation of the Factors based on all variables")
  
  print("FVAR.POOS.r4.h1.04$RMSFE.POOS")
  FVAR.POOS.r4.h1.04$RMSFE.POOS
  print("FVAR.POOS.r4.h4.04$RMSFE.POOS")
  FVAR.POOS.r4.h4.04$RMSFE.POOS
  print("FVAR.POOS.r4.h8.04$RMSFE.POOS")
  FVAR.POOS.r4.h8.04$RMSFE.POOS
  
sink()
