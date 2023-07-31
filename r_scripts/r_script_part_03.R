
#..................................................
# 0) Set-up ----
rm(list = ls())

library(dplyr) # for data manipulation
library(lubridate) # for dealing with dates
# library(tidyverse) # for using tidyverse

library(lmtest) # for robust inference
library(sandwich)

library(dynlm) # for estimation
library(vars)



#..................................................
# 1) Load data and estimated factors ----

# Observable macro variables
macro.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/us_macro_data.txt",
                        header = TRUE,
                        sep = ",",
                        colClasses = c("character", "numeric", "numeric", "numeric"))
# head(macro.dat)

macro.ts <- ts(macro.dat[,-c(1)], frequency = 4, start = c(1955, 2), end = c(2017, 4))
# dim(macro.ts)

# Unobservable/estimated factors
factor.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/factor_all_4.txt",
                        header = TRUE,
                        sep = ",",
                        colClasses = c("character", "numeric", "numeric", "numeric", "numeric"))
head(factor.dat)

factor.ts <- ts(factor.dat[,2:5], frequency = 4, start = c(1959, 3), end = c(2017, 4))
dim(factor.ts)


#..................................................
# 2) Merge time series ----

# use some date sequence as basis for ts
date <- ts(seq.Date(from = as.Date("1950-01-01"), to = as.Date("2025-01-01"), by = "quarter"), frequency = 4, start = c(1950, 1), end = c(2025,1))
TT <- length(date)

# merge
data.all.ts <- cbind(date, macro.ts, factor.ts)
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
                          TSpread = as.numeric(data.all.ts[,4]),
                          F01     = as.numeric(data.all.ts[,5]),
                          F02     = as.numeric(data.all.ts[,6]),
                          F03     = as.numeric(data.all.ts[,7]),
                          F04     = as.numeric(data.all.ts[,8]))
# head(data.all.df)

data.all.df <- data.all.df %>%
  mutate(GDPGR_X = log(GDP / lag(GDP, 1)))
# head(data.all.df)



#..................................................
# 3) POOS FVAR Model ----

FVAR_POOS_function <- function(h, print = FALSE) {
  
  # 1) Choose h here! ----
  # h <- 2
  
  
  
  # 2) Prepare data (depends on model!) ----
  data <- data.all.df %>%
    mutate(GDPGR_h = (400/h) * (log(GDP / lag(GDP, h))),
           # GDPGR_h_L1 = lag(GDPGR, h),
           # GDPGR_h_L2 = lag(GDPGR, h + 1),
           GDPGR_h_L1 = lag(GDPGR_X, h),
           GDPGR_h_L2 = lag(GDPGR_X, h + 1),
           F01_h_L1 = lag(F01, h),
           F02_h_L1 = lag(F02, h),
           F03_h_L1 = lag(F03, h),
           F04_h_L1 = lag(F04, h),
           CONST = 1) %>%
    dplyr::select(date, GDPGR, GDPGR_h, GDPGR_h_L1, GDPGR_h_L2, F01_h_L1, F02_h_L1, F03_h_L1, F04_h_L1, CONST)
  # head(data, 10)

  
  
  # 3) Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575) ----
  
  # complete period
  all.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2017-10-01"), by = "quarter")
  # estimation period
  est.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2002-10-01"), by = "quarter")
  # prediction period
  pre.per <- seq.Date(from = as.Date("2002-10-01"), to = as.Date("2017-10-01"), by = "quarter")
  
  # starting index for estimation/prediction
  est.sta <- which(all.per %in% as.Date("1981-01-01"))
  pre.sta <- which(all.per %in% as.Date("2002-10-01"))
  
  # starting index for s
  s.ii.00 <- which(all.per %in% as.Date("2002-10-01"))
  s.ii.00.h <- s.ii.00 - h
  
  # ending index of s
  s.ii.TT <- which(all.per %in% as.Date("2017-10-01"))
  s.ii.TT.h <- s.ii.TT - h
  
  s.ii.seq.h <- seq(s.ii.00.h, s.ii.TT.h)
  
  pre.sta <- pre.sta - 1
  
  
  
  # 4) Prepare object to collect POOS-analysis-results ----
  
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data for estimation
    est.dat.tmp <- data %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    # estimate model
    coef.tmp <- lm(GDPGR_h ~ GDPGR_h_L1 + GDPGR_h_L2 + F01_h_L1 + F02_h_L1 + F03_h_L1 + F04_h_L1 + CONST - 1,
                   data = est.dat.tmp)$coefficients
    # extract data for prediction Y
    y.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h)
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    X.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h_L1, GDPGR_h_L2, F01_h_L1, F02_h_L1, F03_h_L1, F04_h_L1, CONST)
    # predict y
    y.hat[tt] <- as.numeric(matrix(coef.tmp, nrow = 1)) %*% as.numeric(matrix(X.pre.dat.tmp, ncol = 1))
    # evaluate prediction
    u.til[tt] <- y.act[tt] - y.hat[tt]
    
    if (print == TRUE) {
      
      # print for diagnostics
      cat(paste0("Step ", tt, " from ", length(s.ii.seq.h), "\n"))
      cat(paste0("   Estimation: From ", as.Date(all.per[est.sta]), " to ", as.Date(all.per[s.ii]), "\n"))
      cat(paste0("   Prediction: For ", as.Date(all.per[pre.sta + tt]), "\n"))
      
      cat(paste0("  Y: ", y.pre.dat.tmp, "\n"))
      cat(paste0("  Y (name): ", colnames(y.pre.dat.tmp), "\n"))
      cat(paste0("  b: ", coef.tmp, "\n"))
      cat(paste0("  b (name): ", names(coef.tmp), "\n"))
      cat(paste0("  X: ", X.pre.dat.tmp, "\n"))
      cat(paste0("  X (name): ", colnames(X.pre.dat.tmp), "\n"))
      
    }
    
  }
  
  # MSFE_POOS <- 1/length(u.til) * sum(u.til^2)
  MSFE_POOS <- 1/length(u.til[-seq(1,h),]) * sum(u.til[-seq(1,h),]^2)
  # MSFE_POOS
  
  RMSFE_POOS <- sqrt(MSFE_POOS)
  # RMSFE_POOS
  
  # returns
  ret.lis <- list(MSFE_POOS = MSFE_POOS,
                  RMSFE_POOS = RMSFE_POOS,
                  data = data,
                  y.act = y.act,
                  y.hat = y.hat,
                  u.til = u.til)
  
  return(ret.lis)
  
}

# h=1
FVAR_POOS_h1 <- FVAR_POOS_function(h = 1, print = FALSE)
FVAR_POOS_h1$RMSFE_POOS

# h=4
FVAR_POOS_h4 <- FVAR_POOS_function(h = 4, print = FALSE)
FVAR_POOS_h4$RMSFE_POOS

# h=8
FVAR_POOS_h8 <- FVAR_POOS_function(h = 8, print = FALSE)
FVAR_POOS_h8$RMSFE_POOS



#..................................................
# 4) POOS AR Model ----

AR_POOS_function <- function(h, print = FALSE) {
  
  # 1) Choose h here! ----
  # h <- 2
  
  
  
  # 2) Prepare data (depends on model!) ----
  data <- data.all.df %>%
    mutate(GDPGR_h = (400/h) * (log(GDP / lag(GDP, h))),
           # GDPGR_h_L1 = lag(GDPGR, h),
           # GDPGR_h_L2 = lag(GDPGR, h + 1),
           GDPGR_h_L1 = lag(GDPGR_X, h),
           GDPGR_h_L2 = lag(GDPGR_X, h + 1),
           CONST = 1) %>%
    dplyr::select(date, GDPGR, GDPGR_h, GDPGR_h_L1, GDPGR_h_L2, CONST)
  # head(data, 10)
  
  
  
  # 3) Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575) ----
  
  # complete period
  all.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2017-10-01"), by = "quarter")
  # estimation period
  est.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2002-10-01"), by = "quarter")
  # prediction period
  pre.per <- seq.Date(from = as.Date("2002-10-01"), to = as.Date("2017-10-01"), by = "quarter")
  
  # starting index for estimation/prediction
  est.sta <- which(all.per %in% as.Date("1981-01-01"))
  pre.sta <- which(all.per %in% as.Date("2002-10-01"))
  
  # starting index for s
  s.ii.00 <- which(all.per %in% as.Date("2002-10-01"))
  s.ii.00.h <- s.ii.00 - h
  
  # ending index of s
  s.ii.TT <- which(all.per %in% as.Date("2017-10-01"))
  s.ii.TT.h <- s.ii.TT - h
  
  s.ii.seq.h <- seq(s.ii.00.h, s.ii.TT.h)
  
  pre.sta <- pre.sta - 1
  
  
  
  # 4) Prepare object to collect POOS-analysis-results ----
  
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data for estimation
    est.dat.tmp <- data %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    # estimate model
    coef.tmp <- lm(GDPGR_h ~ GDPGR_h_L1 + GDPGR_h_L2 + CONST - 1,
                   data = est.dat.tmp)$coefficients
    # extract data for prediction Y
    y.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h)
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    X.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h_L1, GDPGR_h_L2, CONST)
    # predict y
    y.hat[tt] <- as.numeric(matrix(coef.tmp, nrow = 1)) %*% as.numeric(matrix(X.pre.dat.tmp, ncol = 1))
    # evaluate prediction
    u.til[tt] <- y.act[tt] - y.hat[tt]
    
    if (print == TRUE) {
      
      # print for diagonstics
      cat(paste0("Step ", tt, " from ", length(s.ii.seq.h), "\n"))
      cat(paste0("   Estimation: From ", as.Date(all.per[est.sta]), " to ", as.Date(all.per[s.ii]), "\n"))
      cat(paste0("   Prediction: For ", as.Date(all.per[pre.sta + tt]), "\n"))
      
      cat(paste0("  Y: ", y.pre.dat.tmp, "\n"))
      cat(paste0("  Y (name): ", colnames(y.pre.dat.tmp), "\n"))
      cat(paste0("  b: ", coef.tmp, "\n"))
      cat(paste0("  b (name): ", names(coef.tmp), "\n"))
      cat(paste0("  X: ", X.pre.dat.tmp, "\n"))
      cat(paste0("  X (name): ", colnames(X.pre.dat.tmp), "\n"))
      
    }
    
  }
  
  # MSFE_POOS <- 1/length(u.til) * sum(u.til^2)
  MSFE_POOS <- 1/length(u.til[-seq(1,h),]) * sum(u.til[-seq(1,h),]^2)
  # MSFE_POOS
  
  RMSFE_POOS <- sqrt(MSFE_POOS)
  # RMSFE_POOS
  
  # returns
  ret.lis <- list(MSFE_POOS = MSFE_POOS,
                  RMSFE_POOS = RMSFE_POOS,
                  data = data,
                  y.act = y.act,
                  y.hat = y.hat,
                  u.til = u.til)
  
  return(ret.lis)
  
}

# h=1
AR_POOS_h1 <- AR_POOS_function(h = 1, print = FALSE)
AR_POOS_h1$RMSFE_POOS

# h=4
AR_POOS_h4 <- AR_POOS_function(h = 4, print = FALSE)
AR_POOS_h4$RMSFE_POOS

# h=8
AR_POOS_h8 <- AR_POOS_function(h = 8, print = FALSE)
AR_POOS_h8$RMSFE_POOS



#..................................................
# 5) POOS ADL Model ----

ADL_POOS_function <- function(h, print = TRUE) {
  
  # 1) Choose h here! ----
  # h <- 2
  
  
  
  # 2) Prepare data (depends on model!) ----
  data <- data.all.df %>%
    mutate(GDPGR_h = (400/h) * (log(GDP / lag(GDP, h))),
           # GDPGR_h_L1 = lag(GDPGR, h),
           # GDPGR_h_L2 = lag(GDPGR, h + 1),
           GDPGR_h_L1 = lag(GDPGR_X, h),
           GDPGR_h_L2 = lag(GDPGR_X, h + 1),
           TSpread_h_L1 = lag(TSpread, h),
           TSpread_h_L2 = lag(TSpread, h + 1),
           CONST = 1) %>%
    dplyr::select(date, GDPGR, GDPGR_h, GDPGR_h_L1, GDPGR_h_L2, TSpread_h_L1, TSpread_h_L2, CONST)
  # head(data, 10)
  
  
  
  # 3) Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575) ----
  
  # complete period
  all.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2017-10-01"), by = "quarter")
  # estimation period
  est.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2002-10-01"), by = "quarter")
  # prediction period
  pre.per <- seq.Date(from = as.Date("2002-10-01"), to = as.Date("2017-10-01"), by = "quarter")
  
  # starting index for estimation/prediction
  est.sta <- which(all.per %in% as.Date("1981-01-01"))
  pre.sta <- which(all.per %in% as.Date("2002-10-01"))
  
  # starting index for s
  s.ii.00 <- which(all.per %in% as.Date("2002-10-01"))
  s.ii.00.h <- s.ii.00 - h
  
  # ending index of s
  s.ii.TT <- which(all.per %in% as.Date("2017-10-01"))
  s.ii.TT.h <- s.ii.TT - h
  
  s.ii.seq.h <- seq(s.ii.00.h, s.ii.TT.h)
  
  pre.sta <- pre.sta - 1
  
  
  
  # 4) Prepare object to collect POOS-analysis-results ----
  
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data for estimation
    est.dat.tmp <- data %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    # estimate model
    coef.tmp <- lm(GDPGR_h ~ GDPGR_h_L1 + GDPGR_h_L2 + TSpread_h_L1 + TSpread_h_L2 + CONST - 1,
                   data = est.dat.tmp)$coefficients
    # extract data for prediction Y
    y.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h)
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    X.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h_L1, GDPGR_h_L2, TSpread_h_L1, TSpread_h_L2, CONST)
    # predict y
    y.hat[tt] <- as.numeric(matrix(coef.tmp, nrow = 1)) %*% as.numeric(matrix(X.pre.dat.tmp, ncol = 1))
    # evaluate prediction
    u.til[tt] <- y.act[tt] - y.hat[tt]
    
    if (print == TRUE) {
      
      # print for diagonstics
      cat(paste0("Step ", tt, " from ", length(s.ii.seq.h), "\n"))
      cat(paste0("   Estimation: From ", as.Date(all.per[est.sta]), " to ", as.Date(all.per[s.ii]), "\n"))
      cat(paste0("   Prediction: For ", as.Date(all.per[pre.sta + tt]), "\n"))
      
      cat(paste0("  Y: ", y.pre.dat.tmp, "\n"))
      cat(paste0("  Y (name): ", colnames(y.pre.dat.tmp), "\n"))
      cat(paste0("  b: ", coef.tmp, "\n"))
      cat(paste0("  b (name): ", names(coef.tmp), "\n"))
      cat(paste0("  X: ", X.pre.dat.tmp, "\n"))
      cat(paste0("  X (name): ", colnames(X.pre.dat.tmp), "\n"))
      
    }
    
  }
  
  # MSFE_POOS <- 1/length(u.til) * sum(u.til^2)
  MSFE_POOS <- 1/length(u.til[-seq(1,h),]) * sum(u.til[-seq(1,h),]^2)
  # MSFE_POOS

  RMSFE_POOS <- sqrt(MSFE_POOS)
  # RMSFE_POOS
  
  # returns
  ret.lis <- list(MSFE_POOS = MSFE_POOS,
                  RMSFE_POOS = RMSFE_POOS,
                  data = data,
                  y.act = y.act,
                  y.hat = y.hat,
                  u.til = u.til)
  
  return(ret.lis)
  
}

# h=1
ADL_POOS_h1 <- ADL_POOS_function(h = 1, print = FALSE)
ADL_POOS_h1$RMSFE_POOS

# h=4
ADL_POOS_h4 <- ADL_POOS_function(h = 4, print = FALSE)
ADL_POOS_h4$RMSFE_POOS

# h=8
ADL_POOS_h8 <- ADL_POOS_function(h = 8, print = FALSE)
ADL_POOS_h8$RMSFE_POOS



#..................................................
# 6) POOS MEAN Model ----

MEAN_POOS_function <- function(h, print = FALSE) {
  
  # 1) Choose h here! ----
  # h <- 1

  
  
  # 2) Prepare data (depends on model!) ----
  data <- data.all.df %>%
    mutate(GDPGR_h = (400/h) * (log(GDP / lag(GDP, h))),
           # GDPGR_h_L1 = lag(GDPGR, h),
           # GDPGR_h_L2 = lag(GDPGR, h + 1),
           GDPGR_h_L1 = lag(GDPGR_X, h),
           GDPGR_h_L2 = lag(GDPGR_X, h + 1),
           CONST = 1) %>%
    dplyr::select(date, GDPGR, GDPGR_h, GDPGR_h_L1, GDPGR_h_L2, CONST)
  # head(data, 10)
  
  
  
  # 3) Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575) ----
  
  # complete period
  all.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2017-10-01"), by = "quarter")
  # estimation period
  est.per <- seq.Date(from = as.Date("1981-01-01"), to = as.Date("2002-10-01"), by = "quarter")
  # prediction period
  pre.per <- seq.Date(from = as.Date("2002-10-01"), to = as.Date("2017-10-01"), by = "quarter")
  
  # starting index for estimation/prediction
  est.sta <- which(all.per %in% as.Date("1981-01-01"))
  pre.sta <- which(all.per %in% as.Date("2002-10-01"))
  
  # starting index for s
  s.ii.00 <- which(all.per %in% as.Date("2002-10-01"))
  s.ii.00.h <- s.ii.00 - h
  
  # ending index of s
  s.ii.TT <- which(all.per %in% as.Date("2017-10-01"))
  s.ii.TT.h <- s.ii.TT - h
  
  s.ii.seq.h <- seq(s.ii.00.h, s.ii.TT.h)
  
  pre.sta <- pre.sta - 1
  
  
  
  # 4) Prepare object to collect POOS-analysis-results ----
  
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data for estimation
    est.dat.tmp <- data %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    # estimate model
    coef.tmp <- lm(GDPGR_h ~ CONST - 1,
                   data = est.dat.tmp)$coefficients
    # extract data for prediction Y
    y.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h)
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    X.pre.dat.tmp <- data %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(CONST)
    # predict y
    y.hat[tt] <- as.numeric(matrix(coef.tmp, nrow = 1)) %*% as.numeric(matrix(X.pre.dat.tmp, ncol = 1))
    # evaluate prediction
    u.til[tt] <- y.act[tt] - y.hat[tt]
    
    if (print == TRUE) {
      
      # print for diagonstics
      cat(paste0("Step ", tt, " from ", length(s.ii.seq.h), "\n"))
      cat(paste0("   Estimation: From ", as.Date(all.per[est.sta]), " to ", as.Date(all.per[s.ii]), "\n"))
      cat(paste0("   Prediction: For ", as.Date(all.per[pre.sta + tt]), "\n"))
      
      cat(paste0("  Y: ", y.pre.dat.tmp, "\n"))
      cat(paste0("  Y (name): ", colnames(y.pre.dat.tmp), "\n"))
      cat(paste0("  b: ", coef.tmp, "\n"))
      cat(paste0("  b (name): ", names(coef.tmp), "\n"))
      cat(paste0("  X: ", X.pre.dat.tmp, "\n"))
      cat(paste0("  X (name): ", colnames(X.pre.dat.tmp), "\n"))
      
    }
    
  }
  
  # MSFE_POOS <- 1/length(u.til) * sum(u.til^2)
  MSFE_POOS <- 1/length(u.til[-seq(1,h),]) * sum(u.til[-seq(1,h),]^2)
  # MSFE_POOS
  
  RMSFE_POOS <- sqrt(MSFE_POOS)
  # RMSFE_POOS
  
  # returns
  ret.lis <- list(MSFE_POOS = MSFE_POOS,
                  RMSFE_POOS = RMSFE_POOS,
                  data = data,
                  y.act = y.act,
                  y.hat = y.hat,
                  u.til = u.til)
  
  return(ret.lis)
  
}

# h=1
MEAN_POOS_h1 <- MEAN_POOS_function(h = 1, print = FALSE)
MEAN_POOS_h1$RMSFE_POOS

# h=4
MEAN_POOS_h4 <- MEAN_POOS_function(h = 4, print = FALSE)
MEAN_POOS_h4$RMSFE_POOS

# h=8
MEAN_POOS_h8 <- MEAN_POOS_function(h = 8, print = FALSE)
MEAN_POOS_h8$RMSFE_POOS



#..................................................
# 7) Combine results ----

res.df <- data.frame(h01 = c(AR_POOS_h1$RMSFE_POOS, ADL_POOS_h1$RMSFE_POOS, FVAR_POOS_h1$RMSFE_POOS, MEAN_POOS_h1$RMSFE_POOS),
                     h04 = c(AR_POOS_h4$RMSFE_POOS, ADL_POOS_h4$RMSFE_POOS, FVAR_POOS_h4$RMSFE_POOS, MEAN_POOS_h4$RMSFE_POOS),
                     h08 = c(AR_POOS_h8$RMSFE_POOS, ADL_POOS_h8$RMSFE_POOS, FVAR_POOS_h8$RMSFE_POOS, MEAN_POOS_h8$RMSFE_POOS))

rownames(res.df) <- c("AR", "ADL", "FVAR", "MEAN")

print("--------------------------------------------------")
print("Compsrison of Direct Forecasts of Cumulative GDP Growth at an Annual Rate (see: Table 17.3, S&W, 2020, p. 680)")
round(res.df, 2)
