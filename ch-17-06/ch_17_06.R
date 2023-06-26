
#..................................................
# 1) Set-up ----
rm(list = ls())

# setwd("./ch-17-06")

# 1.1) helper libraries ----
library(dplyr) # for data manipulation
library(lubridate) # for dealing with dates
# library(tidyverse) # for using tidyverse

library(lmtest) # for robust inference
library(sandwich)

library(dynlm) # for estimation
library(vars)

library(microbenchmark)

# for including matlab results
library(R.matlab)
library(matlab)

# for computing time
library(Matrix)
# see: https://cran.r-project.org/web/packages/Matrix/vignettes/Comparisons.pdf
library(RcppArmadillo)
# see:

# 1.2) helper functions ----

repmat <- function(X, nc) {
  
  # repeat matrix across columns  
  
  X %x% matrix(1, nrow = nc)
  
}

packr <- function(x) {
  
  # delete all columns with an NA
  
  ii <- apply(X = x, FUN = function(x) {ifelse(sum(is.na(x))==0,TRUE,FALSE)}, MARGIN = 2)
  x <- x[,c(ii)]
  
  x
  
}

inv <- function(x) {
  
  solve(x)
  
}

factor_estimation_ls <- function(xdata, nfac_t, nt_min) {
  
  # # 2.1) inputs ----
  # xdata <- xdata
  # nfac_t <- 4# number of total factors
  # nt_min <- 20 # minimum number of observations for any series used to estimate factors
  # 
  # xdata <- fac.est.dat.tmp[,-c(1,2,3,4)]
  # nfac_t <- 4
  # nt_min <- 20
  
  
  
  # dimensions
  nt <- dim(xdata)[1]
  ns <- dim(xdata)[2]
  
  # 2.2) standardization ----
  xmean <- as.matrix(colMeans(xdata, na.rm = TRUE)) # mean across rows for each column % mean (ignoring NaN)
  mult <- sqrt((sum(!is.na(xdata))-1)/sum(!is.na(xdata))) # % num of non-NaN entries for each series
  xstd <- as.matrix(apply(xdata, 2, sd, na.rm = TRUE) * mult) # % std (ignoring NaN)
  xdata_std <- (xdata - (t(xmean) %x% matrix(1, nrow = nt))) / (t(xstd) %x% matrix(1, nrow = nt))
  
  # exclude columns with missing values
  xbal <- packr(xdata_std)
  
  # 2.3) pca ----
  # pca.res <- princomp(xbal)
  # pca.res <- stats::prcomp(xbal)
  k <- min(nrow(xbal), ncol(xbal))
  s <- svd(xbal, nu = k, nv = 0)
  
  diff <- 100
  ssr <- 0
  tol <- 1/10000000000
  X <- xbal
  
  # f <- pca.res$scores[,1:4]
  # f <- pca.res$rotation[,1:4]
  f <- s$u[,1:4]
  fa <- f
  lambda <- NA*matrix(0, nrow = ns, ncol = nfac_t)
  
  # 2.4) iterative ls
  while(diff > tol*(nt*ns)) {
    
    ssr_old <- ssr
    
    # 2.4.1) regression over ns ----
    
    # regress all N (T x 1) columns of X on the (T x r) factor scores (f) to get (N x r) factor loadings (lambda)
    for (i in 1:ns) {
      
      tmp <- t(packr(t(cbind(xdata_std[,i, drop = FALSE], fa))))
      
      if (dim(tmp)[1] >= nt_min) {
        
        y <- tmp[,1, drop = FALSE]
        x <- tmp[,2:ncol(tmp), drop = FALSE]
        
        # lambda[i,] <- t(cbind(lm.fit(x = x, y = y)$coefficients))
        # lambda[i,] <- t(backsolve(chol(crossprod(x)), forwardsolve(chol(crossprod(x)), crossprod(x, y), upper = TRUE, trans = TRUE)))
        lambda[i,] <- t(fastLmPure(x, y)$coefficients)
        
      }
      
    }
    
    # 2.4.2) regression over ns ----
    
    # regress all T (N x 1) rows of X on the (N x r) factor loadings (lambda) to get (N x r) factor scores (f)
    for (t in 1:nt) {
      
      tmp <- t(packr(t(cbind(t(xdata_std[t, , drop = FALSE]), lambda))))
      
      y <- tmp[,1, drop = FALSE]
      x <- tmp[,2:ncol(tmp), drop = FALSE]
      
      # f[t,] <- t(cbind(lm.fit(x = x, y = y)$coefficients))
      # f[t,] <- t(backsolve(chol(crossprod(x)), forwardsolve(chol(crossprod(x)), crossprod(x, y), upper = TRUE, trans = TRUE)))
      f[t,] <- t(fastLmPure(x, y)$coefficients)
      
    }
    
    fa <- f
    
    # 2.4.3) compute residuals ----
    e <- xdata_std - fa %*% t(lambda)
    ssr <- sum(colSums(e^2, na.rm = TRUE))                          
    diff <- abs(ssr_old - ssr)
    
  }
  
  # 2.4) extract results
  f_est <- fa
  
  # fac_est <- matrix(NA, nrow = nrow(est_data), ncol = nfac_t)
  # fac_est[istart:iend,] <- f_est
  fac_est <- f_est
  colnames(fac_est) <- paste("F", seq(1,ncol(fac_est)), sep = "")
  
  check <- FALSE
  if (check) {
    head(fac_est)
    tail(fac_est)
    
    tmp <- readMat(con = "factor_all_4.mat")
    head(tmp$factor.all)
    tail(tmp$factor.all)
    # seems to work
  }
  
  lambda = lambda*t(t(xstd) %x% matrix(1, nrow = nfac_t))
  
  ret.lis <- list(fac_est = fac_est,
                  lambda = lambda)
  
  return(ret.lis)
  
  
}



#..................................................
# 2) load data ----

# observable macro variables
macro.sub.dat <- read.table("https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/data/us_macro_data.txt",
                            header = TRUE,
                            sep = ",",
                            colClasses = c("character", "numeric", "numeric", "numeric"))

macro.sub.ts <- ts(macro.sub.dat[,-c(1)], frequency = 4, start = c(1955, 2), end = c(2017, 4))

# raw data for unobservable factor estimation
tmp <- readMat(con = "https://raw.githubusercontent.com/mmoessler/stock-watson-2020-textbook/main/ch-17-06/datain.mat")

bpdata <- tmp$datain[2,1,1]$bpdata 
bpinclcode <- tmp$datain[10,1,1]$bpinclcode 

est_data <- bpdata[,bpinclcode==1]
n_series <- dim(est_data)[2]

macro.all.ts <- ts(est_data, frequency = 4, start = c(1959, 1), end = c(2017, 4))
colnames(macro.all.ts) <- paste("V_", seq(1,ncol(macro.all.ts)), sep = "")

head(macro.all.ts)[,1:10]

# sample period
istart <- 3 
iend <- 236 

xdata <- est_data[istart:iend,]



#..................................................
# 3) LS factor estimation ----

# microbenchmark(res <- factor_estimation_ls(xdata = xdata, nfac_t = 4, nt_min = 20))
# system.time(res <- factor_estimation_ls(xdata = xdata, nfac_t = 4, nt_min = 20))

# factor.dat <- res$fac_est
# 
# factor.ts <- ts(factor.dat, frequency = 4, start = c(1959, 3), end = c(2017, 4))



#..................................................
# 2) Merge time series ----

# use some date sequence as basis for ts
date.ts <- ts(seq.Date(from = as.Date("1950-01-01"), to = as.Date("2025-01-01"), by = "quarter"), frequency = 4, start = c(1950, 1), end = c(2025,1))
TT <- length(date.ts)

# merge data sets of interest
macro.sub.ts <- cbind(date.ts, macro.sub.ts) # merge gdp and term spread only
macro.all.ts <- cbind(date.ts, macro.all.ts) # merge date and all macro variables
# extract time period of interest
macro.sub.ts <- window(macro.sub.ts, start = c(1959, 3), end = c(2017, 4))
macro.all.ts <- window(macro.all.ts, start = c(1959, 3), end = c(2017, 4))

# transform to data frame
macro.sub.df <- data.frame(date = as.Date(macro.sub.ts[,1]),
                           year = lubridate::year(as.Date(macro.sub.ts[,1])),
                           quarter = lubridate::quarter(as.Date(macro.sub.ts[,1])),
                           GDP     = as.numeric(macro.sub.ts[,2]),
                           GDPGR   = as.numeric(macro.sub.ts[,3]),
                           TSpread = as.numeric(macro.sub.ts[,4]))

macro.all.df <- data.frame(date = as.Date(macro.all.ts[,1]),
                           year = lubridate::year(as.Date(macro.all.ts[,1])),
                           quarter = lubridate::quarter(as.Date(macro.all.ts[,1])))

macro.all.df <- cbind(macro.all.df, macro.all.ts)
dim(macro.all.df)
head(macro.all.df)[,1:10]



#..................................................
# 3) POOS FVAR Model ----

FVAR_POOS_function <- function(h, print = FALSE) {
  
  # h <- 2
  # print <- TRUE
  
  
  

    
  # Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575)
  
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
  
  # prepare (observed) variables of interest
  var.dat.all <- cbind(macro.sub.df) %>%
    mutate(GDPGR_h = (400/h) * log(lead(GDP, h) / GDP),
           GDPGR_L1 = log(GDP / lag(GDP, 1)),
           GDPGR_L2 = log(lag(GDP, 1) / lag(GDP, 2)),
           CONST = 1) %>%
    dplyr::select(date, GDP, GDPGR, GDPGR_h, GDPGR_L1, GDPGR_L2, CONST)
  # Note:
  # -> The last row represents all relevant information to estimate the h-lead of the target variable 
  # -> The h-lead of the target is the lhs variable in the direct forecast regression!

  # prepare object to collect POOS-analysis-results
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  coef.lis <- list() 
  fac.est.tmp.lis <- list() # list to store results of factor estimation
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data to estimate unobserved factors
    fac.est.dat.tmp <- macro.all.df %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    # estimate unobserved factors
    fac.est.tmp <- factor_estimation_ls(xdata = fac.est.dat.tmp[,-c(1,2,3,4)], nfac_t = 4, nt_min = 20)
    
    fac.est.tmp.lis[[tt]] <- fac.est.tmp
    
    # extract (observed) variables of interest
    var.dat.tmp <- var.dat.all %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    
    # merge & select observed variables of interest and estimated factors (depends on model!) ----
    est.dat.tmp <- cbind(var.dat.tmp, fac.est.tmp$fac_est) %>%
      dplyr::select(date, GDPGR, GDPGR_h, GDPGR_L1, GDPGR_L2, F1, F2, F3, F4, CONST)

    # estimate model
    coef.tmp <- lm(GDPGR_h ~ GDPGR_L1 + GDPGR_L2 + F1 + F2 + F3 + F4 + CONST - 1,
                   data = est.dat.tmp)$coefficients
    coef.lis[[tt]] <- coef.tmp
    
    # extract data for prediction Y
    y.pre.dat.tmp <- var.dat.all %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h)
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    X.pre.dat.tmp <- est.dat.tmp %>%
      filter(date == all.per[pre.sta + tt - h]) %>%
      dplyr::select(GDPGR_L1, GDPGR_L2, F1, F2, F3, F4, CONST)
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
  
  MSFE_POOS <- 1/length(u.til[-seq(1,h),]) * sum(u.til[-seq(1,h),]^2)
  RMSFE_POOS <- sqrt(MSFE_POOS)

  # returns
  ret.lis <- list(MSFE_POOS = MSFE_POOS,
                  RMSFE_POOS = RMSFE_POOS,
                  data = data,
                  y.act = y.act,
                  y.hat = y.hat,
                  u.til = u.til,
                  coef.lis = coef.lis)
  
  return(ret.lis)
  
  tmp <- paste("./ch-17-06/fac_est_res_h", h, ".RData", sep="")
  save(fac.est.tmp.lis, file = tmp)
  
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

