
#..................................................
# Set-up ----
rm(list = ls())

# libraries
library(dplyr) # for data manipulation
library(lubridate) # for dealing with dates

# functions
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

#..................................................
# Estimate unobserved factors based on all macro variables ----

factor_estimation_ls <- function(xdata, n.fac, nt.min, init = "svd") {
  
  # dimensions
  xdata <- as.matrix(xdata)
  col.nam <- colnames(xdata)
  nt <- dim(xdata)[1]
  ns <- dim(xdata)[2]

  # 2) standardization ----
  xmean <- as.matrix(colMeans(xdata, na.rm = TRUE)) # mean across rows for each column % mean (ignoring NaN)
  # mult <- sqrt((sum(!is.na(xdata))-1)/sum(!is.na(xdata))) # % num of non-NaN entries for each series
  mult <- sqrt((colSums(!is.na(xdata))-1)/colSums(!is.na(xdata))) # % num of non-NaN entries for each series
  xstd <- as.matrix(apply(xdata, 2, sd, na.rm = TRUE) * mult) # % std (ignoring NaN)
  xdata.std <- (xdata - (t(xmean) %x% matrix(1, nrow = nt))) / (t(xstd) %x% matrix(1, nrow = nt))
  
  # exclude columns with missing values
  xbal <- packr(xdata.std)
  
  # 3) pca ----
  if (init == "princomp") {
    
    pc <- princomp(xbal)
    fac <- pc$scores[,1:n.fac]
    
  } else if (init == "svd") {
    
    sv <- svd(xbal)
    scores <- xbal %*% sv$v
    fac <- scores[,1:n.fac]
    
  }
  
  lam <- NA*matrix(0, nrow = ns, ncol = n.fac)
  
  diff <- 100
  ssr <- 0
  tol <- 1/10000000000
  X <- xbal
  
  # 4) iterative ls
  while(diff > tol*(nt*ns)) {
    
    ssr.old <- ssr
    
    # 4.1) regression over ns ----
    
    # regress all N (T x 1) columns of X on the (T x r) factor scores (f) to get (N x r) factor loading (lambda)
    for (i in 1:ns) {
      
      tmp <- t(packr(t(cbind(xdata.std[,i, drop = FALSE], fac))))
      
      if (dim(tmp)[1] >= nt.min) {
        
        y <- tmp[,1, drop = FALSE]
        x <- tmp[,2:ncol(tmp), drop = FALSE]
        
        lam[i,] <- t(cbind(lm.fit(x = x, y = y)$coefficients))

      }
      
    }
    
    # 4.2) regression over nt ----
    
    # regress all T (N x 1) rows of X on the (N x r) factor loading (lambda) to get (N x r) factor scores (f)
    for (t in 1:nt) {
      
      tmp <- t(packr(t(cbind(t(xdata.std[t, , drop = FALSE]), lam))))
      
      y <- tmp[,1, drop = FALSE]
      x <- tmp[,2:ncol(tmp), drop = FALSE]
      
      fac[t,] <- t(cbind(lm.fit(x = x, y = y)$coefficients))

    }
    
    # 4.3) compute residuals ----
    e <- xdata.std - fac %*% t(lam)
    ssr <- sum(colSums(e^2, na.rm = TRUE))                          
    diff <- abs(ssr.old - ssr)
    
  }
  
  # 5) extract results ----
  colnames(fac) <- paste("FAC_", seq(1,ncol(fac)), sep = "")
  
  # lam = lam*t(t(xstd) %x% matrix(1, nrow = n.fac))
  lam <- lam
  colnames(lam) <- paste("FAC_", seq(1,ncol(fac)), sep = "")
  rownames(lam) <- col.nam
  
  ret.lis <- list(fac = fac,
                  lam = lam)
  
  return(ret.lis)
  
}

# determine the srart and end of period for factor estimation
dat.seq <- seq.Date(from = as.Date("1955-04-01"), to = as.Date("2017-10-01"), by = "quarter")
dat.seq[3]   # start used
dat.seq[236] # end used

dat.tmp <- window(macro.02.ts, start = c(1955, 4), end = c(2014, 1))
xdata <- matrix(as.vector(dat.tmp), nrow = nrow(dat.tmp), ncol = ncol(dat.tmp)) # transform to "matrix"/"array"

factor.est <- factor_estimation_ls(xdata = xdata, n.fac = 4, nt.min = 20, init = "svd")
factor.est.ts <- ts(factor.est$fac, frequency = 4, start = c(1959, 3), end = c(2017, 4))

#..................................................
# Merge time series ----

# use some date sequence as basis for ts
date <- ts(seq.Date(from = as.Date("1950-01-01"), to = as.Date("2025-01-01"), by = "quarter"), frequency = 4, start = c(1950, 1), end = c(2025,1))
TT <- length(date)

# merge all variables/factors
data.all.ts <- cbind(date, macro.01.ts, factor.est.ts, macro.02.ts) # based on estimated factors

dim(data.all.ts)

# extract time period of interest
data.all.ts <- window(data.all.ts, start = c(1959, 3), end = c(2017, 4))
# corresponding date
date <- seq.Date(from = as.Date("1959-07-01"), to = as.Date("2017-10-01"), by = "quarter")

# transform to data frame
data.all.df <- data.frame(date = date,
                          year = lubridate::year(date),
                          quarter = lubridate::quarter(date),
                          GDP     = as.numeric(data.all.ts[,2]),
                          F1      = as.numeric(data.all.ts[,5]),
                          F2      = as.numeric(data.all.ts[,6]),
                          F3      = as.numeric(data.all.ts[,7]),
                          F4      = as.numeric(data.all.ts[,8]))

# .......................................................
# POOS Analysis ----

FVAR_POOS_function <- function(h, print = FALSE) {
  
  data.h.all <-  data.all.df %>%
    mutate(GDPGR_h = (400/h) * (log(GDP / lag(GDP, h))),
           GDPGR_h_L1 = lag(log(GDP / lag(GDP, 1)), h),
           GDPGR_h_L2 = lag(log(GDP / lag(GDP, 1)), h + 1),
           F1_h_L1 = lag(F1, h),
           F2_h_L1 = lag(F2, h),
           F3_h_L1 = lag(F3, h),
           F4_h_L1 = lag(F4, h),
           CONST = 1) %>%
    dplyr::select(date, GDPGR_h, GDPGR_h_L1, GDPGR_h_L2, F1_h_L1, F2_h_L1, F3_h_L1, F4_h_L1, CONST)
  
  #..................................................
  # Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575) ----
  
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
  
  #..................................................
  # Prepare object to collect POOS-analysis-results ----
  
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  coef.lis <- list() 
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data for estimation
    data.h.tmp <- data.h.all
    est.dat.tmp <- data.h.tmp %>%
      filter(date >= all.per[est.sta] & date <= all.per[s.ii])
    # estimate model
    coef.tmp <- lm(GDPGR_h ~ GDPGR_h_L1 + GDPGR_h_L2 + F1_h_L1 + F2_h_L1 + F3_h_L1 + F4_h_L1 + CONST - 1,
                   data = est.dat.tmp)$coefficients
    coef.lis[[tt]] <- coef.tmp
    # extract data for prediction Y
    data.h.tmp <- data.h.all
    y.pre.dat.tmp <- data.h.tmp %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h)
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    data.h.tmp <- data.h.all
    X.pre.dat.tmp <- data.h.tmp %>%
      filter(date == all.per[pre.sta + tt]) %>%
      dplyr::select(GDPGR_h_L1, GDPGR_h_L2, F1_h_L1, F2_h_L1, F3_h_L1, F4_h_L1, CONST)
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
