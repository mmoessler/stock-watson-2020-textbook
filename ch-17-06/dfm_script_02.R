
#..................................................
# Set-up ----
rm(list = ls())

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

# R helper function data preparation 1
lag_fun <- function(x, lag) {
  
  # x <- seq(1,10)
  # lag <- 1
  
  
  if (lag >= 1) {
    c(rep(NA, lag), x[-seq(length(x), length(x)-lag+1)])
  } else {
    c(x[-seq(1, -lag)], rep(NA, -lag))
  }
  
}

# R helper function data preparation 2
dat_man_fun <- function(dat.fra, arg.lis, col.nam.sel, col.nam.out) {
  
  # Note:
  # -> Use margin 1 to apply transformation row-by-row
  # -> Use arg.lis[[1]]=NULL as dummy for dat.fra
  
  # # Checks:
  # dat.fra <- data.all.df
  # arg.lis <- list(X=NULL, MARGIN=2, FUN="L0h.fun", h=1)
  # col.nam.sel <- paste0("V", seq(1,97))
  # col.nam.out <- paste0("V", seq(1,97),"_h_L1")
  
  
  
  # add dat.fra as first element to list
  arg.lis[[1]] <- dat.fra
  # select columns of dat.fra
  dat.fra <- arg.lis[[1]]
  arg.lis[[1]] <- subset(arg.lis[[1]],select=col.nam.sel)

  # do.call(apply,args=arg.lis) # core
  res <- cbind(dat.fra,
               setNames(
                 as.data.frame(
                   do.call(
                     apply,args=arg.lis)),col.nam.out))
  
  return(res)
  
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
  mult <- sqrt((colSums(!is.na(xdata))-1)/colSums(!is.na(xdata))) # num of non-NaN entries for each series
  xstd <- as.matrix(apply(xdata, 2, sd, na.rm = TRUE) * mult) # std (ignoring NaN)
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
    
    # regress each of the N columns of X (T x 1) on the (T x r) factor scores (f) to get (N x r) factor loading (lambda)
    for (i in 1:ns) {
      
      tmp <- t(packr(t(cbind(xdata.std[,i, drop = FALSE], fac)))) # note delete rows/t's with NAs
      
      if (dim(tmp)[1] >= nt.min) {
        
        y <- tmp[,1, drop = FALSE]
        x <- tmp[,2:ncol(tmp), drop = FALSE]
        
        lam[i,] <- t(cbind(lm.fit(x = x, y = y)$coefficients))
        
      }
      
    }
    
    # 4.2) regression over nt ----
    
    # regress each of the T rows of X (N x 1)  on the (N x r) factor loadings (lambda) to get (N x r) factor scores (f)
    for (t in 1:nt) {
      
      tmp <- t(packr(t(cbind(t(xdata.std[t,, drop = FALSE]), lam)))) # note delete columns/i's with NAs
      
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

# determine the start and end of period for factor estimation
dat.seq <- seq.Date(from = as.Date("1955-04-01"), to = as.Date("2017-10-01"), by = "quarter")

print("--------------------------------------------------")
print("Start and end for estimation of unobserved factors")
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

# extract time period of interest
data.all.ts <- window(data.all.ts, start = c(1959, 3), end = c(2017, 4))
# corresponding date
date <- seq.Date(from = as.Date("1959-07-01"), to = as.Date("2017-10-01"), by = "quarter")

#..................................................
# FVAR Model (recursively estimated factor) ----
# (for a fair comparison the factors have to be estimated recursively)

# target variable and dates (data frame)
data.01 <- data.frame(date = date,
                      year = format(date, "%Y"),
                      quarter = quarters(date),
                      GDP     = as.numeric(data.all.ts[,2]))

# predictors for factor analysis (matrix/array)
data.02 <- xdata

# function for fair poos for r=4
FVAR_FAIR_POOS_r04_function <- function(data.01, data.02, h, mis.val = TRUE, print = FALSE) {
  
  # # X) Choose h here! ----
  # data.01 <- data.01
  # data.02 <- data.02
  # h <- 1
  # mis.val <- FALSE
  # print <- FALSE
  

  
  # 1) Prepare data no 1
  
  if (mis.val == TRUE) {
    est.data <- data.02
  } else if (mis.val == FALSE) {
    x <- data.02
    ii <- apply(X = x, FUN = function(x) {ifelse(sum(is.na(x))==0,TRUE,FALSE)}, MARGIN = 2)
    est.data <- x[,c(ii)]
  }
  
  NN <- ncol(est.data)
  TT <- nrow(est.data)
  
  F_ALL <- as.data.frame(est.data)
  dat.fra <- F_ALL
  col.nam.sel <- paste0("V", seq(1,NN))
  col.nam.out <- paste0("V", seq(1,NN),"_h_l1")
  arg.lis <- list(X = NULL, MARGIN = 2, FUN = "lag_fun", lag = 1)
  F_ALL <- dat_man_fun(dat.fra, arg.lis, col.nam.sel, col.nam.out)
  
  
  
  # 2) Prepare data no 2
  
  data <- cbind(data.01, F_ALL,
                data.frame(GDPGR_h = (400/h) * (log(data.01$GDP / lag_fun(data.01$GDP, h))),
                     GDPGR_h_l1 = lag_fun(log(data.01$GDP / lag_fun(data.01$GDP, 1)), h),
                     GDPGR_h_l2 = lag_fun(log(data.01$GDP / lag_fun(data.01$GDP, 1)), h + 1),
                     CONST = 1))
  
  
  tmp <- c("date", "GDPGR", "GDPGR_h", "GDPGR_h_l1", "GDPGR_h_l2", "CONST", paste0("V", seq(1,NN),"_h_l1"))
  ii <- which(colnames(data) %in% tmp)
  data <- data[,ii]
  
  
  
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
    est.dat.tmp <- data[which(data$date %in% seq.Date(from = all.per[est.sta], to = all.per[s.ii], by = "quarter")),]

    # factor analysis
    ii <- which(colnames(est.dat.tmp) %in% paste0("V", seq(1,NN),"_h_l1"))
    
    if (mis.val == TRUE) {
      
      # based on iterative ls
      factor.est <- factor_estimation_ls(xdata = est.dat.tmp[,ii], n.fac = 4, nt.min = 20, init = "svd")
      F.res <- factor.est$fac
      
    } else {
      
      # based on ordinary pca
      est.fac.dat.tmp <- scale(est.dat.tmp[,ii], center = TRUE, scale = FALSE)
      X <- as.matrix(est.fac.dat.tmp)
      
      # # 1) based on princomp (i.e., eigenvalue decomposition)
      # XX <- X %*% t(X)
      # Sig <- dim(X)[2]^(-1) * XX
      # pca.res <- princomp(Sig)
      # F.res <- pca.res$scores[,1:4]
      
      # # 2) based on prcomp (i.e., svd)
      # pca.res <- prcomp(X, retx = TRUE, rank. = 4)
      # F.res <- pca.res$x[,1:4]
      
      # 3) based on svd directly
      sv <- svd(X)
      scores <- X %*% sv$v
      F.res <- scores[,1:4]
      
    }
    
    colnames(F.res) <- c("F1_h_l1","F2_h_l1","F3_h_l1","F4_h_l1")
    est.dat.tmp <- cbind(est.dat.tmp, F.res)
    
    # estimate model
    coef.tmp <- lm(GDPGR_h ~ GDPGR_h_l1 + GDPGR_h_l2 + F1_h_l1 + F2_h_l1 + F3_h_l1 + F4_h_l1 + CONST - 1,
                   data = est.dat.tmp)$coefficients
    # extract data for prediction Y
    y.pre.dat.tmp <- data[which(data$date == all.per[pre.sta + tt]), "GDPGR_h"]
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    # extract data for prediction X
    X.pre.dat.tmp <- data[which(data$date == all.per[pre.sta + tt]), c("GDPGR_h_l1", "GDPGR_h_l2", "CONST")]
    X.pre.dat.tmp <- cbind(matrix(X.pre.dat.tmp[-3], nrow = 1), matrix(F.res[nrow(F.res),], nrow = 1), matrix(1, nrow = 1))
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



# Analysis based on all observations

# h=1
FVAR_FAIR_POOS_r4_h1_01 <- FVAR_FAIR_POOS_r04_function(data.01 = data.01, data.02 = data.02,
                                                       h = 1, mis.val = TRUE, print = FALSE)
FVAR_FAIR_POOS_r4_h1_01$RMSFE_POOS

# h=4
FVAR_FAIR_POOS_r4_h4_01 <- FVAR_FAIR_POOS_r04_function(data.01 = data.01, data.02 = data.02,
                                                       h = 4, mis.val = TRUE, print = FALSE)
FVAR_FAIR_POOS_r4_h4_01$RMSFE_POOS

# h=8
FVAR_FAIR_POOS_r4_h8_01 <- FVAR_FAIR_POOS_r04_function(data.01 = data.01, data.02 = data.02,
                                                       h = 8, mis.val = TRUE, print = FALSE)
FVAR_FAIR_POOS_r4_h8_01$RMSFE_POOS

# Analysis based on variables without missing values

# h=1
FVAR_FAIR_POOS_r4_h1_02 <- FVAR_FAIR_POOS_r04_function(data.01 = data.01, data.02 = data.02,
                                                       h = 1, mis.val = FALSE, print = FALSE)
FVAR_FAIR_POOS_r4_h1_02$RMSFE_POOS

# h=4
FVAR_FAIR_POOS_r4_h4_02 <- FVAR_FAIR_POOS_r04_function(data.01 = data.01, data.02 = data.02,
                                                       h = 4, mis.val = FALSE, print = FALSE)
FVAR_FAIR_POOS_r4_h4_02$RMSFE_POOS

# h=8
FVAR_FAIR_POOS_r4_h8_02 <- FVAR_FAIR_POOS_r04_function(data.01 = data.01, data.02 = data.02,
                                                       h = 8, mis.val = FALSE, print = FALSE)
FVAR_FAIR_POOS_r4_h8_02$RMSFE_POOS

sink(file = "./ch-17-06/results_dfm_02.txt")

  print("FVAR_FAIR_POOS_r4_h1_01$RMSFE_POOS")
  FVAR_FAIR_POOS_r4_h1_01$RMSFE_POOS
  print("FVAR_FAIR_POOS_r4_h4_01$RMSFE_POOS")
  FVAR_FAIR_POOS_r4_h4_01$RMSFE_POOS
  print("FVAR_FAIR_POOS_r4_h8_01$RMSFE_POOS")
  FVAR_FAIR_POOS_r4_h8_01$RMSFE_POOS
  
  print("FVAR_FAIR_POOS_r4_h1_02$RMSFE_POOS")
  FVAR_FAIR_POOS_r4_h1_02$RMSFE_POOS
  print("FVAR_FAIR_POOS_r4_h4_02$RMSFE_POOS")
  FVAR_FAIR_POOS_r4_h4_02$RMSFE_POOS
  print("FVAR_FAIR_POOS_r4_h8_02$RMSFE_POOS")
  FVAR_FAIR_POOS_r4_h8_02$RMSFE_POOS
  
sink()
