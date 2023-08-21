
#..................................................
# Set-up ----
rm(list = ls())

# to save the results into .json
library(jsonlite)

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
  
  # x <- data.01$GDP
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

# Function to estimate factors using iterative ls
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

# function for poos for observed variables and unobserved factor models
FVAR_POOS_function <- function(date.df, data.01, data.02, data.03, h, n.fac, fre,
                               all.per.sta, all.per.end,
                               est.per.sta, est.per.end,
                               pre.per.sta, pre.per.end,
                               mis.val = TRUE, print = FALSE) {
  
  # # X) Choose inputs here!
  # date.df <- date.df # data frame with a date column
  # data.01 <- data.01 # target variable (e.g., h-step GDPGR)
  # data.02 <- data.02 # observed predictors (variables) (data frame)
  # data.03 <- data.03 # observed predictors for (fair) construction of unobserved factors (matrix/array)
  # h <- 1
  # mis.val <- FALSE
  # print <- TRUE
  # 
  # n.fac <- 4
  # 
  # all.per.sta <- as.Date("1981-01-01")
  # all.per.end <- as.Date("2017-10-01")
  # est.per.sta <- as.Date("1981-01-01")
  # est.per.end <- as.Date("2002-10-01")
  # pre.per.sta <- as.Date("2002-10-01")
  # pre.per.end <- as.Date("2017-10-01")
  # fre <- "quarter"
  
  
  
  
  
  # 1) Prepare data ----
  
  if (is.null(data.03)) {
    
    # variable names
    data.01.nam <- colnames(data.01)
    data.02.nam <- colnames(data.02)
    
    # collect all data
    data <- cbind(data.01, data.02)
    tmp <- c(data.01.nam, data.02.nam)
    ii <- which(colnames(data) %in% tmp)
    data <- data[,ii]
    
  } else {
    
    # variable names
    data.01.nam <- colnames(data.01)
    data.02.nam <- colnames(data.02)
    data.03.nam <- colnames(data.03)
    
    # factor estimation
    dat.fra <- as.data.frame(data.03)
    col.nam.sel <- colnames(dat.fra)
    col.nam.out <- paste0("V", seq(1,ncol(dat.fra)),"_h_l1")
    arg.lis <- list(X = NULL, MARGIN = 2, FUN = "lag_fun", lag = 1)
    fac.data <- dat_man_fun(dat.fra, arg.lis, col.nam.sel, col.nam.out)
    
    # collect all data
    data <- cbind(data.01, data.02, fac.data)
    tmp <- c(data.01.nam, data.02.nam, paste0("V", seq(1,ncol(data.03)),"_h_l1"))
    ii <- which(colnames(data) %in% tmp)
    data <- data[,ii]
    
  }
  
  # 2) Prepare periods and POOS-analysis-index-s (see S&W, 2020, p. 575) ----
  
  # complete period
  all.per <- seq.Date(from = all.per.sta, to = all.per.end, by = fre)
  # estimation period
  est.per <- seq.Date(from = est.per.sta, to = est.per.end, by = fre)
  # prediction period
  pre.per <- seq.Date(from = pre.per.sta, to = pre.per.end, by = fre)
  
  # starting index for estimation/prediction
  est.sta <- which(all.per %in% est.per.sta)
  pre.sta <- which(all.per %in% pre.per.sta)
  
  # starting index for s
  s.ii.00 <- which(all.per %in% pre.per.sta)
  s.ii.00.h <- s.ii.00 - h
  
  # ending index of s
  s.ii.TT <- which(all.per %in% pre.per.end)
  s.ii.TT.h <- s.ii.TT - h
  
  s.ii.seq.h <- seq(s.ii.00.h, s.ii.TT.h)
  
  pre.sta <- pre.sta - 1
  
  # 3) POOS-analysis ----
  
  y.hat <- matrix(NA, nrow = length(s.ii.seq.h))
  y.act <- matrix(NA, nrow = length(s.ii.seq.h))
  u.til <- matrix(NA, nrow = length(s.ii.seq.h))
  
  for (tt in 1:length(s.ii.seq.h)) {
    
    # moving POOS-analysis-index-s
    s.ii <- s.ii.seq.h[tt]
    
    # extract data for estimation
    est.dat.tmp <- data[which(date.df$date %in% seq.Date(from = all.per[est.sta], to = all.per[s.ii], by = fre)),]
    
    if (is.null(data.03)) {
      
      # construct formula
      formula <- as.formula(paste(paste0(data.01.nam, " ~ "),
                                  paste(data.02.nam, collapse = " + "), " - 1"))
      
    } else {
      
      # factor analysis
      ii <- which(colnames(est.dat.tmp) %in% paste0("V", seq(1,ncol(data.03)),"_h_l1"))
      
      if (mis.val == TRUE) {
        
        # fa based on iterative ls
        factor.est <- factor_estimation_ls(xdata = est.dat.tmp[,ii], n.fac = n.fac, nt.min = 20, init = "svd")
        F.est <- factor.est$fac
        
      } else {
        
        # fa based on ordinary pca
        est.fac.dat.tmp <- scale(packr(est.dat.tmp[,ii]), center = TRUE, scale = FALSE)
        X <- as.matrix(est.fac.dat.tmp)
        
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
        
      }
      
      colnames(F.est) <- paste0("F", seq(1,n.fac),"_h_l1")
      est.dat.tmp <- cbind(est.dat.tmp, F.est)
      
      # construct formula
      formula <- as.formula(paste(paste0(data.01.nam, " ~ "),
                                  paste(data.02.nam, collapse = " + "),
                                  " + ", paste(paste0("F", seq(1,n.fac),"_h_l1"), collapse="+"), " - 1"))
      
    }
    
    # estimate model
    coef.tmp <- lm(formula = formula, data = est.dat.tmp)$coefficients
    
    # extract data for prediction Y
    y.pre.dat.tmp <- data[which(date.df$date == all.per[pre.sta + tt]), data.01.nam]
    y.act[tt] <- as.numeric(y.pre.dat.tmp)
    
    # extract data for prediction X
    X.pre.dat.tmp <- data[which(date.df$date == all.per[pre.sta + tt]), data.02.nam]
    if (!is.null(data.03)) {
      X.pre.dat.tmp <- cbind(matrix(X.pre.dat.tmp, nrow = 1), matrix(F.est[nrow(F.est),], nrow = 1))
    }
    
    # predict y
    y.hat[tt] <- as.numeric(matrix(coef.tmp, nrow = 1)) %*% as.numeric(matrix(X.pre.dat.tmp, ncol = 1))
    # evaluate prediction
    u.til[tt] <- y.act[tt] - y.hat[tt]
    
    if (print == TRUE) {
      
      # print for diagnostics
      cat("--------------------------------------------------", "\n")
      cat(paste0("Step ", tt, " from ", length(s.ii.seq.h), "\n"))
      cat(paste0("   Estimation: From ", as.Date(all.per[est.sta]), " to ", as.Date(all.per[s.ii]), "\n"))
      cat(paste0("   Prediction: For ", as.Date(all.per[pre.sta + tt]), "\n"))
      
      cat(paste0("  Y: ", y.pre.dat.tmp, "\n"))
      cat(paste0("  b: ", coef.tmp, "\n"))
      cat(paste0("  b (name): ", names(coef.tmp), "\n"))
      cat(paste0("  X: ", X.pre.dat.tmp, "\n"))
      cat(paste0("  Y (hat): ", y.hat[tt], "\n"))
      cat(paste0("  u (til): ", u.til[tt], "\n"))
      
    }
    
  }
  
  # 4) POOS-analysis-results ----
  
  MSFE.POOS <- 1/length(u.til[-seq(1,h),]) * sum(u.til[-seq(1,h),]^2)
  
  RMSFE.POOS <- sqrt(MSFE.POOS)
  
  # returns
  ret.lis <- list(MSFE.POOS = MSFE.POOS,
                  RMSFE.POOS = RMSFE.POOS,
                  data = data,
                  y.act = y.act,
                  y.hat = y.hat,
                  u.til = u.til)
  
  return(ret.lis)
  
}
