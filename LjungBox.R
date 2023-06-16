"LjungBox" <-
function(obj,lags=seq(5,30,5),order=0,season=1,squared.residuals=FALSE){
     class.obj = class(obj)[1]
     TestType <- "0"
    if (class.obj == "ts" || class.obj == "numeric" || class.obj == 
        "matrix" || class.obj == "mts") 
        TestType <- "1"
    if (class.obj == "ar" || class.obj == "arima0" || class.obj == 
        "Arima" || class.obj == "ARIMA" || class.obj == "varest" || class.obj == "lm"
        || class.obj == "glm" || class.obj == "list") 
        TestType<-"2"
    if (TestType == "0") 
        stop("obj must be class ar, arima0, Arima, (ARIMA forecast_ARIMA Arima), varest, lm, (glm lm), ts, numeric, matrix, (mts ts), or list")
     Maxlag <- max(lags)
     if (TestType=="1")
       res <- as.ts(obj)
     else{ 
          GetResid <- GetResiduals(obj)
          res <- GetResid$res
          order <- GetResid$order
     }
     if (squared.residuals) 
         res <- res^2
       n <- NROW(res)
       k <- NCOL(res)
       if (Maxlag*season >= n)
         stop("Maximum value of arguments lags * season can't exceed n!")
        df <- k^2*(lags-order)
        NegativeDF <- which(df<0)
        df[NegativeDF] <- 0
	     Accmat <- stats::acf(res, lag.max = Maxlag*season, plot = FALSE, type = "correlation")$acf
	     inveseR0 <- solve(Accmat[1,,])
             prodvec <- numeric(Maxlag*season)
	  for(l in 1:Maxlag){
            tvecR <- t(as.vector(Accmat[l*season+1,,]))
	        prodvec[l] <- 1/(n-l)*crossprod(t(tvecR),crossprod(t(kronecker(inveseR0,inveseR0)),t(tvecR)))
	  }
       Q <- n*(n+2)*cumsum(prodvec)
       STATISTIC <- Q[lags]
       PVAL <- 1 - stats::pchisq(STATISTIC,df)
       PVAL[NegativeDF] <- NA
       summary <- matrix(c(lags,STATISTIC,df,PVAL),ncol=4)
       dimnames(summary) <- list(rep("", length(STATISTIC)),c("lags","statistic","df","p-value"))
    return(summary)
}
