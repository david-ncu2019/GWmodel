# Renamed as gwr.multiscale, i.e. multiscale GWR
# Geographically Weighted Regression with Parameter-Specific Distance Metrics  
# With Flexiable bandwidth selection
# Method (Yang, p44, 2014): select an optimum bandwidth for each coefficient(intercept inclusive),
# within each step of back-fitting, employing an existing criteria for a basic GWR

# Improvements: Set the corresponding threshold for each selected bandwidth, i.e. fix the bandwidth if the change of the selected bandwidths is less than threshold
gw.fitted <- function(X, beta) {
  fitted(X, beta)
}

gwr.multiscale <- function(formula, data, kernel="bisquare", adaptive=FALSE, criterion="dCVR", max.iterations=2000,threshold=0.00001, dMats, var.dMat.indx, p.vals, theta.vals,longlat=FALSE,
                           bws0=NULL, bw.seled, approach = "AIC", bws.thresholds, bws.reOpts=5, verbose=F, hatmatrix=T, 
                           predictor.centered=rep(T, length(bws0)-1),nlower = 10, parallel.method=F,parallel.arg=NULL, force.armadillo=F)
{
  ##Record the start time
  timings <- list()
  timings[["start"]] <- Sys.time()
  ###################################macth the variables
  this.call <- match.call()
  p4s <- as.character(NA)
  ##Data points{
  if (inherits(data, "Spatial"))
  {
    p4s <- proj4string(data)
    dp.locat<-coordinates(data)
    regression.points <- data
    data <- as(data, "data.frame")
  }
  else if(inherits(data, "sf")) {
    if(any((st_geometry_type(data)=="POLYGON")) | any(st_geometry_type(data)=="MULTIPOLYGON"))
       dp.locat <- st_coordinates(st_centroid(st_geometry(data)))
    else
       dp.locat <- st_coordinates(st_geometry(data))
    regression.points <- data
  }
  else
  {
    stop("Given regression data must be Spatial*DataFrame")
  }
  #########Distance matrix is given or not
  dp.n <- nrow(dp.locat)
  ######Extract the data frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  var.n <- ncol(x)
  #x2 <- x
  idx1 <- match("(Intercept)", colnames(x))
  ##################################################
  #####Linear regression
  lms <- lm(formula,data=data)
  #lms <-fastLm(formula,data=data) 
  lms$x <- x
  lms$y <- y
  #colnames(x)[1]<-"Intercept"
  #################################################
  if(missing(bw.seled))
  {
    bw.seled <- rep(F, var.n)
  }
  else {
     bw.seled <- rep_len(bw.seled, var.n)
  }
  if(missing(bws.thresholds))
  {
    bws.thresholds <- rep(0.1, var.n)
  }
  else {
     bws.thresholds<- rep_len(bws.thresholds, var.n)
  }
  ##########Centering and scaling predictors
  if(missing(predictor.centered))
  {
    predictor.centered <- rep(T, length.out=var.n-1)
  }
  else {
     if(!is.na(idx1))
     {
       if(length(predictor.centered)!=var.n-1)
       {
         predictor.centered <- rep(T, length.out=var.n-1)
         warning("All the predictors will be centered, please check the parameter predictor.centered")
       }
     }
     else
     {
       if(length(predictor.centered)!=var.n)
       {
         predictor.centered <- rep(T, length.out=var.n)
         warning("All the predictors will be centered, please check the parameter predictor.centered")
       }
     }
  }
  
  
  n.cent <- length(which(predictor.centered))
  if(!is.na(idx1))
  {
    colnames(x)[idx1]<-"Intercept"
    x1 <- x[,-idx1]
  }
  else
  {
    x1 <- x
  }
  #n.scal <- length(which(predictor.centered))
  if (is.null(nrow(x1))) x1 <- matrix(x1, nrow=length(x1))
  if(n.cent>1)
    predictors.centered.means <- colMeans(x1[,predictor.centered])
  else 
    predictors.centered.means <- mean(x1[,predictor.centered])
  if(n.cent>1)
  {
    #meancent <- partial(scale, scale=FALSE)
    x1[,predictor.centered] <- scale(x1[,predictor.centered], scale=FALSE)
    
  }
  if(!is.na(idx1))
  {
    x1 <- cbind(1, x1)
  }
  colnames(x1) <-  colnames(x)
  
  #####Centering predictors
  #####Centered predictors x1
  #################################################
  
  ##Calculate the initial bandwidths: find an optimum bandwidth for each univariate model (dependent variable~each independent variable)
  allvars <- all.vars(formula)
  DeVar <- allvars[1]
  InDevars <- colnames(x)
  ##########################Check the distance matrices
  ###check the distance matrices and bandwidths
  if(missing(dMats))
  {
    dMats <- list()
    if(missing(p.vals))
    {
      dMat <- gw.dist(dp.locat,longlat=longlat)
      dMats[[1]] <- dMat
	  var.dMat.indx <- rep(1, var.n)
    }
    else
    {
      p.vals <- rep_len(p.vals, var.n)
      if(missing(theta.vals))
	    {
	      theta.vals <- rep_len(0, var.n)
	    } 
      else
        theta.vals<- rep_len(theta.vals, var.n)
      for(i in 1:var.n)
      {
        dMats[[i]] <- gw.dist(dp.locat, p=p.vals[i], theta=theta.vals[i],longlat=longlat)
      }
	    var.dMat.indx <- 1:var.n
    }
  }
  else if(is.list(dMats))
  {
    if(length(dMats) < max(var.dMat.indx))
      stop("Please specify a distance matrix index via the var.dMat.indx for each independent variable!")
    else
    {
      for(i in 1:length(dMats))
      {
        dim.dMat<-dim(dMats[[i]])
        if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
          stop("Dimensions of dMat are not correct")
      }
    }
  }
  else
  {
    stop("dMats are not correct")
  }  
  ############################### Intialize the bandwidth
  if(is.null(bws0))
  {
    bws0 <- numeric(var.n)
    cat("------ Calculate the initial bandwidths for each independent variable ------\n")
    for(i in 1:var.n)
    {
      if(InDevars[i]=="Intercept")
        fml <- Generate.formula(DeVar,c("1"))
      else
      {
        fml <- Generate.formula(DeVar,c(InDevars[i]))
        #fml <- paste(fml, "1", sep="-")
      }
      cat("Now select an optimum bandwidth for the model: ", fml,"\n")
      part1<-paste("bws0[i]<-bw.gwr(",fml,sep="")
      part2<-"data=regression.points,kernel=kernel,approach=approach,adaptive=adaptive,dMat=dMats[[var.dMat.indx[i]]], parallel.method=parallel.method,parallel.arg=parallel.arg)"
      expression<-paste(part1,part2,sep=",")
      print(expression)
      eval(parse(text=expression))
    }
    cat("------            The end for the initial selections              ------\n")
  }
  else
  {
    bws0 <- rep(bws0, length.out=var.n)
  }
  #################
  if(length(bw.seled)!=var.n)
    bw.seled <- rep(F, length.out=var.n)
  
  #################
  # Eigen code have only been tested for parallel.method = "omp" and parallel.method = F.
  # If another parallelization framework is requested, default to original implementation.
  if ((parallel.method == "omp" || parallel.method == F) && !force.armadillo) {
    kerneln <- NA
    if (kernel == "gaussian") kerneln <- 0
    if (kernel == "exponential") kerneln <- 1
    if (kernel == "bisquare") kerneln <- 2
    if (kernel == "tricube") kerneln <- 3
    if (kernel == "boxcar") kerneln <- 4
    if (is.na(kerneln)) kerneln <- 2
    
    approachn <- NA
    if (approach == "AIC" || approach == "aic" || approach == "AICc") approachn <- 0
    if (approach == "BIC" || approach == "bic") approachn <- 1
    if (approach == "CV" || approach == "cv") approachn <- 2
    if (is.na(approachn)) approachn <- 0
    
    crin <- NA
    if (criterion == "dCVR") crin <- 0
    if (criterion == "CVR") crin <- 1
    if (is.na(crin)) crin <- 0
    
    num.cores <- parallel.arg
    if (is.null(num.cores)) {
      num.cores <- 0
    }
    else if (!is.numeric(num.cores) || is.nan(num.cores) || parallel.method == F) {
      num.cores <- 1
    }
    
    new_ret <- new_multiscale(x, x1, dMats, dp.locat, y, bws0, var.dMat.indx,
                              adaptive, verbose, nlower, hatmatrix, 
                              max.iterations, threshold, num.cores, InDevars,
                              kerneln, approachn, crin, bws.reOpts)
    #####Output
    yhat <- as.matrix(new_ret[[1]])
    residual <- as.matrix(new_ret[[2]])
    betas <- new_ret[[3]]
    Shat <- new_ret[[4]]
    S.arrays <- array(new_ret[[5]], c(var.n, dp.n, dp.n))
    Beta_SE <- new_ret[[6]]
    Beta_TV <- new_ret[[7]]
    bws0 <- as.numeric(new_ret[[8]])
    bws.vars <- new_ret[[9]]
    mgwr.diag <- new_ret[[10]]
    GW.diagnostic <- NA
    ###############
    if(hatmatrix)
    {
      AIC <- mgwr.diag[[1]]
      AICc <- mgwr.diag[[2]]
      edf <- mgwr.diag[[3]]
      enp <- mgwr.diag[[4]]
      RSS.gw <- mgwr.diag[[5]]
      R2.val <- mgwr.diag[[6]]
      R2adj <- mgwr.diag[[7]]
      BIC <- mgwr.diag[[10]]
      tr.Shat <- mgwr.diag[[8]]
      
      GW.diagnostic<-list(RSS.gw=RSS.gw,AICc=AICc,AIC=AIC,BIC=BIC,R2.val=R2.val, R2adj = R2adj, edf=edf, enp=enp)
    }
  } else {
    cat("------   Calculate the initial beta0 from the above bandwidths    ------\n")
    #dMat <- gw.dist(dp.locat=dp.locat, p=2, theta=0, longlat=longlat)
    dMat <- dMats[[1]]
    bw.int0 <- bw.gwr2(x1, y, dp.locat,approach=approach,kernel=kernel,adaptive=adaptive,dMat,verbose=verbose, nlower = nlower, parallel.method=parallel.method,parallel.arg=parallel.arg)
    #betas <- gwr.q(x, y, dp.locat, adaptive=adaptive, bw=bw.int0, kernel=kernel,dMat=dMat)
    #######Re-initialize the betas by a simple back-fitting process
    #betas <- gwr.backfit(x, y, betas,dp.locat,dp.locat, FALSE,criterion, adaptive, bws0, kernel,dMats, max.iterations,threshold*100)
    #####Hatmatrix for the whole process
    if(hatmatrix)
    {
      Shat <- matrix(nrow=dp.n, ncol=dp.n)
      S.arrays <- array(dim=c(var.n, dp.n, dp.n))
      C <- array(dim=c(dp.n,var.n,dp.n))
      ####SEs and t-values
      Beta_SE <- matrix(nrow=dp.n, ncol=var.n)
      Beta_TV <- matrix(nrow=dp.n, ncol=var.n)
    }
    
    res <- gwr.q2(x1, y, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix,bw=bw.int0, kernel=kernel,dMat=dMat)
    betas <- res[[1]]
    if(hatmatrix)
    {
      Shat <- res[[2]]
      C <- res[[3]]
      idm <- diag(var.n)
      for(i in 1:var.n)
        for(j in 1:dp.n)
        {
          S.arrays[i,j,] <- x1[j,i]*(idm[i,]%*%C[j,,])
        }
    }
    #betas <- gwr.backfit(x, y, betas,dp.locat,dp.locat, FALSE,criterion, adaptive, bws0, kernel,dMats, max.iterations,0.0001)
    cat("------            The end for calculating the initial beta0              ------\n") 
    cat("------ Select the optimum bandwidths for each independent variable via ", approach, " aproach ------\n")
    ieration <-0 
    #bws1 <- bws0
    bws.vars <- bws0
    bws.change.NO <- numeric(var.n)
    criterion.val <- 10000000  
    #yhat.i <- betas*x
    #print(yhat.i)
    resid.i <- y - gw_fitted(x1, betas)
    RSS0 <- sum(resid.i^2)
    RSS1 <- 0
    RSS.vals <- c(RSS0, RSS1, criterion.val)
    cat("*****  The back-fitting process for model calibration with bandwiths selected *****\n")
    while((ieration < max.iterations) && criterion.val > threshold)
    { 
      #AICcs <- numeric(var.n)
      cat("    Iteration ", ieration+1, ":\n")   
      for(i in 1:var.n)
      {
        dMat <- dMats[[var.dMat.indx[i]]]  
        f.i <- betas[,i]*x1[,i]
        y.i <- resid.i + f.i
        if(bw.seled[i])
          bw.i <- bws0[i]
        else
        {
          cat("Now select an optimum bandwidth for the variable: ", InDevars[i], "\n")
          bw.i <- bw.gwr2(matrix(x1[, i], ncol = 1), y.i, dp.locat, approach = approach, kernel = kernel, adaptive = adaptive, dMat, verbose = verbose, nlower = nlower,parallel.method=parallel.method,parallel.arg=parallel.arg)
          cat("The newly selected bandwidth for variable ", InDevars[i])
          cat(" is: ", bw.i, "\n")
          cat("The bandwidth used in the last ieration is: ", bws0[i])
          cat(" and the difference between these two bandwidths is: ", abs(bw.i - bws0[i]), "\n")
          if (abs(bw.i - bws0[i]) > bws.thresholds[i]) {
            cat("The bandwidth for variable ", InDevars[i])
            cat(" will be continually selected in the next ieration.\n")
            bws.change.NO[i] <- 0
          }
          else {
            bws.change.NO[i] <- bws.change.NO[i] + 1
            if(bws.change.NO[i] < bws.reOpts)
            {
              cat("The bandwidth for variable ", InDevars[i])
              cat(" seems to be converged for ", bws.change.NO[i], " times.")
              cat("It will be continually optimized in the next ", bws.reOpts-bws.change.NO[i], " times\n")
            }
            else
            {
              cat("The bandwidth for variable ", InDevars[i])
              cat(" seems to be converged and will be kept the same in the following ierations.\n")
              bw.seled[i] <- T
            }
          }
        }
        bws0[i] <- bw.i       
        res <- gwr.q2(matrix(x1[,i], ncol=1), y.i, dp.locat, adaptive=adaptive, hatmatrix = hatmatrix,bw=bw.i, kernel=kernel,dMat=dMat)
        betai <- res[[1]]
        if(hatmatrix)
        {
          Si <- res[[2]]
          ###See Yu et al. 2018, eq. 18 on P8
          S.arrayi <- S.arrays[i,,]
          S.arrays[i,,] <- Si%*%S.arrayi + Si - Si%*%Shat
          Shat <- Shat- S.arrayi + S.arrays[i,,] 
        }
        
        #betai <- gwr.q(matrix(x[,i], ncol=1), y.i, loc=dp.locat, adaptive=adaptive, bw=bw.i, kernel=kernel,dMat=dMat)
        #AICcs[i] <- gwr.aic(bw.i, matrix(x[,i], ncol=1), y.i, kernel, adaptive, dp.locat, dMat=dMat, verbose=F)
        #yhat.i[,i] <- betai*x[,i]
        betas[,i] <- betai
        resid.i <- y - gw_fitted(x1, betas)
        #resid.i <- ehat(y.i, matrix(x[,i], ncol=1), matrix(betas[,i], ncol=1))
        #betas[,i] <- betai
      }
      bws.vars <- rbind(bws.vars, bws0)
      #AICc.vals <- rbind(AICc.vals, AICcs) 
      RSS1 <- sum((y - gw_fitted(x1, betas))^2)   
      if(criterion=="CVR")
      {
        criterion.val <- abs(RSS1-RSS0)
        cat("    Ieration ", ieration, "the change value of RSS (CVR) is: ", criterion.val,"\n")
      }
      else
      {
        criterion.val <- sqrt(abs(RSS1-RSS0)/RSS1)
        cat("    Ieration ", ieration+1, "the differential change value of RSS (dCVR) is: ", criterion.val,"\n") 
      }
      RSS0 <- RSS1  
      cat("----------End of    Iteration ", ieration+1, "----------\n") 
      RSS.vals <- rbind(RSS.vals, c(RSS0, RSS1, criterion.val)) 
      ieration <- ieration+1            
    }
    #####Output
    yhat <- gw_fitted(x1, betas)
    residual <- y - yhat
    GW.diagnostic <- NA
    ###############
    # RSS.gw <- RSS1
    # sigma.hat21 <- RSS.gw/dp.n
    # yss.g <- sum((y - mean(y))^2)
    # R2.val <-  1-RSS.gw/yss.g
    # R2adj <- NA
    # AICc <- NA
    # if(hatmatrix)
    # {
      # tr.Shat <- sum(diag(Shat))
      # tr.StShat<-sum(Shat^2)
      # edf<- dp.n - 2*tr.Shat + tr.StShat
      # AICc <- dp.n*log(sigma.hat21) + dp.n*log(2*pi) + dp.n *((dp.n + tr.Shat) / (dp.n - 2 - tr.Shat))
      # R2adj <- 1-(1-R2.val)*(dp.n-1)/(edf-1)
    # }
    if(hatmatrix)
    {
      mgwr.diag <- gwr_diag(y, x1, betas, Shat)
    	AIC <- mgwr.diag[[1]]
    	AICc <- mgwr.diag[[2]]
    	edf <- mgwr.diag[[3]]
    	enp <- mgwr.diag[[4]]
    	RSS.gw <- mgwr.diag[[5]]
    	R2.val <- mgwr.diag[[6]]
    	R2adj <- mgwr.diag[[7]]
    	BIC <- mgwr.diag[[10]]
      tr.Shat <- mgwr.diag[[8]]
      ####Calculate the SEs and t-values
      sigma.hat11 <- RSS.gw/(dp.n-tr.Shat)
      for(i in 1:var.n)
      {
        Ci <- diag(1/x[,i])%*%S.arrays[i,,]
        Beta_SE[,i] <- sqrt(diag(Ci%*%t(Ci)*sigma.hat11))
        Beta_TV[,i] <- betas[,i]/Beta_SE[,i]
      }
  	  GW.diagnostic<-list(RSS.gw=RSS.gw,AICc=AICc,AIC=AIC,BIC=BIC,R2.val=R2.val, R2adj = R2adj, edf=edf, enp=enp)
    }
  }  
  #sigma.hat11<-RSS.gw/(dp.n-2*tr.Shat+tr.StShat)
  
  #########################Retrive the intercept
  if(!is.na(idx1) && n.cent>=1)
  {
    idx.cent <- which(predictor.centered)
    betas.cent <- c()
    for(i in 1:n.cent)
      betas.cent <- cbind(betas.cent, betas[,idx.cent[i]+1]*predictors.centered.means[i])
    beta0 <- betas[,idx1] - apply(betas.cent,1, sum)
    betas[,idx1] <- beta0
  }
  ############################
  if(hatmatrix)
  {
    vdgwr.df <- data.frame(betas, yhat, residual, Beta_SE, Beta_TV)
    colnames(vdgwr.df) <- c(colnames(x), "yhat", "residual", paste(colnames(x), "SE", sep="_"),paste(colnames(x), "TV", sep="_"))
  }
  else
  {
    vdgwr.df <- data.frame(betas, yhat, residual)
    #colnames(vdgwr.df) <- c(colnames(x), "yhat", "residual",paste(colnames(betas), "SE", sep="_"), paste(colnames(betas), "TV", sep="_"))
    colnames(vdgwr.df) <- c(colnames(x), "yhat", "residual")
  }
  griddedObj <- F
  if(inherits(regression.points, "Spatial")) 
  { 
    if (is(regression.points, "SpatialPolygonsDataFrame"))
    {
      polygons<-polygons(regression.points)
      #SpatialPolygons(regression.points)
      #rownames(gwres.df) <- sapply(slot(polygons, "polygons"),
      #  function(i) slot(i, "ID"))
      SDF <-SpatialPolygonsDataFrame(Sr=polygons, data=vdgwr.df,match.ID=F)
    }
    else
    {
      griddedObj <- gridded(regression.points)
      SDF <- SpatialPointsDataFrame(coords=dp.locat, data=vdgwr.df, proj4string=CRS(p4s), match.ID=F)
      gridded(SDF) <- griddedObj 
    }
  }
  else if(inherits(regression.points, "sf"))
  {
     SDF <- st_sf(vdgwr.df, geometry = st_geometry(regression.points))
  }
  else
    SDF <- SpatialPointsDataFrame(coords=dp.locat, data=vdgwr.df, proj4string=CRS(p4s), match.ID=F)  
  #GW.arguments<-list(formula=formula,rp.given=T,hatmatrix=T,criterion="RSS",bws=bws0, kernel=kernel,adaptive=adaptive)  
  timings[["stop"]] <- Sys.time()
  GW.arguments<-list(formula=formula,criterion=criterion,bws=bws0, kernel=kernel,adaptive=adaptive, hatmatrix=hatmatrix)  
  #res <- list(SDF=SDF, GW.arguments=GW.arguments,lm=lms,timings=timings,this.call=this.call)
  #class(res) <- "psdmgwr"
  #invisible(res)
  res <- list(SDF=SDF,GW.arguments=GW.arguments,GW.diagnostic=GW.diagnostic,lm=lms, bws.vars,timings=timings,this.call=this.call)
  class(res) <- "multiscalegwr"
  invisible(res) 
}

############################Layout function for outputing the Multiscale(PSDM) GWR results
##Author: BL
print.multiscalegwr<-function(x, ...)
{
  if(!inherits(x, "multiscalegwr")) stop("It's not a multi-scale gwr object")
  cat("   ***********************************************************************\n")
  cat("   *                       Package   GWmodel                             *\n")
  cat("   ***********************************************************************\n")
  cat("   Program starts at:", as.character(x$timings$start), "\n")
  cat("   Call:\n")
  cat("   ")
  print(x$this.call)
  vars<-all.vars(x$GW.arguments$formula)
  var.n<-length(x$lm$coefficients)
	cat("\n   Dependent (y) variable: ",vars[1])
	cat("\n   Independent variables: ",vars[-1])
	dp.n<-length(x$lm$residuals)
	cat("\n   Number of data points:",dp.n)
	#########################################################################
	cat("\n   ***********************************************************************\n")
    cat("   *                       Multiscale (PSDM) GWR                          *\n")
	cat("   ***********************************************************************\n")
	cat("\n   *********************Model calibration information*********************\n")
	cat("   Kernel function:", x$GW.arguments$kernel, "\n")
	if(x$GW.arguments$adaptive)
	   cat("   Adaptive bandwidths for each coefficient(number of nearest neighbours): \n") 
  else
     cat("   Fixed bandwidths for each coefficient: \n")
	bws <- matrix(x$GW.arguments$bws,nrow=1)
  rownames(bws) <- c("   Bandwidth ")
  colnames(bws) <- names(x$lm$coefficients)
  printCoefmat(bws)
	cat("\n   *************Summary of multiscale GWR coefficient estimates:***************\n")       
		if(inherits(x$SDF, "Spatial"))
       df0 <- as(x$SDF, "data.frame")[,1:var.n, drop=FALSE]
    else
       df0 <- st_drop_geometry(x$SDF)[,1:var.n, drop=FALSE]
        if (any(is.na(df0))) {
            df0 <- na.omit(df0)
            warning("NAs in coefficients dropped")
        }
	CM <- t(apply(df0, 2, summary))[,c(1:3,5,6)]
	if(var.n==1) 
    { 
      CM <- matrix(CM, nrow=1)
      colnames(CM) <- c("Min.", "1st Qu.", "Median", "3rd Qu.", "Max.")
      rownames(CM) <- names(x$SDF)[1]
    }
	rnames<-rownames(CM)
		for (i in 1:length(rnames))
			 rnames[i]<-paste("   ",rnames[i],sep="")
	rownames(CM) <-rnames 
	printCoefmat(CM)
	cat("   ************************Diagnostic information*************************\n")
	#rownames(diag.mat) <- c("   Residual sum of squares", "   AICc value:","   R-square value","   Adjusted R-square value")
	options("scipen"=999)
	if(x$GW.arguments$hatmatrix)
	{
	  cat("   Number of data points:", dp.n, "\n")
	  cat("   Effective number of parameters (2trace(S) - trace(S'S)):", x$GW.diagnostic$enp, "\n")
	  cat("   Effective degrees of freedom (n-2trace(S) + trace(S'S)):", x$GW.diagnostic$edf, "\n")
	  cat("   AICc value: ", x$GW.diagnostic$AICc, "\n")
	  cat("   AIC value: ", x$GW.diagnostic$AIC, "\n")
	  cat("   BIC value: ", x$GW.diagnostic$BIC, "\n")
	  cat("   Residual sum of squares: ", x$GW.diagnostic$RSS.gw, "\n")
	  cat("   R-square value: ", x$GW.diagnostic$R2.val, "\n")
	  cat("   Adjusted R-square value: ", x$GW.diagnostic$R2adj, "\n")
	}
	cat("\n   ***********************************************************************\n")
	cat("   Program stops at:", as.character(x$timings$stop), "\n")
	invisible(x)
}

##Back-fitting algorithm, see Wenbai's thesis
gwr.backfit <- function(x, y, betas,dp.locat,rp.locat,hatmatrix, criterion="CVR", adaptive=F, bws, kernel="bisquare",dMats, max.iterations=100,threshold=0.0001)
{
   ieration <-0
   y.i <- y
   var.n <- ncol(x)
   dp.n <- nrow(x)
   rp.n <- nrow(rp.locat)
   criterion.val <- 10000000
   resid.i <- y - gw_fitted(x, betas)
   #resid.i <- ehat(y, x, betas)
   RSS0 <- sum(resid.i^2)
   RSS1 <- 0
   #if(criterion=="RSS")
   #{
     #
     #
   
   #  RSS0 <- rss(y, x, betas)
  #   
   #}
   #else
   #{
    #  betas0 <- betas  
   #}  
   cat("*************The back-fitting process*************\n")
   while((ieration < max.iterations) && criterion.val > threshold)
   {    
     for(i in 1:var.n)
     {
       bw <- bws[i]
       dMat <- dMats[[i]] 
       f.i <- betas[,i]*x[,i]
       y.i <- f.i+resid.i
       betai <- gwr.q(matrix(x[,i], ncol=1), y.i, loc=dp.locat, adaptive=adaptive, bw=bw, kernel=kernel,dMat=dMat)
       
       resid.i <- y.i - betai*x[,i]
       #resid.i <- ehat(y.i, matrix(x[,i], ncol=1), matrix(betas[,i], ncol=1))
       betas[,i] <- betai
     }
     #RSS1 <- rss(y,x,betas)
     RSS1 <- sum((y - gw_fitted(x, betas))^2)   
     if(criterion=="CVR")
     {
         #RSS1 <- sum((y - gwr.fitted(x, betas))^2)
         #cat("RSS1: ", RSS1,"\n")
         #cat("RSS0: ", RSS0,"\n")
         criterion.val <- abs(RSS1-RSS0)
         cat("    Ieration ", ieration, "the change value of RSS (CVR) is: ", criterion.val,"\n")
     }
     else
     {
          criterion.val <- sqrt(abs(RSS1-RSS0)/RSS1)
          cat("    Ieration ", ieration, "the differential change value of RSS (dCVR) is: ", criterion.val,"\n") 
     }
     RSS0 <- RSS1  
     ieration <- ieration+1            
   }
   #calculate the hatmatrix
   betas
}
###For bandwidth selection
bw.gwr2<-function(x, y, dp.locat,approach="AIC",kernel="gaussian",adaptive=FALSE,dMat, verbose=F, parallel.method=F,parallel.arg=NULL, nlower = 10)
{
  dp.n <-  nrow(dp.locat)
  if(adaptive)
  {
    upper <- dp.n
    lower <- nlower
  }
  else
  {
      upper<-range(dMat)[2]
      lower<-upper/5000
  }
  ########################## Now the problem for the golden selection is too computationally heavy
    #Select the bandwidth by golden selection
    bw<-NA
    # ## make cluseter
  if (parallel.method == "cluster") {
    if (missing(parallel.arg)) {
      cl.n <- max(detectCores() - 4, 2)
      parallel.arg <- makeCluster(cl.n)
    } else cl.n <- length(parallel.arg)
    clusterCall(parallel.arg, function() { library(GWmodel) })
  }
  # ## call for functions
  if(approach == "bic" || approach == "BIC")
      bw <- gold(gwr.bic, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p=2, theta=0, longlat=F, dMat, verbose, parallel.method, parallel.arg)
  else if(approach == "aic" || approach == "AIC" || approach == "AICc")
      bw <- gold(gwr.aic, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p=2, theta=0, longlat=F, dMat, verbose, parallel.method, parallel.arg)    
  else 
      bw <- gold(gwr.cv, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p=2, theta=0, longlat=F, dMat, verbose, parallel.method, parallel.arg)
  # ## stop cluster
  if (parallel.method == "cluster") {
    if (missing(parallel.arg)) stopCluster(parallel.arg)
  }
  bw
}
gwr.cv2 <- function(bw, X, Y, kernel,adaptive, dp.locat,dMat)
{
   dp.n<-length(dp.locat[,1])
   var.n <- ncol(X)
   betas <- matrix(nrow=dp.n,ncol=var.n) 
  CV<-numeric(dp.n)
  for (i in 1:dp.n)
  {
    dist.vi<-dMat[,i]
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    ##lm.i <- try(lm.wfit(y = y, x = x, w = w.i))
    fun1<-function(X,Y,W.i) {betai<- solve(t(X*W.i)%*%X)%*%{t(X*W.i)%*%Y}}
    betas[i,] <- try(fun1(X,Y,W.i))
    #gw.resi <- try(lm.wfit(y = Y, x = X, w = W.i))
    W.i[i]<-0
    gw.resi<-try(fun1(X,Y,W.i))
    if(!inherits(gw.resi, "try-error"))
    {
      #b <- coefficients(gw.resi)
      yhat.noi<-X[i,]%*%gw.resi
      #CV[i] <- Y[i] - (t(b) %*% X[i,])
      CV[i]<-Y[i]-yhat.noi
      
    }
    else
    {
      CV[i]<-Inf
      break
    }
  }
  CV.score<-t(CV) %*% CV
  res<- list(CV.score, betas)
}
gwr.q2 <- function(x, y, loc, adaptive=F, hatmatrix,bw=sqrt(var(loc[,1])+var(loc[,2])),
                  kernel, p, theta, longlat,dMat, wt2=rep(1,nrow(loc)))
{
  if (missing(dMat))
     DM.given <- F
  else
     DM.given <- T
  dp.n <- nrow(loc)
  var.n <- ncol(x)
  betas <- matrix(nrow=dp.n, ncol=var.n)
  S <- matrix(nrow=dp.n, ncol=dp.n)
  C <- array(dim=c(dp.n, var.n, dp.n))
  for (i in 1:dp.n)
  {
    if(DM.given)
       dist.vi <- dMat[,i]
    else
       dist.vi <- gw.dist(loc, focus=i, p, theta, longlat)
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    gw.resi<-gw_reg(x,y,as.vector(W.i*wt2),hatmatrix=hatmatrix,i)
    betas[i,]<-gw.resi[[1]]
    if(hatmatrix)
    {
      S[i,]<-gw.resi[[2]]
      C[i,,]<-gw.resi[[3]]
    }
  }
  colnames(betas) <- colnames(x)
  res <- list(betas, S,C)
  res
}

####Calculate the AICc with a given bandwidth
##Author: Binbin Lu
gwr.aic1<-function(bw, X, Y, kernel="bisquare",adaptive=FALSE, dp.locat, p=2, theta=0, longlat=F,dMat, verbose=T)
{
   dp.n<-length(dp.locat[,1])
   var.n <- ncol(X)
   #########Distance matrix is given or not

  if (is.null(dMat))
      DM.given<-F
  else
  {
    DM.given<-T
    dim.dMat<-dim(dMat)
    if (dim.dMat[1]!=dp.n||dim.dMat[2]!=dp.n)
    stop ("Dimensions of dMat are not correct")
  }
  ############################################AIC
  ###In this function, the whole hatmatrix is not fully calculated and only the diagonal elements are computed
  S<-matrix(nrow=dp.n,ncol=dp.n)
  betas <-matrix(nrow=dp.n, ncol=var.n)
  for (i in 1:dp.n)
  {
    if (DM.given)
         dist.vi<-dMat[,i]
    else
    {
       dist.vi<-gw.dist(dp.locat=dp.locat, focus=i, p=p, theta=theta, longlat=longlat)
    }
    W.i<-gw.weight(dist.vi,bw,kernel,adaptive)
    res<- try(gw_reg(X,Y,W.i,TRUE,i))
    #Ci=solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
    #fun2<-function(X,W.i) {Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}}
    #Ci<-try(fun2(X,W.i))
    
    #Ci<-solve(t(X*W.i)%*%X)%*%{t(X*W.i)}
   # gw.resi<-gw.reg(X,Y,W.i,hatmatrix=T,focus=i)
    #betas[i,]<-gw.resi[[1]] ######See function by IG
    #S[i,]<-gw.resi[[2]]
    if(!inherits(res, "try-error"))
    {
      S[i,]<-res[[2]]   
      betas[i,] <- res[[1]]
    }
    else
    {
      S[i,]<-Inf
      break
    }  
  }
  
  if (!any(is.infinite(S)))
  {
    #tr.S<-sum(diag(S))
#    RSS.gw<-t(Y)%*%t(diag(dp.n)-S)%*%(diag(dp.n)-S)%*%Y
#    sigma.hat2 <- RSS.gw/dp.n
#    AICc<-dp.n*log(sigma.hat2) + dp.n*log(2*pi) + dp.n *((dp.n + tr.S) / (dp.n - 2 - tr.S))
     AICc<-AICc(Y,X,betas, S)
  }
  else
    AICc<-Inf
  if(verbose)
  {     
    if(adaptive)
      cat("Adaptive bandwidth (number of nearest neighbours):", bw, "AICc value:", AICc, "\n")
    else
      cat("Fixed bandwidth:", bw, "AICc value:", AICc, "\n")
  }
  if(is.nan(AICc))
      AICc <- Inf
  AICc
}




















