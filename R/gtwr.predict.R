################################################################################
# gtwr.predict: Enhanced spatiotemporal prediction with GTWR models
#
# Extends standard GTWR prediction capabilities with timestamp preservation,
# statistical metrics, and original value tracking for validation.
################################################################################

gtwr.predict <- function(formula, data, obs.tv, predictdata, pred.tv, st.bw,
                         kernel = "bisquare", adaptive = FALSE, p = 2, theta = 0, 
                         longlat = FALSE, lamda = 0.05, t.units = "auto", ksi = 0, 
                         st.dMat1, st.dMat2, parallel_processing = FALSE, 
                         chunk_size = 1000, calculate_variance = TRUE) {
  
  ## Record execution time and setup
  timings <- list()
  timings[["start"]] <- Sys.time()
  this.call <- match.call()
  p4s <- as.character(NA)
  predict.SPDF <- NULL
  
  #=============================================================================
  # Track prediction mode (calibration vs. new locations)
  #=============================================================================
  predicting_at_calib_points <- FALSE
  
  #=============================================================================
  # Input validation
  #=============================================================================
  
  # Validate time variables
  validate_time_format <- function(time_var, var_name) {
    if (!inherits(time_var, c("Date", "POSIXlt", "POSIXct", "numeric", "yearmon", "yearqtr"))) {
      stop(paste0(var_name, " must be in an accepted time format: Date, POSIXlt, POSIXct, numeric, yearmon, or yearqtr"))
    }
  }
  
  if (missing(obs.tv)) stop("Time stamps 'obs.tv' for calibration data required")
  validate_time_format(obs.tv, "obs.tv")
  
  # Validation of predictdata and pred.tv will happen in their respective sections
  
  #=============================================================================
  # Process calibration data
  #=============================================================================
  
  # Extract spatial information from calibration data
  if (inherits(data, "Spatial")) {
    p4s <- proj4string(data)
    fd.locat <- coordinates(data)
    data.df <- as(data, "data.frame")
  } else if (inherits(data, "sf")) {
    p4s <- st_crs(data)
    if (any((st_geometry_type(data) == "POLYGON")) | any(st_geometry_type(data) == "MULTIPOLYGON")) {
      fd.locat <- st_coordinates(st_centroid(st_geometry(data)))
    } else {
      fd.locat <- st_coordinates(st_geometry(data))
    }
    data.df <- st_drop_geometry(data)
  } else {
    stop("Calibration data 'data' must be a Spatial* or sf object")
  }
  fd.n <- nrow(fd.locat)
  
  # Validate observation timestamps
  if (length(obs.tv) != fd.n) {
    stop("Number of timestamps in 'obs.tv' must match number of observations in 'data'")
  }
  
  # Extract variables from formula
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$data <- data.df  # Use data.df instead of data to avoid spatial object issues
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.extract(mf, "response")
  x <- model.matrix(mt, mf)
  var.n <- ncol(x)
  
  # Store original response values for validation metrics
  original_values <- y
  
  # Standardize column names
  idx1 <- match("(Intercept)", colnames(x))
  if (!is.na(idx1)) {
    colnames(x)[idx1] <- "Intercept"
  }
  inde_vars <- colnames(x)
  if ("Intercept" %in% inde_vars) {
    inde_vars <- inde_vars[inde_vars != "Intercept"]
  }
  
  #=============================================================================
  # Process prediction data
  #=============================================================================
  
  pd.given <- !missing(predictdata)
  
  if (!pd.given) {
    # Use calibration data if prediction data not provided
    warning("Prediction data not provided. Using calibration data for prediction.")
    predictdata <- data
    pred.tv <- obs.tv
    pd.locat <- fd.locat
    predict.SPDF <- data
    predictdata.df <- data.df
    predicting_at_calib_points <- TRUE
  } else {
    # Validate and extract prediction data
    if (missing(pred.tv)) stop("Time stamps 'pred.tv' for prediction data required")
    validate_time_format(pred.tv, "pred.tv")
    
    # Check time format consistency between calibration and prediction
    if (class(obs.tv)[1] != class(pred.tv)[1]) {
      warning(paste0("Different time formats used: obs.tv (", class(obs.tv)[1], 
                     ") vs pred.tv (", class(pred.tv)[1], "). Consider standardizing."))
    }
    
    # Extract spatial information from prediction data
    if (inherits(predictdata, "Spatial")) {
      if (is.na(p4s)) p4s <- proj4string(predictdata)
      pd.locat <- coordinates(predictdata)
      predict.SPDF <- predictdata
      predictdata.df <- as(predictdata, "data.frame")
    } else if (inherits(predictdata, "sf")) {
      if (is.na(p4s)) p4s <- st_crs(predictdata)
      if (any((st_geometry_type(predictdata) == "POLYGON")) | any(st_geometry_type(predictdata) == "MULTIPOLYGON")) {
        pd.locat <- st_coordinates(st_centroid(st_geometry(predictdata)))
      } else {
        pd.locat <- st_coordinates(st_geometry(predictdata))
      }
      predict.SPDF <- predictdata
      predictdata.df <- st_drop_geometry(predictdata)
    } else {
      stop("Prediction data 'predictdata' must be a Spatial* or sf object")
    }
    
    # Validate prediction timestamps
    if (length(pred.tv) != nrow(pd.locat)) {
      stop("Number of timestamps in 'pred.tv' must match number of observations in 'predictdata'")
    }
    
    # Check if all required variables exist in prediction data
    if (!all(inde_vars %in% names(predictdata.df))) {
      missing_vars <- inde_vars[!inde_vars %in% names(predictdata.df)]
      stop(paste("Missing variables in prediction data:", paste(missing_vars, collapse=", ")))
    }
    
    # Check if prediction points are identical to calibration points
    if (!is.null(predict.SPDF) && identical(coordinates(data), coordinates(predictdata))) {
      predicting_at_calib_points <- TRUE
    }
  }
  
  pd.n <- nrow(pd.locat)
  
  # Create prediction design matrix
  if (length(inde_vars) > 0) {
    x.p <- as.matrix(predictdata.df[, inde_vars, drop = FALSE])
    x.p <- cbind(Intercept = rep(1, pd.n), x.p)
  } else {
    # Intercept-only model
    x.p <- matrix(1, nrow = pd.n, ncol = 1, dimnames = list(NULL, "Intercept"))
  }
  
  # Ensure column order matches calibration data
  if (!identical(colnames(x.p), colnames(x))) {
    x.p <- x.p[, colnames(x), drop = FALSE]
  }
  
  #=============================================================================
  # Calculate spatiotemporal distance matrices
  #=============================================================================
  
  st.DM1.given <- !missing(st.dMat1)
  st.DM2.given <- !missing(st.dMat2)
  
  # Calibration-to-prediction points distance matrix
  if (!st.DM1.given) {
    cat("Calculating spatiotemporal distances between calibration and prediction points...\n")
    
    # Use chunking for large matrices
    if (fd.n * pd.n > 10000) {
      cat("Large matrices detected - using chunked calculation...\n")
      st.dMat1 <- matrix(0, fd.n, pd.n)
      n_chunks <- ceiling(fd.n / chunk_size)
      
      for (i in 1:n_chunks) {
        start_idx <- (i - 1) * chunk_size + 1
        end_idx <- min(i * chunk_size, fd.n)
        chunk_indices <- start_idx:end_idx
        
        cat(sprintf("  Processing chunk %d of %d (rows %d to %d)...\n", 
                   i, n_chunks, start_idx, end_idx))
        
        st.dMat1[chunk_indices, ] <- st.dist(
          dp.locat = fd.locat[chunk_indices, , drop = FALSE],
          rp.locat = pd.locat,
          obs.tv = obs.tv[chunk_indices],
          reg.tv = pred.tv,
          p = p, theta = theta, longlat = longlat,
          lamda = lamda, t.units = t.units, ksi = ksi
        )
        
        # Force garbage collection
        gc()
      }
    } else {
      st.dMat1 <- st.dist(
        dp.locat = fd.locat,
        rp.locat = pd.locat,
        obs.tv = obs.tv,
        reg.tv = pred.tv,
        p = p, theta = theta, longlat = longlat,
        lamda = lamda, t.units = t.units, ksi = ksi
      )
    }
    
    st.DM1.given <- TRUE
  } else {
    # Validate provided distance matrix
    if (!is.matrix(st.dMat1) || nrow(st.dMat1) != fd.n || ncol(st.dMat1) != pd.n) {
      stop("Provided st.dMat1 dimensions incorrect (expected ", fd.n, "×", pd.n, ")")
    }
    cat("Using provided spatiotemporal distance matrix for prediction.\n")
  }
  
  # Calibration-to-calibration points distance matrix (for variance)
  if (calculate_variance && !st.DM2.given) {
    cat("Calculating spatiotemporal distances between calibration points (for variance)...\n")
    
    # Use chunking for large matrices
    if (fd.n * fd.n > 10000) {
      cat("Large matrix detected - using chunked calculation...\n")
      st.dMat2 <- matrix(0, fd.n, fd.n)
      n_chunks <- ceiling(fd.n / chunk_size)
      
      for (i in 1:n_chunks) {
        start_idx <- (i - 1) * chunk_size + 1
        end_idx <- min(i * chunk_size, fd.n)
        chunk_indices <- start_idx:end_idx
        
        cat(sprintf("  Processing chunk %d of %d (rows %d to %d)...\n", 
                   i, n_chunks, start_idx, end_idx))
        
        st.dMat2[chunk_indices, ] <- st.dist(
          dp.locat = fd.locat[chunk_indices, , drop = FALSE],
          rp.locat = fd.locat,
          obs.tv = obs.tv[chunk_indices],
          reg.tv = obs.tv,
          p = p, theta = theta, longlat = longlat,
          lamda = lamda, t.units = t.units, ksi = ksi
        )
        
        # Force garbage collection
        gc()
      }
    } else {
      st.dMat2 <- st.dist(
        dp.locat = fd.locat,
        rp.locat = fd.locat,
        obs.tv = obs.tv,
        reg.tv = obs.tv,
        p = p, theta = theta, longlat = longlat,
        lamda = lamda, t.units = t.units, ksi = ksi
      )
    }
    
    st.DM2.given <- TRUE
  } else if (calculate_variance && st.DM2.given) {
    # Validate provided distance matrix
    if (!is.matrix(st.dMat2) || nrow(st.dMat2) != fd.n || ncol(st.dMat2) != fd.n) {
      stop("Provided st.dMat2 dimensions incorrect (expected ", fd.n, "×", fd.n, ")")
    }
    cat("Using provided spatiotemporal distance matrix for calibration points.\n")
  } else if (!calculate_variance) {
    cat("Skipping calibration distance matrix calculation (variance not requested).\n")
  }
  
  #=============================================================================
  # Estimate coefficients at prediction locations
  #=============================================================================
  
  cat("Estimating GTWR coefficients at prediction locations...\n")
  betas1 <- matrix(nrow = pd.n, ncol = var.n)
  colnames(betas1) <- colnames(x)
  xtxinv <- array(0, dim = c(pd.n, var.n, var.n))
  wt <- matrix(nrow = fd.n, ncol = pd.n)
  
  # Tracking for errors during estimation
  estimation_errors <- 0
  max_failures <- round(pd.n * 0.1)  # 10% threshold
  
  # Parallel processing implementation
  if (parallel_processing && requireNamespace("parallel", quietly = TRUE)) {
    cat("Using parallel processing for coefficient estimation...\n")
    
    # Determine optimal number of cores
    n_cores <- min(parallel::detectCores() - 1, pd.n, 8)  # Cap at 8 cores for safety
    cat("  Using", n_cores, "cores\n")
    
    # Create cluster and export necessary objects
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterExport(cl, c("st.dMat1", "x", "y", "st.bw", 
                                 "kernel", "adaptive", "pd.given"))
    
    # Load required functions in cluster
    parallel::clusterEvalQ(cl, {
      # Load required functions if not already in environment
      if (!exists("gw.weight", mode = "function")) {
        # Assumes these functions are available when sourced
        source("path/to/gw.weight.R")  # Adjust path as needed
      }
      if (!exists("gw_reg_1", mode = "function")) {
        # These functions typically come from a package
        # If using GWmodel package, load it
        if (requireNamespace("GWmodel", quietly = TRUE)) {
          library(GWmodel)
        }
      }
    })
    
    # Process each prediction point in parallel
    cat("  Starting parallel computations...\n")
    results <- parallel::parLapply(cl, 1:pd.n, function(i) {
      dist.vi <- st.dMat1[, i]
      W.i <- gw.weight(dist.vi, st.bw, kernel, adaptive)
      
      # Handle self-prediction case
      if (!pd.given) W.i[i] <- 0
      
      # Perform weighted regression
      gw.resi <- try(gw_reg_1(x, y, W.i), silent = TRUE)
      
      if (!inherits(gw.resi, "try-error")) {
        return(list(
          beta = gw.resi[[1]], 
          xtxinv = gw.resi[[2]], 
          wt = W.i, 
          error = FALSE
        ))
      } else {
        return(list(
          beta = rep(NA, ncol(x)), 
          xtxinv = matrix(NA, ncol(x), ncol(x)), 
          wt = W.i, 
          error = TRUE,
          error_msg = conditionMessage(attr(gw.resi, "condition"))
        ))
      }
    })
    
    # Close cluster
    parallel::stopCluster(cl)
    
    # Extract results and count errors
    for (i in 1:pd.n) {
      betas1[i, ] <- results[[i]]$beta
      xtxinv[i, , ] <- results[[i]]$xtxinv
      wt[, i] <- results[[i]]$wt
      
      if (results[[i]]$error) {
        estimation_errors <- estimation_errors + 1
        if (estimation_errors <= 5) { # Limit warning messages
          warning(paste("Coefficient estimation failed for prediction point", i, 
                       "- Error:", results[[i]]$error_msg))
        }
      }
    }
    
    if (estimation_errors > 5) {
      warning(paste("Additional", estimation_errors - 5, "estimation failures occurred"))
    }
    
  } else {
    # Sequential processing
    cat("Using sequential processing for coefficient estimation...\n")
    
    pb <- txtProgressBar(min = 0, max = pd.n, style = 3)
    for (i in 1:pd.n) {
      dist.vi <- st.dMat1[, i]
      W.i <- gw.weight(dist.vi, st.bw, kernel, adaptive)
      
      # Handle self-prediction case
      if (!pd.given) W.i[i] <- 0
      
      wt[, i] <- W.i
      
      # Perform weighted regression
      gw.resi <- try(gw_reg_1(x, y, W.i), silent = TRUE)
      
      if (!inherits(gw.resi, "try-error")) {
        betas1[i, ] <- gw.resi[[1]]
        xtxinv[i, , ] <- gw.resi[[2]]
      } else {
        estimation_errors <- estimation_errors + 1
        betas1[i, ] <- NA
        xtxinv[i, , ] <- NA
        
        if (estimation_errors <= 5) { # Limit warning messages
          warning(paste("Coefficient estimation failed for prediction point", i, 
                       "- Error:", conditionMessage(attr(gw.resi, "condition"))))
        }
      }
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    if (estimation_errors > 5) {
      warning(paste("Additional", estimation_errors - 5, "estimation failures occurred"))
    }
  }
  
  # Check if estimation was largely successful
  if (estimation_errors > max_failures) {
    stop(paste("Excessive coefficient estimation failures (>", max_failures, 
               "points). Check model specification and data."))
  }
  
  valid_betas <- !is.na(betas1[, 1])  # Check NAs based on first coefficient
  cat("Coefficient estimation completed with", sum(!valid_betas), "failures out of", pd.n, "points.\n")
  
  #=============================================================================
  # Calculate predictions
  #=============================================================================
  
  cat("Calculating predictions...\n")
  gw.predict <- rep(NA, pd.n)
  
  if (any(valid_betas)) {
    # Use gw_fitted for matrix multiplication: yhat_i = x_i * beta_i
    gw.predict <- gw_fitted(x.p, betas1)
  } else {
    warning("No valid coefficients estimated. All predictions will be NA.")
  }
  
  #=============================================================================
  # Calculate prediction variance (if requested)
  #=============================================================================
  
  calculate_variance <- calculate_variance && st.DM2.given
  predict.var <- rep(NA, pd.n)
  sigma.hat <- NA
  
  if (calculate_variance) {
    cat("\nCalculating prediction variance (Note: Method adapted from GWR)...\n")
    cat("  This adaptation to GTWR has not been fully validated in statistical literature.\n")
    cat("  Interpret variance estimates with appropriate caution.\n\n")
    
    # First, estimate sigma.hat based on fitting the model to calibration data
    cat("Estimating model error variance from calibration data fit...\n")
    S <- matrix(0, nrow = fd.n, ncol = fd.n)
    betas2 <- matrix(NA, nrow = fd.n, ncol = var.n)
    
    # Calculate hat matrix and coefficients at calibration points
    estimation_ok <- TRUE
    variance_errors <- 0
    max_var_failures <- round(fd.n * 0.1)  # 10% threshold
    
    pb <- txtProgressBar(min = 0, max = fd.n, style = 3)
    for (j in 1:fd.n) {
      dist.vj <- st.dMat2[, j]
      W.j <- gw.weight(dist.vj, st.bw, kernel, adaptive)
      
      # Perform weighted regression with hat matrix
      gw.resi.calib <- try(gw_reg(x, y, W.j, TRUE, j), silent = TRUE)
      
      if (!inherits(gw.resi.calib, "try-error")) {
        betas2[j, ] <- gw.resi.calib[[1]]
        S[j, ] <- gw.resi.calib[[2]]  # Store row j of hat matrix
      } else {
        variance_errors <- variance_errors + 1
        betas2[j, ] <- NA
        S[j, ] <- NA
        
        if (variance_errors <= 3) {  # Limit warnings
          warning(paste("Variance estimation failed for calibration point", j, 
                       "- Error:", conditionMessage(attr(gw.resi.calib, "condition"))))
        }
        
        if (variance_errors > max_var_failures) {
          warning("Excessive failures in variance estimation. Aborting variance calculation.")
          estimation_ok <- FALSE
          break
        }
      }
      
      setTxtProgressBar(pb, j)
    }
    close(pb)
    
    if (estimation_ok) {
      # Check for too many NA rows in hat matrix
      na_rows <- rowSums(is.na(S)) > 0
      if (sum(na_rows) > max_var_failures) {
        warning("Too many NA rows in hat matrix. Variance calculation may be unreliable.")
      }
      
      # Handle NA rows in hat matrix - replace with zeros to allow calculation to proceed
      # This is a compromise solution - better would be to re-estimate or use a more sophisticated approach
      if (any(na_rows)) {
        S[na_rows, ] <- 0
        warning(paste(sum(na_rows), "rows in hat matrix had missing values and were replaced with zeros."))
      }
      
      # Calculate necessary traces
      tr.S <- sum(diag(S))  # Trace of S
      
      # For tr.StS, we need a more accurate computation than sum(S^2)
      # We need tr(S'S) = sum(diag(S'S)) = sum of squared elements of S
      tr.StS <- sum(S^2)
      
      # Check if traces are valid
      if (!is.finite(tr.S) || !is.finite(tr.StS)) {
        warning("Trace calculations resulted in non-finite values. Skipping variance calculation.")
        estimation_ok <- FALSE
      } else {
        # Calculate Residual Sum of Squares (RSS)
        # Directly calculate residuals using estimated coefficients
        yhat_calib <- numeric(fd.n)
        for (j in 1:fd.n) {
          if (!any(is.na(betas2[j, ]))) {
            yhat_calib[j] <- sum(x[j, ] * betas2[j, ])
          } else {
            yhat_calib[j] <- NA
          }
        }
        
        # Handle potential NA values in yhat_calib
        na_yhat <- is.na(yhat_calib)
        if (any(na_yhat)) {
          warning(paste(sum(na_yhat), "fitted values in calibration had missing values."))
          if (sum(!na_yhat) < 0.5 * fd.n) {
            warning("Too few valid fitted values for reliable variance estimation.")
            estimation_ok <- FALSE
          }
        }
        
        if (estimation_ok) {
          residual_calib <- y - yhat_calib
          RSS.gw <- sum(residual_calib[!na_yhat]^2)
          
          # Calculate effective degrees of freedom
          edf <- fd.n - 2 * tr.S + tr.StS
          if (edf <= 0) {
            warning("Effective degrees of freedom non-positive (", edf, 
                   "). Variance calculation may be unreliable.")
            edf <- max(1, edf)  # Force minimum of 1 to proceed
          }
          
          # Estimate sigma_hat^2 (model error variance)
          sigma.hat <- RSS.gw / edf
          
          if (!is.finite(sigma.hat) || sigma.hat < 0) {
            warning("Estimated error variance is ", 
                   ifelse(!is.finite(sigma.hat), "non-finite", "negative"), 
                   " (", sigma.hat, "). Skipping variance calculation.")
            estimation_ok <- FALSE
          } else {
            cat("  Estimated sigma.hat (model error variance):", sigma.hat, "\n")
            
            # Calculate variance for each prediction point
            cat("  Calculating variance for each prediction point...\n")
            pb <- txtProgressBar(min = 0, max = pd.n, style = 3)
            
            for (i in 1:pd.n) {
              setTxtProgressBar(pb, i)
              if (!valid_betas[i]) next  # Skip if coefficients were NA
              
              # Use stored xtxinv and weights wt[,i]
              w2 <- wt[, i] * wt[, i]
              w2x <- sweep(x, 1, w2, "*")  # Weight rows of x by w2
              xtw2x <- t(x) %*% w2x
              
              # Extract the stored (X'W_iX)^-1 for prediction point i
              xtxinvp <- xtxinv[i, , ]
              if (any(is.na(xtxinvp))) next
              
              # Calculate s0 = (X'W_iX)^-1 (X'W_i^2 X) (X'W_iX)^-1
              s0 <- try(xtxinvp %*% xtw2x %*% xtxinvp, silent = TRUE)
              if (inherits(s0, "try-error") || any(!is.finite(s0))) {
                next
              }
              
              # Calculate s1 = x_pred_i * s0 * x_pred_i'
              x.pi <- matrix(x.p[i, ], nrow = 1)
              s1 <- try(x.pi %*% s0 %*% t(x.pi), silent = TRUE)
              if (inherits(s1, "try-error") || !is.finite(s1[1,1])) {
                next
              }
              
              # Prediction variance
              s1_value <- as.numeric(s1)
              pse_squared <- sigma.hat * (1 + s1_value)
              
              if (is.finite(pse_squared) && pse_squared >= 0) {
                predict.var[i] <- pse_squared
              }
            }
            close(pb)
          }
        }
      }
    }
  } else {
    cat("Skipping prediction variance calculation (requires st.dMat2).\n")
  }
  
  #=============================================================================
  # Calculate evaluation metrics (for calibration points prediction)
  #=============================================================================
  
  metrics <- list()
  if (predicting_at_calib_points) {
    cat("Calculating model evaluation metrics...\n")
    valid_indices <- !is.na(gw.predict) & !is.na(original_values)
    
    if (sum(valid_indices) > 0) {
      residuals <- original_values[valid_indices] - gw.predict[valid_indices]
      
      # Calculate RMSE (Root Mean Square Error)
      metrics$rmse <- sqrt(mean(residuals^2))
      
      # Calculate MAE (Mean Absolute Error)
      metrics$mae <- mean(abs(residuals))
      
      # Calculate MAPE (Mean Absolute Percentage Error)
      if (all(original_values[valid_indices] != 0)) {
        metrics$mape <- mean(abs(residuals/original_values[valid_indices])) * 100
      } else {
        metrics$mape <- NA
        warning("MAPE calculation skipped due to zero values in original data")
      }
      
      # Calculate R-squared
      ss_total <- sum((original_values[valid_indices] - mean(original_values[valid_indices]))^2)
      ss_residual <- sum(residuals^2)
      metrics$r_squared <- 1 - (ss_residual / ss_total)
      
      # Calculate adjusted R-squared
      n_valid <- sum(valid_indices)
      p_vars <- ncol(x)
      metrics$adj_r_squared <- 1 - (1 - metrics$r_squared) * ((n_valid - 1) / (n_valid - p_vars - 1))
      
      cat("  RMSE:", metrics$rmse, "\n")
      cat("  MAE:", metrics$mae, "\n")
      if (!is.na(metrics$mape)) cat("  MAPE:", metrics$mape, "%\n")
      cat("  R²:", metrics$r_squared, "\n")
      cat("  Adjusted R²:", metrics$adj_r_squared, "\n")
    } else {
      warning("Cannot calculate metrics - no valid predictions or original values")
    }
  }
  
  #=============================================================================
  # Package and return results
  #=============================================================================
  
  cat("Packaging results...\n")
  
  # Create enhanced output dataframe with timestamps
  result_df <- data.frame(matrix(nrow = pd.n, ncol = 0))
  
  # Add coefficients
  for (i in 1:var.n) {
    col_name <- paste0(colnames(betas1)[i], "_coef")
    result_df[[col_name]] <- betas1[, i]
  }
  
  # Add predictions
  result_df[["prediction"]] <- gw.predict
  
  # Add prediction variance if calculated
  if (calculate_variance && !is.na(sigma.hat)) {
    result_df[["prediction_var"]] <- predict.var
    result_df[["prediction_se"]] <- sqrt(predict.var)
  }
  
  # Add original values for validation when predicting at calibration points
  if (predicting_at_calib_points) {
    result_df[["original_value"]] <- original_values
    result_df[["residual"]] <- original_values - gw.predict
  }
  
  # Add timestamps for temporal reference
  # Convert different time formats to character for consistency
  time_column <- "time_value"
  if (inherits(pred.tv, "Date")) {
    result_df[[time_column]] <- as.character(pred.tv)
  } else if (inherits(pred.tv, "POSIXct") || inherits(pred.tv, "POSIXlt")) {
    result_df[[time_column]] <- format(pred.tv, "%Y-%m-%d %H:%M:%S")
  } else {
    # For numeric or other formats, convert to character
    result_df[[time_column]] <- as.character(pred.tv)
  }
  
  # Create spatial object with predictions
  if (inherits(predict.SPDF, "Spatial")) {
    # Create SpatialPointsDataFrame or SpatialPolygonsDataFrame
    if (inherits(predict.SPDF, "SpatialPolygonsDataFrame")) {
      SDF <- SpatialPolygonsDataFrame(
        Sr = polygons(predict.SPDF),
        data = result_df,
        match.ID = FALSE
      )
    } else {
      SDF <- SpatialPointsDataFrame(
        coords = pd.locat,
        data = result_df,
        proj4string = CRS(as.character(p4s)),
        match.ID = FALSE
      )
      
      # Preserve gridded structure if input was gridded
      if (inherits(predict.SPDF, "SpatialPixelsDataFrame") || 
         inherits(predict.SPDF, "SpatialGridDataFrame")) {
        gridded(SDF) <- TRUE
      }
    }
  } else if (inherits(predict.SPDF, "sf")) {
    # Create sf object
    SDF <- st_sf(result_df, geometry = st_geometry(predict.SPDF))
  } else {
    # Just return the data frame if not spatial
    SDF <- result_df
  }
  
  # Record end time
  timings[["stop"]] <- Sys.time()
  
  # Store arguments for reference
  GTW.arguments <- list(
    formula = formula, 
    st.bw = st.bw, 
    kernel = kernel, 
    adaptive = adaptive,
    p = p, 
    theta = theta, 
    longlat = longlat, 
    lamda = lamda, 
    ksi = ksi,
    t.units = t.units, 
    fd.n = fd.n, 
    pd.n = pd.n
  )
  
  # Create result object with enhanced properties
  res <- list(
    GTW.arguments = GTW.arguments, 
    SDF = SDF, 
    metrics = metrics,
    predicting_at_calib_points = predicting_at_calib_points,
    sigma.hat = sigma.hat,
    timings = timings, 
    this.call = this.call,
    valid_predictions = sum(!is.na(gw.predict)),
    has_variance = calculate_variance && !is.na(sigma.hat)
  )
  
  class(res) <- "gtwrm.pred"
  
  # Calculate and show execution time
  execution_time <- as.numeric(difftime(timings[["stop"]], timings[["start"]], units = "secs"))
  cat("\nPrediction completed in", round(execution_time, 2), "seconds\n")
  cat("Successfully predicted", sum(!is.na(gw.predict)), "out of", pd.n, "locations\n")
  
  invisible(res)
}

#===============================================================================
# Enhanced print method for GTWR prediction results
#===============================================================================

print.gtwrm.pred <- function(x, ...) {
  if (!inherits(x, "gtwrm.pred")) 
    stop("Object is not a gtwr prediction result ('gtwrm.pred')")
  
  # Header information
  cat("===========================================================================\n")
  cat(" Geographically and Temporally Weighted Regression Prediction Results\n")
  cat("===========================================================================\n")
  cat("Start time:", as.character(x$timings$start), "\n")
  cat("End time:  ", as.character(x$timings$stop), "\n")
  cat("Execution time:", round(difftime(x$timings$stop, x$timings$start, units="secs"), 2), "seconds\n\n")
  
  # Call information  
  cat("Call:\n")
  print(x$this.call)
  
  # Formula information
  vars <- all.vars(x$GTW.arguments$formula)
  cat("\nDependent variable predicted:", vars[1])
  cat("\nIndependent variables:", paste(vars[-1], collapse=", "))
  
  # Data summary
  cat("\n\nData summary:")
  cat("\n  Calibration points:", x$GTW.arguments$fd.n)
  cat("\n  Prediction points:", x$GTW.arguments$pd.n)
  cat("\n  Valid predictions:", x$valid_predictions, 
      sprintf("(%.1f%%)", 100*x$valid_predictions/x$GTW.arguments$pd.n))
  
  # Model information
  cat("\n\nModel parameters:")
  cat("\n  Kernel function:", x$GTW.arguments$kernel)
  
  if (x$GTW.arguments$adaptive) {
    cat("\n  Adaptive ST bandwidth:", x$GTW.arguments$st.bw, "(number of nearest neighbours)")
  } else {
    cat("\n  Fixed ST bandwidth:", x$GTW.arguments$st.bw)
  }
  
  cat("\n  Spatial distance: ", ifelse(x$GTW.arguments$longlat, 
                                     "Great Circle", 
                                     paste0("Minkowski p=", x$GTW.arguments$p)))
  cat("\n  Temporal units:", x$GTW.arguments$t.units)
  cat("\n  ST combination: lambda =", x$GTW.arguments$lamda, ", ksi =", x$GTW.arguments$ksi)
  
  # Display evaluation metrics if available
  if (length(x$metrics) > 0) {
    cat("\n\nModel evaluation metrics (predicting at calibration points):\n")
    if (!is.null(x$metrics$rmse)) cat("  RMSE:", round(x$metrics$rmse, 4), "\n")
    if (!is.null(x$metrics$mae)) cat("  MAE:", round(x$metrics$mae, 4), "\n")
    if (!is.null(x$metrics$mape) && !is.na(x$metrics$mape)) 
      cat("  MAPE:", round(x$metrics$mape, 2), "%\n")
    if (!is.null(x$metrics$r_squared)) 
      cat("  R²:", round(x$metrics$r_squared, 4), "\n")
    if (!is.null(x$metrics$adj_r_squared)) 
      cat("  Adjusted R²:", round(x$metrics$adj_r_squared, 4), "\n")
  }
  
  # Prepare summary statistics
  if (inherits(x$SDF, "Spatial")) {
    SDF.df <- as(x$SDF, "data.frame")
  } else if (inherits(x$SDF, "sf")) {
    SDF.df <- st_drop_geometry(x$SDF)
  } else {
    SDF.df <- x$SDF
  }
  
  # Identify coefficient columns
  coef_cols <- grep("_coef$", colnames(SDF.df), value = TRUE)
  
  # Print coefficient summary if available
  if (length(coef_cols) > 0) {
    cat("\n\nSummary of coefficient estimates at prediction points:\n")
    
    df0 <- SDF.df[, coef_cols, drop = FALSE]
    CM <- t(apply(df0, 2, summary, na.rm = TRUE))
    
    # Format the coefficient summary with proper alignment
    cm_formatted <- format(round(CM, 4), digits = 4, nsmall = 4)
    rnames <- rownames(cm_formatted)
    
    # Print with proper indentation
    for (i in 1:nrow(cm_formatted)) {
      cat("  ", rnames[i], ":", sep="")
      cat(paste(colnames(CM), "=", cm_formatted[i,]), sep="  ", collapse="\n                 ")
      cat("\n")
    }
  }
  
  # Print prediction summary
  cat("\nSummary of predictions:\n")
  if ("prediction" %in% names(SDF.df)) {
    pred_summary <- summary(SDF.df$prediction, na.rm = TRUE)
    cat("  Min:", format(round(pred_summary[1], 4), nsmall=4), 
        "  1st Qu:", format(round(pred_summary[2], 4), nsmall=4), "\n",
        "  Median:", format(round(pred_summary[3], 4), nsmall=4), 
        "  Mean:", format(round(pred_summary[4], 4), nsmall=4), "\n",
        "  3rd Qu:", format(round(pred_summary[5], 4), nsmall=4), 
        "  Max:", format(round(pred_summary[6], 4), nsmall=4), "\n", sep="")
    
    if (sum(is.na(SDF.df$prediction)) > 0) {
      cat("  NA count:", sum(is.na(SDF.df$prediction)), "\n")
    }
  } else {
    cat("  No prediction results found.\n")
  }
  
  # Print time value summary if available
  if ("time_value" %in% names(SDF.df)) {
    unique_times <- unique(SDF.df$time_value)
    cat("\nTime periods in prediction results:", length(unique_times), "\n")
    if (length(unique_times) <= 10) {
      cat("  Time values:", paste(unique_times, collapse=", "), "\n")
    } else {
      cat("  Time range:", min(unique_times), "to", max(unique_times), 
          "(showing 2 of", length(unique_times), "values)\n")
    }
    
    # Calculate predictions by time period
    cat("\nPredictions by time period:\n")
    time_groups <- split(SDF.df$prediction, SDF.df$time_value)
    time_stats <- data.frame(
      Time = names(time_groups),
      Count = sapply(time_groups, length),
      Mean = sapply(time_groups, mean, na.rm = TRUE),
      Min = sapply(time_groups, min, na.rm = TRUE),
      Max = sapply(time_groups, max, na.rm = TRUE)
    )
    
    if (nrow(time_stats) <= 10) {
      # Print full summary for few time periods
      print(time_stats, row.names = FALSE)
    } else {
      # Print abbreviated summary for many time periods
      print(time_stats[1:5,], row.names = FALSE)
      cat("  ... and", nrow(time_stats) - 5, "more time periods\n")
    }
  } else {
    cat("\nNo time information available in prediction results.\n")
  }
  
  # Print residual summary if original values are available
  if (all(c("original_value", "residual") %in% names(SDF.df))) {
    cat("\nResidual summary (predictions vs. original values):\n")
    res_summary <- summary(SDF.df$residual, na.rm = TRUE)
    cat("  Min:", format(round(res_summary[1], 4), nsmall=4), 
        "  1st Qu:", format(round(res_summary[2], 4), nsmall=4), "\n",
        "  Median:", format(round(res_summary[3], 4), nsmall=4), 
        "  Mean:", format(round(res_summary[4], 4), nsmall=4), "\n",
        "  3rd Qu:", format(round(res_summary[5], 4), nsmall=4), 
        "  Max:", format(round(res_summary[6], 4), nsmall=4), "\n", sep="")
    
    # Calculate absolute percentage error
    valid_idx <- !is.na(SDF.df$original_value) & !is.na(SDF.df$prediction) & 
                 SDF.df$original_value != 0
    if (sum(valid_idx) > 0) {
      ape <- abs(SDF.df$residual[valid_idx] / SDF.df$original_value[valid_idx]) * 100
      cat("  Absolute Percentage Error (APE):\n")
      ape_summary <- summary(ape)
      cat("    Min:", format(round(ape_summary[1], 2), nsmall=2), "% ", 
          "  Median:", format(round(ape_summary[3], 2), nsmall=2), "% ", 
          "  Mean:", format(round(mean(ape), 2), nsmall=2), "% ", 
          "  Max:", format(round(ape_summary[6], 2), nsmall=2), "%\n", sep="")
    }
  }
  
  # Print variance summary if available
  if ("prediction_var" %in% names(SDF.df) && x$has_variance) {
    cat("\nPrediction variance information:\n")
    cat("  Estimated sigma.hat (model error variance):", format(round(x$sigma.hat, 6), nsmall=6), "\n")
    
    var_summary <- summary(sqrt(SDF.df$prediction_var), na.rm = TRUE)
    cat("  Prediction std.error summary:\n")
    cat("  Min:", format(round(var_summary[1], 4), nsmall=4), 
        "  1st Qu:", format(round(var_summary[2], 4), nsmall=4), "\n",
        "  Median:", format(round(var_summary[3], 4), nsmall=4), 
        "  Mean:", format(round(var_summary[4], 4), nsmall=4), "\n",
        "  3rd Qu:", format(round(var_summary[5], 4), nsmall=4), 
        "  Max:", format(round(var_summary[6], 4), nsmall=4), "\n", sep="")
    
    if (sum(is.na(SDF.df$prediction_var)) > 0) {
      cat("  NA count:", sum(is.na(SDF.df$prediction_var)), "\n")
    }
    
    cat("\n  Note: Variance estimates adapted from GWR methodology.\n")
    cat("        Interpret with caution as GTWR variance theory differs.\n")
  } else if ("prediction_var" %in% names(SDF.df) && !x$has_variance) {
    cat("\nPrediction variance calculation was attempted but failed.\n")
  } else {
    cat("\nPrediction variance was not calculated.\n")
  }
  
  cat("\n===========================================================================\n")
  
  invisible(x)
}