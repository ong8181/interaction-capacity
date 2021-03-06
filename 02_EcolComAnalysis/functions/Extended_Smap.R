####
#### Spatial S-map to incorporate spatial dependence into EDM
####

#---------- Extended smap ----------#
extended_smap <- function(vectors,
                          target,
                          lib_indices,
                          pred_indices,
                          theta,
                          regularized = FALSE,
                          lambda = 0.1,
                          alpha = 0, # default is the ridge regression. If alpha = 1, then do lasso regression
                          glmnet_parallel = FALSE, # require(doParallel)
                          save_smap_coefficients = FALSE)
{
  #require(glmnet)
  
  # set E here
  E <- NCOL(vectors)
  
  # setup output
  pred_vals <- rep.int(NaN, times = length(target))
  smap_coefficient_vals <- matrix(NaN, ncol = NCOL(vectors) + 1, nrow = length(target))
  
  # exclude libs that contains NaN in target
  lib_nona <- !apply(cbind(vectors, target), 1, function(x) any(is.na(x)))
  lib_indices <- lib_nona & lib_indices
  
  # Add NaN if vectors contains NaN #<--- should be deleted?
  pred_vals[apply(vectors, 1, function(x) any(is.na(x))) & pred_indices] <- NaN
  smap_coefficient_vals[apply(vectors, 1, function(x) any(is.na(x))) & pred_indices] <- NaN
  
  # Make new pred_indices (exclude indices that contains NaN)
  pred_indices <- !apply(vectors, 1, function(x) any(is.na(x))) & pred_indices
  
  #-------------------- Main loop to make predictions --------------------#
  for(p in which(pred_indices))
  {
    temp_lib <- lib_indices[p]
    lib_indices[p] <- FALSE
    libs <- which(lib_indices)
    
    # compute distances of temporal information
    q <- matrix(rep(vectors[p,], length(libs)), nrow = length(libs), byrow = T)
    distances <- sqrt(rowSums((vectors[libs,] - q)^2))
    
    # compute temporal weights
    d_bar <- mean(distances, na.rm = TRUE)
    c_ws <- exp(- theta * distances / d_bar)
    
    if(regularized){
      # do regularized S-map
      # Currently rely on "glmnet" package of R
      if(E == 1){
        A <- cbind(vectors[libs,], 1)# * c_ws
        B <- cbind(target[libs])# * c_ws
        if(is.null(lambda)){
          # make prediction
          fit <- cv.glmnet(A, B, weights = c_ws, type.measure = "mae", alpha = alpha, family = "gaussian", nfolds = 10,
                           parallel = glmnet_parallel, intercept = TRUE)
          pred_vals[p] <- predict(fit, s = fit$lambda.1se, newx = matrix(c(vectors[p,], 1), nrow = 1))
          smap_coefficient_vals[p,] <- matrix(t(coef(fit, s = fit$lambda.1se)), nrow = 1)[c(2,1)]
        }else{
          # make prediction
          fit <- glmnet(A, B, weights = c_ws, alpha = alpha, family = "gaussian", lambda = lambda, intercept = TRUE)
          pred_vals[p] <- predict(fit, s = fit$lambda, newx = matrix(c(vectors[p,], 1), nrow = 1))
          smap_coefficient_vals[p,] <- matrix(t(coef(fit, s = fit$lambda)), nrow = 1)[c(2,1)]
        }
      }else{
        A <- cbind(vectors[libs,])# * c_ws
        B <- cbind(target[libs])# * c_ws
        if(is.null(lambda)){
          # make prediction
          fit <- cv.glmnet(A, B, weights = c_ws, type.measure = "mae", alpha = alpha, family = "gaussian", nfolds = 10,
                           parallel = glmnet_parallel, intercept = TRUE)
          pred_vals[p] <- predict(fit, s = fit$lambda.1se, newx = matrix(c(vectors[p,]), nrow = 1))
          smap_coefficient_vals[p,] <- matrix(t(coef(fit, s = fit$lambda.1se)), nrow = 1)[c(2:(NCOL(vectors)+1),1)]
        }else{
          # make prediction
          fit <- glmnet(A, B, weights = c_ws, alpha = alpha, family = "gaussian", lambda = lambda, intercept = TRUE)
          pred_vals[p] <- predict(fit, s = fit$lambda, newx = matrix(c(vectors[p,]), nrow = 1))
          smap_coefficient_vals[p,] <- matrix(t(coef(fit, s = fit$lambda)), nrow = 1)[c(2:(NCOL(vectors)+1),1)]
        }
      }
      
      lib_indices[p] <- temp_lib 
    }else{
      
      # do singular-value decomposition
      A <- cbind(vectors[lib_indices,], 1) * c_ws
      A_svd <- svd(A)
      
      # remove singular values that are too small
      s <- A_svd$d
      s_inv <- matrix(0, nrow = E+1, ncol = E+1)
      for(i in seq_along(s))
      {
        if(s[i] >= max(s) * 1e-5)
          s_inv[i,i] <- 1/s[i]
      }
      
      # perform back-substitute to solve        
      map <- A_svd$v %*% s_inv %*% t(A_svd$u) %*% (c_ws * target[lib_indices])
      
      # make prediction
      pred_vals[p] <- sum(map * c(vectors[p,], 1))
      smap_coefficient_vals[p,] <- t(map)
      
      lib_indices[p] <- temp_lib
    }
  }
  #-------------------- Main loop to finished --------------------#
  
  # Add colnames for smap_coefficient_vals
  colnames(smap_coefficient_vals) <- c(paste0("c_", 1:NCOL(vectors)), "c_0")
  
  # return output & stats
  if(save_smap_coefficients){
    return(list(pred = pred_vals, stats = compute_stats_SSR(target[pred_indices], pred_vals[pred_indices]),
                smap_coefficients = smap_coefficient_vals))
  }else{
    return(list(pred = pred_vals, stats = compute_stats_SSR(target[pred_indices], pred_vals[pred_indices])))
  }
}

