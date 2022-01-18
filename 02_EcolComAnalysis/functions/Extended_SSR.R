####
#### Extended SSR to incorporate spatial dependence into EDM
####

#---------- Main extended lnlp function ----------#
extended_lnlp <- function(block_time,
                          lib = c(1, NROW(block_time)),
                          pred = lib,
                          tp = 1,
                          target_column = 1,
                          lib_column = 1:NCOL(block_time),
                          num_neighbors = NCOL(block_time) + 1,
                          theta = 0,
                          method = "simplex",
                          regularized = FALSE,
                          lambda = NULL,
                          alpha = 0, # default is the ridge regression. If alpha = 1, then do lasso regression
                          glmnet_parallel = FALSE,
                          save_smap_coefficients = FALSE)
{
  # do multivariate prediction using simplex projection or s-map
  # theta = relative weighting of neighbors based on Euclidean distance in a state space

  if(!is.matrix(block_time)) block_time <- as.matrix(block_time)
  n <- NROW(block_time)
  lib <- matrix(lib, ncol = 2)
  pred <- matrix(pred, ncol = 2)
  
  # setup vectors
  vectors <- matrix(block_time[,lib_column], ncol = length(lib_column))
  
  # setup lib_indices
  lib_indices <- rep.int(FALSE, times = n)
  for(i in 1:NROW(lib))
  {
    row_start <- lib[i, 1]
    row_end <- lib[i, 2] - tp
    if(row_end > row_start)
      lib_indices[row_start:row_end] <- TRUE
  }
  
  # setup pred_indices
  pred_indices <- rep.int(FALSE, times = n)
  for(i in 1:NROW(pred))
  {
    row_start <- pred[i, 1]
    row_end <- pred[i, 2] - tp
    if(row_end > row_start)
      pred_indices[row_start:row_end] <- TRUE
  }
  
  # setup target
  target <- rep.int(NaN, times = n)
  target[1:(n-tp)] <- block_time[(1+tp):n, target_column]
  
  if(method == "simplex")
  {
    # perform simplex projection and summarize results
    simplex_out <- extended_simplex(vectors, target, lib_indices, pred_indices, num_neighbors)
    pred_df <- data.frame(time = 1:nrow(vectors),
                          obs = c(rep(NaN, tp), head(target, n = n - tp)),
                          pred = c(rep(NaN, tp), head(simplex_out$pred, n = n - tp)))
    stats <- simplex_out$stats
    
    # Return results
    return(list(model_output = pred_df, stats = stats))
    
  }else if(method == "s-map" |method == "s_map" |method == "smap"){
    smap_out <- extended_smap(vectors, target,
                              lib_indices, pred_indices,
                              theta,
                              regularized = regularized,
                              lambda = lambda,
                              alpha = alpha,
                              glmnet_parallel = glmnet_parallel,
                              save_smap_coefficients = save_smap_coefficients)
    pred_df <- data.frame(time = 1:nrow(vectors),
                          obs = c(rep(NaN, tp), head(target, n = n - tp)),
                          pred = c(rep(NaN, tp), head(smap_out$pred, n = n - tp)))
    
    stats <- smap_out$stats
    
    # Return results
    if(save_smap_coefficients){
      smap_coef_df <- cbind(data.frame(time = 1:nrow(vectors)),
                            matrix(NaN, ncol = dim(smap_out$smap_coefficients)[2],
                                   nrow = dim(smap_out$smap_coefficients)[1]))
      smap_coef_df[,2:(dim(smap_out$smap_coefficients)[2]+1)] <- smap_out$smap_coefficients
      colnames(smap_coef_df)[2:(dim(smap_out$smap_coefficients)[2]+1)] <- colnames(smap_out$smap_coefficients)
      
      return(list(model_output = pred_df, stats = stats,
                  smap_coefficients = smap_coef_df))
    }else{
      return(list(model_output = pred_df, stats = stats))
    }
  }
}

#---------- Optional: extended_lnlp function for univariate time series ----------#
extended_univariate_lnlp <- function(block_time,
                                     E,
                                     lib = c(1, length(block_time)),
                                     pred = lib,
                                     tau = 1,
                                     tp = 1,
                                     num_neighbors = E + 1,
                                     theta = 0,
                                     method = "simplex",
                                     regularized = FALSE,
                                     lambda = NULL,
                                     alpha = 0, # default is the ridge regression. If alpha = 1, then do lasso regression
                                     glmnet_parallel = FALSE,
                                     save_smap_coefficients = FALSE)
{
  # do univariate prediction using simplex projection or s-map
  # theta = relative weighting of neighbors based on Euclidean distance in a state space

  #if(!is.matrix(block_time)) block_time <- as.matrix(block_time)
  n <- length(block_time)
  lib <- matrix(lib, ncol = 2)
  pred <- matrix(pred, ncol = 2)
  
  # setup vectors
  vectors <- matrix(NaN, nrow = n, ncol = E)
  lag <- 0
  for (i in 1:E){
    vectors[(lag+1):n,i] <- block_time[1:(n-lag)]
    lag <- lag + tau
  }
  
  # setup lib_indices
  lib_indices <- rep.int(FALSE, times = n)
  for(i in 1:NROW(lib))
  {
    row_start <- lib[i, 1] + (E-1)*tau
    row_end <- lib[i, 2] - tp
    if(row_end > row_start)
      lib_indices[row_start:row_end] <- TRUE
  }
  
  # setup pred_indices
  pred_indices <- rep.int(FALSE, times = n)
  for(i in 1:NROW(pred))
  {
    row_start <- pred[i, 1] + (E-1)*tau
    row_end <- pred[i, 2] - tp
    if(row_end > row_start)
      pred_indices[row_start:row_end] <- TRUE
  }
  
  # setup target
  target <- rep.int(NaN, times = n)
  target[1:(n-tp)] <- block_time[(1+tp):n]
  
  if(method == "simplex")
  {
    # perform simplex projection and summarize results
    simplex_out <- extended_simplex(vectors, target, lib_indices, pred_indices, num_neighbors)
    pred_df <- data.frame(time = 1:nrow(vectors),
                          obs = c(rep(NaN, tp), head(target, n = n - tp)),
                          pred = c(rep(NaN, tp), head(simplex_out$pred, n = n - tp)))
    stats <- simplex_out$stats
    
    # Return results
    return(list(model_output = pred_df, stats = stats))

  }else if(method == "s-map" |method == "s_map" |method == "smap"){
    smap_out <- extended_smap(vectors, target,
                              lib_indices, pred_indices,
                              theta,
                              regularized = regularized,
                              lambda = lambda,
                              alpha = alpha,
                              glmnet_parallel = glmnet_parallel,
                              save_smap_coefficients = save_smap_coefficients)
    pred_df <- data.frame(time = 1:nrow(vectors),
                          obs = c(rep(NaN, tp), head(target, n = n - tp)),
                          pred = c(rep(NaN, tp), head(smap_out$pred, n = n - tp)))
    stats <- smap_out$stats
    
    # Return results
    if(save_smap_coefficients){
      smap_coef_df <- cbind(data.frame(time = 1:nrow(vectors)),
                            matrix(NaN, ncol = dim(smap_out$smap_coefficients)[2],
                                   nrow = dim(smap_out$smap_coefficients)[1]))
      smap_coef_df[,2:(dim(smap_out$smap_coefficients)[2]+1)] <- smap_out$smap_coefficients
      colnames(smap_coef_df)[2:(dim(smap_out$smap_coefficients)[2]+1)] <- colnames(smap_out$smap_coefficients)
      
      return(list(model_output = pred_df, stats = stats, smap_coefficients = smap_coef_df))
    }else{
      return(list(model_output = pred_df, stats = stats))
    }
  }
  
}

#---------- Function to calcluate stats ----------#
compute_stats_SSR <- function(obs, pred)
{
  # computes performance metrics for how well predictions match observations
  # obs = vector of observations
  # pred = vector of prediction
  
  N = sum(is.finite(obs) & is.finite(pred))
  rho = cor(obs, pred, use = "pairwise.complete.obs")
  mae = mean(abs(obs-pred), na.rm = TRUE)
  rmse = sqrt(mean((obs-pred)^2, na.rm = TRUE))
  return(data.frame(N = N, rho = rho, mae = mae, rmse = rmse))
}
