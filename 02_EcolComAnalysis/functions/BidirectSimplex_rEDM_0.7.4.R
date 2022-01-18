####
#### Function: Bidirectional simplex projection
####


#require(rEDM)

# Define function
bidirect_simplex <- function(time_series, E,
                             complete_case_only = FALSE, # if TRUE, averaged predictions will be calculated only for complete pairs
                             lib = c(1, length(time_series)),
                             pred = lib,
                             tau = 1,
                             tp_forward = 1,
                             tp_backward = -E * tau - tp_forward + 1,
                             num_neighbors = "e+1",
                             stats_only = TRUE,
                             silent = TRUE,
                             exclusion_radius = NULL,
                             epsilon = NULL){
  # Perform forward and backward simplex projection
  simplex_forward <- simplex(time_series, E = E, tp = tp_forward,
                             lib = lib, pred = pred,
                             tau = tau, num_neighbors = num_neighbors,
                             stats_only = F, silent = silent)
  simplex_backward <- simplex(time_series, E = E, tp = tp_backward,
                              lib = lib, pred = pred,
                              tau = tau, num_neighbors = num_neighbors,
                              stats_only = F, silent = silent)
  
  if(complete_case_only){
    # Match time indices of the model output
    time_shared <- as.numeric(na.omit(intersect(simplex_forward$model_output[[1]]$time,
                                                simplex_backward$model_output[[1]]$time)))
    # Extract model output
    time_forward <- simplex_forward$model_output[[1]][match(time_shared, simplex_forward$model_output[[1]]$time),]
    time_backward <- simplex_backward$model_output[[1]][match(time_shared, simplex_backward$model_output[[1]]$time),]
    
    # Averaging forward and backward predictions
    averaged_prediction <- (time_forward$pred + time_backward$pred)/2
    
    simplex_bidirect <- data.frame(time = time_forward$time,
                                   obs = time_forward$obs,
                                   pred = averaged_prediction)
  }else{
    # Match time indices of the model output
    simplex_forward_nona <- simplex_forward$model_output[[1]][!is.na(simplex_forward$model_output[[1]]$time),]
    simplex_backward_nona <- simplex_backward$model_output[[1]][!is.na(simplex_backward$model_output[[1]]$time),]
    time_covered <- sort(unique(c(simplex_forward_nona$time, simplex_backward_nona$time)))
    
    # Identify time index of simplex_forward and simplex_backward
    time_forward_id <- match(simplex_forward_nona$time, time_covered)
    time_backward_id <- match(simplex_backward_nona$time, time_covered)
    
    # Create simplex_bidirect object
    simplex_bidirect0 <- data.frame(time = time_covered)
    simplex_bidirect0$forward_pred <- simplex_bidirect0$forward_obs <- simplex_bidirect0$forward_time <- NaN
    simplex_bidirect0$backward_pred <- simplex_bidirect0$backward_obs <- simplex_bidirect0$backward_time <- NaN
    simplex_bidirect0[time_forward_id,]$forward_time <- simplex_forward_nona$time
    simplex_bidirect0[time_forward_id,]$forward_obs <- simplex_forward_nona$obs
    simplex_bidirect0[time_forward_id,]$forward_pred <- simplex_forward_nona$pred
    simplex_bidirect0[time_backward_id,]$backward_time <- simplex_backward_nona$time
    simplex_bidirect0[time_backward_id,]$backward_obs <- simplex_backward_nona$obs
    simplex_bidirect0[time_backward_id,]$backward_pred <- simplex_backward_nona$pred
    
    # Sanity check
    if(F){
      sum(simplex_bidirect0$forward_time - simplex_bidirect0$backward_time, na.rm = T) == 0
      sum(simplex_bidirect0$time - simplex_bidirect0$backward_time, na.rm = T) == 0
      sum(simplex_bidirect0$forward_obs - simplex_bidirect0$backward_obs, na.rm = T) == 0
    }

    # Calculate average prediction
    obs_combined <- apply(cbind(simplex_bidirect0$forward_obs, simplex_bidirect0$backward_obs), 1, function(x) mean(x, na.rm = T))
    averaged_prediction <- apply(cbind(simplex_bidirect0$forward_pred, simplex_bidirect0$backward_pred), 1, function(x) mean(x, na.rm = T))
    
    # Summary
    simplex_bidirect <- data.frame(time = simplex_bidirect0$time,
                                   obs = obs_combined,
                                   pred = averaged_prediction)
  }
  
  # Collecting parameters
  params_forward <- simplex_forward[1:16]
  params_backward <- simplex_backward[1:16]
  params_all <- cbind(data.frame(method = c("forward_simplex","backward_simplex")),
                     rbind(params_forward, params_backward))
  
  # Calculate prediction accuracy
  stats_forward <- simplex_forward[,c("num_pred", "E", "tp", "rho", "mae", "rmse")]
  stats_backward <- simplex_backward[,c("num_pred", "E", "tp", "rho", "mae", "rmse")]
  stats_bidirection <- compute_stats2(simplex_bidirect$obs, simplex_bidirect$pred)
  stats_bidirection$E <- simplex_forward$E
  stats_bidirection$tp <- sprintf("%s, %s", simplex_forward$tp, simplex_backward$tp)
  stats_bidirection <- stats_bidirection[,c("num_pred", "E", "tp", "rho", "mae", "rmse")]
  
  stats_all <- cbind(data.frame(method = c("bidirect_simplex","forward_simplex","backward_simplex"),
                     rbind(stats_bidirection, stats_forward, stats_backward)))

  # Return output
  if(stats_only){
    simplex_bidirect_all <- data.frame(stats_all)
  }else{
    simplex_bidirect_all <- list(params = params_all, model_output = simplex_bidirect, stats = stats_all)
  }
  return(simplex_bidirect_all)
}

compute_stats2 <- function(obs, pred){
  # computes performance metrics for how well predictions match observations
  # obs = vector of observations
  # pred = vector of prediction
  
  num_pred <- sum(is.finite(obs) & is.finite(pred))
  rho <- cor(obs, pred, use = "pairwise.complete.obs")
  mae <- mean(abs(obs-pred), na.rm = TRUE)
  rmse <- sqrt(mean((obs-pred)^2, na.rm = TRUE))
  return(data.frame(num_pred = num_pred, rho = rho, mae = mae, rmse = rmse))
}

# Calculate best E using bidirectional simplex projection
bestE_bidirect <- function(time_series,
                           E_range,
                           lib = c(1, length(time_series)),
                           complete_case_only = F,
                           criteria = "rmse",
                           show_fig = F,
                           save_stats = F){
  simp_res <- bidirect_simplex(time_series, E = 1, lib = lib, stats_only = T, complete_case_only = complete_case_only)[1,]
  for(i in 2:length(E_range)){
    simp_res[i,] <- bidirect_simplex(time_series, E = E_range[i],
                                     lib = lib,
                                     stats_only = T,
                                     complete_case_only = complete_case_only)[1,]
  }
  best_E <- simp_res[which.min(simp_res[,criteria]), "E"]
  
  if(show_fig) plot(simp_res$E, simp_res[,criteria], type = "b", ylab = criteria, xlab = "E") 
  if(save_stats){
    all_res <- list(E = best_E, stats = simp_res)
  }else{
    all_res <- best_E
  }
  return(all_res)
}