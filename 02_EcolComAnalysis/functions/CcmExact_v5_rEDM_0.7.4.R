#----------------------------------------------------------------------------------------------------#
# Exact CCM function
#----------------------------------------------------------------------------------------------------#

ccm_exact = function (
  block, lib = c(1, NROW(block)), pred = lib, norm = 2, E = 1, tau = 1, tp = 0,
  num_neighbors = "e+1", lib_column = 1, target_column = 2, first_column_time = FALSE,
  RNGseed = NULL, exclusion_radius = NULL, epsilon = NULL)
{
  require(rEDM)
  lib  = rbind(lib)
  pred = rbind(pred)
  lib_size = sum(lib[,2] - lib[,1] + 1)
  
  ccm_wrap = function (block, E, tp, nn = "e+1") {
    ccm(
      block, lib = lib, pred = pred, norm = norm, E = E, tau = tau, tp = tp,
      num_neighbors = nn, lib_sizes = lib_size, num_samples = 1, replace = FALSE,
      lib_column = lib_column, target_column = target_column, first_column_time = first_column_time,
      RNGseed = RNGseed, exclusion_radius = exclusion_radius, epsilon = epsilon,
      stats_only = TRUE, silent = TRUE)    
  }
  
  xmapH = ccm_wrap(block, E = E, tp = tp)
  pred_id = unique(unlist(lapply(1:nrow(pred), function(i) pred[i,1]:pred[i,2])))
  err = abs(block[pred_id,target_column] - mean(block[pred_id,target_column], na.rm = TRUE))
  nnH = xmapH$nn
  xmap0 = data.frame(
    E = E, tau = tau, tp = 0, nn = nnH, lib_column = lib_column, target_column = target_column,
    lib_size = 0, num_pred = sum(!is.na(err)), rho = 0,
    mae = (1+1/nnH)*mean(err, na.rm = TRUE), rmse = sqrt((1+1/nnH)*mean(err^2, na.rm = TRUE)))
  if(E == 1)
    xmapL = xmap0
  else {
    xmapL = ccm_wrap(block, E = E - 1, tp = tp + tau, nn = xmapH$nn)
    xmapL$E   = xmapH$E
    xmapL$tau = xmapH$tau
    xmapL$lib_size = 0
  }
  op = rbind(xmap0 = xmap0, xmapL = xmapL, xmapH = xmapH)
  op$r2 = 1 - op$rmse^2 / op$rmse[1]^2
  op
}

#----------------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------------#
# Exact CCM nonaive function
#----------------------------------------------------------------------------------------------------#
ccm_exact_nonaive = function (
  block, lib = c(1, NROW(block)), pred = lib, norm = 2, E = 1, tau = 1, tp = 0,
  num_neighbors = "e+1", lib_column = 1, target_column = 2, first_column_time = FALSE,
  RNGseed = NULL, exclusion_radius = NULL, epsilon = NULL)
{
  require(rEDM)
  lib  = rbind(lib)
  pred = rbind(pred)
  lib_size = sum(lib[,2] - lib[,1] + 1)
  
  ccm_wrap = function (block, E, tp, nn = "e+1") {
    ccm(
      block, lib = lib, pred = pred, norm = norm, E = E, tau = tau, tp = tp,
      num_neighbors = nn, lib_sizes = lib_size, num_samples = 1, replace = FALSE,
      lib_column = lib_column, target_column = target_column, first_column_time = first_column_time,
      RNGseed = RNGseed, exclusion_radius = exclusion_radius, epsilon = epsilon,
      stats_only = TRUE, silent = TRUE)    
  }
  
  xmapH = ccm_wrap(block, E = E, tp = tp)
  pred_id = unique(unlist(lapply(1:nrow(pred), function(i) pred[i,1]:pred[i,2])))
  err = abs(block[pred_id,target_column] - mean(block[pred_id,target_column], na.rm = TRUE))
  nnH = xmapH$nn
  if(E == 1){
    xmapL = data.frame(
      E = E, tau = tau, tp = 0, nn = nnH, lib_column = lib_column, target_column = target_column,
      lib_size = 0, num_pred = sum(!is.na(err)), rho = 0,
      mae = (1+1/nnH)*mean(err, na.rm = TRUE), rmse = sqrt((1+1/nnH)*mean(err^2, na.rm = TRUE)))
    #xmapL = xmap0
  } else {
    xmapL = ccm_wrap(block, E = E - 1, tp = tp + tau, nn = xmapH$nn)
    xmapL$E   = xmapH$E
    xmapL$tau = xmapH$tau
    xmapL$lib_size = 0
  }
  op = rbind(xmapL = xmapL, xmapH = xmapH)
  op$r2 = 1 - op$rmse^2 / op$rmse[1]^2
  op
}

#----------------------------------------------------------------------------------------------------#
