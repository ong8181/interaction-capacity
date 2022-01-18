####
#### Calculate 95% CI for surrogate data
####
####

#require(pforeach)

ccm_exact_surrogate <- function(effect_ts, cause_ts, surrogate_ts,
                                surrogate = "effect",
                                E_range = 1:15, # only valid when surrogate == "effect"
                                E_fix = 2, # only valid when surrogate == "cause"
                                lib = c(1, length(effect_ts)),
                                pred = lib,
                                tp = 0,
                                simplex_criteria = "rmse",
                                simplex_method = "multivariate_simplex")
{
  #require(pforeach)
  # do CCM for the surrogate data
  if(surrogate == "effect"){
    # Main loop
    surrogate_all <- pforeach(i = 1:ncol(surrogate_ts), .c=rbind)({
      # Set time series
      effect_sur <- surrogate_ts[,i]
      block <- cbind(effect_sur, cause_ts)
      if(simplex_method == "unidirection"){
        E_effect <- bestE(effect_sur, E = E_range, criteria = simplex_criteria, lib = lib, show_fig = F)
      }else if(simplex_method == "bidirection"){
        E_effect <- bestE_bidirect(effect_sur, E = E_range, criteria = simplex_criteria, lib = lib, show_fig = F)
      }else if(simplex_method == "multivariate_simplex"){
        E_effect <- bestE_bidirect_block_lnlp(block, E_range = E_range, lib = lib, criteria = simplex_criteria, show_fig = F) + 1
      }
      
      # Do CCM
      ccm_exact_res <- ccm_exact_nonaive(block, E = E_effect, tp = tp, lib = lib, pred = pred)
      fore_skill <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")]["xmapH",]
      fore_improvement <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")]["xmapH",] - ccm_exact_res[,c("rho", "mae", "rmse", "r2")]["xmapL",]
      colnames(fore_improvement) <- c("d_rho", "d_mae", "d_rmse", "d_r2")
      
      # Summarize results
      # Naive predictions are excluded because they are not used to calculate p-values
      # and because to speed up the calculations
      fore_skill_all <- data.frame(cbind(fore_skill, fore_improvement))
    })
    
  }else if(surrogate == "cause"){
    # Main loop
    surrogate_all <- pforeach(i = 1:ncol(surrogate_ts), .c=rbind)({
      # Set time series
      cause_sur <-  surrogate_ts[,i]
      block <- cbind(effect_ts, cause_sur)

      # Do CCM
      ccm_exact_res <- ccm_exact_nonaive(block, E = E_fix, tp = tp, lib = lib, pred = pred)
      fore_skill <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")]["xmapH",]
      fore_improvement <- ccm_exact_res[,c("rho", "mae", "rmse", "r2")]["xmapH",] - ccm_exact_res[,c("rho", "mae", "rmse", "r2")]["xmapL",]
      colnames(fore_improvement) <- c("d_rho", "d_mae", "d_rmse", "d_r2")

      # Summarize results
      # Naive predictions are excluded because they are not used to calculate p-values
      # and because to speed up the calculations
      fore_skill_all <- data.frame(cbind(fore_skill, fore_improvement))
    })
  }
  
  # Add rownames
  rownames(surrogate_all) <- sprintf("sur_%05d", 1:ncol(surrogate_ts))

  return(surrogate_all)
}
 

# Function to summarize stats of surrogate CCM
summarize_sur <- function(surrogate_res,
                          quantile_vals = c(0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.1, 0.5,
                                            0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9995)){
  return(apply(surrogate_res, 2, function(x) quantile(x, quantile_vals, na.rm = TRUE)))
}


# Function to plot surrogate results
plot_exact_ccm <- function(original_result, surrogate_result, surrogate_stat, criteria = "rmse", plot_title = NA){
  xmax <- max(surrogate_result[,criteria], original_result[,criteria]) + 0.02
  xmin <- min(surrogate_result[,criteria], original_result[,criteria]) - 0.02
  ymax <- max(surrogate_result[,sprintf("d_%s", criteria)], original_result[,sprintf("d_%s", criteria)]) + 0.02
  ymin <- min(surrogate_result[,sprintf("d_%s", criteria)], original_result[,sprintf("d_%s", criteria)]) - 0.02
  
  plot(surrogate_result[,criteria], surrogate_result[,sprintf("d_%s", criteria)],
       xlab = criteria, ylab = sprintf("d_%s", criteria),
       main = plot_title, xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  points(original_result[,criteria], original_result[,sprintf("d_%s", criteria)], cex = 2, pch = 23, bg = "red3")
  abline(v=original_result[,criteria], lty = 1)
  abline(h=original_result[,sprintf("d_%s", criteria)], lty = 1)
  abline(v=surrogate_stat[c("2.5%", "97.5%"),criteria], lty = 2)
  abline(h=surrogate_stat[c("2.5%", "97.5%"),sprintf("d_%s", criteria)], lty = 2)
}