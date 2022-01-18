####
#### Helper functions for CERrice2017 time series
#### Functions for meta analysis
####

# ggplot function
geom_obspred <- function(df, plot_title){
  ggobj <- ggplot(df, aes(y = exp(observed), x = exp(predicted)))
  ggobj <- ggobj + geom_point(alpha = 0.5, size = 2) + geom_smooth(method = "lm", color = "red3")
  ggobj <- ggobj + geom_abline(intercept = 0, slope = 1, linetype = 2)
  ggobj <- ggobj + xlab("Predicted diversity") + ylab("Observed diversity")
  ggobj <- ggobj + ggtitle(plot_title)
  return(ggobj)
}

# Data frame for gam plot
make_df_gam <- function(gam_result, dataset_name, scale_x = TRUE){
  if(scale_x){
    df_gam_effect <- data.frame(log_temp_coef = plot(gam_result, select = 1, seWithMean = T)[[1]]$fit,
                                log_abun_coef = plot(gam_result, select = 2, seWithMean = T)[[2]]$fit,
                                log_temp = scale(plot(gam_result, select = 1, seWithMean = T)[[1]]$x),
                                log_abun = scale(plot(gam_result, select = 2, seWithMean = T)[[2]]$x),
                                dataset = dataset_name)
  }else{
    df_gam_effect <- data.frame(log_temp_coef = plot(gam_result, select = 1, seWithMean = T)[[1]]$fit,
                                log_abun_coef = plot(gam_result, select = 2, seWithMean = T)[[2]]$fit,
                                log_temp = exp(plot(gam_result, select = 1, seWithMean = T)[[1]]$x) - 273.15,
                                log_abun = plot(gam_result, select = 2, seWithMean = T)[[2]]$x,
                                dataset = dataset_name)
  }
  return(df_gam_effect)
}