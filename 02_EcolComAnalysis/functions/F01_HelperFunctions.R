####
#### F01. Collection of helper functions
#### For data analysis
####

# Compile rice data
adjust_edna_df_format <- function(edna_samdf, edna_otudf, plot_number, edna_data = edna_all0, date_data = date_for_df){
  # Require phyloseq package
  #require(phyloseq)
  #require(lubridate)
  
  # Check date indices
  match_date_id <- match(mdy(edna_samdf[edna_samdf$plot == plot_number, "date"]), date_data)
  no_growth_id1 <- 1:(match_date_id[1]-1)
  no_growth_id2 <- (rev(match_date_id)[1]+1):length(date_data)
  
  # Format non-numeric data
  edna_samplename <- as.factor(c(rep(NA, length(no_growth_id1)), as.character(edna_samdf[edna_samdf$plot == plot_number,][,"Sample_Name2"]), rep(NA, length(no_growth_id2))))
  edna_description <- as.factor(c(rep(NA, length(no_growth_id1)), as.character(edna_samdf[edna_samdf$plot == plot_number,][,"Description"]), rep(NA, length(no_growth_id2))))
  edna_plot <- as.character(c(rep(NA, length(no_growth_id1)), edna_samdf[edna_samdf$plot == plot_number,][,"plot"], rep(NA, length(no_growth_id2))))
  edna_samplenc <- as.factor(c(rep(NA, length(no_growth_id1)), as.character(edna_samdf[edna_samdf$plot == plot_number,][,"sample_nc"]), rep(NA, length(no_growth_id2))))
  edna_filt045ml <- as.numeric(c(rep(NA, length(no_growth_id1)), edna_samdf[edna_samdf$plot == plot_number,][,"filt045_ml"], rep(NA, length(no_growth_id2))))
  edna_filt022ml <- as.numeric(c(rep(NA, length(no_growth_id1)), edna_samdf[edna_samdf$plot == plot_number,][,"filt022_ml"], rep(NA, length(no_growth_id2))))
  edna_validity <- as.factor(c(rep(NA, length(no_growth_id1)), as.character(edna_samdf[edna_samdf$plot == plot_number,][,"std_validity"]), rep(NA, length(no_growth_id2))))
  
  # Format numeric meta data
  numeric_cols <- colnames(edna_samdf)[9:ncol(edna_samdf)]
  first_meta_nas <- data.frame(matrix(rep(NA, length(no_growth_id1)*length(numeric_cols)), ncol = length(numeric_cols)))
  last_meta_nas <- data.frame(matrix(rep(NA, length(no_growth_id2)*length(numeric_cols)), ncol = length(numeric_cols)))
  colnames(first_meta_nas) <- colnames(last_meta_nas) <- numeric_cols
  edna_numeric_meta_vars <- rbind(first_meta_nas, edna_samdf[edna_samdf$plot == plot_number,][,numeric_cols], last_meta_nas)

  # Format numeric otu (ASV) data
  first_otu_nas <- data.frame(matrix(rep(NA, length(no_growth_id1)*ncol(edna_otudf)), ncol = ncol(edna_otudf)))
  last_otu_nas <- data.frame(matrix(rep(NA, length(no_growth_id2)*ncol(edna_otudf)), ncol = ncol(edna_otudf)))
  colnames(first_otu_nas) <- colnames(last_otu_nas) <- colnames(edna_otudf)
  edna_numeric_otu_vars <- rbind(first_otu_nas, edna_otudf, last_otu_nas)
  
  # Combine all meat data including date
  edna_formatted_meta_df <- cbind(date_data, edna_samplename, edna_description, edna_plot, edna_samplenc, edna_filt045ml, edna_filt045ml, edna_validity, edna_numeric_meta_vars)
  colnames(edna_formatted_meta_df)[1:8] <- c("date", "Sample_Name2", "Description", "plot", "sample_nc", "filt045_ml", "filt022_ml", "std_validity")
  
  # Combine all DNA data and meta data
  edna_formatted_df <- cbind(edna_formatted_meta_df, edna_numeric_otu_vars)
  edna_formatted_df$plot <- as.character(edna_formatted_df$plot)
  edna_formatted_df$plot[!is.na(edna_formatted_df$plot)] <- as.character(plot_number)

  return(edna_formatted_df)
}

# Functions to calculate cumulative values
calculate_cum_vals <- function(ts_df, var_name, cum_range = 2:14, cum_func = "sum"){
  if(cum_func == "sum"){
    for(cum_i in cum_range){
      cum_col_name <- sprintf("%s_cum%03ddays", var_name, cum_i)
      ts_df[,cum_col_name] <- apply(embed(c(rep(NaN, cum_i - 1), ts_df[,var_name]), dimension = cum_i), 1, sum)
    }
  }else if(cum_func == "mean"){
    for(cum_i in cum_range){
      cum_col_name <- sprintf("%s_mean%03ddays", var_name, cum_i)
      ts_df[,cum_col_name] <- apply(embed(c(rep(NaN, cum_i - 1), ts_df[,var_name]), dimension = cum_i), 1, mean)
    }
  }
  return(ts_df)
}
