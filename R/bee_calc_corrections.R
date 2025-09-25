#' Compute the impact of the selected definition of extreme events on the number
#' of days of extreme event
#' 
#' @description It computes the number of modifications between the raw 
#' binerised Spatraster and the Spartraster produced using the 
#' BEE.calc.true_event. 
#' This indicates of the differences in the number of days above the
#' threshold/baseline and the number of days that are considered to be
#' extreme events, according to the criteria set in BEE.calc.true_event.
#' @param Events_corrected The second element of the BEE.calc.true_event ouput,
#' which is a list of data.table containing information 
#' about the value of each pixel before and after definition criteria are 
#' applied to distinguish isolated extreme days from extreme events).
#'
#' @returns A text summarising the modifications induced.
#' @export
#'
#' @examples 
#' TO BE ADD
#' 
#'

# La fonction analyse change était au tout début du script mais je l'ai déplacée
# après pour pas avoir de problème avec le Roxygen et l'aide.

#Apply to the list of dataframe
# results = c()


BEE.calc.corrections <- function(Events_corrected) {
  # Number of cores available
  num_cores <- parallel::detectCores() - 2
  # Calculate changes
  results <- parallel::mclapply(Events_corrected, analyze_changes, mc.cores = num_cores)
  results_df <- do.call(rbind, results)
  rm(results)
  
  #Calculate percentages (data for the whole area and the whole time period)
  abs_one_to_zero <- as.numeric(sum(results_df$one_to_zero, na.rm = TRUE))
  abs_1 <-  as.numeric(sum(results_df$total_ones, na.rm = TRUE))
  abs_zero_to_one <-  as.numeric(sum(results_df$zero_to_one, na.rm = TRUE))
  abs_0 <-  as.numeric(sum(results_df$total_zeros, na.rm = TRUE))
  #percentage of 1 converted to 0 :
  percentage_one_to_zero <- abs_one_to_zero * 100 / abs_1 
  percentage_zero_to_one <- abs_zero_to_one * 100 / abs_0
  
  #Calculate mean per year and pixel
  total_nb_info = sum(results_df$total_not_NA)
  abs_one_to_zero_Ymean <- abs_one_to_zero * 100 / total_nb_info 
  abs_1_Ymean <- abs_1 / total_nb_info
  # sd_1_Ymean <- as.numeric(sd(results_df$total_ones, na.rm=TRUE)) / 
  # total_nb_info # variability of being 1 among the pixel
  # corrections among the pixel :
  sd_1to0_Ymean <- as.numeric(sd(results_df$one_to_zero, na.rm = TRUE))  
  abs_zero_to_one_Ymean <- abs_zero_to_one / total_nb_info
  abs_0_Ymean <- abs_0 / total_nb_info
  # corrections among the pixel :
  sd_1to0_Ymean <- as.numeric(sd(results_df$zero_to_one, na.rm = TRUE))  
  cat(
    "Absolut values for the whole area and whole time period :",
    "\n",
    "Original total number of days above threshold (1):",
    "\n",
    abs_1,
    "\n",
    "Percentage of 1 converted to 0:",
    "\n",
    percentage_one_to_zero,
    "%",
    ", (percentage among all initial '1' values)",
    "\n",
    "Original total number of days below threshold (0):",
    "\n",
    abs_0,
    "\n",
    "Percentage of 0 converted to 1:",
    "\n",
    percentage_zero_to_one,
    "%",
    "\n",
    "\n",
    "Mean values per pixel and date :",
    "\n",
    "Mean percentage of pixel value above threshold :",
    "\n",
    abs_1_Ymean * 100,
    "%",
    "\n",
    "Mean percentage of corrections from 1 to 0:",
    "\n",
    abs_one_to_zero_Ymean,
    "%",
    "(percentage among all values (1 & 0))",
    "\n",
    "Standard deviation of corrections among pixel(1 to 0):",
    "\n",
    sd_1to0_Ymean,
    "%",
    "\n",
    "Mean percentage of pixel value bellow thershold :",
    "\n",
    abs_0_Ymean * 100,
    "%",
    "\n",
    "Mean percentage of corrections to 1 ",
    "\n",
    abs_zero_to_one_Ymean,
    "%",
    "\n",
    "Standard deviation of corrections among pixel (0 to 1)",
    "\n",
    sd_1to0_Ymean,
    "%",
    "\n",
    "\n",
    "NA values in your dataset :",
    "\n",
    "Total number of NA values :",
    sum(results_df$total_NA),
    "\n",
    "Total number of pixel that are always NA :",
    sum(sapply(Events_corrected, function(df)
      all(is.na(df$Original_value)))),
    "\n",
    "Total number of NA values excluding those from pixels that are 
    always NA :",
    # total NA value - number of pixels that alwayes are NA * number of days 
    # in the time serie :
    sum(results_df$total_NA) - sum(sapply(Events_corrected, function(df)
all(is.na(df$Original_value)))) * length(Events_corrected[[1]]$Original_value),
    "\n"  
  )
  
}

analyze_changes <- function(dataframe) {
  # Count modifications
  one_to_zero <- sum(dataframe$Original_value == 1 &
                       dataframe$Cleanned_value == 0,
                     na.rm = TRUE)
  zero_to_one <- sum(dataframe$Original_value == 0 &
                       dataframe$Cleanned_value == 1,
                     na.rm = TRUE)
  
  # Count orginal total value
  total_ones <- sum(dataframe$Original_value == 1, na.rm = TRUE)
  total_zeros <- sum(dataframe$Original_value == 0, na.rm = TRUE)
  
  # Count NA
  total_NA <- sum(is.na(dataframe$Original_value))
  total_not_NA <- sum(!is.na(dataframe$Original_value))
  tot_value <- total_not_NA + total_NA
  
  return(
    data.frame(
      one_to_zero = one_to_zero,
      zero_to_one = zero_to_one,
      total_ones = total_ones,
      total_zeros = total_zeros,
      total_NA = total_NA,
      total_not_NA = total_not_NA,
      tot_value = tot_value
    )
  )
}
