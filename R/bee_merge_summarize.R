#' Merge the outputs of the metrics functions and summarize them trought time.
#'
#'@description
#'
#'@param data_metrics_point the output of the BEE.calc.metrics_point computed 
#' using the argument group_by_event = FALSE.
#' 
#'@param data_metrics_morpho the output of the BEE.calc.metrics_morpho 
#' 
#'@param data_escape the output of the BEE.calc.escape computed 
#' using the argument group_by_event = FALSE.
#' 
#'@param summarize_by takes the followings options, into "" :
#' - "extrem_event" for each metrics, a mean, median, variance, min and max will
#' be computed for each extrem event
#' - "week" for each metrics, a mean, median, variance, min and max will
#' be computed weekly. Note that this imply to provide daily output 
#' (BEE.calc.escape(only_days_EE = FALSE))
#' - "two_weeks" for each metrics, a mean, median, variance, min and max will
#' be computed every two weeks. Note that this imply to provide daily output 
#' (BEE.calc.escape(only_days_EE = FALSE))
#' - "month" for each metrics, a mean, median, variance, min and max will
#' be computed monthly. Note that this imply to provide daily output 
#' (BEE.calc.escape(only_days_EE = FALSE))
#' - "year" for each metrics, a mean, median, variance, min and max will
#' be computed yearly. Note that this imply to provide daily output 
#' (BEE.calc.escape(only_days_EE = FALSE))
#'@return 
#' 
#'@examples
#'
#'@export
#-------------------------------------------------------------------------------

# data_metrics_point <- points_metrics ; 
# data_metrics_morpho <- list_morpho_metrics ; 
# data_escape <- dist_to_escape

BEE.merge_summarize <- function(
  data_metrics_point = NULL,
  data_metrics_morpho = NULL,
  data_escape = NULL,
  summarize_by = "extrem_event"
) {
# Identify which datasets have been provided
  ## Is there at least two datasets provided ?
  check_1 <- paste0(is.null(data_metrics_point),
                   is.null(data_metrics_morpho),
                   is.null(data_escape))

  if(check_1 == "TRUETRUETRUE"){
    warnings("You haven't provided any data to be merged and summarised. Please
    provide at least two datasets.")
  }
  if(check_1 == "FALSETRUETRUE"){
    message("As you have only provided the metrics_point dataset, 
    summarising through time is all that will be done.")
  }
  if(check_1 == "TRUEFALSETRUE"){
    message("As you have only provided the metrics_morpho dataset, 
    summarising through time is all that will be done.")
  }
  if(check_1 == "TRUETRUEFALSE"){
    message("As you have only provided the escape dataset, 
    summarising through time is all that will be done.")
  }
  
  ## Is it a daily resolution ?
  if(!is.null(data_metrics_point)){
    dates <- data_metrics_point[[1]]$Date
    if(length(dates)!=dates[length(dates)]-dates[1]+1){
      warnings("Some days are missing in data_metrics_point or the time 
      resolution is not daily. Please, make sure to use 'group_by_event = FALSE'
      when running BEE.calc.metrics_point")
    }
  }
  if(!is.null(data_metrics_morpho)){
    data_metrics_morpho <- data_metrics_morpho[[2]]
    dates <- data_metrics_morpho[ vapply(data_metrics_morpho, ncol, integer(1)) == 6 ][[1]]$date
    dates <- as.Date(dates)
    if(length(dates)!=dates[length(dates)]-dates[1]+1){
      warnings("Some days are missing in data_metrics_point or the time 
      resolution is not daily. Please, make sure to use 'group_by_event = FALSE'
      when running BEE.calc.metrics_point")
    }
  }
  if(!is.null(data_escape)){
    # TO COMPLET
  }
  # tester homogénéité des périodes fournies
  # tester homogénéité spatiale
}
