#'@param summarize_by takes the followings options, into "" :
#' - "extreme_event" for each metrics, a mean, median, variance, min and max will
#' be computed for each extrem event
#' - "day" metrics are not summarized through time, the function keeps a daily
#' resolution and just merge the datasets.
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

BEE.data.summarise <- function(
  data = data,
  variable = variable,
  summarise_by = "extreme_event"
) {
  # Check that the summarize_by is a valid option
  # ("extreme_event"/"weak"/"month"/"year") :
  if (!(summarise_by %in% c("extreme_event", "day", "weak", "month", "year"))) {
    warnings(
      "The summarise_by argument is not a suitable option, please 
    provide one from the following list : extreme_event, day, weak, month, year"
    )
  }

  ### Summarize by time
  if (summarise_by == "day") {
    merged_df <- split(merged_df, merged_df$pixel_id)
    return(merged_df)
  }
  if (summarize_by == "extreme_event") {
    #

    # Delete daily value:

    merged_df <- split(merged_df, merged_df$pixel_id)
    return(merged_df)
  }
}


#' To compute mean median variance max min and day of maximum of every metrics
#' to be summarized
#' @noRd
#'
summarize_ID <- function(data = df, variable = "var_name") {
  #data <- merged_df
  data[[variable]] = as.numeric(data[[variable]])
  funs <- list(
    mean = function(x) mean(x, na.rm = TRUE),
    var = function(x) stats::var(x, na.rm = TRUE),
    min = function(x) min(x, na.rm = TRUE),
    max = function(x) max(x, na.rm = TRUE)
  )

  res_list <- lapply(
    funs,
    function(f) {
      tapply(
        X = data[[variable]],
        INDEX = data$ID_df_point,
        FUN = f
      )
    }
  )

  output <- data.frame(
    ID_df_point = names(res_list[[1]]),
    do.call(
      cbind,
      lapply(
        names(res_list),
        function(nm) as.numeric(res_list[[nm]])
      )
    ),
    row.names = NULL
  )

  names(output)[-1] <- paste0(variable, "_", names(funs))
  return(output)
}
