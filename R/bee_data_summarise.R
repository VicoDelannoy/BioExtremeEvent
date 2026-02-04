#' Merge the outputs of the metrics functions.
#'
#'@description To merge the daily outputs of at least two datasets of the
#' following fonctions from BEE package : bee.calc.metrics_point() ;
#' bee.calc.metrics_morpho() ; bee.calc.escape(). This function can also
#' summarise metrics over different time periods: extreme events, weak periods,
#' monthly periods and yearly periods.
#'
#'@param data
#'@param variable description
#'@param summarise_by takes the followings options:
#' - "extreme_event" for each metric, the following will be computed for each
#' extreme event: mean, median, variance, minimum and maximum, and the dates of
#' the maximum and minimum values.
#' - a numeric value indicating the time step for which you want to compute the
#' mean, median, variance, minimum and maximum values, as well as the dates of
#' the maximum and minimum values.
#'
#' #'@return
#'
#'@examples
#' # TO BE ADDED
#'
#'@export
#'
#-------------------------------------------------------------------------------
# data = merged_output ; variable = "perim_pixel_ratio" ; summarise_by = "extreme_event"

BEE.data.summarise <- function(
  data = data,
  variable = variable,
  summarise_by = "extreme_event"
) {
  ########################## CHECKS ############################################
  # Check that the summarize_by is a valid option
  # ("extreme_event"/"weak"/"month"/"year") :
  if (summarise_by != "extreme_event" | !is.numeric(summarise_by)) {
    warnings(
      "The summarise_by argument is not a suitable option, please 
    provide one from the following list : extreme_event, day, weak, month, year"
    )
  }
  if (is.list(data) == TRUE) {
    data <- do.call(rbind, data)
  }
  variable <- as.character(variable)
  # Is it a variable from BEE ?
  BEE_col_names <- c(
    "pixel_id",
    "date",
    "original_value",
    "cleanned_value",
    "event_id",
    "duration",
    "ID_df_point",
    "value",
    "evolution_rate_lag_1",
    "variance_value_lag_1",
    "x_df_point",
    "y_df_point",
    "first_date",
    "last_date",
    "mean_value",
    "median_value",
    "max_value",
    "date_max_value",
    "days_onset_abs",
    "raw_onset_rate_abs",
    "days_offset_abs",
    "raw_offset_rate_abs",
    "daily_rates",
    "baseline_qt",
    "baseline_mean",
    "anomaly_qt",
    "anomaly_mean",
    "anomaly_unit",
    "daily_category",
    "cumulative_anomaly_qt",
    "sum_anomaly_qt",
    "max_category",
    "mean_anomaly_qt",
    "sd_anomaly_qt",
    "max_anomaly_qt",
    "mean_anomaly_mean",
    "sd_anomaly_mean",
    "max_anomaly_mean",
    "most_frequent_category",
    "x_df_morpho",
    "y_df_morpho",
    "patch_id",
    "bordering",
    "patch_area",
    "core_area",
    "core_area_index",
    "n_pixel",
    "pixel_EE_tot",
    "cover_percent",
    "core_pixel",
    "ID_df_morpho",
    "perimeter",
    "centroid_x",
    "centroid_y",
    "perim_pixel_ratio",
    "perim_area_ratio_m",
    "shape_index",
    "fractal_cor_index",
    "circle_ratio_index",
    "contiguity_index",
    "x_df_escape",
    "y_df_escape",
    "to_x",
    "to_y",
    "pixel_to_id",
    "distance",
    "azimut",
    "ID_df_escape"
  )
  if (!(variable %in% BEE_col_names)) {
    warnings(paste0(
      "'",
      variable,
      "'",
      "is not a column name from the BEE package. 
    We recommend not changing the name of the outputs before using this function
    because some variables require special functions. For instance, 'azimut' is
    a circular variable so a special function is used to compute mean, sd etc.
    These variables are recognised by this function, if their names are changed,
    they will be treated as regular variables which may result in silent
    errors."
    ))
  }
  ############################ CODE #############################################
  ### Identify a column ID
  id <- colnames(data)[which(
    colnames(data) %in% c("ID_df_point", "ID_df_morpho", "ID_df_escape")
  )][1]
  if (summarise_by == "extreme_event") {
    ### Particular cases:
    if (variable == "azimut") {
      data$azimut_num <- as.numeric(data$azimut)
      data$azimut_circ <- circular::circular(
        data$azimut_num,
        units = "degrees",
        template = "geographics"
      )
      azimut_mean <- tapply(
        data$azimut_circ,
        data[, id],
        circular::mean.circular,
        na.rm = FALSE
      )
      data$azimut_mean <- (azimut_mean[data[, id]] + 360) %% 360

      azimut_median <- tapply(
        data$azimut_circ,
        data[, id],
        median_if_no_na
      )
      data$azimut_median <- azimut_median[data[, id]]

      azimut_sd <- tapply(
        data$azimut_circ,
        data[, id],
        sd_if_no_na
      )
      data$azimut_sd <- azimut_sd[data[, id]]
    } else {
      output <- summarise_ID(data = data, variable = variable)
    }
    groups <- split(seq_len(nrow(data)), data[[id]])
    pixel_id <- lapply(groups, function(idx) {
      unique(sub("_.*", "", data[idx, id]))
    })
    output$pixel_id <- as.numeric(unlist(pixel_id))
    output <- split(output, output$pixel_id)
    return(output)
  }
  ### regular situation:
  if (is.numeric(summarise_by)) {
    ### Particular cases:
    if (variable == "azimut") {
      # conversion of angle format to "circular"
      data$azimut_num <- as.numeric(data$azimut)
      data$azimut_circ <- circular::circular(
        data$azimut_num,
        units = "degrees",
        template = "geographics"
      )

      data$azimut_mean <- NA
      data$azimut_median <- NA
      data$azimut_sd <- NA
      for (i in seq(summarise_by, length(data[, variable]), 1)) {
        start <- i - summarise_by + 1
        index <- start:i # indices à utiliser
        data$azimut_mean[i] <- (circular::mean.circular(
          data[index, "azimut_circ"],
          na.rm = FALSE
        ) +
          360) %%
          360

        if (any(is.na(data[index, "azimut_circ"]))) {
          data$azimut_median[i] <- NA
        } else {
          data$azimut_median[i] <- (circular::median.circular(
            data[index, "azimut_circ"],
            na.rm = FALSE
          ) +
            360) %%
            360
        }
        data$azimut_sd[i] <- (circular::sd.circular(
          data[index, "azimut_circ"],
          na.rm = FALSE
        ) +
          360) %%
          360
      }
    } else {
      output <- summarise_lag(data = data, variable = variable)
    }
  }
}


#' To compute mean median variance max min and day of maximum of every metrics
#' to be summarized PER ID
#' @noRd
#'
summarise_ID <- function(data = data, variable = "var_name") {
  data[[variable]] = as.numeric(data[[variable]])
  funs <- list(
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE),
    var = function(x) stats::var(x, na.rm = TRUE),
    min = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
    max = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
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

  List <- split(seq_len(nrow(data)), data$ID_df_point)
  dates_max_min <- lapply(
    List,
    function(L_id) {
      data_mini <- data[L_id, ]
      x <- data_mini[[as.character(variable)]] # <- ici
      d <- data_mini$date

      if (all(is.na(x))) {
        list(NA, NA)
      } else {
        d_max <- d[which.max(x)]
        d_min <- d[which.min(x)]
        list(d_max, d_min)
      }
    }
  )
  date_max <- as.Date(sapply(date_max, function(x) x[[1]]))
  date_min <- as.Date(sapply(dates_max_min, function(x) x[[2]]))
  output$date_max <- date_max[output[, id]]
  output$date_min <- date_min[output[, id]]

  return(output)
}

#' To compute azimut median taking in account NA
#' @noRd

median_if_no_na <- function(x) {
  if (any(is.na(x))) {
    return(NA) # retourne NA si au moins un NA
  } else {
    return((circular::median.circular(x, na.rm = FALSE) + 360) %% 360)
  }
}

#' To compute azimut sd taking in account NA
#' @noRd

sd_if_no_na <- function(x) {
  if (any(is.na(x))) {
    return(NA) # retourne NA si au moins un NA
  } else {
    return((circular::sd.circular(x, na.rm = FALSE)) * 180 / pi)
  }
}

#' To compute mean median variance max min and day of maximum of every metrics
#' to be summarized PER TIME LAG
#' @noRd
#'
summarise_lag <- function(data = data, variable = "var_name") {
  data[[variable]] = as.numeric(data[[variable]])
  funs <- list(
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE),
    var = function(x) stats::var(x, na.rm = TRUE),
    min = function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE),
    max = function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE),
    sum = function(x) sum(x, na.rm = TRUE)
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

  List <- split(seq_len(nrow(data)), data$ID_df_point)
  dates_max_min <- lapply(
    List,
    function(L_id) {
      data_mini <- data[L_id, ]
      x <- data_mini[[as.character(variable)]] # <- ici
      d <- data_mini$date

      if (all(is.na(x))) {
        list(NA, NA)
      } else {
        d_max <- d[which.max(x)]
        d_min <- d[which.min(x)]
        list(d_max, d_min)
      }
    }
  )
  date_max <- as.Date(sapply(date_max, function(x) x[[1]]))
  date_min <- as.Date(sapply(dates_max_min, function(x) x[[2]]))
  output$date_max <- date_max[output[, id]]
  output$date_min <- date_min[output[, id]]

  return(output)
}
