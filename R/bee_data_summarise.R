#' Summarise the outputs of the 'merge' function.
#'
#'@description
#' To summarise one of the daily outputs of the
#' following fonctions from BEE package : bee.calc.metrics_point() ;
#' bee.calc.metrics_morpho() ; bee.calc.escape() ; bee.data.merge(). Metrics
#' values can be summarised per extreme event or over a sliding time window.
#'
#' @details
#' In addition to the mean value, the function provides the median, standard
#' deviation, minimum and maximum values, sum, date of minimum and maximum
#' values (within the time window), and date of minimum value (also within the
#' time window).
#'
#'@param data the dataset containing the column you want to process (e.g. the
#' ouput of BEE.data.merge(), as a list of data frame).
#'@param variable the name of the column for which you want to compute mean,
#' median etc. over a time window.
#'@param summarise_by takes the following options:
#' - "extreme_event" for each metric, the following will be computed for each
#' extreme event: mean, median, variance or sd, minimum and maximum, sum, and
#' the dates of the maximum and minimum values.
#' - a numeric value indicating the time step for which you want to compute the
#' mean, median, etc.
#'
#' #'@return a list of dataframes, with one dataframe per pixel, with the dates
#' in the rows and the values summarised.
#'
#'@examples
#' # TO BE ADDED
#'
#'@export
#'
#-------------------------------------------------------------------------------
# data = merged_output ; variable = "perim_pixel_ratio" ;
# summarise_by = "extreme_event" ; summarise_by = 3 ; variable = "azimut"

BEE.data.summarise <- function(
  data = data,
  variable = variable,
  summarise_by = "extreme_event"
) {
  ########################## CHECKS ############################################
  # Check that the summarize_by is a valid option
  # ("extreme_event"/"weak"/"month"/"year") :
  if (summarise_by != "extreme_event" & !is.numeric(summarise_by)) {
    warnings(
      "The summarise_by argument is not a suitable option, please 
    provide 'extreme_event' or a number of time unit."
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
  ############################ CODE ############################################
  ### Identify a column ID
  id <- colnames(data)[which(
    colnames(data) %in% c("ID_df_point", "ID_df_morpho", "ID_df_escape")
  )][1]

  if (summarise_by == "extreme_event") {
    ### Particular cases:
    if (variable == "azimut") {
      data$azimut_num <- as.numeric(data$azimut)
      n_rows <- length(unique(data[, id]))
      output <- data.frame(
        ID = unique(data[, id]),
        azimut_mean = circular::circular(rep(NA, n_rows)),
        azimut_median = circular::circular(rep(NA, n_rows)),
        azimut_sd = circular::circular(rep(NA, n_rows))
      )

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
      azimut_mean <- (azimut_mean + 360) %% 360
      output$azimut_mean <- azimut_mean[match(output$ID, names(azimut_mean))]

      azimut_median <- tapply(
        data$azimut_circ,
        data[, id],
        median_if_no_na
      )
      azimut_median <- (azimut_median + 360) %% 360
      output$azimut_median <- azimut_median[match(
        output$ID,
        names(azimut_median)
      )]

      azimut_sd <- tapply(
        data$azimut_circ,
        data[, id],
        sd_if_no_na
      )
      azimut_sd <- (azimut_sd + 360) %% 360
      output$azimut_sd <- azimut_sd[match(output$ID, names(azimut_sd))]
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
      data$azimut_num <- as.numeric(data$azimut)
      n_rows <- length(data[, id])
      output <- data.frame(
        ID = data[, id],
        date = data[, "date"],
        azimut = circular::circular(
          data$azimut_num,
          units = "degrees",
          template = "geographics"
        ),
        data[, "azimut_num"],
        azimut_mean = circular::circular(rep(NA, n_rows)),
        azimut_median = circular::circular(rep(NA, n_rows)),
        azimut_sd = circular::circular(rep(NA, n_rows)),
        pixel_id = data[, "pixel_id"]
      )

      for (i in seq(summarise_by, length(data[, variable]), 1)) {
        start <- i - summarise_by + 1
        index <- start:i

        output$azimut_mean[i] <- ifelse(
          is.na(output[i, "azimut"]),
          NA,
          circular::mean.circular(
            output[index, "azimut"],
            na.rm = TRUE
          )
        )

        output$azimut_median[i] <- ifelse(
          is.na(output[i, "azimut"]),
          NA,
          circular::median.circular(
            output[index, "azimut_circ"],
            na.rm = TRUE
          )
        )

        output$azimut_sd[i] <- ifelse(
          is.na(output[i, "azimut"]),
          NA,
          circular::sd.circular(
            output[index, "azimut_circ"],
            na.rm = TRUE
          )
        )
      }
      output$azimut_mean <- (circular::conversion.circular(
        output$azimut_mean,
        units = "degrees"
      ) +
        360) %%
        360
      output$azimut_median <- (circular::conversion.circular(
        output$azimut_median,
        units = "degrees"
      ) +
        360) %%
        360
      output$azimut_sd <- (circular::conversion.circular(
        output$azimut_sd,
        units = "degrees"
      ) +
        360) %%
        360
      output <- split(output, output$pixel_id)
      return(output)
    } else {
      output <- summarise_lag(
        data = data,
        variable = variable,
        lag_time = summarise_by
      )
      output <- split(output, output$pixel_id)
      return(output)
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
        INDEX = data[, id],
        FUN = f
      )
    }
  )

  output <- data.frame(
    ID = names(res_list[[1]]),
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

  List <- split(seq_len(nrow(data)), data[, id])
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
  date_max <- as.Date(sapply(dates_max_min, function(x) x[[1]]))
  date_min <- as.Date(sapply(dates_max_min, function(x) x[[2]]))
  output[["date_max"]] <- date_max[match(output$ID, names(date_max))]
  output[, "date_min"] <- date_min[match(output$ID, names(date_min))]

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
summarise_lag <- function(
  data = data,
  variable = "var_name",
  lag_time = summarise_by
) {
  data[[variable]] = as.numeric(data[[variable]])
  output <- as.data.frame(
    stats::setNames(
      list(
        data[id],
        data$date,
        data[, as.character(variable)],
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        NA_real_,
        as.Date(NA),
        as.Date(NA),
        data[["pixel_id"]]
      ),
      c(
        "ID",
        "date",
        as.character(variable),
        paste0(variable, "_mean"),
        paste0(variable, "_median"),
        paste0(variable, "_variance"),
        paste0(variable, "_min"),
        paste0(variable, "_max"),
        paste0(variable, "_sum"),
        paste0(variable, "_date_min"),
        paste0(variable, "_date_max"),
        "pixel_id"
      )
    )
  )

  for (i in seq(lag_time, length(data[, variable]), 1)) {
    start <- i - lag_time + 1
    index <- start:i

    output[i, paste0(as.character(variable), "_", "mean")] <- ifelse(
      is.na(data[i, variable]),
      NA,
      mean(
        data[index, variable],
        na.rm = TRUE
      )
    )
    output[i, paste0(as.character(variable), "_", "median")] <- ifelse(
      is.na(data[i, variable]),
      NA,
      median(
        data[index, variable],
        na.rm = TRUE
      )
    )
    output[i, paste0(as.character(variable), "_", "variance")] <- ifelse(
      is.na(data[i, variable]),
      NA,
      stats::var(
        data[index, variable],
        na.rm = TRUE
      )
    )

    output[i, paste0(as.character(variable), "_min")] <- ifelse(
      is.na(data[i, variable]),
      NA,
      min(data[index, variable], na.rm = TRUE)
    )

    output[i, paste0(as.character(variable), "_max")] <- ifelse(
      is.na(data[i, variable]),
      NA,
      max(data[index, variable], na.rm = TRUE)
    )

    output[i, paste0(as.character(variable), "_sum")] <- ifelse(
      is.na(data[i, variable]),
      NA,
      sum(data[index, variable], na.rm = TRUE)
    )

    output[i, paste0(as.character(variable), "_", "date_min")] <- ifelse(
      is.na(data[i, variable]),
      as.Date(NA),
      as.Date(data$date[
        i -
          which.min(
            data[index, paste0(as.character(variable))]
          ) +
          1
      ])
    )

    output[i, paste0(as.character(variable), "_", "date_max")] <- ifelse(
      is.na(data[i, variable]),
      as.Date(NA),
      as.Date(data$date[
        i -
          which.max(
            data[index, paste0(as.character(variable))]
          ) +
          1
      ])
    )
  }

  return(output)
}
