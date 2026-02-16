#' Compute metrics through time for specific locations.
#'
#' @param events_corrected :
#'  Is the list of data tables produced by:
#'  bee_calc_true_event (the second element of the output). For each pixel, it
#'  contains a data.table with dates in the rows and a column indicating whether
#'  it is a day belonging to a heatwave (1) or not (0).
#' @param yourspatraster :
#'  Is the spatraster with the values of the studded parameter, through time and
#'  space.
#' @param gps :
#'  Is a data frame containing the positions for which you want to compute
#'  metrics. It must contain columns labelled 'x' and 'y', which should contain
#'  longitudes and latitudes, respectively.
#' @param start_date :
#'  First day on which you want to start computing metrics. Optional argument.
#' @param end_date :
#'  Last day on which you want to start computing metrics. Optional argument.
#' @param baseline_qt :
#'  Spatraster of the 90th percentile baseline (or the 10th percentile
#'  baseline). Optional argument.
#' @param baseline_fixed :
#'  A numerical value representing a fixed threshold for which you want to
#'  calculate the distance between the observed value and the threshold.
#'  Optional argument.
#' @param baseline_mean :
#'  Spatraster of the mean value baseline. Optional argument.
#' @param time_lapse_vector :
#'  A vector of time laps on which to compute mean evolution rate and variance
#'  of the studdied parameter. Optional argument.
#' @param group_by_event :
#'  Whether you want an output summarise by extreme event or not. If not, you
#'  just get daily metrics. Optional argument.
#'
#' @note
#'  *BEE.calc.metrics_point()* is not designed to work on 4D data
#'  (time + spatial 3D).
#'
#' @returns
#'  A list of dataframe (one per gps point), each dataframe contains
#'  informations on the date of the extrem events, mean, median, max and min
#'  values, peak day, onset-rate, off-set rate, mean anomaly, maximum category
#'  *etc.* Categories are defined in Hobday *et al.* 2018.
#'
#' @examples
#' # Prepare function arguments:
#' file_name_1 <- system.file(file.path("extdata",
#'                                      "copernicus_data_celsius.tiff"),
#'                                      package = "BioExtremeEvent")
#' copernicus_data_celsius <- terra::rast(file_name_1)
#' file_name_2 <- system.file(file.path("extdata",
#'                                      "binarized_corrected_df.rds"),
#'                                      package = "BioExtremeEvent")
#' binarized_corrected_df <- readRDS(file_name_2)
#' file_name_3 <- system.file(file.path("extdata",
#'                                      "baseline_qt90_smth_15.tiff"),
#'                                      package = "BioExtremeEvent")
#' baseline_qt90_smth_15 <- terra::rast(file_name_3)
#' file_name_4 <- system.file(file.path("extdata",
#'                                      "baseline_mean_smth_15.tiff"),
#'                                      package = "BioExtremeEvent")
#' baseline_mean_smth_15 <- terra::rast(file_name_4)
#' gps <- data.frame(x = c(3.6883,4.4268), # Sète, Saintes-Marie-de-la-Mer
#'                   y = c(43.3786,43.4279))
#'
#' # Get daily value (required for BEE.merge() and BEE.summarise()):
#' metrics_points_day <- BEE.calc.metrics_point(
#'  events_corrected = binarized_corrected_df,
#'  yourspatraster = copernicus_data_celsius,
#'  gps = gps,
#'  start_date = NULL,
#'  end_date = NULL,
#'  time_lapse_vector = NULL,
#' # number of time unit on which to compute the warming rates and cooling rates,
#' #  NULL means it will not be computated
#'  baseline_qt = baseline_qt90_smth_15,
#'  baseline_mean = baseline_mean_smth_15,
#'  baseline_fixed = 12.5,
#'  group_by_event = FALSE
#')
#'
#' # Get mean, min, max, sd and median per event (extreme event and btw extreme
#' # events):
#' metrics_points_ee <- BEE.calc.metrics_point(
#'  events_corrected = binarized_corrected_df,
#'  yourspatraster = copernicus_data_celsius,
#'  gps = gps,
#'  start_date = NULL,
#'  end_date = NULL,
#'  time_lapse_vector = NULL,
#' # number of time unit on which to compute the warming rates and cooling rates,
#' #  NULL means it will not be computated
#'  baseline_qt = baseline_qt90_smth_15,
#'  baseline_mean = baseline_mean_smth_15,
#'  group_by_event = TRUE
#')
#'
#' @export
#'

#-------------------------------------------------------------------------------
# start_date <- "2024-06-01" ; end_date <- "2024-10-31" ; yourspatraster <- ds ; p <- 1;
# gps <- data.frame(x = c(3.7659, 5.386, 3.146), y = c(43.4287, 43.183, 42.781));
# group_by_event = TRUE; time_lapse_vector = c(1,3,5,7,14,21) ;
# baseline_qt = baseline_qt90
BEE.calc.metrics_point <- function(
  events_corrected,
  yourspatraster,
  gps,
  start_date = NULL,
  end_date = NULL,
  time_lapse_vector = NULL, # number of time unit on which to compute the warming rates and cooling rates, NULL means it will not be computated
  baseline_qt = NULL,
  baseline_mean = NULL,
  baseline_fixed = NULL,
  group_by_event = FALSE
) {
  terra::set.names(yourspatraster, as.Date(terra::time(yourspatraster)))

  ############################### WARNINGS #####################################
  if (is.null(start_date) | is.null(end_date)) {
    message(
      "You didn't specify a begining date and a ending date (see argument 
      'start_date' and 'end_date'), the first date and last date in your time 
      yourspatraster will be used."
    )
    start_date <- min(as.Date.character(terra::names(yourspatraster)))
    end_date <- max(as.Date.character(terra::names(yourspatraster)))
  }
  #Check that start date and end_date are within the SpatRasters provided
  if (
    !(start_date %in% terra::names(yourspatraster)) |
      !(end_date %in% terra::names(yourspatraster))
  ) {
    warning(
      "One or both of the specified layers is/are not present in the SpatRaster 
      containing the corrected binarised extreme event. The provided dates may
      be outside the SpatRaster timeframe or in the wrong format. Please check 
      terra::time(values) to see in which format start_date and end_date must be
      provided in."
    )
  }
  yourspatraster_extent <- terra::ext(yourspatraster)
  #Check that all gps points are within yourspatraster and period_of_interest extent
  if (
    !all(
      abs(gps[1]) <= abs(yourspatraster_extent$xmax) &
        abs(gps[1]) >= abs(yourspatraster_extent$xmin)
    )
  ) {
    warning(
      "At least one longitude coordinate falls outside the extent of the 
      SpatRaster containing the binarized corrected events. Ensure the first 
      column of the dataframe contains valid longitude (x) values for analysis."
    )
  }
  if (
    !all(
      abs(gps[2]) <= abs(yourspatraster_extent$ymax) &
        abs(gps[2]) >= abs(yourspatraster_extent$ymin)
    )
  ) {
    warning(
      "At least one latitude coordinate falls outside the extent of the 
      SpatRaster containing the binarized corrected events. Ensure the second
      column of the dataframe contains valid latitude (y) values for analysis."
    )
  }
  # Check that one of the gps position is not in a pixel that is always an NA
  # (it may indicates that it falls in an area that is not interesting for the
  # study)
  ## List of pixel that are always NA:
  NA_pixels <- which(vapply(events_corrected, is.null, logical(1)))
  ## Identify the pixels corresponding to the gps position provided:
  gps$pixel <- terra::cellFromXY(yourspatraster, gps)
  if (any(gps$pixel %in% NA_pixels)) {
    wrong_position <- which(gps$pixel %in% NA_pixels)
    warning(
      paste(
        "Problematic gps positions:\n",
        paste(
          utils::capture.output(print(gps[wrong_position, ])),
          collapse = "\n"
        ),
        "These pixels fall on pixels that are always marked as NA, no matter the 
      date.  Please check the accuracy of the gps position(s) provided and 
      ensure that it is in the correct format."
      ),
      call. = FALSE
    )
  }

  ############################### CODE #########################################
  #Extract yourspatraster for the given gps position
  df_list <- lapply(gps$pixel, function(p) {
    events_corrected[[p]]
  }) # on df per points/pixel
  yourspatraster <- t(terra::extract(yourspatraster, gps[, 3]))

  # Subset both dataset so they match the timeframe provided with 'start_date'
  # and 'end_date'.
  Date <- rownames(yourspatraster)
  yourspatraster <- yourspatraster[
    which(
      as.Date(rownames(yourspatraster)) >= as.Date(start_date) &
        as.Date(rownames(yourspatraster)) <= as.Date(end_date)
    ),
  ]
  yourspatraster <- as.matrix(yourspatraster)
  df_list <- Map(
    function(df, col_idx) {
      df <- df[df$date >= start_date & df$date <= end_date, ]
      df$value <- yourspatraster[, col_idx] #Merge df_list and yourspatraster  # Ad to each
      # dataframe the corresponding column of pixel value
      return(df)
    },
    df_list,
    seq_along(df_list)
  )

  # For each event I want : duration, maximum intensity, mean and median
  # intensity, category, sum of anomalies, date of maximum intensity, position
  # of the day of maximum intensity, mean increasing slop, mean decreasing
  # slope, start_date, end_date

  #Create a list of dataframe to store the information using one element (df)
  # per pixel and one row per event
  colnames(gps) <- c("x", "y", names(gps[3]))

  #Get daily anomaly to baseline_qt and to baseline_mean
  if (!is.null(baseline_qt)) {
    qt <- as.data.frame(t(terra::extract(baseline_qt, gps[, 3])))
    qt$dates <- format(
      seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day"),
      "%m-%d"
    )
  }
  if (!is.null(baseline_mean)) {
    mean <- as.data.frame(t(terra::extract(baseline_mean, gps[, 3])))
    mean$dates <- format(
      seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day"),
      "%m-%d"
    )
  }

  metrics <- lapply(1:length(df_list), function(p) {
    # Go through each pixel
    df <- df_list[[p]]
    df$date <- as.Date(df$date)
    df$dates <- format(as.Date(df$date), "%m-%d")

    # Check that the first/last extreme event didn't start before/after the
    # selected timeframe.
    if (sub(".*?_(\\d{4}-\\d{2}-\\d{2}).*", "\\1", df$ID[1]) < df$date[1]) {
      warnings(
        "The first event starts before the start of the selected time 
               frame. Please note that metrics will only be computed using data
               within the selected time frame. Thus, metrics computed for the 
               first event will be based on a sub-selection of values belonging 
               to that event, using only values that are in the selected time
               frame."
      )
    }
    if (
      sub(".*_(\\d{4}-\\d{2}-\\d{2})_.*", "\\1", df$ID[nrow(df)]) >
        df$date[nrow(df)]
    ) {
      warnings(
        "The last event end after the last day of the selected time 
               frame. Please note that metrics will only be computed using data
               within the selected time frame. Thus, metrics computed for the 
               last event will be based on a sub-selection of values belonging 
               to that event, using only values that are in the selected time
               frame."
      )
    }

    ## Evolution rate on the given time laps
    if (!is.null(time_lapse_vector)) {
      for (lag in time_lapse_vector) {
        col_name <- paste0("evolution_rate_lag_", lag)
        df[[col_name]] <- NA_real_
        for (r in seq(lag + 1, nrow(df))) {
          df[[col_name]][r] <- (df$value[r] - df$value[r - lag]) / lag
        }
      }
    }

    ## Variance during the given time laps
    if (!is.null(time_lapse_vector)) {
      for (lag in time_lapse_vector) {
        col_name <- paste0("variance_value_lag_", lag)
        df[[col_name]] <- NA_real_
        for (r in seq(lag + 1, nrow(df))) {
          df[[col_name]][r] <- stats::var(df$value[(r - lag):r])
        }
      }
    }

    df <- df |>
      dplyr::arrange(ID, date) |> # put data in pixel order and chronological order
      dplyr::group_by(ID) |> # Create columns with info on the previous group
      dplyr::mutate(row_num = dplyr::row_number()) |> # useful later
      dplyr::ungroup() |>
      dplyr::mutate(
        prev_value = dplyr::lag(value),
        prev_ID = dplyr::lag(ID),
        last_value_prev_group = ifelse(
          row_num == 1 &
            ID != prev_ID,
          prev_value,
          NA
        )
      ) |>
      tidyr::fill(last_value_prev_group, .direction = "down") |>
      dplyr::group_by(ID) |>
      dplyr::mutate(
        # Get gps position of each point
        x = gps$x[p],
        y = gps$y[p],

        # Duration of each event (assuming duration is constant per ID)
        duration = data.table::first(duration),

        # First and last day of each event
        first_date = as.Date(min(date)),
        last_date = as.Date(max(date)),

        # Mean, median, max temperature per event
        mean_value = mean(value),
        median_value = stats::median(value),
        max_value = max(value),

        # Day of absolute maximum temperature
        date_max_value = as.Date(date[which.max(value)]),

        # Onset rate (up to the absolut maximum value during the event)
        days_onset_abs = as.numeric(date_max_value - first_date),
        ## raw
        raw_onset_rate_abs = ifelse(
          days_onset_abs > 0,
          (max_value - value[which(date == first_date)]) / days_onset_abs,
          value[1] - data.table::first(last_value_prev_group)
        ),

        # Offset rate (from the absolut maximum value during the event)
        days_offset_abs = as.numeric(last_date - date_max_value),
        ## raw
        raw_offset_rate_abs = ifelse(
          days_offset_abs > 0,
          (value[which(date == last_date)] - max_value) / days_offset_abs,
          utils::tail(value, 1) - utils::tail(value, 2)[1]
        )
      ) |>
      dplyr::ungroup()

    df <- df |>
      dplyr::mutate(
        #Daily rates
        daily_rates = (value - dplyr::lag(value)) /
          as.numeric(date - dplyr::lag(date)),
      )

    # Add daily baseline_qt and/or daily baseline_mean and/or baseline_fixed
    if (!is.null(baseline_qt)) {
      df <- df |>
        dplyr::left_join(
          qt[, c("dates", paste0("V", as.character(p)))],
          by = "dates"
        ) |>
        dplyr::rename(baseline_qt = paste0("V", as.character(p)))
    }
    if (!is.null(baseline_mean)) {
      df <- df |>
        dplyr::left_join(
          mean[, c("dates", paste0("V", as.character(p)))],
          by = "dates"
        ) |>
        dplyr::rename(baseline_mean = paste0("V", as.character(p)))
    }
    if (!is.null(baseline_fixed)) {
      df$baseline_fixed <- rep(baseline_fixed, nrow(df))
    }
    # Calculate anomalies
    if (!is.null(baseline_qt)) {
      df <- df |>
        dplyr::mutate(anomaly_qt = value - baseline_qt) |>
        dplyr::group_by(ID) |>
        dplyr::mutate(
          cumulative_anomaly_qt = cumsum(anomaly_qt),
          sum_anomaly_qt = sum(anomaly_qt, na.rm = TRUE)
        ) |>
        dplyr::ungroup()
      # Add mean value and standard deviation of anomaly_qt and anomaly_mean to
      # the ouputs + add maximal category of each event
      df <- df |>
        dplyr::group_by(ID) |>
        dplyr::mutate(
          mean_anomaly_qt = mean(anomaly_qt, na.rm = TRUE),
          sd_anomaly_qt = stats::sd(anomaly_qt, na.rm = TRUE),
          max_anomaly_qt = max(anomaly_qt, na.rm = TRUE)
        )
    }

    if (!is.null(baseline_mean)) {
      df <- df |>
        dplyr::mutate(
          anomaly_mean = value - baseline_mean
        )
      # Add mean value and standard deviation of anomaly_qt and anomaly_mean to
      # the ouputs + add maximal category of each event
      df <- df |>
        dplyr::group_by(ID) |>
        dplyr::mutate(
          mean_anomaly_mean = mean(anomaly_mean, na.rm = TRUE),
          sd_anomaly_mean = stats::sd(anomaly_mean, na.rm = TRUE),
          max_anomaly_mean = max(anomaly_mean, na.rm = TRUE)
        )
    }
    if (!is.null(baseline_fixed)) {
      df <- df |>
        dplyr::mutate(
          anomaly_fixed = value - baseline_fixed
        )
      # Add mean value and standard deviation of anomaly_qt and anomaly_mean to
      # the ouputs + add maximal category of each event
      df <- df |>
        dplyr::group_by(ID) |>
        dplyr::mutate(
          mean_anomaly_fixed = mean(anomaly_fixed, na.rm = TRUE),
          sd_anomaly_fixed = stats::sd(anomaly_fixed, na.rm = TRUE),
          max_anomaly_fixed = max(anomaly_fixed, na.rm = TRUE)
        )
    }
    if (!is.null(baseline_mean) & !is.null(baseline_qt)) {
      df <- df |>
        dplyr::mutate(
          anomaly_unit = baseline_qt - baseline_mean,
          daily_category = dplyr::case_when(
            cleaned_value == 0 ~ "No extreme event",
            anomaly_qt < anomaly_unit ~ "Category I",
            anomaly_qt < 2 * anomaly_unit & anomaly_qt >= 0 ~ "Category II",
            anomaly_qt < 3 * anomaly_unit & anomaly_qt >= 0 ~ "Category III",
            anomaly_qt >= 3 * anomaly_unit & anomaly_qt >= 0 ~ "Category IV",
            TRUE ~ NA_character_
          )
        )
      # Keep the maximum category reached by each event
      df <- df |>
        dplyr::mutate(
          daily_category = factor(
            daily_category,
            levels = c(
              "No extreme event",
              "Category I",
              "Category II",
              "Category III",
              "Category IV"
            ),
            ordered = TRUE
          )
        )

      max_category <- df |>
        dplyr::filter(!is.na(daily_category)) |>
        dplyr::group_by(ID) |>
        dplyr::summarise(max_category = max(daily_category), .groups = 'drop')
      df <- df |>
        dplyr::left_join(max_category, by = "ID")

      # Add mean value and standard deviation of anomaly_qt and anomaly_mean to
      # the ouputs + add maximal category of each event
      df <- df |>
        dplyr::group_by(ID) |>
        dplyr::mutate(
          max_category = names(sort(table(daily_category), decreasing = TRUE))[
            1
          ],
          # Catégorie la plus fréquente
          .groups = 'drop'
        )
    }

    if (group_by_event) {
      # Delete daily values that are not usefull to describe the full event
      columns_to_delete <- c(
        "baseline_qt",
        "baseline_mean",
        "anomaly_mean", # anomaly to the baseline_mean
        "anomaly_qt",
        "anomaly_unit",
        "cumulative_anomaly_qt",
        "daily_category",
        "daily_rates",
        "date",
        "daily_category",
        "row_num",
        "prev_value",
        "prev_ID",
        "last_value_prev_group",
        ".groups"
      )
      df <- df |>
        dplyr::select(
          -dplyr::starts_with("evolution_rate_lag_"),
          -dplyr::starts_with("variance_value_lag_")
        ) |>
        dplyr::distinct(ID, .keep_all = TRUE) # One row per ID
      df <- df[, !(names(df) %in% columns_to_delete)]
    } else {
      df <- df |>
        dplyr::select(
          -dates,
          -row_num,
          -prev_value,
          -prev_ID,
          -last_value_prev_group,
          -.groups
        )
    }
    return(df)
  })
}

# OTHER PART TO DEVELOP :
# For each pixels :
# - anomalie cumulée par événements
