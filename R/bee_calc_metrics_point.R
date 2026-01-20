# BEE.calc.metrics_point is not designed to work on 4D data (time + spatial 3D).

#' BEE.calc.metrics_point compute metrics through time for specific locations.
#'
#' @param Events_corrected Is the list of data tables produced by:
#' bee_calc_true_event (the second element of the output). For each pixel, it
#' contains a data.table with dates in the rows and a column indicating whether
#' it is a day belonging to a heatwave (1) or not (0).
#' @param Values is the spatraster with the values of the studded parameter,
#'  through time and space.
#' @param GPS is a data frame containing the positions for which you want to
#' compute metrics. It must contain columns labelled 'x' and 'y', which should
#' contain longitudes and latitudes, respectively.
#' @param start_date First day on which you want to start computing metrics.
#' @param end_date Last day on which you want to start computing metrics.
#' @param baseline_qt Spatraster of the 90th percentile baseline (or the 10th
#'  percentile baseline)
#' @param baseline_mean Spatraster of the 90mean value baseline.
#' @param time_lapse_vector a vector of time laps on which to compute mean
#' evolution rate and variance of the studdied parameter.
#' @param group_by_event Whether you want an output summarise by extreme event
#' or not. If not, you just get daily metrics.
#'
#' @returns A list of dataframe (one per GPS point), each dataframe contains
#' informations on the date of the extrem events, mean, median, max and min
#' values, peak day, onset-rate, off-set rate, mean anomaly, maximum category
#' etc. Categories are defined in Hobday et al. 2018.
#'
#' @examples
#' # To be added
#'
#' @export
#'

#-------------------------------------------------------------------------------
# start_date <- "2024-06-01" ; end_date <- "2024-10-31" ; Values <- ds ; p <- 1
# GPS <- data.frame(x = c(3.7659, 5.386, 3.146), y = c(43.4287, 43.183, 42.781))
# group_by_event = TRUE, time_lapse_vector = c(3,5,7,14,21)
BEE.calc.metrics_point <- function(
  Events_corrected,
  Values,
  GPS,
  start_date = NULL,
  end_date = NULL,
  time_lapse_vector = NULL, # number of time unit on which to compute the warming rates and cooling rates, NULL means it will not be computated
  baseline_qt = baseline_qt90,
  baseline_mean = baseline_mean,
  group_by_event = TRUE
) {
  terra::set.names(Values, as.Date(terra::time(Values)))
  # WARNINGS
  if (is.null(start_date) | is.null(end_date)) {
    warning(
      "You didn't specify a begining date and a ending date (see argument 
      'start_date' and 'end_date'), the first date and last date in your time 
      Values SpatRaster will be used."
    )
    start_date <- min(as.Date.character(terra::names(Values)))
    end_date <- max(as.Date.character(terra::names(period_of_interest)))
  }
  if (
    class(start_date) != class(terra::names(Values[[1]])) |
      class(end_date) != class(terra::names(Values[[1]])) |
      class(start_date) != class(end_date)
  ) {
    warning(
      "The date formats are inconsistent between start_date, end_date and 
      Values. Please ensure that the layer names follow a consistent date 
      format. Use class(YourObject) to verify the current format."
    )
  }
  #Check that start date and end_date are within the SpatRasters provided
  if (
    !(start_date %in% terra::names(Values)) |
      !(end_date %in% terra::names(Values))
  ) {
    warning(
      "One or both the specified layers are not present in the SpatRaster 
      containning corrected binarized extreme event."
    )
  }
  Values_extent <- terra::ext(Values)
  #Check that all GPS points are within Values and period_of_interest extent
  if (
    !all(
      abs(GPS[1]) <= abs(Values_extent$xmax) &
        abs(GPS[1]) >= abs(Values_extent$xmin)
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
      abs(GPS[2]) <= abs(Values_extent$ymax) &
        abs(GPS[2]) >= abs(Values_extent$ymin)
    )
  ) {
    warning(
      "At least one latitude coordinate falls outside the extent of the 
      SpatRaster containing the binarized corrected events. Ensure the second
      column of the dataframe contains valid latitude (y) values for analysis."
    )
  }
  # Check that one of the GPS position is not in a pixel that is always an NA
  # (it may indicates that it falls in an area that is not interesting for the
  # study)
  NA_pixels <- which(sapply(Events_corrected, function(df) {
    all(is.na(
      df$Original_value
    ))
  })) # List of pixel that are always NA
  GPS$pixel <- terra::cellFromXY(Values, GPS) # List of the pixels corresponding
  # to the GPS position provided
  if (any(GPS$pixel %in% NA_pixels)) {
    wrong_position <- which(GPS$pixel %in% NA_pixels)
    message("Problematic GPS positions:\n")
    message(utils::capture.output(print(GPS[wrong_position, ])))
  }

  # CODE
  #Extract values for the given GPS position
  df_list <- lapply(GPS$pixel, function(p) {
    Events_corrected[[p]]
  }) # on df per points/pixel
  Values <- t(terra::extract(Values, GPS[, 3]))

  # Subset both dataset so they match the timeframe provided with 'start_date'
  # and 'end_date'.
  Date <- rownames(Values)
  Values <- Values[
    which(
      as.Date(rownames(Values)) >= as.Date(start_date) &
        as.Date(rownames(Values)) <= as.Date(end_date)
    ),
  ]
  Values <- as.matrix(Values)
  df_list <- Map(
    function(df, col_idx) {
      # col_idx is the id
      df$Date <- Date
      df <- df[df$Date >= start_date & df$Date <= end_date, ]
      df$value <- Values[, col_idx] #Merge df_list and Values  # Ad to each
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
  ##Contrary to Events_corrected, Values didn't kept this information of which
  # extreme event were merge together because they were separated by d days or
  # less. Thus, we need to use Events corrected if we want to calculate a metric
  # that describe an event.

  #Create a list of dataframe to store the information using one element (df)
  # per pixel and one row per event
  colnames(GPS) <- c("x", "y", names(GPS[3]))

  #Get daily anomaly to baseline_qt90 and to baseline_mean
  qt90 <- as.data.frame(t(terra::extract(baseline_qt90, GPS[, 3])))
  qt90$dates <- format(
    seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day"),
    "%m-%d"
  )
  mean <- as.data.frame(t(terra::extract(baseline_mean, GPS[, 3])))
  mean$dates <- format(
    seq(as.Date("2024-01-01"), as.Date("2024-12-31"), by = "day"),
    "%m-%d"
  )

  metrics <- lapply(1:length(df_list), function(p) {
    # Go through each pixel
    df <- df_list[[p]]
    df$Date <- as.Date(df$Date)
    df$dates <- format(as.Date(df$Date), "%m-%d")

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
      sub(".*_(\\d{4}-\\d{2}-\\d{2})$", "\\1", df$ID[nrow(df)]) >
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
      dplyr::arrange(ID, Date) |> # put data in pixel order and chronological order
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
        # Get GPS position of each point
        x = GPS$x[p],
        y = GPS$y[p],

        # Get the ID of each event
        event_ID = unique(ID),

        # Duration of each event (assuming Nb_days is constant per ID)
        Nb_days = data.table::first(Nb_days),

        # First and last day of each event
        first_date = as.Date(min(Date)),
        last_date = as.Date(max(Date)),

        # Mean, median, max temperature per event
        mean_value = mean(value),
        median_value = stats::median(value),
        max_value = max(value),

        # Day of absolute maximum temperature
        date_max_value = as.Date(Date[which.max(value)]),

        #Daily rates
        daily_rates = (dplyr::lead(value) - value) /
          as.numeric(dplyr::lead(Date) - Date),

        # Onset rate (up to the absolut maximum value during the event)
        days_onset_abs = as.numeric(date_max_value - first_date),
        ## raw
        raw_onset_rate_abs = ifelse(
          days_onset_abs > 0,
          (max_value - value[which(Date == first_date)]) / days_onset_abs,
          value[1] - data.table::first(last_value_prev_group)
        ),
        ## mean
        mean_onset_rate_abs = ifelse(
          days_onset_abs > 0,
          mean(
            daily_rates[
              Date >= first_date &
                Date <= date_max_value
            ],
            na.rm = TRUE
          ),
          raw_onset_rate_abs
        ),
        ## Median onset rate
        median_onset_rate_abs = ifelse(
          days_onset_abs > 0,
          stats::median(
            daily_rates[
              Date >= first_date &
                Date <= date_max_value
            ],
            na.rm = TRUE
          ),
          raw_onset_rate_abs
        ),
        ## Standard deviation
        sd_onset_rate_abs = ifelse(
          days_onset_abs > 0,
          stats::sd(
            daily_rates[
              Date >= first_date &
                Date <= date_max_value
            ],
            na.rm = TRUE
          ),
          NA
        ),

        # Offset rate (from the absolut maximum value during the event)
        days_offset_abs = as.numeric(last_date - date_max_value),
        ## raw
        raw_offset_rate_abs = ifelse(
          days_offset_abs > 0,
          (value[which(Date == last_date)] - max_value) / days_offset_abs,
          utils::tail(value, 1) - utils::tail(value, 2)[1]
        ),
        ## mean
        mean_offset_rate_abs = ifelse(
          days_offset_abs > 0,
          mean(
            daily_rates[
              Date >= date_max_value &
                Date <= last_date
            ],
            na.rm = TRUE
          ),
          raw_offset_rate_abs
        ),

        median_offset_rate_abs = ifelse(
          days_offset_abs > 0,
          stats::median(
            daily_rates[
              Date >= date_max_value &
                Date <= last_date
            ],
            na.rm = TRUE
          ),
          raw_offset_rate_abs
        ),
        ## Standard deviation
        sd_offset_rate_abs = ifelse(
          days_offset_abs > 0,
          stats::sd(
            daily_rates[
              Date >= date_max_value &
                Date <= last_date
            ],
            na.rm = TRUE
          ),
          NA
        )
      ) |>
      dplyr::ungroup()
    # Add daily baseline_qt90 and daily baseline_mean
    df <- df |>
      dplyr::left_join(
        qt90[, c("dates", paste0("V", as.character(p)))],
        by = "dates"
      ) |>
      dplyr::rename(baseline_qt90 = paste0("V", as.character(p))) |>
      dplyr::left_join(
        mean[, c("dates", paste0("V", as.character(p)))],
        by = "dates"
      ) |>
      dplyr::rename(baseline_mean = paste0("V", as.character(p)))
    # Calculate anomalies
    df <- df |>
      dplyr::mutate(
        anomaly_qt90 = value - baseline_qt90,
        anomaly_mean = value - baseline_mean,
        anomaly_unit = baseline_qt90 - baseline_mean,
        daily_category = dplyr::case_when(
          Cleanned_value == 0 ~ "No extreme event",
          anomaly_qt90 < anomaly_unit ~ "Category I",
          anomaly_qt90 < 2 * anomaly_unit & anomaly_qt90 >= 0 ~ "Category II",
          anomaly_qt90 < 3 * anomaly_unit & anomaly_qt90 >= 0 ~ "Category III",
          anomaly_qt90 >= 3 * anomaly_unit & anomaly_qt90 >= 0 ~ "Category IV",
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
    df <- df |> dplyr::mutate(event_ID = ID)
    max_category <- df |>
      dplyr::filter(!is.na(daily_category)) |>
      dplyr::group_by(event_ID) |>
      dplyr::summarise(max_category = max(daily_category), .groups = 'drop')
    df <- df |>
      dplyr::left_join(max_category, by = "event_ID")

    # Add mean value and standard deviation of anomaly_qt90 and anomaly_mean to
    # the ouputs + add maximal category of each event
    summary_stats <- df |>
      dplyr::group_by(ID) |>
      dplyr::summarise(
        mean_anomaly_qt90 = mean(anomaly_qt90, na.rm = TRUE),
        sd_anomaly_qt90 = stats::sd(anomaly_qt90, na.rm = TRUE),
        max_anomaly_qt90 = max(anomaly_qt90, na.rm = TRUE),
        mean_anomaly_mean = mean(anomaly_mean, na.rm = TRUE),
        sd_anomaly_mean = stats::sd(anomaly_mean, na.rm = TRUE),
        max_anomaly_mean = max(anomaly_mean, na.rm = TRUE),
        max_category = names(sort(table(daily_category), decreasing = TRUE))[1],
        # Catégorie la plus fréquente
        .groups = 'drop'
      )
    df <- df |>
      dplyr::left_join(summary_stats, by = "ID")
    data.table::setnames(
      df,
      old = "max_category.y",
      new = "most_frequent_category"
    )
    data.table::setnames(df, old = "max_category.x", new = "max_category")
    if (group_by_event) {
      df <- df |>
        # Delete daily values that are not usefull to describe the full event
        dplyr::select(
          -baseline_qt90,
          -baseline_mean,
          -anomaly_mean,
          -anomaly_qt90,
          -anomaly_unit,
          -daily_category,
          -date,
          -ID,
          -daily_category,
          -row_num,
          -prev_value,
          -prev_ID,
          -last_value_prev_group,
          -dplyr::starts_with("evolution_rate_lag_"),
          -dplyr::starts_with("variance_value_lag_")
        ) |>
        dplyr::distinct(event_ID, .keep_all = TRUE) # Une seule ligne par ID
    } else {
      df <- df |>
        dplyr::select(
          -date,
          -ID,
          -row_num,
          -prev_value,
          -prev_ID,
          -last_value_prev_group
        )
    }
    return(df)
  })
}

# OTHER PART TO DEVELOP :
# For each pixels :
# - anomalie cumulée par événements
# - catégories
# - onset rate, offset_rate (dans l'événement et depuis une période donnée
# avant et après dont la durée est fixée par l'utilisateur)
#For the all period I want : Frequency, maximum intensity, mean and median
# intensity, sum of anomalies, date of maximum intensity, first date 1, last
# date 1

#The categorie of each event are determined according according Hobday et al.
# 2018 definition

#Category I : btw the 90th percentile and twice the value of the anomalies
# between mean value end 90th percentile. Category II : btw twice and 3 times
# the anomalie Category III : btw 3 times and four time the anomaly
# Category IV : observed value are higher than four times the anomaly btw mean
# and 90th percentile + the mean value.
