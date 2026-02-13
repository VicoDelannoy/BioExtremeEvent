#' Calculate baseline value for each date and pixel.
#'
#' @description
#'  Calculate the daily baseline value using the mean or a quantile of a
#'  timeserie.
#'
#' @param yourspatraster :
#'  The SpatRaster that contains the values to be used to calculate the baseline
#'  (your reference time serie). Each layer must have a date and there must be
#'  **no duplicates**.
#' @param start_date :
#'  The date to use as the first day of your reference period. It must be in the
#'  format **YYYY-MM-DD**. (Note that the provided dataset can cover a longer
#'  time period than the one you want to use as a reference.)
#' @param end_date :
#'  The date to use as the first day of your reference period. It must be in the
#'  format **YYYY-MM-DD**.
#' @param threshold :
#'  Specifies whether you want to use a percentile or the mean of the observed
#'  values as the threshold.
#'  The accepted arguments are *"qt"*, to specify that you want to compute a
#'  percentile of the observed values in the reference timeframe, or *"mean"*,
#'  if you want to use the mean of the observed values as threshold.
#'  If you want to use a fixed value as the threshold (e.g. biological optimum),
#'  you can skip this step and use BEE.id.extreme_day() directly.
#' @param quantile_value :
#'  Indicates the desired percentile value. This must be between 0 and 1.
#' @param time_window :
#'  Number of days on either side of day 'd' that are used to calculate the
#'  threshold value associated for day 'd'.
#'  For example, if *"time_window = 5"*, the thershold value for 'd' day will
#'  be calculated using data from five days before day d, day d itself,
#'  and the five days after day d, from all years between *"start_date"* and
#'  *"end_date"*.
#' @param smooth_window :
#'  Number of days on either side of day 'd' that are used to calculate the mean
#'  value of the threshold assigned to the 'd' day. This allows your threshold
#'  value to be smoothed across days.
#'  For example, if *"smooth_window"* is set to 10 days, the final value on
#'  day 'd' will be the mean of the baseline/thershold values calculated for
#'  days: d - 10 to d + 10 (eleven values).
#'
#' @return
#'  A SpatRaster with one day of the year per layer (366 layers) that has
#'  the same extent, pixel resolution and CRS as the provided SpatRaster. Each
#'  pixel contains the threshold value at which an extreme event is identified
#'  (the baseline value).
#'
#' @examples
#' ### Load the example dataset in R environement :
#' file_name <- system.file(file.path("extdata", "copernicus_data_celsius.tiff"),
#'                                   package = "BioExtremeEvent")
#' copernicus_data_celsius <- terra::rast(file_name)
#'
#' ### **90th percentile threshold:** 90 % of the values observed at a given day
#' # accross all years of the reference time period are bellow the baseline
#' # value.
#' #### No smoothing: # 16s to run
#' baseline_qt90_smth_0 <- BEE.calc.baseline(yourspatraster = copernicus_data_celsius,
#'                                    start_date = "2023-01-01",
#'                                    end_date = "2025-12-31",
#'                                    threshold = "qt",
#'                                    quantile_value = 0.9,
#'                                    time_window = 5,
#'                                    smooth_window = 0)
#' #### Smoothing over 15 days: # 18s to run
#' baseline_qt90_smth_15 <- BEE.calc.baseline(yourspatraster = copernicus_data_celsius,
#'                                    start_date = "2023-01-01",
#'                                    end_date = "2025-12-31",
#'                                    threshold = "qt",
#'                                    quantile_value = 0.9,
#'                                    time_window = 5,
#'                                    smooth_window = 7) #7+1+7=15
#' ###**Mean value threshold:**
#' baseline_mean <- BEE.calc.baseline(yourspatraster = copernicus_data_celsius,
#'                                    start_date = "2023-01-01",
#'                                    end_date = "2025-12-31",
#'                                    threshold = "mean",
#'                                    time_window = 5,
#'                                    smooth_window = 7)
#' #if threshold = "mean", any quantile_value provided will be ignored.
#'
#' @export
#'
#-------------------------------------------------------------------------------
# For tests : yourspatraster = ds; start_date = "1982-01-01"
# end_date = "2010-12-31"; threshold = "qt"; quantile_value = 0.9
# time_window = 5; smooth_window = 7

BEE.calc.baseline <- function(
  yourspatraster,
  start_date = NULL,
  end_date = NULL,
  threshold,
  quantile_value = NULL,
  time_window = 5,
  smooth_window = 10
) {
  # Subset the dataset to the chosen start_date and end-date if necessary
  if (is.null(start_date)) {
    start_date <- terra::names(yourspatraster)[1]
    message(
      "You have provided an 'end_date' but no 'start_date', first ",
      "date of the dataset (",
      start_date,
      ") have been used as first ",
      "day of the baseline."
    )
  }
  if (is.null(end_date)) {
    end_date <- terra::names(yourspatraster)[terra::nlyr(yourspatraster)]
    message(
      "You have provided a 'start_date' but no 'end_date', last date of ",
      "the dataset (",
      end_date,
      ") have been used as last day of the baseline."
    )
  }
  # Cut yourspatraster to only keep the part between start_date and end_date
  dates <- as.Date(terra::time(yourspatraster))
  if (!is.null(start_date) & !is.null(end_date)) {
    if (!(start_date %in% dates) | !(end_date %in% dates)) {
      warnings(
        "start_date and/or end_date are not among yourspatraster dates.
    The provided dates may not be within the specified time series or may be in
    a different format. Please ensure that the start_date and end_date fields on
    your Sparster account are in the same format and time period."
      )
    }
  }

  to_keep <- which(
    dates >= as.Date(start_date) &
      dates <= as.Date(end_date)
  )
  yourspatraster <- yourspatraster[[to_keep]]
  # Function to detect leap years
  is_leap_year <- function(years) {
    return((years %% 4 == 0 & years %% 100 != 0) | (years %% 400 == 0))
  }
  dates <- as.Date(terra::time(yourspatraster)) # Dates
  years <- as.integer(format(dates, "%Y")) # Years
  # Build a dataframe with the layer number, the doy and the corrected doy
  # (366, 60 don't exist in non leap year)
  df <- data.frame(
    date = dates,
    N_layer = seq(1, terra::nlyr(yourspatraster), 1),
    leap_year = is_leap_year((years)),
    doy = as.integer(format(dates, "%j"))
  ) #day of the year (with correction taking in account leapyear)
  df$doy2 <- ifelse(
    df$doy >= 60 &
      !df$leap_year,
    df$doy + 1,
    df$doy
  )
  # For each day of the year from 0 to 366, get a list of the layers to use to
  # calculate baseline
  doy_indices <- vector("list", 366)
  for (j in 1:366) {
    days_in_window <- ((j - time_window):(j + time_window)) %% 366 # %% 366 is
    #a 'modulo' which mean that when a number reach 366 it is converted to 0 and
    # the series that over from 0, 1, 2...
    days_in_window[days_in_window == 0] <- 366 # to avoid 0
    doy_indices[[j]] <- df$N_layer[which(df$doy2 %in% days_in_window)]
  }

  # Calculate baseline according the chosen methodology
  if (threshold == "qt") {
    baseline <- lapply(1:366, function(d) {
      terra::app(yourspatraster[[doy_indices[[d]]]], fun = function(v) {
        stats::quantile(v, probs = quantile_value, na.rm = TRUE)
      })
    })
  } else if (threshold == "mean") {
    # if the baseline chosen is a mean value
    baseline <- lapply(1:366, function(d) {
      terra::app(
        yourspatraster[[doy_indices[[d]]]],
        fun = function(v) {
          mean(v, na.rm = TRUE)
        }
      )
    })
  }
  baseline <- terra::rast(baseline)
  # Smooth the value with an 11 days moving average
  smooth_indices <- vector("list", 366)
  for (s in 1:366) {
    days_in_window <- ((s - smooth_window):(s + smooth_window)) %% 366
    days_in_window[days_in_window == 0] <- 366
    smooth_indices[[s]] <- days_in_window
  }
  baseline <- lapply(1:366, function(i) {
    terra::app(baseline[[smooth_indices[[i]]]], fun = mean, na.rm = TRUE)
  })
  baseline <- terra::rast(baseline)
  return(baseline)
}
