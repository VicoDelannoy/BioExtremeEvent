#' Identify the extreme event
#'
#' @description 
#'  Identify the extreme event according to the set of constraints.
#'
#' @param binarized_spatraster :
#'  A Spatraster built using the BEE.id.extreme_days()
#'  function. It only contains 1 and 0.
#' @param n :
#'  The minimum number of consecutive days above the baseline/threshold
#'  required for a series of days more extrem than the baseline to be considered
#'  as an extreme event (≥).
#' @param d :
#'  The maximum number of days for which two series of values above the
#'  threshold/baseline, separated by some days below the threshold/baseline, can
#'  be considered a single extreme event (the maximum distance in days between
#'  the two events to be merged <=).
#' @param nbis :
#'  The minimum number of days in a series of 1 to allow merging with
#'  serie of 1 longer than *'n'* and distant from *'d'* days or fewer.
#'  /!\ *'nbis'* must be inferior or equal to *'n'*.
#'  Example :
#'  If you set *'nbis'* to 3, a short series of 1 (<*'n'*) can be merged with a
#'  long series of 1 (≥*'n'*) if the two series are separated by no more than
#'  *'d'* days, and that the shorter of the two series contains at least *'nbis'*
#'  days above the threshold.
#'  Default: NULL
#'
#' @param w :
#'  NOT AVAILBALE YET Window w between to event in which you
#'  allow a certain proportion p of days (within the window) to be bellow the
#'  threshold, windows is bordered by series of consecutive days above
#'  threshold. w and p replace d. If w is an odd number, the minimum number of
#'  days above threshold within w will be the value rounded to the superior
#'  round number. For instance if w = 7 and p = 0.3 a window of 7 days between
#'  to event of n consecutive days will not create two distinct event if at
#'  least 4 days are above threshold (wp = 7*0.5 = 3.5 rounded to 4).
#'  Example : n = 5; w = 10; p = 0.5 -> w * p = 5 . Default: NULL.
#'
#' @param p :
#'  NOT AVAILBALE YET, see "w" explanations. Default: NULL
#'
#' @return 
#'  Returns a list with 2 elements. The first one is called
#'  **"binarized_spatraster_corrected"**. This is a Spatraster containing 
#'  binarised values that have been corrected according to the definition used
#'  in the function. It is used in the next funciton of the pipeline.
#'  The second element is called **"Event_corrected"**. This is a list with one
#'  data.table per pixel.
#'  Each data.table has a date for each row (and as many dates as there are 
#'  layers in yourspatraster). The "original_value" column gives the value of
#'  the raw binarisation before applying this function. "cleaned_value" provides
#'  the corrected values for each date. "ID" is the identifier of the event,
#'  built using the pixel number, the first day of the event and the last day of
#'  the event. It is proper to each pixel. "Duration" indicates the duration (in
#'  days) of the event. "Date" indicates the date corresponding to the row.
#'
#' @details
#'  Setting *'nbis'* = 1, has the same effect as setting nbis= NULL. When
#'  *'nbis'* = NULL any series of 1 distant from less than *'d'* to a series of
#'  1 at least as long as *'n'* days will be merged into one big series event,
#'  even if the serie last only one day. If you provide an argument for 
#'  *'nbis'*, only series at least as long as *'nbis'* can be merge to serie as
#'  long as *'n'* (if the two are distant from less than *'d'* days).
#'  Days above threshold that don't qualify has "extreme event" (e.g. to short
#'  series or to isolated) are transformed to 0. The second output of the
#'  function provided for each pixel a dataframe containning daily values before
#'  and after corrections.
#'  *"w"* and *"p"* settings are not available yet. BEE.id.extreme_events is not
#'  designed to work on 4D data (time + spatial 3D).
#'
#' @examples
#' # This function apply a filter to withdraw isolated 1 and 0. For instance,
#' # is you consider that an event that last only 2 days should not be accounted
#' # for in future analysis, you can correct those value by 0. Here are the
#' # different filter you can apply:
#'
#' ## Low complexity:
#' # - Minimum number of consecutive day above threshold to
#' #  consider that there is an extreme event: n  (>=)
#' #
#' # -   Maximum number of days bellow thershold to consider two series of
#' # value above threshold as a single extreme event (distance max in days btw
#' # 2 events to merge them): d (<=)
#' #
#' # Examples: n >= 5  d>= 2
#' ### Load the example dataset in R environement:
#' file_name <- system.file(file.path("extdata", "copernicus_binarized.tiff"),
#'                                   package = "BioExtremeEvent")
#' binarized <- terra::rast(file_name)
#' ### Apply filter:
#' binarized_EE <- BEE.id.extreme_events(binarized_spatraster = binarized,
#' n = 5,
#' d = 2)
#'
#' ## Medium complexity :
#' # -   Minimum number of consecutive day above threshold to
#' #  consider that there is an extreme event: n  (>=)
#' # -   Maximum number of days bellow thershold to consider two series of
#' # value above threshold as a single extreme event (distance max in days btw
#' # 2 events to merge them): d (<=)
#' # -   Minimum number of consecutive day above threshold to allow merging to
#' # an event at least as long as n. (nbis must be inferior or equal to n)
#' #
#' # Examples: n >= 5  d>= 2 nbis>=3
#' ### Load the example dataset in R environement:
#' file_name <- system.file(file.path("extdata", "copernicus_binarized.tiff"),
#'                                   package = "BioExtremeEvent")
#' binarized <- terra::rast(file_name)
#' ### Apply filter:
#' binarized_EE <- BEE.id.extreme_events(binarized_spatraster = binarized,
#' n = 5,
#' d = 2,
#' nbis = 3)
#'
#' @export
#'
#-------------------------------------------------------------------------------

# n = 5; d= 3 ; nbis= 2 # for medium complexity
BEE.id.extreme_events <- function(
  binarized_spatraster,
  n,
  d,
  nbis = NULL,
  w = NULL,
  p = NULL
) {
  pixel_time_series <- terra::extract(
    binarized_spatraster,
    #one row per pixel, 1 column per time step
    1:terra::ncell(binarized_spatraster)
  )
  # We need each line of pixel_time_series in a numeric format but the
  # convertion inside the fonction correct_lowcomplexity_n_d consumes a
  # lot of RAM, thus we are converting the data frame to a numeric matrix,
  # conversion everything at one is more efficient than row by row in the
  # function ;) :
  pixel_time_series <- as.matrix(pixel_time_series)
  storage.mode(pixel_time_series) <- "numeric" # NB : if this still to much
  # memory we could swich to sparseMatrix but then, the manipulation on
  # that matrix may be slower.

  all_dates <- as.Date(terra::time(binarized_spatraster))

  if (
    !is.null(n) &
      !is.null(d) &
      is.null(nbis) &
      is.null(w) &
      is.null(p)
  ) {
    Event_corrected <- list()
    #We want to avoid processing a pixel that is always NA :
    indices_all_na <- which(rowSums(!is.na(pixel_time_series)) == 0)
    indices_to_do <- which(rowSums(!is.na(pixel_time_series)) != 0)

    Event_corrected <- vector("list", length(indices_to_do))
    Event_corrected <- lapply(indices_to_do, function(c) {
      correct_lowcomplexity_n_d(
        pixel_time_series,
        c,
        n,
        d,
        all_dates
      )
    })
  }

  if (
    !is.null(n) &
      !is.null(d) &
      !is.null(nbis) &
      is.null(w) &
      is.null(p)
  ) {
    Event_corrected <- list()
    #We want to avoid processing a pixel that is always NA :
    indices_all_na <- which(rowSums(!is.na(pixel_time_series)) == 0)
    indices_to_do <- which(rowSums(!is.na(pixel_time_series)) != 0)

    Event_corrected <- vector("list", length(indices_to_do))
    Event_corrected <- lapply(indices_to_do, function(c) {
      correct_mediumcomplexity_n_d_nbis(
        pixel_time_series,
        c,
        n,
        d,
        nbis,
        all_dates
      )
    })
  }
  if (
    !is.null(n) &
      is.null(d) &
      is.null(nbis) &
      !is.null(w) &
      !is.null(p)
  ) {
    #here add code for high complexity
  }
  ### Modify the spatraster :
  # Build a matrice pixel x time using cleanned_value column from Event_corrected
  n_pixels <- terra::ncell(binarized_spatraster)
  n_dates <- length(all_dates)
  # empty matrice :
  mat <- matrix(NA_real_, nrow = n_pixels, ncol = n_dates)
  # add corrected values from Event_correct (column cleanned value of each dt)
  for (p in seq(1, length(indices_to_do), 1)) {
    #1.4 s 419 MB
    # no correction needed for pixels that are always NA
    mat[indices_to_do[p], ] <- Event_corrected[[p]]$cleanned_value
  }
  binarized_spatraster_corrected <- binarized_spatraster
  terra::values(binarized_spatraster_corrected) <- mat
  rm(mat)
  gc()

  ### Convert data.table list to a dataframe list with empty position for NA
  # pixels:
  list_df <- vector("list", n_pixels)

  list_df[indices_to_do] <- lapply(
    seq(1, length(indices_to_do), 1),
    function(i) {
      dt <- Event_corrected[[i]]
      as.data.frame(dt)
    }
  )

  return(
    list(
      binarized_spatraster_corrected = binarized_spatraster_corrected,
      Event_corrected = list_df
    )
  )
}


#' ###############  Correct the dt list according "low complexity" framework
#'
#' @noRd
correct_lowcomplexity_n_d <- function(pixel_time_series, c, n, d, all_dates) {
  # (go through each pixels of pixel_time_series)
  # c=675
  pixel_values <- pixel_time_series[c, ]

  rle_series <- rle(pixel_values)
  # Clean isolated '1'
  One_to_0 <- which(
    rle_series$values == 1 &
      rle_series$lengths < n
  ) ##list isolated '1'
  ## (series < n ) by 0
  rle_series$values[One_to_0] <- 0 # replace isolated 1 by 0
  pixel_values_cleanned <- inverse.rle(rle_series) # values of each
  ## pixel trough time
  rle_series_cleanned <- rle(pixel_values_cleanned) # re-identify
  ## the series of 0, this way series of zeros are re-defined to include 0
  ## and 1 transformed to 0 in one same serie when they are next to each
  ## other.
  # Clean isolated '0'
  Zero_to_1 <- which(
    rle_series_cleanned$values == 0 &
      rle_series_cleanned$lengths <= d
  ) ##list isolated '1'
  ## (series < n ) by 1
  rle_series_cleanned$values[Zero_to_1] <- 1 # replace isolated 0 by 1
  pixel_values_cleanned <- inverse.rle(rle_series_cleanned) # values of each
  ## pixel trough time
  rle_series_cleanned <- rle(pixel_values_cleanned) # re-identify
  ## the series of 1, this way series of ones are re-defined to include 1
  ## and 0 transformed to 1 in one same serie when they are next to each
  ## other.

  pixel <- local({
    dt <- data.table::data.table(
      pixel_id = rep(c, length(pixel_values)),
      original_value = as.vector(pixel_values),
      #need to go through a dataframe because I don't know how to modify the
      # event_id directly in rle object
      cleanned_value = unlist(pixel_values_cleanned),
      event_id = rep(
        seq_along(rle_series_cleanned$lengths),
        rle_series_cleanned$lengths
      ),
      duration = rep(
        rle_series_cleanned$lengths,
        rle_series_cleanned$lengths
      )
    )
    dt
  })

  # add dates:
  pixel$date <- sort(all_dates)

  # Get the first and last date of each extreme event
  first_date <- stats::ave(pixel$date, pixel$event_id, FUN = function(x) {
    min(x, na.rm = TRUE)
  })
  last_date <- stats::ave(pixel$date, pixel$event_id, FUN = function(x) {
    max(x, na.rm = TRUE)
  })

  # Create new extreme event id has pixel_first_day_last_day
  pixel$ID <- paste0(
    pixel$pixel_id,
    "_",
    first_date,
    "_",
    last_date,
    "_",
    pixel$event_id
  )

  return(pixel)
}

#' #########  Correct the dt list according "medium complexity" framework
#'
#' @noRd
correct_mediumcomplexity_n_d_nbis <- function(
  pixel_time_series,
  c,
  n,
  d,
  nbis,
  all_dates
) {
  # (go through each pixels of pixel_time_series)
  # c=675
  pixel_values <- pixel_time_series[c, ]

  rle_series <- rle(pixel_values)

  ## 1) Identify series of 1 with a length superior or equal to n
  series_1_n <- which(rle_series$values == 1 & rle_series$lengths >= n)
  ## 2) Identify series of 1 with a length superior or equal to nbis and inferior
  #  to n that are distant from a serie identified in 1) from d or less day
  series_1_nbis_potentials <- which(
    rle_series$values == 1 & rle_series$lengths < n & rle_series$lengths >= nbis
  )
  series_1_nbis <- series_1_nbis_potentials[which(
    rle_series$lengths[as.integer(series_1_nbis_potentials) - 1] <= d &
      as.integer(series_1_nbis_potentials - 2) %in% as.integer(series_1_n) |
      rle_series$lengths[as.integer(series_1_nbis_potentials) + 1] <= d &
        as.integer(series_1_nbis_potentials + 2) %in% as.integer(series_1_n)
  )] #before a serie of 1 there is always a serie of 0

  # Clean isolated '1' (series < n), with "lowcomplexity" all of those should
  # disappear
  One_to_0 <- which(
    rle_series$values == 1 &
      rle_series$lengths < n
  )
  # 3) With "mediumcomplexity" we want to save the series of ones with a length
  # btw nbis and n if they are close (<d) to a serie of ones longer or equal to
  # n:
  One_to_0 <- One_to_0[which(!(One_to_0 %in% series_1_nbis))] #save the series
  # identified above.
  # 4) correct the VERY isolated ones:
  rle_series$values[One_to_0] <- 0 # replace isolated 1 by 0
  pixel_values_cleanned <- inverse.rle(rle_series) # values of each
  ## pixel trough time
  rle_series_cleanned <- rle(pixel_values_cleanned) # re-identify
  ## the series of 0, this way series of zeros are re-defined to include 0
  ## and 1 transformed to 0 in one same serie when they are next to each
  ## other.

  # Clean isolated '0'
  Zero_to_1 <- which(
    rle_series_cleanned$values == 0 &
      rle_series_cleanned$lengths <= d
  )

  rle_series_cleanned$values[Zero_to_1] <- 1 # replace isolated 0 by 1
  pixel_values_cleanned <- inverse.rle(rle_series_cleanned) # values of each
  ## pixel trough time
  rle_series_cleanned <- rle(pixel_values_cleanned) # re-identify
  ## the series of 1, this way series of ones are re-defined to include 1
  ## and 0 transformed to 1 in one same serie when they are next to each
  ## other.

  pixel <- local({
    dt <- data.table::data.table(
      pixel_id = rep(c, length(pixel_values)),
      original_value = as.vector(pixel_values),
      #need to go through a dataframe because I don't know how to modify the
      # event_id directly in rle object
      cleanned_value = unlist(pixel_values_cleanned),
      event_id = rep(
        seq_along(rle_series_cleanned$lengths),
        rle_series_cleanned$lengths
      ),
      duration = rep(
        rle_series_cleanned$lengths,
        rle_series_cleanned$lengths
      )
    )
    dt
  })

  # add dates:
  pixel$date <- sort(all_dates)

  # Get the first and last date of each extreme event
  first_date <- stats::ave(pixel$date, pixel$event_id, FUN = function(x) {
    min(x, na.rm = TRUE)
  })
  last_date <- stats::ave(pixel$date, pixel$event_id, FUN = function(x) {
    max(x, na.rm = TRUE)
  })

  # Create new extreme event id has pixel_first_day_last_day
  pixel$ID <- paste0(
    pixel$pixel_id,
    "_",
    first_date,
    "_",
    last_date,
    "_",
    pixel$event_id
  )

  return(pixel)
}

#' ################# Apply corrections to the binarize spatraster
#'
#' @noRd

correct_raster <- function(raster, t, Event_corrected, indices_to_do) {
  # t=12625 ; raster <- binarized_spatraster[[t]]

  # get values of each layers (taking in account the NA pixels that are no
  # longer in Event_corrected)
  ## Create a vector full of NA
  n <- terra::ncell(raster)
  full_values <- rep(NA_real_, n)
  ## Get corrected values of non NA pixels
  Matching_cleaning_values <- sapply(Event_corrected, function(event) {
    event[t]$cleanned_value
  })
  ## Add non NA pixels to the NA vector
  full_values[indices_to_do] <- Matching_cleaning_values

  # replace raster values by corrected values
  terra::values(raster) <- full_values
  return(raster)
}

#  High complexity : NOT AVAILABLE YET
# - Unstable window * w * between to event in which you allow a certain
#  proportion p of days of w to be bellow threshold if this windows is
# bordered by series of consecutive days above threshold on each size.
# w  and  p  replace d. If w is an odd number, the minimum number of days
# above threshold within w will be the value rounded to the superior round
#  number. For instance if w = 7 and p = 0.3 a window of 7 days between to
#  event of * n * consecutive days will not create two distinct event if at
#  least 4 days are above threshold.
#(w * \ * * p * = 7 \ * 0.5 = 3.5 rounded to 4)
#
# Examples : n = 5
# w = 10
# p = 0.5  ->  w*p = 5
# case that lead to a single big event of 15 days:
# 0 0 0 0 0 1 1 1 1 1 0 1 1 1 0 0 0 1 1 0 1 1 1 1 1 0 0 0 0 0
# 0 0 0 0 0 1 1 1 1 1 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
# (you can also use low complexity n = 5 and d = 3)
# 0 0 0 0 0 1 1 1 1 1 0 1 0 1 0 1 0 1 0 1 1 1 1 1 1 0 0 0 0 0
# case that lead to 2 distinct events:
# 0 0 0 0 0 1 1 1 1 1 0 1 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0 0 0 0
# 0 0 0 0 0 1 1 1 1 1 0 0 0 1 1 1 1 0 0 0 1 1 1 1 1 0 0 0 0 0
# 0 0 0 0 0 1 1 1 1 1 0 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 (
# if you want to detect 1 event of 10 days then a distinct event of 5 days
# instead of 2 events of 5 days only you can use :
# low complexity with n = 4 and d = 1)
