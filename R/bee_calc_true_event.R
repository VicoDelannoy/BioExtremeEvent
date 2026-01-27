#' Identify the extreme event
#'
#' @description Identify the extreme event according to the set of constraints.
#'
#' @param binarized_EE A list of Spatraster built using the BEE.calc.binarize
#'  function.
#' @param n The minimum number of consecutive days above the baseline/threshold
#' required for a series of "1" to be considered as an extreme event (≥).
#' @param d The maximum number of days below the baseline/thershold for
#' considering two series of values above the threshold/baseline but separated
#' by some days below threshold/baseline as a single extreme event (maximum
#' distance in days between the two events to be merged, ≤).
#' @param nbis NOT AVAILBALE YET The minimum number of days in each series of
#' 1 to allow merging if there are d days or fewer between the to series.
#' Example :
#'   if you set 'nbis' to 3, short series of 1 (n) can be merged with long
#'   series of 1 (≥n) if the two series are separated by no more than d days and
#'   there are at least nbis days between them (with at least one of the series
#'   having to be longer than n).
#'   Default: NULL
#' @param w NOT AVAILBALE YET Window w between to event in which you
#'  allow a certain proportion p of days (within the window) to be bellow the
#'  threshold, windows is bordered by series of consecutive days above
#'  threshold. w and p replace d. If w is an odd number, the minimum number of
#'  days above threshold within w will be the value rounded to the superior
#'  round number. For instance if w = 7 and p = 0.3 a window of 7 days between
#'  to event of n consecutive days will not create two distinct event if at
#'  least 4 days are above threshold (wp = 7*0.5 = 3.5 rounded to 4).
#'  Example : n = 5; w = 10; p = 0.5 -> w * p = 5 . Default: NULL.
#'
#' @param p NOT AVAILBALE YET, see "w" explanations. Default: NULL
#'
#' @return Returns a list with 2 elements, first one is called
#' "stacked_rasters_corrected". This is a Spatraster containing binarised values
#'  corrected according the definition used in the function. The second element
#'   is called "Event_corrected". This is a list with one data.table per pixel.
#'   Each data.table as a date per rows (and as many date as there are layers in
#'    YourSpatraster). The columns "original_value" gives the value of the raw
#'     binarisation, before appling the definition. "cleanned_value" gives for
#'      each dates the values after correction. "event_ID" is an event_ID of the event,
#'      built using the pixel event_ID, first day of the event and last day of the
#'      event. "duration" indicates the duration (in unit of time) of the event.
#'      "date" indicates the date corresponding to the row.
#' @details w and p settings are not available yet. BEE.calc.true_event is not
#'  designed to work on 4D data (time + spatial 3D).
#' @examples
#' # This function apply a filter to withdraw isolated events. For instance,
#' # is you consider that an event that last only 2 days should not be accounted
#' # for in future analysis, you can correct those value by 0. Here are the
#' # different filter you can apply:
#' #
#' # Low complexity:
#' # - Minimum number of consecutive day above threshold to
#' #  consider that there is an extreme event: n  (>=)
#' #
#' # -   Maximum number of days bellow thershold to consider two series of
#' # value above threshold as a single extreme event (distance max in days btw
#' # 2 events to merge them): d (<=)
#' #
#' # Examples: n >= 5  d>= 2
#' #       0 0 0 1 1 1 1 ** 1 ** 0 0 0 = 1 event,
#' #   but 0 0 0 1 1 1 1 0 0 0 = no event
#' #       0 0 0 1 1 1 1 1 0 0 1 1 1 0 0 0 = 1 event of 5 days,
#' #   but 0 0 0 1 1 1 1 1 0 0 1 1 1 **1 1 ** 0 0 0 = 1 event of 12 days,
#' #   but 0 0 0 1 1 1 1 1 0 0 ** 0 ** 1 1 1 1 1 0 0 0 = 2 instincts
#' # event of 5 days each.
#' #
#' # Medium complexity (Hobday et al. 2016, 2018):
#' #- Minimum number of days distant from d days or less to an event
#' #  of n days : nbis .
#' #
#' # Example: n = 5, d = 2, nbis = 3
#' # This is for the case you want to consider
#' #     0 0 0 1 1 1 1 1 0 0 1 1 ** 1 ** 0 0 0 as one event of 10 days
#' # but 0 0 0 1 1 1 1 1 0 0 1 1 0 0 0 as one event of 5 days.
#' # In the first case, if there are not five consecutive days above threshold
#' # right before or after (d = n still fixed to 5).
#' #
#' #  High complexity : NOT AVAILABLE YET
#' # - Unstable window * w * between to event in which you allow a certain
#' #  proportion p of days of w to be bellow threshold if this windows is
#' # bordered by series of consecutive days above threshold on each size.
#' # w  and  p  replace d. If w is an odd number, the minimum number of days
#' # above threshold within w will be the value rounded to the superior round
#' #  number. For instance if w = 7 and p = 0.3 a window of 7 days between to
#' #  event of * n * consecutive days will not create two distinct event if at
#' #  least 4 days are above threshold.
#' #(w * \ * * p * = 7 \ * 0.5 = 3.5 rounded to 4)
#' #
#' # Examples : n = 5
#' # w = 10
#' # p = 0.5  ->  w*p = 5
#' # case that lead to a single big event of 15 days:
#' # 0 0 0 0 0 1 1 1 1 1 0 1 1 1 0 0 0 1 1 0 1 1 1 1 1 0 0 0 0 0
#' # 0 0 0 0 0 1 1 1 1 1 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0
#' # (you can also use low complexity n = 5 and d = 3)
#' # 0 0 0 0 0 1 1 1 1 1 0 1 0 1 0 1 0 1 0 1 1 1 1 1 1 0 0 0 0 0
#' # case that lead to 2 distinct events:
#' # 0 0 0 0 0 1 1 1 1 1 0 1 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0 0 0 0
#' # 0 0 0 0 0 1 1 1 1 1 0 0 0 1 1 1 1 0 0 0 1 1 1 1 1 0 0 0 0 0
#' # 0 0 0 0 0 1 1 1 1 1 0 1 1 1 1 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 (
#' # if you want to detect 1 event of 10 days then a distinct event of 5 days
#' # instead of 2 events of 5 days only you can use :
#' # low complexity with n = 4 and d = 1)
#'
#' @export
#'
#-------------------------------------------------------------------------------

# n = 5; d= 3 ; nbis= 2 # for medium complexity
BEE.calc.true_event <- function(
  binarized_EE,
  n,
  d,
  nbis = NULL,
  w = NULL,
  p = NULL
) {
  # Extract all rasters (If we have to think to a faster way we could filter so
  # pixel that are always NA (land) are not processed and we could create
  # "all_dates" starting from scratch, only knowing first day and last day and
  # assuming there is no missing day)
  all_data <- lapply(seq_along(binarized_EE), function(i) {
    j_indices <- seq(1, length(names(binarized_EE[[i]])), 1)
    list(
      dates = names(binarized_EE[[i]][j_indices]),
      rasters = binarized_EE[[i]][[j_indices]]
    )
  })
  # Merge data to save memory
  all_dates <- unlist(lapply(all_data, `[[`, "dates"), recursive = FALSE)
  all_rasters <- terra::rast(lapply(all_data, `[[`, "rasters"))
  all_rasters_names <- unlist(
    lapply(all_data, `[[`, "dates"),
    recursive = FALSE
  )
  # Clean memory
  rm(all_data)
  gc()

  # Reorder informations in chronological order :
  ## get a list of the raster layers
  all_rasters <- lapply(1:terra::nlyr(all_rasters), function(i) {
    all_rasters[[i]]
  })
  ## Transforme SpatRaster to list of rasters so it is easier to modify the
  # order.
  sorted_indices <- order(all_rasters_names)
  ## get a list of the current position of the rasters : the first number tell
  # you where is currently the raster that should be first in chronological
  # order
  sorted_rasters <- all_rasters[sorted_indices] #sort the rasters
  # Transform the list of sorted raster into a multilayers raster to be able to
  # use terra:extract on cells
  stacked_rasters <- terra::rast(sorted_rasters)
  rm(all_rasters, all_rasters_names, sorted_indices, sorted_rasters)
  gc()
  pixel_time_series <- terra::extract(
    stacked_rasters,
    #one row per pixel, 1 column per time step
    1:terra::ncell(stacked_rasters)
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

  if (
    !is.null(n) &
      !is.null(d) &
      is.null(nbis) &
      is.null(w) &
      is.null(p)
  ) {
    gc()
    Event_corrected <- list()
    #We want to avoid processing a pixel that is always NA :
    indices_all_na <- which(rowSums(!is.na(pixel_time_series)) == 0)
    indices_to_do <- which(rowSums(!is.na(pixel_time_series)) != 0)

    Event_corrected <- vector("list", length(indices_to_do))
    Event_corrected <- lapply(indices_to_do, function(c) {
      correct_lowcomplexity_n_d(pixel_time_series, c, n, d, all_dates)
    })
  }

  if (
    !is.null(n) &
      !is.null(d) &
      !is.null(nbis) &
      is.null(w) &
      is.null(p)
  ) {
    gc()
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
  # Aply corrections to raster of 1 and 0
  stacked_rasters_corrected <- lapply(
    1:terra::nlyr(stacked_rasters),
    function(t) {
      #stacked raster are rasters in the same order than in Event_corrected but
      # before modifications
      correct_raster(stacked_rasters[[t]], t, Event_corrected, indices_to_do)
    }
  )

  stacked_rasters_corrected <- terra::rast(stacked_rasters_corrected)
  names(stacked_rasters_corrected) <- sort(all_dates)

  return(
    list(
      stacked_rasters_corrected = stacked_rasters_corrected,
      Event_corrected = Event_corrected
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
  rle_series_cleanned <- rle(unlist(pixel_values_cleanned)) # re-identify
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
  rle_series_cleanned <- rle(unlist(pixel_values_cleanned)) # re-identify
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
  rle_series_cleanned <- rle(unlist(pixel_values_cleanned)) # re-identify
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
  rle_series_cleanned <- rle(unlist(pixel_values_cleanned)) # re-identify
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
  # t=1 ; raster <- stacked_rasters[[t]]

  # get values of each layers (taking in account the NA pixels that are no
  # longer in Event_corrected)
  ## Create a vector full of NA
  n <- terra::ncell(raster)
  full_values <- rep(NA_real_, n)
  ## Get corrected values of non NA pixels
  Matching_cleanning_values <- sapply(Event_corrected, function(event) {
    event[t]$cleanned_value
  })
  Matching_cleaning_values <- vapply(
    Event_corrected,
    function(event) event$cleanned_value[t],
    numeric(1)
  )
  ## Add non NA pixels to the NA vector
  full_values[indices_to_do] <- Matching_cleaning_values

  # replace raster values by corrected values
  terra::values(raster) <- full_values
  return(raster)
}
