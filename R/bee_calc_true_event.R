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
#'    YourSpatraster). The columns "Original_value" gives the value of the raw
#'     binarisation, before appling the definition. "Cleanned_value" gives for
#'      each dates the values after correction. "ID" is an ID of the event,
#'      built using the pixel ID, first day of the event and last day of the
#'      event. "Nb_days" indicates the duration (in unit of time) of the event.
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

# n = 5, d= 3
BEE.calc.true_event <- function(
  binarized_EE,
  n,
  d = NULL,
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

  # start_date <- min(unlist(all_dates), na.rm = TRUE) # for later when I'll
  # give the possibility to only work on a subset of the time période
  # end_date <- max(all_dates, na.rm = TRUE)

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
  # Transform to character the series of value for each pixels

  if (
    !is.null(n) &
      !is.null(d) &
      is.null(nbis) &
      is.null(w) &
      is.null(p)
  ) {
    # code for low complexity

    correct_lowcomplexity_n_d <- function(pixel_time_series, c, n, d, pix) {
      # pixel_values <- pixel_time_series[1,] and pix=1; "pix" is a dynamic
      # argument that is automatically indexed when calling the function,
      # it identify each pixels (~ each rows of pixel_time_series)

      pixel_values <- pixel_time_series[c, ]
      rle_series <- rle(pixel_values)
      # Clean isolated '1'
      One_to_0 <- which(
        rle_series$values == 1 &
          rle_series$lengths < n
      ) #list isolated '1'
      # (series < n ) by 0
      rle_series$values[One_to_0] <- 0 # replace isolated 1 by 0
      pixel_values_cleanned <- inverse.rle(rle_series) # values of each
      # pixel trough time
      rle_series_cleanned <- rle(unlist(pixel_values_cleanned)) # re-identify
      # the series of 0, this way series of zeros are re-defined to include 0
      # and 1 transformed to 0 in one same serie when they are next to each
      # other.
      pixel <- local({
        dt <- data.table::data.table(
          Original_value = as.vector(pixel_values),
          #need to go through a dataframe because I don't know how to modify the
          # ID directly in rle object
          Cleanned_value = unlist(pixel_values_cleanned),
          ID = rep(
            seq_along(rle_series_cleanned$lengths),
            rle_series_cleanned$lengths
          ),
          Nb_days = rep(
            rle_series_cleanned$lengths,
            rle_series_cleanned$lengths
          )
        )
        dt
      })

      # Rename event that should gathered because they are separated by d days
      # (or less) bellow threshold
      zero_series <- unique(pixel$ID[
        (is.na(pixel$Cleanned_value) |
          pixel$Cleanned_value == 0) &
          pixel$Nb_days <= d
      ]) #identify series of
      # 0 or NA shorter than d (or equal to d) This imply that a sequence of NA
      # shorter than d will not be considered as a limit of event
      zero_series <- setdiff(zero_series, c(1, max(pixel$ID))) # withdraw the
      # beginning and end of the timeserie as it cannot be surrounded by
      # '1' values
      to_modify <- sort(c(zero_series, zero_series + 1)) #list of futur EE
      # index to rename zero series and one series under one same name when
      # necessary, using the first one serie's name (only when 0 serie length
      # <=d)
      modifier <- zero_series - 1 #sequence of one before the short sequence of
      # 0, <-> the sequence of one that will be merge with the next sequence of
      # 1
      # associate the ID value to replace (firt line) with the value of
      # replacement (second line)
      replacement_map <- stats::setNames(rep(modifier, each = 2), to_modify)
      # Here above, if you have something like 111110111110111110 you will get
      # aaabb as ID instead of aaaaa to adress this issue, I've added the lines
      # bellow :
      pixel$ID <- ifelse(
        pixel$ID %in% names(replacement_map),
        as.numeric(replacement_map[as.character(pixel$ID)]),
        pixel$ID
      )

      pixel$Nb_days <- rep(rle(pixel$ID)$lengths, rle(pixel$ID)$lengths)

      pixel$ID <- paste0(c, "_", pixel$ID)

      pixel$date <- sort(all_dates)

      pixel$id_part <- data.table::tstrsplit(pixel$ID, "_", keep = 2)[[1]]

      # Get the first and last date of each extreme event
      first_date <- stats::ave(pixel$date, pixel$id_part, FUN = function(x) {
        min(x, na.rm = TRUE)
      })
      last_date <- stats::ave(pixel$date, pixel$id_part, FUN = function(x) {
        max(x, na.rm = TRUE)
      })

      # Create new extreme event id has pixel_first_day_last_day
      pixel$ID <- paste0(pixel$id_part, "_", first_date, "_", last_date)
      pixel$id <- NULL # Supprimer la colonne temporaire "id"

      return(pixel)
    }
    #pixel_time_series <- as.data.table(pixel_time_series)
    gc()

    Event_corrected <- list()
    #We want to avoid processing a pixel that is always NA :
    indices_all_na <- which(rowSums(!is.na(pixel_time_series)) == 0)
    indices_to_do <- which(rowSums(!is.na(pixel_time_series)) != 0)

    Event_corrected <- vector("list", nrow(pixel_time_series))
    Event_corrected[indices_to_do] <- lapply(indices_to_do, function(c) {
      correct_lowcomplexity_n_d(pixel_time_series, c, n, d, c)
    })
    n_columns <- dim(pixel_time_series)[2]
    if (length(indices_all_na) > 0) {
      Event_corrected[indices_all_na] <- lapply(indices_all_na, function(c) {
        data.table::data.table(
          Original_value = rep(NA, n_columns),
          Cleanned_value = rep(NA, n_columns),
          ID = rep(paste0(c, "_"), n_columns),
          Nb_days = rep(0, n_columns)
        )
      })
    }
  }

  if (
    !is.null(n) &
      !is.null(d) &
      !is.null(nbis) &
      is.null(w) &
      is.null(p)
  ) {
    #here add code for medium complexity
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
      correct_raster(stacked_rasters[[t]], t, Event_corrected)
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


#function used in the big function :

# Apply corrections
correct_raster <- function(raster, t, Event_corrected) {
  # get values of each layers
  corrected_values <- terra::values(raster)
  Matching_cleanning_values <- sapply(Event_corrected, function(event) {
    event$Cleanned_value[t]
  })
  # Identify where corrections are needed /!\
  correction_indices <- which(
    !is.na(corrected_values) &
      !is.na(Matching_cleanning_values) &
      corrected_values != Matching_cleanning_values
  )
  # Apply corrections
  corrected_values[correction_indices] <-
    Matching_cleanning_values[correction_indices]
  # Modify raster layer values
  terra::values(raster) <- corrected_values
  return(raster)
}
