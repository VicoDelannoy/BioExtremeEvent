#' Binarise the spatraster
#'
#' @description
#'  The function binarise a Spatraster by comparing the observed value with the 
#' baseline. If "direction" argument is "above", the values above the baseline 
#' will be converted to 1, and the value bellow baseline will be converted to 0.
#' The process works in reverse when "direction" = below.
#'
#' @param YourSpatraster The Spatraster containing the values of the studied 
#' parameter, with one layer per timestep.
#' @param baseline A Spatraster built using the BEE.calc.baseline function that
#'  contains the baseline value for each pixel on a given day of the year (with
#'   366 layers).
#' @param direction The accepted values for this argument are "above" or
#'  "below". This tells the function whether it values above baseline or below
#'   baseline are to be considered as extreme events.
#'
#' @return The function returns a list containing a Spatraster for each day of
#'a typical year (366 Spatraster, with one layer per year provided in 
#'YourSpatraster). These rasters have the same properties as "YourSpatraster" 
#'but the values have been replace with 1 if the pixel value was more extreme 
#'than the baseline and with 0 if it was not.
#'
#' @details For computational purposes in the next functions of the pipeline, 
#' the order of the layers differs from that of "YourSpatraster". 
#' BEE.calc.binarize is not design to work on 4D data (3D in space + time). 
#' To do so, please refer to BEE.calc.4thD.
#' #---------------------------------------------------------------------------
#' 
# baseline <- baseline_qt90 ; direction <- "above"; YourSpatraster <- ds

BEE.calc.binarize <- function(YourSpatraster, baseline, direction) {
  if (class(baseline)[1] == "SpatRaster") {
    # Retrieve the extent of YourSpatraster
    new_extent <- terra::ext(YourSpatraster)
    # Change the extent to match that of YourSpatraster
    baseline <- terra::resample(baseline, YourSpatraster, method = "near")
    # Crop to obtain the same area as YourSpatraster
    #baseline <- crop(baseline, new_extent, snap = "in") 
    # Retrieve the dates
    time_list_exp <- terra::time(YourSpatraster, format = "days") 
    # Create dates in the format %m.%d for a leap year
    time_list_base <- generate_month_day(2024) 
    names(YourSpatraster) <- time_list_exp # Rename the layers with the dates
    names(baseline) <- time_list_base # Rename the layers with the dates
    # Check if there is a fourth hidden dimension. For instance there could be
    #several depths or several altitudes which would result in several layer 
    #having the same dates and thus the same length.
    if (any(table(terra::time(YourSpatraster)) != 1)) {
      #at least one date has several layers
      rep <- sum(table(terra::time(YourSpatraster))) / 
        length(table(terra::time(YourSpatraster)))
      if (rep == round(rep)) {
        # all dates are repeted the same number of time
        rast_list <- list()
        n <- terra::nlyr(YourSpatraster)
        for (l in seq(1, rep, 1)) {
          # create a list with all days from the same depths
          layers <- seq(l, n, rep)
          rast_list[[l]] <- terra::rast(YourSpatraster[[layers]])
        }
        delta <- list()
        for (d in seq(1, length(rast_list))) {
          delta[[d]] <- binarize_spat(rast_list[[d]], baseline, direction)
        }
        # Binarize the pixels to 1 & 0 (function from the packacge itself, 
        # see bellow) /!\ this part doesn't check if the value is above for 
        # five consecutive days or any duration
        delta <- binarize_spat(YourSpatraster, baseline, direction) 
        # Create a list of SpatRaster by date in the format %d.%m
        delta_list <- RastToList(delta) 
        return(delta_list) # Output the list
      }
    }
    # Binarize the pixels to 1 & 0 (function from the packacge itself, 
    # see bellow) /!\ this part doesn't check if the value is above for five 
    # consecutive days or any duration
    delta <- binarize_spat(YourSpatraster, baseline, direction) 
    # Create a list of SpatRaster by date in the format %d.%m
    delta_list <- RastToList(delta) 
    return(delta_list) # Output the list
  }
  
  if (is.numeric(baseline)) {
    # Retrieve the dates
    time_list_exp <- terra::time(YourSpatraster, format = "days") 
    # Rename the layers with the dates
    names(YourSpatraster) <- time_list_exp 
    if (direction == "above") {
      # Binarize the cells based on the condition
      delta <- terra::ifel(YourSpatraster < baseline, 0, 1) 
    }
    if (direction == "bellow") {
      # Binarize the cells based on the condition
      delta <- terra::ifel(YourSpatraster < baseline, 1, 0) 
    }
    # Create a list of SpatRaster by date in the format %d.%m
    delta_list <- RastToList(delta) 
    return(delta_list) # Output the list
  }
}


#' Generate dates that take in account lead year
#'
#' @noRd

generate_month_day <- function(year) {
  days_per_month <- c(31, 28 + (year %% 4 == 0 &
                      (year %% 100 != 0 |
                      year %% 400 == 0)),
                      31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  months <- sprintf("%02d", rep(1:12, days_per_month))
  days <- sprintf("%02d", unlist(sapply(days_per_month, seq_len)))
  return(paste(months, days, sep = "."))
}


#' Compare pixel value and baseline value
#'
#' @noRd

binarize_spat <- function(obs, base, direction = "above") {
  #adjust dimension of 'base' so it matches those of 'obs'
  dates <- as.Date(terra::time(obs))  # Dates
  years <- as.integer(format(dates, "%Y")) # Years
  is_leap_year <- function(years) {
    return((years %% 4 == 0 & years %% 100 != 0) | (years %% 400 == 0))
  }
  leap <- is_leap_year(unique(years))
  base_list <- lapply(seq_along(leap), function(y) {
    if (!leap[y]) {
      # Pour une année non bissextile (exclure le jour 60)
      return(c(base[[1:59]], base[[61:366]]))
    } else {
      # Pour une année bissextile (inclus toutes les couches)
      return(base)
    }
  })
  base_extended <- terra::rast(base_list)
  # In case data dont start a first of January and:or don't end the 31th of
  # december, I need to adapt base_extended
  first_year <- min(years)
  last_year <- max(years)
  first_date <- min(dates[format(dates, "%Y") == first_year])
  last_date <- max(dates[format(dates, "%Y") == last_year])
  doy_first_year <- as.integer(format(first_date, "%j"))
  nb_last_day <- which(dates == last_date)
  last_layer <- doy_first_year + nb_last_day - 1
  base_extended <- base_extended[[doy_first_year:last_layer]]
  terra::time(base_extended) <- dates
  #calculate anomalia
  delta <- obs - base_extended
  
  if (direction == "above") {
    delta <- terra::ifel(delta <= 0, 0, 1)
  }
  if (direction == "bellow") {
    delta <- terra::ifel(delta <= 0, 1, 0)
  }
  return(delta)
}


#' Convert from raster to list
#'
#' @noRd

RastToList <- function(x) {
  md <- generate_month_day(2024)
  indices_strata <- list()
  noms_strata <- names(x)
  for (i in unique(md)) {
    indices <- which(grepl(paste0(i, "$"), noms_strata))
    if (length(indices) == 0) {
      next
    }
    indices_strata[[i]] <- indices
  }
  return(lapply(
    indices_strata,
    FUN = function(y) {
      indices <- x[[y]]
    }
  ))
}
