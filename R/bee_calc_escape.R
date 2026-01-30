#' Computes metrics relative to the distance to the border of the extreme event
#' (in surface, 2D)
#'
#'@description
#'This function calculate the median distance, mean distance and standard
#'deviation of distance to escape form an extreme event 'as the crow files'
#'(but avoiding NA pixels, which are suppose to represent a place that cannot be
#'crossed) through time for a given GPS position or for all pixels.
#'
#'@param true_event_output is the ouput of function BEE.calc.true_event().
#'@param start_date defines the beginning of the timeframe during which
#' you want to analyse distance to escape. If no dates are provided, computation
#' will be done on all days_provided.
#' @param end_date defines the endding of the timeframe during which
#' you want to analyse distance to escape. If no dates are provided, computation
#' will be done on all days_provided
#'@param only_days_EE = 'TRUE' if you want to perform calculations
#' using all the days. Use account = 'FALSE' when you want to calculate only
#' for the days when your reference point is experiencing an extreme event.
#' By default, account_all_days is set on FALSE
#'@param group_by_event allows to specify if you want your result by extreme
#'  event or for the whole time periode specified. "group_by_event" = TRUE : compute
#'  one mean distance per event. "group_by_event" = FALSE : mean distance will be
#'  based on distances from all days in the time frame (belonging to an extreme
#'  event or not),  including all days belonging to this event event if they
#'  are a '0' day, TRUE : compute one mean distance accounting only for days
#'  that belong to an EE.
#'@param pixel is a vector with x,y coordinates or a df of x,y coordinates or
#' "all" if you want to compute distance metrics for all pixels in the raster
#' provided through true_event_output
#'
#'@note all distance metrics are limited by the size of the Spatraster your
#' are providing.
#'
#' @export
#'
#-------------------------------------------------------------------------------

# BEE.calc.escape is not designed to work on 4D data (time + spatial 3D).
# true_event_output <- result ; start_date = "2024-08-01";
# end_date = "2024-08-20" ; pixel = GPS ; only_days_EE = FALSE ;
# group_by_event = FALSE

BEE.calc.escape <- function(
  true_event_output,
  start_date = NULL,
  end_date = NULL,
  pixel,
  only_days_EE = FALSE,
  group_by_event = FALSE
) {
  ### Recreate a Spatraster using Events_corrected. In this list, there are one
  # df per pixel and one raw per dates
  ## Get data and shell
  data <- true_event_output[[2]]
  data <- Filter(Negate(is.null), data)
  shell <- true_event_output[[1]]
  ## Select timeframe of interest to save computation time later
  # Check that dates from the two products of result match
  if (length(data[[1]]$original_value) != terra::nlyr(shell)) {
    warning(
      "The two elements you provided as 'true_event_output don't cover the same
      number of days, please use the output of BEE.calc.true_event."
    )
  }
  date_indices <- match(c(start_date, end_date), names(shell))
  dates <- names(shell)[date_indices[1]:date_indices[2]]
  data <- lapply(data, \(df) df[date_indices[1]:date_indices[2], ])

  # Subset the layers in the timeframe of interest
  rasters <- true_event_output[[1]][[
    names(true_event_output[[1]]) >= start_date &
      names(true_event_output[[1]]) <= end_date
  ]]

  # coordinates of every pixel in the dataset :
  coords_all <- terra::crds(
    rasters[[1]],
    df = TRUE,
    na.rm = FALSE,
    na.all = TRUE
  )

  # Compute by data/layer
  dist_dir <- lapply(rasters, function(x_r) {
    # x_r <- rasters[[4]] # 1 : no MHW , 700 : 1 MHW # <-> iterate through dates
    values_x <- terra::values(x_r)
    if (only_days_EE == TRUE) {
      pixels_from <- which(values_x == 1)
    }
    if (only_days_EE == FALSE) {
      pixels_from <- which(!is.na(values_x))
    }
    pixels_to <- which(values_x == 0)

    if (is.vector(pixel) | is.data.frame(pixel)) {
      if (is.data.frame(pixel)) {
        pixels_to_do <- terra::cellFromXY(
          true_event_output$stacked_rasters_corrected,
          pixel
        )
      } else {
        if (!is.character(pixel)) {
          pixels_to_do <- terra::cellFromXY(
            true_event_output$stacked_rasters_corrected,
            t(matrix(pixel))
          )
        }
        if (pixel == "all") {
          pixels_to_do <- which(!is.na(terra::values(x_r)))
        } else {
          warnings(
            "The current 'pixel' argument format is not accepted, 
          please try with a dataframe with only an x and y colum or using 'pixel = all'."
          )
        }
      }
      pixels_from <- pixels_from[which(pixels_from %in% pixels_to_do)]
    }

    # If there are some pixels to flee but no pixel where to escape :
    if (
      length(pixels_from) != 0 &
        length(pixels_to) == 0
    ) {
      # the whole area is an EE
      coords_from <- coords_all[pixels_from, ] # Coord to flee
      if (as.character(pixel)[1] == "all") {
        points <- data.table::data.table(
          date = terra::time(x_r),
          x = coords_from[, 1],
          y = coords_from[, 2],
          to_x = "unkown",
          to_y = "unkown",
          pixel_id = as.integer(pixels_from),
          pixel_to_id = "no escape",
          distance = "outside of the raster",
          azimut = "not possible to compute"
        )
      }
      if (is.data.frame(pixel)) {
        points <- data.table::data.table(
          date = rep(terra::time(x_r), nrow(pixel)),
          x = pixel[, 1],
          y = pixel[, 2],
          to_x = "unkown",
          to_y = "unkown",
          pixel_id = as.integer(pixels_from),
          pixel_to_id = "no escape",
          distance = "outside of the raster",
          azimut = "not possible to compute"
        )
      }
      if (is.vector(pixel) & as.character(pixel)[1] != "all") {
        points <- data.table::data.table(
          date = rep(
            terra::time(x_r),
            nrow(t(
              matrix(pixel)
            ))
          ),
          x = t(matrix(pixel))[, 1],
          y = t(matrix(pixel))[, 2],
          to_x = "unkown",
          to_y = "unkown",
          pixel_id = as.integer(pixels_from),
          pixel_to_id = "no escape",
          distance = "outside of the raster",
          azimut = "not possible to compute"
        )
      }
    }
    if (length(pixels_from) == 0) {
      # nothing to escape from AND only_days_EE==TRUE (because when
      # only_days_EE==FALSE, all pixel are considered as pixel to escape
      # because we want to compute the distance to escape that are = 0.)
      if (as.character(pixel)[1] == "all") {
        points <- data.table::data.table(
          date = terra::time(x_r),
          x = NA,
          y = NA,
          to_x = NA,
          to_y = NA,
          pixel_id = NA,
          pixel_to_id = NA,
          distance = 0,
          azimut = NA
        )
      }
      if (is.data.frame(pixel)) {
        # the pixel of intrest is not in an EE
        points <- data.table::data.table(
          date = rep(terra::time(x_r), nrow(pixel)),
          x = pixel[, 1],
          y = pixel[, 2],
          to_x = NA,
          to_y = NA,
          pixel_id = pixels_to_do,
          pixel_to_id = NA,
          distance = rep(0, length(pixels_to_do)),
          azimut = NA
        )
      }
      if (is.vector(pixel) & as.character(pixel)[1] != "all") {
        points <- data.table::data.table(
          date = rep(
            terra::time(x_r),
            nrow(t(
              matrix(pixel)
            ))
          ),
          x = t(matrix(pixel))[, 1],
          y = t(matrix(pixel))[, 2],
          to_x = NA,
          to_y = NA,
          pixel_id = pixels_to_do,
          pixel_to_id = NA,
          distance = rep(0, length(pixels_to_do)),
          azimut = NA
        )
      }
    }

    if (
      length(pixels_from) != 0 &
        length(pixels_to) != 0
    ) {
      #Basic situation
      # Get coordinates of pixel to flee and pixel where to escape
      coords_from <- coords_all[pixels_from, ] # Coord to flee
      coords_to <- coords_all[pixels_to, ] # Coord refuge
      # Create a data.table
      points <- data.table::data.table(
        date = rep(terra::time(x_r), nrow(coords_from) * nrow(coords_to)),
        x = rep(coords_from[[1]], each = nrow(coords_to)),
        y = rep(coords_from[[2]], each = nrow(coords_to)),
        to_x = rep(coords_to[[1]], times = nrow(coords_from)),
        to_y = rep(coords_to[[2]], times = nrow(coords_from)),
        pixel_id = rep(pixels_from, each = nrow(coords_to)),
        pixel_to_id = rep(pixels_to, times = nrow(coords_from))
      )
      # Compute distances btw each points
      points$distance <- geosphere::distHaversine(
        cbind(points$x, points$y),
        cbind(points$to_x, points$to_y)
      )
      data.table::setorder(points, pixel_id, distance) # sort point by id (from_id)
      # and then by shortest distance, when ex aequo, it keeps the same order
      # as in 'points'
      points <- points[!duplicated(points$pixel_id), ] # keep only the first
      # occurrence of each 'from_id' <-> the occurence with the shortest

      # distance
      points$azimut <- (geosphere::bearing(
        cbind(points$x, points$y),
        cbind(points$to_x, points$to_y)
      ) +
        360) %%
        360
    }
    if (length(pixels_from) == 0 & length(pixels_to) == 0) {
      warning(
        "There are no pixel of value 1 AND there are no pixel of value 0, 
        the raster are probably not binarized (please see BEE.calc.binarized_EE
        and BEE.calc.true_event) or the raster is fully covered by NA."
      )
    }
    # Make sure that when they are no reason to leave a pixel (distance = 0),
    # the azimut is NA
    points$azimut <- ifelse(points$distance == 0, NA, points$azimut)
    return(points)
  })

  warnings(
    "When several pixels are the 'closest pixel', the one with the smallest 
    number/id is kept to compute shortest distance and azimut. Moreover, this 
    function uses 'geosphere::distHaversine' to compute distances, it accounts 
    for earth rotondity, which is better than euclidian distance, but it 
    consider the earth as a sphere, thus near the poles, this methods is less 
    percise than geosphere::distVincentyEllipsoid. This is a compromise between
    computation speed and precision. The estimate error using distHaversine 
    compare to distVincentyEllipsoid is between 1 km and 5 km for a distance of 
    300 km in polar area."
  )

  dist_dir <- data.table::rbindlist(dist_dir)

  no_event <- dist_dir[which(dist_dir$pixel_to_id == "no escape"), ] # saving the lines
  # where there are no distances to compute for the case only_days_EE == FALSE

  if (only_days_EE == TRUE) {
    dist_dir <- dist_dir[dist_dir$distance != 0, ]
    if (nrow(dist_dir) == 0) {
      message("There were no extreme event for the given pixels and timeframe.")
    }
  }

  if (any(dist_dir$pixel_to_id == "no escape", na.rm = TRUE)) {
    dangerous_date <- dist_dir$date[which(dist_dir$pixel_to_id == "no escape")]
    message(
      "During the following dates, the raster was fully covered by an EE and it
      was no possible to compute a distance to escape or an azimut.",
      paste(unique(dangerous_date), collapse = ", ")
    )
  }

  if (group_by_event == FALSE) {
    names(data) <- sapply(data, function(df) df$pixel_id[1])
    data <- data[names(data) %in% as.character(unique(dist_dir$pixel_id))]
    data <- data.table::rbindlist(data)
    data$date <- as.Date(data$date)
    dist_dir$date <- as.Date(dist_dir$date)
    dist_dir <- merge(
      dist_dir,
      data[, c("pixel_id", "date", "ID")],
      by = c("pixel_id", "date"),
      all.x = TRUE
    )
    return(dist_dir)
  }
  if (group_by_event == TRUE) {
    # Here we want to give 'summary' type value compute across all day of a
    # same event for each pixel (and event)
    # First, we need to re-identify all the days that belong to a same event :
    ## Create a data.table from  true_event_output[[2]] (data) with a column to
    # easly identify the pixel represented in each row :
    names(data) <- sapply(data, function(df) df$pixel_id[1])
    data <- data[names(data) %in% as.character(unique(dist_dir$pixel_id))]
    data <- data.table::rbindlist(data)
    data$date <- as.Date(data$date)
    dist_dir$date <- as.Date(dist_dir$date)
    dist_dir <- merge(
      dist_dir,
      data[, c("pixel_id", "date", "ID")],
      by = c("pixel_id", "date"),
      all.x = TRUE
    )

    ## Distance
    dist_dir$distance <- as.numeric(dist_dir$distance)
    tmp_mean <- aggregate(
      distance ~ ID,
      data = dist_dir,
      FUN = function(df) mean(df, na.rm = TRUE)
    )
    dist_dir$distance_mean <- tmp_mean[
      match(dist_dir$ID, tmp_mean$ID),
      2
    ]

    distance_sd <- tapply(
      as.numeric(dist_dir$distance),
      dist_dir$ID,
      stats::sd,
      na.rm = TRUE
    )
    dist_dir$distance_sd <- distance_sd[dist_dir$ID]

    distance_median <- tapply(
      as.numeric(dist_dir$distance),
      dist_dir$ID,
      stats::median,
      na.rm = TRUE
    )
    dist_dir$distance_median <- distance_median[dist_dir$ID]

    distance_min <- tapply(
      as.numeric(dist_dir$distance),
      dist_dir$ID,
      min,
      na.rm = TRUE
    )
    dist_dir$distance_min <- distance_min[dist_dir$ID]

    distance_max <- tapply(
      as.numeric(dist_dir$distance),
      dist_dir$ID,
      max,
      na.rm = TRUE
    )
    dist_dir$distance_max <- distance_max[dist_dir$ID]

    ## Since degree are a circular unit (after 360=0 comes 1,2...) we need a
    # special way to compute it :
    # Conversion to circular angle :
    dist_dir$azimut_num <- as.numeric(dist_dir$azimut)
    dist_dir$azimut_circ <- ifelse(
      !is.na(dist_dir$azimut_num),
      circular::circular(
        dist_dir$azimut_num,
        units = "degrees",
        template = "geographics"
      ),
      NA
    )
    azimut_mean <- tapply(
      as.numeric(dist_dir$azimut_circ),
      dist_dir$ID,
      mean,
      na.rm = TRUE
    )
    dist_dir$azimut_mean <- azimut_mean[dist_dir$ID]
    azimut_med <- tapply(
      as.numeric(dist_dir$azimut_circ),
      dist_dir$ID,
      stats::median,
      na.rm = TRUE
    )
    dist_dir$azimut_med <- azimut_med[dist_dir$ID]

    tmp_sd <- stats::aggregate(
      azimut_circ ~ ID,
      data = dist_dir,
      FUN = azimut_sd_fun
    )
    dist_dir$azimut_sd <- tmp_sd[match(dist_dir$ID, tmp_sd$ID), 2]

    azimut_min <- tapply(
      as.numeric(dist_dir$azimut_circ),
      dist_dir$ID,
      min,
      na.rm = TRUE
    )
    dist_dir$azimut_min <- azimut_min[dist_dir$ID]

    azimut_max <- tapply(
      as.numeric(dist_dir$azimut_circ),
      dist_dir$ID,
      max,
      na.rm = TRUE
    )
    dist_dir$azimut_max <- azimut_max[dist_dir$ID]

    # Delete column that refer to daily value and not to value compute on all
    # the event :
    dist_dir[,
      c(
        'date',
        'distance',
        'azimut',
        'azimut_num',
        'azimut_circ',
        'to_x',
        'to_y',
        'pixel_to_id'
      )
    ] <- NULL

    # Keep only one row per EE (<-> per value of ID column) :
    dist_dir <- dist_dir[!duplicated(dist_dir$ID), ] # keep the first line of every
    # group of same value of ID
    dist_dir <- dist_dir[!is.na(dist_dir$pixel_id), ] # withrdraw the line of NA

    return(dist_dir)
  }
  message("this combination of argument is not endle by the function")
}


#' To compute azimut sd taking in account NA
#' @noRd

azimut_sd_fun <- function(x_a) {
  az_group <- x_a[!is.na(x_a)]
  # when there are to many
  # decimals in azimut_circ, sin and cos used inside circular::sd induce
  # little imprecisions that leads to NA value when azimut_circ is constant,
  # to deal with this, when azimut-circ is constant, a value of 0 is forced
  # into the dt.
  if (
    length(az_group) <= 1 ||
      max(az_group) - min(az_group) < 1e-6
  ) {
    return(0)
  }

  val <- suppressWarnings(as.numeric(circular::sd(az_group)))

  if (is.nan(val)) 0 else val
}
