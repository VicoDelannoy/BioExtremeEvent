#' Computes metrics relative to the distance to the border of the extreme event
#' (in surface, 2D)
#'
#'@description
#'This function calculate the mediane distance, mean distance and standard
#'deviation of distance to escape form an extreme event 'as the crow files'
#'(but avoiding NA pixels, which are suppose to represent a place that cannot be
#'crossed) through time for a given GPS position or for all pixels.
#'
#'@param true_event_output is Corrected_rasters, a binarized Spatraster
#' obtained using function
#' BEE.calc.true_event()[["stacked_rasters_corrected"]]).
#' (0 outside of event, 1 extreme event, NA land or missing data) OR
#' Events_corrected, a list of df (obtained using
#' BEE.calc.true_event()[["Event_corrected"]] ).
#'@param start_date and @param end_date defines the timeframe in during which
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

#-------------------------------------------------------------------------------

# BEE.calc.escape is not designed to work on 4D data (time + spatial 3D).
# true_event_output <- result ; start_date = "2022-08-01";
# end_date = "2022-08-12" ; pixel = GPS ; only_days_EE = FALSE ;
# group_by_event = FALSE

BEE.calc.escape <- function(
  true_event_output,
  start_date = NULL,
  end_date = NULL,
  pixel,
  only_days_EE = TRUE,
  group_by_event = TRUE
) {
  ### Recreate a Spatraster using Events_corrected. In this list, there are one
  # df per pixel and one raw per dates
  ## Get data and shell
  data <- true_event_output[[2]]
  shell <- true_event_output[[1]]
  ## Select timeframe of interest to save computation time later
  # Check that dates from the two products of result match
  if (length(data[[1]]$Original_value) != terra::nlyr(shell)) {
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
  print("preparation data ok")

  # Compute by data/layer
  dist_dir <- lapply(rasters, function(x) {
    # x <- rasters[[70]] # 1 : no MHX , 700 : 1 MHW # <-> iterate through dates
    print(x)
    values_x <- terra::values(x)
    if (only_days_EE == TRUE) {
      pixels_from <- which(values_x == 1)
    }
    if (only_days_EE == FALSE) {
      pixels_from <- which(!is.na(values_x))
    }
    pixels_to <- which(values_x == 0)

    if (is.vector(pixel) | class(pixel) == "data.frame") {
      if (class(pixel) == "data.frame") {
        pixels_to_do <- terra::cellFromXY(
          true_event_output$stacked_rasters_corrected,
          pixel
        )
      } else {
        pixels_to_do <- terra::cellFromXY(
          true_event_output$stacked_rasters_corrected,
          t(matrix(pixel))
        )
      }
      pixels_from <- pixels_from[which(pixels_from %in% pixels_to_do)]
    }

    # If there are some pixels to flee but no pixel where to escape :
    if (length(pixels_from) != 0 &
        length(pixels_to) == 0
    ) {
      # the whole area is an EE
      coords_from <- coords_all[pixels_from, ] # Coord to flee
      if (as.character(pixel)[1] == "all") {
        points <- data.table::data.table(
          date = terra::time(x),
          from_x = coords_from[, 1],
          from_y = coords_from[, 2],
          to_x = "unkown",
          to_y = "unkown",
          pixel_from_id = as.integer(pixels_from),
          pixel_to_id = "no escape",
          distance = "outside of the raster",
          azimut = "not possible to compute"
        )
      }
      if (class(pixel) == "data.frame") {
        points <- data.table::data.table(
          date = rep(terra::time(x), nrow(pixel)),
          from_x = pixel[, 1],
          from_y = pixel[, 2],
          to_x = "unkown",
          to_y = "unkown",
          pixel_from_id = as.integer(pixels_from),
          pixel_to_id = "no escape",
          distance = "outside of the raster",
          azimut = "not possible to compute"
        )
      }
      if (is.vector(pixel)) {
        points <- data.table::data.table(
          date = rep(
            terra::time(x),
            nrow(t(
              matrix(pixel)
            ))
          ),
          from_x = t(matrix(pixel))[, 1],
          from_y = t(matrix(pixel))[, 2],
          to_x = "unkown",
          to_y = "unkown",
          pixel_from_id = as.integer(pixels_from),
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
          date = terra::time(x),
          from_x = NA,
          from_y = NA,
          to_x = NA,
          to_y = NA,
          pixel_from_id = NA,
          pixel_to_id = NA,
          distance = 0,
          azimut = NA
        )
      }
      if (class(pixel) == "data.frame") {
        # the pixel of intrest is not in an EE
        points <- data.table::data.table(
          date = rep(terra::time(x), nrow(pixel)),
          from_x = pixel[, 1],
          from_y = pixel[, 2],
          to_x = NA,
          to_y = NA,
          pixel_from_id = pixels_to_do,
          pixel_to_id = NA,
          distance = rep(0, length(pixels_to_do)),
          azimut = NA
        )
      }
        if (is.vector(pixel)) {
          points <- data.table::data.table(
            date = rep(
              terra::time(x),
              nrow(t(
                matrix(pixel)
              ))
            ),
            from_x = t(matrix(pixel))[, 1],
            from_y = t(matrix(pixel))[, 2],
            to_x = NA,
            to_y = NA,
            pixel_from_id = pixels_to_do,
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
        date = rep(terra::time(x), nrow(coords_from) * nrow(coords_to)),
        from_x = rep(coords_from[[1]], each = nrow(coords_to)),
        from_y = rep(coords_from[[2]], each = nrow(coords_to)),
        to_x = rep(coords_to[[1]], times = nrow(coords_from)),
        to_y = rep(coords_to[[2]], times = nrow(coords_from)),
        pixel_from_id = rep(pixels_from, each = nrow(coords_to)),
        pixel_to_id = rep(pixels_to, times = nrow(coords_from))
      )
      # Compute distances btw each points
      points$distance <- geosphere::distHaversine(
        cbind(points$from_x, points$from_y),
        cbind(points$to_x, points$to_y)
      )
      data.table::setorder(points, pixel_from_id, distance) # sort point by id (from_id)
      # and then by shortest distance, when ex aequo, it keeps the same order
      # as in 'points'
      print("avant modif")
      points <- points[!base::duplicated(points$pixel_from_id), ] # keep only the first
      # occurrence of each 'from_id' <-> the occurence with the shortest
      print("modif bien runnee")
      # distance
      points$azimut <- (geosphere::bearing(
        cbind(points$from_x, points$from_y),
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
  print("boucle points_data terminee")
  warnings(
    "When several pixels are the 'closest pixel', the one with the smallest 
    number/id is kept to compute shortest distance and azimut. Moreover, this 
    function uses 'geosphere::distHaversine' to compute distances, it accounts 
    for earth rotondity, which is better than euclidian distance, but it 
    consider the earth as a sphere, thus near the poles, this methods is less 
    percise than geosphere::distVincentyEllipsoid. This is a compromise between
    computation speed and precision. The estimate error using distHaversine 
    compare to distVincentyEllipsoid is between 1 km and 5 km for a distance of 
    300 km."
  )
  print("warning passe")
  dist_dir <- data.table::rbindlist(dist_dir)
  print("dist_dir cree")
  no_event <- dist_dir[which(dist_dir$pixel_to_id == "no escape"), ] # saving the lines
  # where there are no distances to compute for the case only_days_EE == FALSE
  print("no event identifie")
  if (only_days_EE == TRUE) {
    print("only_days_EE = T ")
    dist_dir <- dist_dir[dist_dir$distance != 0, ]
    if (nrow(dist_dir) == 0) {
      message("There were no extreme event for the given pixels and timeframe.")
    }
  }
  print("si only_days_EE=TRUE days 0 supprimes")
  if (any(dist_dir$pixel_to_id == "no escape", na.rm = TRUE)) {
    dangerous_date <- dist_dir$date[which(dist_dir$pixel_to_id == "no escape")]
    message(
      "During the following dates, the raster was fully covered by an EE and it
      was no possible to compute a distance to escape or an azimut.",
      paste(dangerous_date, collapse = ", ")
    )
  }

  if (group_by_event == FALSE) {
    print("return en cours")
    # dist_dir <- points
    return(dist_dir)
  }
  if (group_by_event == TRUE) {
    print("groupement par event commence")
    # Here we want to give 'summary' type value compute across all day of a
    # same event for each pixel (and event)
    # First, we need to re-identify all the days that belong to a same event :
    ## Create a data.table from  true_event_output[[2]] (data) with a column to
    # easly identify the pixel represented in each row :
    data <- data[unique(dist_dir$pixel_from_id)]
    names(data) <- unique(dist_dir$pixel_from_id)
    data <- data.table::rbindlist(data, idcol = "pixel_from_id")
    ## Recreate a column for the time interval of the event using the
    # informaiton in the ID, this will be necessay to use foverlaps, the fastest
    # function to identify if a date belongs to a specific timeframe or not :) :
    data$start_date <- as.Date(sapply(strsplit(data$ID, "_"), `[`, 2))
    data$end_date <- as.Date(sapply(strsplit(data$ID, "_"), `[`, 3))
    data$start_date <- as.Date(data$start_date)
    data$end_date <- as.Date(data$end_date)
    dist_dir$date_start <- dist_dir$date
    dist_dir$date_end <- dist_dir$date
    print("date correctement gerees")
    ## Convert the pixel_id to the same format in both data.table :
    dist_dir$pixel_from_id <- as.integer(dist_dir$pixel_from_id)
    data$pixel_from_id <- as.integer(data$pixel_from_id)
    ## We need to delete 'event' that are bellow threshold (<-> distance == 0)
    print("pixel_id ok")
    ## Define the 'key' <-> the combination of colum used to create sub-group
    # in a data.table format :
    data.table::setkey(dist_dir, pixel_from_id, date_start, date_end) 
    data.table::setkey(data, pixel_from_id, start_date, end_date)
print("sous groupes ok")
    dist_dir <- data.table::foverlaps(
      # magic function that found if the line from dt x is inside the time frame
      # of a line from dt y and identify which one AND joint the 2 dt respecting
      # correspondanies btw several column.
      dist_dir,
      data,
      by.x = c("pixel_from_id", "date_start", "date_end"),
      # dt for which we need to know if it is inside the time frames of the
      # other dt
      by.y = c("pixel_from_id", "start_date", "end_date"),
      # the other dt, that contains all the reference timeframe
      nomatch = NA
    )
    print("boucle foverlaps ok")
    ## Clean the columns date_start and date_end that are redondant with date
    dist_dir$date_start <- NULL
    dist_dir$date_end <- NULL

    # Now we can compute some metrics of the metrics distance and azimut !!
    old_warn <- options("warn") # Save curent warning settings
    options(warn = -1) # unactivate warnings
    ## Distance
    tmp_mean <- aggregate(
      as.numeric(distance) ~ ID,
      data = dist_dir,
      FUN = mean,
      na.rm = TRUE
    )
    dist_dir$distance_mean <- tmp_mean[match(dist_dir$ID, tmp_mean$ID), 2]
    print("distance ok")
    distance_sd <- tapply(
      as.numeric(dist_dir$distance),
      dist_dir$ID,
      sd,
      na.rm = TRUE
    )
    dist_dir$distance_sd <- distance_sd[dist_dir$ID]
    distance_median <- tapply(
      as.numeric(dist_dir$distance),
      dist_dir$ID,
      median,
      na.rm = TRUE
    )
    print("distance sd ok")
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
print("distance median min et max ok")
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
print("azimut circ ok")
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
      median,
      na.rm = TRUE
    )
    dist_dir$azimut_med <- azimut_med[dist_dir$ID]
print("azimut mean et median ok")
    tmp_sd <- aggregate(
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
    options(warn = old_warn$warn) # reactivate warnings
    print("azimut sd min et max ok")
    # Delete column that refer to daily value and not to value compute on all
    # the event :
    dist_dir[,
      c('date', 'distance', 'azimut', 'azimut_num', 'azimut_circ')
    ] <- NULL
    print("column with daily value delated")
    # Keep only one row per EE (<-> per value of ID column) :
    dist_dir <- dist_dir[!duplicated(dist_dir$ID), ] # keep the first line of every
    # group of same value of ID
    print("duplicated lines deleted")
    # delete the column x_to and y_to as the represent the pixel where to
    # escape on first day and not a characteristic of the all event
    dist_dir <- dist_dir |> dplyr::select(-to_x, -to_y, -pixel_to_id)
    dist_dir <- dist_dir[!is.na(dist_dir$pixel_from_id)] # withrdraw the line of NA
    print("fin code")
    return(dist_dir)
  }
  message("this combination of argument is not endle by the function")
}



#' To compute azimut sd taking in account NA
#' @noRd

azimut_sd_fun <- function(x) {
  az_group <- x[!is.na(x)]
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

  val <- suppressWarnings(as.numeric(stats::sd(az_group)))

  if (is.nan(val)) 0 else val
}
