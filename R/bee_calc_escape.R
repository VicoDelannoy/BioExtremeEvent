#' Computes metrics relative to the distance to the border of the extreme event
#' (in surface, 2D)
#'
#' @description
#'  This function calculate the median distance, mean distance and standard
#'  deviation of distance to escape form an extreme event. Two methods are used:
#'  First is the "strait" distance, measured 'as the crow files'this method
#'  doesn't avoid the NA pixels and consider the earth as a sphere for distance
#'  computation.
#'  The second method, "bio", also account for the sphericity of the earth and
#'  compute distances over trajectory that avoid NA pixel. To do so the
#'  trajectory is built step by step which induce more imprecision than the
#'  "strait" method.
#'  If the difference between the two distances is small than the spatraster
#'  resolution, it is very likely that the difference is due to the lower
#'  precision of the "bio" methods. If the differences is close or superior to
#'  the spatraster resolution, it is very likely that it is because the "bio"
#'  distance take in account the extra distance caused by going round an
#'  obstacle (NA pixel).
#'  (see the package vignette for more detail).
#'  The function also provided the azimut between the pixel to flee and the
#'  closest refuge pixel according each method.
#'
#' @param extreme_events_output :
#'  Is the ouput of function BEE.id.extreme_events().
#' @param start_date :
#'  Defines the beginning of the timeframe during which you want to analyse
#'  distance to escape. If no dates are provided, computation will be done
#'  on all days_provided.
#' @param end_date :
#'  Defines the endding of the timeframe during which you want to analyse
#'  distance to escape. If no dates are provided, computation will be done
#'  on all days_provided
#' @param only_days_EE :
#'  = 'TRUE' if you want to perform calculations using all the days.
#'  Use account = 'FALSE' when you want to calculate only for the days when your
#'  reference point is experiencing an extreme event.
#'  By default, only_days_EE is set on FALSE
#' @param group_by_event :
#'  Allows to specify if you want your result by extreme event or for the whole
#'  time periode specified. "group_by_event" = TRUE : compute one mean distance
#'  per event. "group_by_event" = FALSE : mean distance will be based on
#'  distances from all days in the time frame (belonging to an extreme event or
#'  not), including all days belonging to this event event if they are a '0'
#'  day, TRUE : compute one mean distance accounting only for days that belong
#'  to an EE.
#' @param pixel :
#'  Is a vector with x,y coordinates or a df of x,y coordinates or "all" if you
#'  want to compute distance metrics for all pixels in the raster provided
#'  through extreme_events_output.
#'
#' @note
#'  All distance metrics are limited by the size of the Spatraster your
#'  are providing. See function is not developped to work on 4D data (time + 3D).
#'
#' @return
#'  A list with one dataframe per pixel. Please make sure to call element of the
#' list by their names and not by their positions.
#'
#' @examples
#' # Load data for example: library(BioExtremeEvent)
#'# reconstruct the output of BEE.id.extreme_event() from data saved in local:
#'file_name_1 <- system.file(file.path("extdata",
#'                                      "binarized_corrected_spatraster.tiff"),
#'                                      package = "BioExtremeEvent")
#'extreme_event_spatraster <- terra::rast(file_name_1)
#'file_name_2 <- system.file(file.path("extdata",
#'                                      "binarized_corrected_df.rds"),
#'                                      package = "BioExtremeEvent")
#'extreme_event_list <- readRDS(file_name_2)
#'
#'extreme_event <- list(extreme_event_spatraster,extreme_event_list)
#'
#'# create a dataframe with the positions of the studdied sites:
#'GPS <- data.frame(x = c(3.7659),
#'                  y = c(43.4287))
#'
#' #' # Get distance, azimut and pixel coordinates to escape:
#' library(BioExtremeEvent)
#'
#'escape <- BioExtremeEvent::BEE.calc.escape(
#'  extreme_events_output = extreme_event,
#'  start_date = "2024-05-01",
#'  end_date = "2024-11-30",
#'  pixel= GPS,
#'  only_days_EE = FALSE,
#'  group_by_event = FALSE
#')
#'
#' @export
#'
#-------------------------------------------------------------------------------

# BEE.calc.escape is not designed to work on 4D data (time + spatial 3D).
# extreme_events_output <- result ; start_date = "2024-08-01";
# end_date = "2024-08-20" ; pixel = gps ; only_days_EE = FALSE ;
# group_by_event = FALSE

BEE.calc.escape <- function(
  extreme_events_output,
  start_date = NULL,
  end_date = NULL,
  pixel,
  only_days_EE = FALSE,
  group_by_event = FALSE
) {
  ### Recreate a Spatraster using events_corrected. In this list, there are one
  # df per pixel and one raw per dates
  ## Get data and shell
  data <- extreme_events_output[[2]]
  data <- Filter(Negate(is.null), data)
  shell <- extreme_events_output[[1]]
  ## Select timeframe of interest to save computation time later
  # Check that dates from the two products of result match
  if (length(data[[1]]$original_value) != terra::nlyr(shell)) {
    warning(
      "The two elements you provided as 'extreme_events_output' don't cover the same
      number of days, please use the output of BEE.id.extreme_events."
    )
  }
  date_indices <- terra::match(c(start_date, end_date), terra::names(shell))
  dates <- terra::names(shell)[date_indices[1]:date_indices[2]]
  data <- lapply(data, \(df) df[date_indices[1]:date_indices[2], ])

  # Subset the layers in the timeframe of interest
  rasters <- extreme_events_output[[1]][[
    terra::names(extreme_events_output[[1]]) >= start_date &
      terra::names(extreme_events_output[[1]]) <= end_date
  ]]

  # coordinates of every pixel in the dataset :
  coords_all <- terra::crds(
    rasters[[1]],
    df = TRUE,
    na.rm = FALSE,
    na.all = TRUE
  )

  ###--- Create the cost raster (NA pixels can not be crossed) for distance "bio"
  # (avoinding NA pixels)
  ## step1 : create connexion between pixels neighbours_matrix
  neighbours_matrix <- matrix(
    nrow = 7,
    ncol = 7,
    data = c(
      1,
      1,
      1,
      0,
      1,
      1,
      1,
      1,
      1,
      1,
      0,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      0,
      0,
      1,
      1,
      1,
      0,
      0,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      1,
      0,
      1,
      1,
      1,
      1,
      1,
      1,
      0,
      1,
      1,
      1
    )
  )
  non_NA_pixels <- which(!is.na(terra::values(rasters[[1]])))
  neighbours <- terra::adjacent(
    rasters[[1]],
    cells = non_NA_pixels, #not intrested in the
    # neighbour of NA pixels
    directions = neighbours_matrix,
    pairs = TRUE
  ) #give 2 colonnes instead fo a matrix

  #step2a, withdraw the connexions starting from an NA pixel
  neighbours <- neighbours[which(neighbours[, "from"] %in% non_NA_pixels), ]
  #step2b, withdraw the connexions arriving in an NA pixel
  neighbours <- neighbours[which(neighbours[, "to"] %in% non_NA_pixels), ]
  ## step3 : get the coordinates of all combinaison departure-arrival
  x_coords <- coords_all$x
  y_coords <- coords_all$y
  from_all <- neighbours[, 1]
  to_all <- neighbours[, 2]
  p1 <- cbind(x_coords[from_all], y_coords[from_all]) # all starting pts coords
  p2 <- cbind(x_coords[to_all], y_coords[to_all]) #all associated arrival points

  ## step 3 : identify the paths that fly over NA pixels and discard them
  vals <- terra::values(rasters[[1]])
  n <- 50
  m <- nrow(p1)
  all_cells <- vector("list", m)
  valid <- logical(m)
  for (i in seq_len(m)) {
    # almost 6 minute for 434
    #i = 5896 diagonale from 17 to 113
    pts <- geosphere::gcIntermediate(
      p1[i, ],
      p2[i, ],
      n = n,
      addStartEnd = TRUE
    )
    cells <- terra::cellFromXY(rasters[[1]], pts)
    valid[i] <- !any(is.na(vals[cells]))
  }
  #withdraw path going over NA pixel
  p1 <- p1[valid, ]
  p2 <- p2[valid, ]

  cost_bio <- geosphere::distGeo(p1, p2)
  graph_bio <- cppRouting::makegraph(data.frame(
    from_vertex = neighbours[valid, 1],
    to_vertex = neighbours[valid, 2],
    cost = cost_bio
  ))

  #4 compute shortest path between each pixel
  nodes <- graph_bio$dict$ref
  shortest_paths_bio <- cppRouting::get_distance_matrix(
    graph_bio,
    from = nodes,
    to = nodes
  )

  ###--- Distance between every pixel "as the crows fly"

  p1_strait <- cbind(
    x_coords[sort(rep(non_NA_pixels, length(non_NA_pixels)))],
    y_coords[sort(rep(non_NA_pixels, length(non_NA_pixels)))]
  ) # all starting pts coords
  p <- cbind(x_coords[non_NA_pixels], y_coords[non_NA_pixels])
  p2_strait <- do.call(
    rbind,
    replicate(length(non_NA_pixels), p, simplify = FALSE)
  )
  cost_strait <- geosphere::distGeo(p1_strait, p2_strait)
  graph_strait <- cppRouting::makegraph(data.frame(
    from_vertex = sort(rep(non_NA_pixels, length(non_NA_pixels))),
    to_vertex = rep(non_NA_pixels, length(non_NA_pixels)),
    cost = cost_strait
  ))
  # compute shortest path between each pixel
  nodes_strait <- graph_strait$dict$ref
  shortest_paths_strait <- cppRouting::get_distance_matrix(
    graph_strait,
    from = nodes_strait,
    to = nodes_strait
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
          shell,
          pixel
        )
      } else {
        if (!is.character(pixel)) {
          pixels_to_do <- terra::cellFromXY(
            shell,
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
          distance_bio = "outside of the raster",
          distance_strait = "outside of the raster",
          azimut_bio = "not possible to compute",
          azimut_strait = "not possible to compute"
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
          distance_bio = "outside of the raster",
          distance_strait = "outside of the raster",
          azimut_bio = "not possible to compute",
          azimut_strait = "not possible to compute"
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
          distance_bio = "outside of the raster",
          distance_strait = "outside of the raster",
          azimut_bio = "not possible to compute",
          azimut_strait = "not possible to compute"
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
          distance_bio = 0,
          distance_strait = 0,
          azimut_bio = NA,
          azimut_strait = NA
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
          distance_bio = rep(0, length(pixels_to_do)),
          distance_strait = rep(0, length(pixels_to_do)),
          azimut_bio = NA_real_,
          azimut_strait = NA_real_
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
          distance_bio = rep(0, length(pixels_to_do)),
          distance_strait = rep(0, length(pixels_to_do)),
          azimut_bio = NA_real_,
          azimut_strait = NA_real_
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
        pixel_to_id = rep(pixels_to, times = nrow(coords_from)),
        distance_bio = NA_real_,
        distance_strait = NA_real_,
        azimut_bio = NA_real_,
        azimut_strait = NA_real_
      )
      # Compute distances btw each points
      id_to_pos <- stats::setNames(
        seq_len(nrow(shortest_paths_bio)),
        rownames(shortest_paths_bio)
      )
      rows <- unname(id_to_pos[as.character(points$pixel_to_id)])
      cols <- unname(id_to_pos[as.character(points$pixel_id)])

      points$distance_bio <- shortest_paths_bio[cbind(rows, cols)]
      points$distance_strait <- shortest_paths_strait[cbind(rows, cols)]

      points_bio <- points[
        order(points$pixel_id, points$distance_bio),
      ]
      points_strait <- points[
        order(points$pixel_id, points$distance_bio),
      ]
      res_bio <- points_bio[!duplicated(points_bio$pixel_id), ]
      res_strait <- points_strait[
        !duplicated(points_strait$pixel_id),
      ]
      points <- rbind(res_bio, res_strait)
      points <- unique(points)

      points$azimut_bio <- (geosphere::bearing(
        cbind(points$x, points$y),
        cbind(points$to_x, points$to_y)
      ) +
        360) %%
        360
      points$azimut_strait <- (geosphere::bearing(
        cbind(points$x, points$y),
        cbind(points$to_x, points$to_y)
      ) +
        360) %%
        360
    }
    if (length(pixels_from) == 0 & length(pixels_to) == 0) {
      warning(
        "There are no pixel of value 1 AND there are no pixel of value 0, 
        the raster are probably not binarized (please see BEE.id.extreme_day
        and BEE.id.extreme_events) or the raster is fully covered by NA."
      )
    }
    # Make sure that when they are no reason to leave a pixel (distance = 0),
    # the azimut is NA
    points$azimut_bio <- ifelse(points$distance_bio == 0, NA, points$azimut_bio)
    points$azimut_strait <- ifelse(
      points$distance_strait == 0,
      NA,
      points$azimut_strait
    )
    return(points)
  })

  message(
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
    dist_dir <- dist_dir[
      dist_dir$distance_bio != 0 & dist_dir$distance_strait != 0,
    ]
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
    dist_dir <- data.table::merge(
      dist_dir,
      data[, c("pixel_id", "date", "ID")],
      by = c("pixel_id", "date"),
      all.x = TRUE
    )
    dist_dir <- data.table::split(dist_dir, dist_dir$pixel_id)
    return(dist_dir)
  }
  if (group_by_event == TRUE) {
    # Here we want to give 'summary' type value compute across all day of a
    # same event for each pixel (and event)
    # First, we need to re-identify all the days that belong to a same event :
    ## Create a data.table from  extreme_events_output[[2]] (data) with a column to
    # easly identify the pixel represented in each row :
    names(data) <- sapply(data, function(df) df$pixel_id[1])
    data <- data[names(data) %in% as.character(unique(dist_dir$pixel_id))]
    data <- data.table::rbindlist(data)
    data$date <- as.Date(data$date)
    dist_dir$date <- as.Date(dist_dir$date)
    dist_dir <- data.table::merge(
      dist_dir,
      data[, c("pixel_id", "date", "ID")],
      by = c("pixel_id", "date"),
      all.x = TRUE
    )

    ## Distance BIO
    dist_dir$distance_bio <- as.numeric(dist_dir$distance_bio)
    tmp_mean <- stats::aggregate(
      distance_bio ~ ID,
      data = dist_dir,
      FUN = function(df) mean(df, na.rm = TRUE)
    )
    dist_dir$distance_bio_mean <- tmp_mean[
      match(dist_dir$ID, tmp_mean$ID),
      2
    ]

    distance_bio_sd <- tapply(
      as.numeric(dist_dir$distance_bio),
      dist_dir$ID,
      stats::sd,
      na.rm = TRUE
    )
    dist_dir$distance_bio_sd <- distance_bio_sd[dist_dir$ID]

    distance_bio_median <- tapply(
      as.numeric(dist_dir$distance_bio),
      dist_dir$ID,
      stats::median,
      na.rm = TRUE
    )
    dist_dir$distance_bio_median <- distance_bio_median[dist_dir$ID]

    distance_bio_min <- tapply(
      as.numeric(dist_dir$distance_bio),
      dist_dir$ID,
      min,
      na.rm = TRUE
    )
    dist_dir$distance_bio_min <- distance_bio_min[dist_dir$ID]

    distance_bio_max <- tapply(
      as.numeric(dist_dir$distance_bio),
      dist_dir$ID,
      max,
      na.rm = TRUE
    )
    dist_dir$distance_bio_max <- distance_bio_max[dist_dir$ID]

    ## Since degree are a circular unit (after 360=0 comes 1,2...) we need a
    # special way to compute it :
    # Conversion to circular angle :

    dist_dir$azimut_bio_num <- as.numeric(dist_dir$azimut_bio)
    dist_dir$azimut_bio_circ <- ifelse(
      !is.na(dist_dir$azimut_bio_num),
      circular::circular(
        dist_dir$azimut_bio_num,
        units = "degrees",
        template = "geographics"
      ),
      NA
    )

    azimut_bio_mean <- tapply(
      dist_dir$azimut_bio_circ,
      dist_dir$ID,
      circular::mean.circular,
      na.rm = TRUE
    )
    dist_dir$azimut_bio_mean <- (azimut_bio_mean[dist_dir$ID] + 360) %% 360

    azimut_bio_med <- tapply(
      as.numeric(dist_dir$azimut_bio_circ),
      dist_dir$ID,
      circular::median.circular,
      na.rm = TRUE
    )
    dist_dir$azimut_bio_med <- (azimut_bio_med[dist_dir$ID] + 360) %% 360

    tmp_sd <- stats::aggregate(
      azimut_bio_circ ~ ID,
      data = dist_dir,
      FUN = azimut_sd_fun
    )
    dist_dir$azimut_bio_sd <- (tmp_sd[match(dist_dir$ID, tmp_sd$ID), 2] +
      360) %%
      360

    ## Distance STRAIT
    dist_dir$distance_strait <- as.numeric(dist_dir$distance_strait)
    tmp_mean <- stats::aggregate(
      distance_strait ~ ID,
      data = dist_dir,
      FUN = function(df) mean(df, na.rm = TRUE)
    )
    dist_dir$distance_strait_mean <- tmp_mean[
      match(dist_dir$ID, tmp_mean$ID),
      2
    ]

    distance_strait_sd <- tapply(
      as.numeric(dist_dir$distance_strait),
      dist_dir$ID,
      stats::sd,
      na.rm = TRUE
    )
    dist_dir$distance_strait_sd <- distance_strait_sd[dist_dir$ID]

    distance_strait_median <- tapply(
      as.numeric(dist_dir$distance_strait),
      dist_dir$ID,
      stats::median,
      na.rm = TRUE
    )
    dist_dir$distance_strait_median <- distance_strait_median[dist_dir$ID]

    distance_strait_min <- tapply(
      as.numeric(dist_dir$distance_strait),
      dist_dir$ID,
      min,
      na.rm = TRUE
    )
    dist_dir$distance_strait_min <- distance_strait_min[dist_dir$ID]

    distance_strait_max <- tapply(
      as.numeric(dist_dir$distance_strait),
      dist_dir$ID,
      max,
      na.rm = TRUE
    )
    dist_dir$distance_strait_max <- distance_strait_max[dist_dir$ID]

    ## Since degree are a circular unit (after 360=0 comes 1,2...) we need a
    # special way to compute it :
    # Conversion to circular angle :

    dist_dir$azimut_strait_num <- as.numeric(dist_dir$azimut_strait)
    dist_dir$azimut_strait_circ <- ifelse(
      !is.na(dist_dir$azimut_strait_num),
      circular::circular(
        dist_dir$azimut_strait_num,
        units = "degrees",
        template = "geographics"
      ),
      NA
    )

    azimut_strait_mean <- tapply(
      dist_dir$azimut_strait_circ,
      dist_dir$ID,
      circular::mean.circular,
      na.rm = TRUE
    )
    dist_dir$azimut_strait_mean <- (azimut_strait_mean[dist_dir$ID] + 360) %%
      360

    azimut_strait_med <- tapply(
      as.numeric(dist_dir$azimut_strait_circ),
      dist_dir$ID,
      circular::median.circular,
      na.rm = TRUE
    )
    dist_dir$azimut_strait_med <- (azimut_strait_med[dist_dir$ID] + 360) %% 360

    tmp_sd <- stats::aggregate(
      azimut_strait_circ ~ ID,
      data = dist_dir,
      FUN = azimut_sd_fun
    )
    dist_dir$azimut_strait_sd <- (tmp_sd[match(dist_dir$ID, tmp_sd$ID), 2] +
      360) %%
      360

    # Delete column that refer to daily value and not to value compute on all
    # the event :
    dist_dir[,
      c(
        'date',
        'distance_bio',
        'azimut_bio',
        'azimut_bio_num',
        'azimut_bio_circ',
        'distance_strait',
        'azimut_strait',
        'azimut_strait_num',
        'azimut_strait_circ',
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
  message("this combination of argument is not handle by the function")
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

  val <- suppressWarnings(as.numeric(circular::sd.circular(az_group)))

  if (is.nan(val)) 0 else val
}
