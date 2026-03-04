#'Compute metrics spatial features of the extreme event (only in surface, 2D)
#'
#' @description
#'  This function give the morphological metrics for an extreme event,
#'  extreme_event_spatraster is the SpatRaster you want to analyze (it needs to be
#'  binarise with the mhw function)
#'
#' @param extreme_event_spatraster :
#'  The SpatRaster you want to analyse (first output of
#'  *BEE.id.extreme_events()*).
#' @param start_date :
#'  Allows to perform the analysis on a specific  time_frame, this allows to
#'  save computation time.
#' @param end_date :
#'  Allows to perform the analysis on a specific  time_frame, this allows to
#'  save computation time.
#' @param per_pix :
#'  Use TRUE if you want a list with one dt per pixel as an output. Use TRUE to
#'  be able to use the output in BEE.data.merge().
#'
#' @return Units related to area are in m².
#'
#' @note
#'  BEE.calc.metrics_morpho is not designed to work on 4D data
#'  (time+spatial 3D).
#'
#' @examples
#' # Load example data:
#' library(BioExtremeEvent)
#'file_name_1 <- system.file(file.path("extdata",
#'                                      "binarized_corrected_spatraster.tiff"),
#'                                      package = "BioExtremeEvent")
#'extreme_event_spatraster <- terra::rast(file_name_1)
#'
#' # Get metrics value per pixel x day:
#' list_morpho_metrics <- BEE.calc.metrics_morpho(extreme_event_spatraster,
#' start_date = "2024-05-01",
#' end_date = "2024-11-30",
#' per_pix=TRUE)
#'
#' dataframe_list <- list_morpho_metrics[[1]] # One df per pixel
#' patch_spatraster <- list_morpho_metrics[[2]]
#'
#'
#' @export
#'
#-------------------------------------------------------------------------------

# start_date <- "2024-07-01" ; end_date <- "2024-08-31" ;per_pix=TRUE ;
# extreme_event_spatraster <- Corrected_rasters

BEE.calc.metrics_morpho <- function(
  extreme_event_spatraster,
  start_date = NULL,
  end_date = NULL,
  per_pix = FALSE
) {
  # Retrive the dataframe if required
  if (!is.null(start_date) | !is.null(end_date)) {
    ## Check that date format matches
    format_start_date <- lubridate::guess_formats(
      start_date,
      orders = c("dmy", "ymd", "mdy")
    )
    format_end_date <- lubridate::guess_formats(
      end_date,
      orders = c("dmy", "ymd", "mdy")
    )
    format_layers_names <- lubridate::guess_formats(
      names(extreme_event_spatraster[[1]]),
      orders = c("dmy", "ymd", "mdy")
    )
    same_format <- all(format_start_date %in% format_layers_names)
    same_format <- c(same_format, all(format_end_date %in% format_layers_names))
    if (any(same_format == FALSE)) {
      warning(
        "The format of start_date or end_date does not match the format of the
        dates in Corrected_raster. Please run names(Corrected_raster[[1]]) to 
        check the required format."
      )
      stop()
    } else {
      ## Retrive
      start_date <- ifelse(
        is.null(start_date),
        names(extreme_event_spatraster[[1]]),
        start_date
      )
      end_date <- ifelse(
        is.null(end_date),
        names(extreme_event_spatraster[[terra::nlyr(
          extreme_event_spatraster
        )]]),
        end_date
      )
      rasters <- extreme_event_spatraster[[which(
        names(extreme_event_spatraster) >= start_date &
          names(extreme_event_spatraster) <= end_date
      )]]
    }
  }
  if (is.null(start_date) & is.null(end_date)) {
    ## Adjust variable name to match the rest of the function
    rasters <- extreme_event_spatraster
  }

  nb_NA <- unique(terra::global(
    # number of NA in each raster
    rasters,
    fun = function(x) {
      sum(is.na(x))
    }
  ))
  if (length(nb_NA) > 1) {
    message(
      "The number of NA pixels varies through the layers, this may leads to 
      unexpected results in the output of that funciton."
    )
  }

  ################################ CODE  #######################################

  patch_list <- lapply(
    rasters,
    FUN = function(rasters) {
      patch <- terra::patches(
        rasters,
        directions = 8,
        zeroAsNA = T,
        allowGaps = T
      )
      return(patch)
    }
  )

  non_NA_pixels <- which(
    terra::app(rasters, fun = function(x) all(!is.na(x)))[] == 1
  )
  cell_size <- terra::values(terra::cellSize(rasters, unit = "m"))
  area_studied <- sum(cell_size[non_NA_pixels, 1]) # m

  dist_list <- lapply(patch_list, function(x) {
    #  # 1 only NA # 41 one patch et 45 2 patches,
    # 100 4 patch dont un gros en premier #  x <- patch_list[[8]]
    # Create a df with on row per pixel
    vals <- terra::values(x)

    # previous line : brakets to get vector instead of matrix
    boundary <- terra::boundaries(x, directions = 8)
    coordonates <- terra::xyFromCell(x, 1:terra::ncell(x)) #each pixel coordinates
    d <- terra::time(x)

    data <- data.table::data.table(
      pixel_id = seq(1, length(vals), 1),
      x = coordonates[, 1],
      y = coordonates[, 2],
      values = vals,
      # value of the patch
      cell_size = cell_size[, 1],
      boundary = terra::values(terra::boundaries(x, directions = 8))
    )
    data.table::setnames(
      data,
      old = c("values.patches", "boundary.patches"),
      new = c("patch_id", "bordering")
    )

    polygon <- terra::as.polygons(x, aggregate = T, round = F)
    patch_area <- terra::expanse(polygon, unit = "m")
    patch_table <- data.frame(
      patch_id = polygon$patches,
      patch_area = patch_area
    )
    data$patch_area <- patch_table$patch_area[
      match(data$patch_id, patch_table$patch_id)
    ]
    # Split the 'data' by patch id to process faster :
    sp <- split(data, data$patch_id)

    core_area <- vapply(
      sp,
      function(df) {
        if (any(is.na(df$patch_id))) {
          NA_real_
        } else {
          sum(as.numeric(df$cell_size) * as.numeric(df$bordering == 0))
        }
      },
      numeric(1)
    )

    # Remapping to data
    data$core_area <- core_area[match(data$patch_id, names(core_area))]

    sp <- split(data, data$patch_id)
    # add metrics :
    core_area_index <- vapply(
      sp,
      function(df) {
        if (length(df$core_area) == 0 || is.nan(df$core_area[1])) {
          NA_real_
        } else {
          as.numeric(df$core_area[1]) / as.numeric(df$patch_area[1])
        }
      },
      #only one  value per patch so corea_area[1] is ok and it prevent bugs linked to the
      # linked to the numeric(1) bellow
      numeric(1)
    )
    n_pixel <- vapply(
      sp,
      function(df) {
        if (all(is.nan(df$patch_id))) NA_real_ else nrow(df)
      },
      numeric(1)
    )

    # Remapping to data
    data$core_area_index <- core_area_index[match(
      data$patch_id,
      names(core_area_index)
    )]
    data$n_pixel <- n_pixel[match(data$patch_id, names(n_pixel))]
    data$pixel_EE_tot <- sum(n_pixel, na.rm = TRUE)

    # add metrics
    sp <- split(data, data$patch_id)
    cover_percent <- vapply(
      sp,
      function(df) as.numeric(df$patch_area[1]) / as.numeric(area_studied),
      numeric(1)
    )

    n_core_pixel <- vapply(
      sp,
      function(df) as.numeric(df$n_pixel[1]) - sum(df$bordering),
      numeric(1)
    )

    data$cover_percent <- cover_percent[match(
      data$patch_id,
      names(cover_percent)
    )]
    data$n_core_pixel <- n_core_pixel[match(
      data$patch_id,
      names(n_core_pixel)
    )]

    data$date <- rep(d, nrow(data))

    data$ID <- ifelse(
      is.nan(data$patch_id),
      NA,
      paste0(seq_len(nrow(data)), "_", d, "_", data$patch_id)
    )

    perim_m <- data.table::data.table(perimeter = terra::perim(polygon))

    perim_m$patch_id <- unique(data$patch_id[!is.na(data$patch_id)]) # in meters

    data$perimeter <- perim_m$perimeter[match(data$patch_id, perim_m$patch_id)]

    if (length(polygon) == 0) {
      data$centroid_x <- rep(NA_real_, nrow(data))
      data$centroid_y <- rep(NA_real_, nrow(data))
      data$perim_pixel_ratio <- rep(NA_real_, nrow(data))
      data$perim_area_ratio <- rep(NA_real_, nrow(data))
      data$shape_index <- rep(NA_real_, nrow(data))
      data$fractal_cor_index <- rep(NA_real_, nrow(data))
      data$contiguity_index <- rep(NA_real_, nrow(data))
      data$ell_dir <- rep(NA_real_, nrow(data))
      data$ell_crop_area_m2 <- rep(NA_real_, nrow(data))
      data$ell_tot_area_m2 <- rep(NA_real_, nrow(data))
      data$ell_long_axis <- rep(NA_real_, nrow(data))
      data$ell_short_axis <- rep(NA_real_, nrow(data))
      data$patch_ell_ratio <- rep(NA_real_, nrow(data))
    } else {
      centro <- terra::centroids(polygon)
      centro <- terra::extract(x, centro, xy = TRUE)
      centro <- centro |>
        dplyr::select(patches, x, y) |>
        dplyr::rename(
          centroid_x = x,
          centroid_y = y,
          patch_id = patches
        )
      data <- merge(data, centro, by = "patch_id", all.x = TRUE)
      #Perimeter ratio area in pixels :
      perim_pixel_ratio <- vapply(
        sp, # each df contains all info about a patch, thus noudary_length and
        # n_pixel have the same value for a given patch / for all rows of df
        function(df) sum(df$bordering) / unique(df$n_pixel),
        numeric(1)
      )
      data$perim_pixel_ratio <- perim_pixel_ratio[match(
        data$patch_id,
        names(perim_pixel_ratio)
      )]
      #Perimeter ratio area in meters :
      sp <- split(data, data$patch_id)
      perim_area_ratios_m <- vapply(
        sp,
        function(df) sum(df$perimeter) / unique(df$patch_area),
        numeric(1)
      )
      data$perim_area_ratio <- perim_area_ratios_m[match(
        data$patch_id,
        names(perim_area_ratios_m)
      )]
      #Shape Index :
      shape_index <- vapply(
        sp,
        function(df) 0.25 * df$perimeter[1] / sqrt(df$patch_area[1]),
        numeric(1)
      )
      data$shape_index <- shape_index[match(
        data$patch_id,
        names(shape_index)
      )]
      #Fractal correlation :
      fractal_cor_indexes <- vapply(
        sp,
        function(df) 2 * log(0.25 * df$perimeter[1]) / log(df$patch_area[1]),
        numeric(1)
      )
      data$fractal_cor_index <- fractal_cor_indexes[match(
        data$patch_id,
        names(fractal_cor_indexes)
      )]

      #Ratio btw patch shape and ellipsoide
      ellipses <- min_ellipse_from_polygon(x, data)
      data <- merge(data, ellipses, by = "patch_id", all.x = TRUE)

      #Contiguity :
      Contiguity_indexes <- landscapemetrics::lsm_p_contig(x)
      data$contiguity_index <- Contiguity_indexes$value[match(
        data$patch_id,
        Contiguity_indexes$id
      )]
    }
    return(data)
  })

  names(dist_list) <- terra::time(rasters)

  # Add a list of patch_ID per pixel (/!\ there are several patch_ID per EE
  # because we cannot do spatio-temporal analysis here)
  if (per_pix == TRUE) {
    dist_tab <- dplyr::bind_rows(dist_list)
    # Withdraw lines of pixels that are always NA :
    na_cells <- which(is.na(terra::values(rasters[[1]])))
    dist_tab <- dist_tab[!(dist_tab$pixel_id %in% na_cells), ]
    # split with one df per pixel
    data.table::setDF(dist_tab)
    dist_tab <- split(dist_tab, dist_tab$pixel_id)
    # prepare outputs:
    patch_list <- terra::rast(patch_list)
    names(patch_list) <- terra::time(rasters)
    output <- list(dist_tab, patch_list)
    return(output)
  }
  # Summarise by patch
  dist_list <- lapply(dist_list, data.table::setDF)
  data_summarised <- lapply(dist_list, function(x) {
    x <- x[!is.na(x[["patch_id"]]), ] # enlever NA / NaN
    x[, c(
      "x",
      "y",
      "pixel_id",
      "bordering",
      "centroid_x",
      "centroid_y"
    )] <- NULL
    x <- x[!duplicated(x[["patch_id"]]), ] # garder la 1re occurrence
    x
  })

  patch_list <- terra::rast(patch_list)
  names(patch_list) <- terra::time(rasters)

  data_summarised <- data.table::rbindlist(data_summarised, use.names = TRUE)

  output <- list(data_summarised, patch_list)
  return(output)
}


#' Compute de smallest ellipsoide in which a patch fits and compute area
#'
#' @noRd
#'

#Minimum ellipse:

min_ellipse_from_polygon <- function(x, data2) {
  patch_ids <- unique(data2$patch_id[
    !is.na(data2$patch_id) & !is.nan(data2$patch_id)
  ])

  ell_dir <- vector()
  ell_crop_area_m2 <- vector()
  ell_tot_area_m2 <- vector()
  patch_ell_ratio <- vector()
  ell_long_axis <- vector()
  ell_short_axis <- vector()
  i = 1

  for (p in patch_ids) {
    #p <- patch_ids[1]
    # pixels coordinates from patch p
    data_p <- data2[
      !is.na(data2$patch_id) &
        data2$patch_id == p &
        !is.na(data2$x) &
        !is.na(data2$y),
    ]

    #matrix of points: (pixels coordinates that must goes into the ellipse)
    coords_m <- as.matrix(
      data_p[, c("x", "y")]
    )
    if (nrow(coords_m) == 1) {
      ell_dir[i] <- NA
      ell_crop_area_m2[i] <- 0
      ell_tot_area_m2 <- 0
      patch_ell_ratio[i] <- NA
      ell_long_axis[i] <- 0
      ell_short_axis[i] <- 0
      i <- 1 + i
    } else if (
      nrow(coords_m) >= 2 &
        (length(unique(coords_m[, 1])) == 1 |
          length(unique(coords_m[, 2])) == 1)
    ) {
      # the patch formes a line
      if (length(unique(coords_m[, 1])) == 1) {
        #fixed longitude
        ell_dir[i] <- 0
        patch_ell_ratio[i] <- NA
        ell_crop_area_m2[i] <- 0
        ell_tot_area_m2 <- 0
        #distance
        data_p_deg <- data2[data2$patch_id == p, c("x", "y")]
        data_p_deg_v <- terra::vect(data_p_deg, crs = terra::crs(x))

        ell_long_axis[i] <- max(terra::distance(data_p_deg_v)) / 2
        ell_short_axis[i] <- 0
        i <- 1 + i
      }
      if (length(unique(coords_m[, 2])) == 1) {
        #fixed longitude
        ell_dir[i] <- 90
        patch_ell_ratio[i] <- NA
        ell_crop_area_m2[i] <- 0
        ell_tot_area_m2 <- 0
        #distance
        data_p_deg <- data2[data2$patch_id == p, c("x", "y")]
        data_p_deg_v <- terra::vect(data_p_deg, crs = terra::crs(x))

        ell_long_axis[i] <- max(terra::distance(data_p_deg_v)) / 2
        ell_short_axis[i] <- 0
        i <- 1 + i
      }
    } else if (
      nrow(coords_m) >= 2 &
        (length(unique(coords_m[, 1])) == nrow(coords_m) &
          length(unique(coords_m[, 2])) == nrow(coords_m))
    ) {
      #diagonal
      ell_dir[i] <- geosphere::bearing(coords_m[1, ], coords_m[2, ])
      patch_ell_ratio[i] <- NA
      ell_crop_area_m2[i] <- 0
      ell_tot_area_m2 <- 0
      #distance
      data_p_deg <- data2[data2$patch_id == p, c("x", "y")]
      data_p_deg_v <- terra::vect(data_p_deg, crs = terra::crs(x))

      ell_long_axis[i] <- max(terra::distance(data_p_deg_v)) / 2
      ell_short_axis[i] <- 0
      i <- 1 + i
    } else if (
      nrow(coords_m) >= 2 &
        length(unique(coords_m[, 1])) > 1 &
        length(unique(coords_m[, 2])) > 1
    ) {
      # /!\ In between those line, every output is in degrees and distances and
      # angle are wrong because cluster::ellipsoidhull assume that coordinates
      # are planer, in our case we just want the coordinates of the points that
      # formes the ellipse, thoses coordinates will be in degres, once transformed
      # transformed into a polygon, terra will deal with the lon/lat to get the
      # most precise estimation of the covered are. This will be more precise
      # than convert to metric crs and then compute ellipse and its area.
      #-------------------------------------------------------------------------
      #make ellipse:
      ell <- cluster::ellipsoidhull(coords_m)

      #extract ellipse caracteristics:
      ## center coordinates:
      center_m <- ell$loc
      ## ellipse matrix (covariance):
      cov_m <- ell$cov

      ##PCA ellispe: (to get axis)
      eigen_res <- eigen(cov_m)
      eigenvalues <- eigen_res$values
      eigenvectors <- eigen_res$vectors
      ##axis
      semi_major_m <- sqrt(ell$d2 * eigenvalues[1])
      semi_minor_m <- sqrt(ell$d2 * eigenvalues[2])

      ## Total area of the ellipse
      ellipse_area_m2 <- pi * semi_major_m * semi_minor_m

      ##ellispe orientation relatively to first axis in the matrix provided as data_p
      #angle btw longest axis and a latitude
      principal_vector <- eigenvectors[, 1]
      angle_rad <- atan2(
        principal_vector[2],
        principal_vector[1]
      )
      angle_deg <- angle_rad * 180 / pi # par rapport à l'axe des x dans mon repère

      #To compute ratio btw ellispe area and patch area we need to withdraw from
      #the ellipse the portions that fall outside of the raster:
      #convert the ellipse into a polygon:
      ##draw ellipse:
      t <- seq(0, 2 * pi, length.out = 300)
      x0 <- semi_major_m * cos(t)
      y0 <- semi_minor_m * sin(t)
      x_rot <- x0 * cos(angle_rad) - y0 * sin(angle_rad)
      y_rot <- x0 * sin(angle_rad) + y0 * cos(angle_rad)
      x_ellipse <- x_rot + center_m[1]
      y_ellipse <- y_rot + center_m[2]
      ##convert to polygon
      ellipse_mat <- cbind(x_ellipse, y_ellipse)
      ellipse_mat <- rbind(
        ellipse_mat,
        ellipse_mat[1, ]
      )
      #-------------------------------------------------------------------------
      ellipse_poly <- terra::vect(
        ellipse_mat,
        type = "polygons",
        crs = terra::crs(x)
      )

      extent_poly <- terra::as.polygons(terra::ext(x), crs = terra::crs(x))

      #get ellipse only inside raster:
      ellipse_area_crop <- terra::intersect(ellipse_poly, extent_poly)
      terra::crs(ellipse_area_crop) <- terra::crs(x)
      ellipse_area_in_raster_m <- terra::expanse(ellipse_area_crop, unit = "m")

      #get patch/ellipse ratio
      patch_area <- sum(data_p$cell_size)
      patch_ell_ratio_i <- patch_area / ellipse_area_in_raster_m

      # Now proprely deduce axis:
      #in raster (long axis can outside of raster but not short one)

      major_axis <- max(terra::distance(
        terra::crds(ellipse_area_crop),
        lonlat = T,
        unit = "m"
      )) /
        2
      # area=pi*major_axis_REAL*minor_axis :
      pts_ell_poly <- terra::crds(ellipse_poly)
      major_axis_real <- max(terra::distance(
        pts_ell_poly,
        lonlat = T,
        unit = "m"
      )) /
        2
      minor_axis <- terra::expanse(ellipse_poly, unit = "m") /
        (pi * major_axis_real)

      #get long axis azimut

      dist_mat <- as.matrix(terra::distance(
        pts_ell_poly,
        lonlat = T,
        unit = "m"
      ))
      max_idx <- which(dist_mat == max(dist_mat, na.rm = TRUE), arr.ind = TRUE)
      pt1 <- pts_ell_poly[max_idx[1, 1], ]
      pt2 <- pts_ell_poly[max_idx[1, 2], ]
      ell_azimut <- geosphere::bearing(pt1, pt2) # computed on ellipsoid :)

      #save results:
      ell_crop_area_m2[i] <- ellipse_area_in_raster_m
      ell_tot_area_m2[i] <- terra::expanse(ellipse_poly, unit = "m")
      patch_ell_ratio[i] <- patch_ell_ratio_i
      ell_long_axis[i] <- major_axis
      ell_short_axis[i] <- minor_axis
      ell_dir[i] <- ell_azimut
      i <- i + 1
    }
  }
  ellipse <- data.frame(
    patch_id = patch_ids,
    ell_dir = ell_dir,
    ell_crop_area_m2 = ell_crop_area_m2,
    ell_tot_area_m2 = ell_tot_area_m2,
    patch_ell_ratio = patch_ell_ratio,
    ell_long_axis = ell_long_axis,
    ell_short_axis = ell_short_axis
  )

  return(ellipse)
}
