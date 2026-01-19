#'Compute metrics spatial features of the extreme event (only in surface, 2D)
#'
#'@description This function give the morphological metrics for an extreme
#' event, Corrected_rasters is the SpatRaster you want to analyze (it needs to be
#' binarise with the mhw function)
#' @param Corrected_rasters the SpatRaster you want to analyse
#' @param start_date allows to perform the analysis on a
#' specific  time_frame, this allows to save computation time.
#' @param end_date allows to perform the analysis on a
#' specific  time_frame, this allows to save computation time.
#' @param per_pix use TRUE if you want a list with one dt per pixel as an output
#' @param crs to get accurate length and area, data in longitude latitude must
#' be converted into meters data that take in account earth sphericity, the
#' conversion depend on studdied region, thus you must specify the one the most
#' appropriated for your data. See help here to find a suitable EPSG you can
#' check on https://epsg.io/
#'@note Your SpatRaster need to be binarised with the mhw function before you
#' put it in this function
#' BEE.calc.metrics_morpho is not designed to work on 4D data (time+spatial 3D).
#'
#-------------------------------------------------------------------------------

#start_date <- "2024-06-01" ; end_date <- "2024-09-31" ;per_pix=TRUE ;
#crs = "EPSG:2154" # for europe

BEE.calc.metrics_morpho <- function(
  Corrected_rasters,
  start_date = NULL,
  end_date = NULL,
  per_pix = FALSE,
  crs = NULL
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
      names(Corrected_rasters[[1]]),
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
        names(Corrected_rasters[[1]]),
        start_date
      )
      end_date <- ifelse(
        is.null(end_date),
        names(Corrected_rasters[[terra::nlyr(Corrected_rasters)]]),
        end_date
      )
      rasters <- Corrected_rasters[[which(
        names(Corrected_rasters) >= start_date &
          names(Corrected_rasters) <= end_date
      )]]
    }
  }
  if (is.null(start_date) & is.null(end_date)) {
    ## Adjust variable name to match the rest of the function
    rasters <- Corrected_rasters
  }

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

  nb_NA <- unique(terra::global(
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
  nb_pixel_studied <- terra::ncell(rasters[[1]]) - nb_NA

  dist_list <- lapply(patch_list, function(x) {
    #  # 1 only NA # 41 one patch et 45 2 patches,
    # 100 4 patch dont un gros en premier #  x <- patch_list[[100]]
    # Create a df with on row per pixel
    vals <- terra::values(x)
    cell_size <- terra::values(terra::cellSize(x))[, 1]
    boundary <- terra::boundaries(x, directions = 8)
    coordonates <- terra::xyFromCell(x, 1:terra::ncell(x)) #each pixel coordinates
    d <- terra::time(x)

    data <- data.table::data.table(
      x = coordonates[, 1],
      y = coordonates[, 2],
      values = vals,
      # value of the patch
      cell_size = cell_size,
      boundary = terra::values(terra::boundaries(x, directions = 8))
    )
    data.table::setnames(
      data,
      old = c("values.patches", "boundary.patches"),
      new = c("patch_id", "bordering")
    )

    # Split the 'data' by patch id to process faster :
    sp <- split(data, data$patch_id)

    # add metrics :
    patch_area <- vapply(
      sp,
      function(df) if (any(is.na(df$patch_id))) NA_real_ else sum(df$cell_size),
      numeric(1)
    )

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
    data$patch_area <- patch_area[match(data$patch_id, names(patch_area))]
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
      function(df) as.numeric(df$n_pixel[1]) / as.numeric(nb_pixel_studied),
      numeric(1)
    )

    data$cover_percent <- cover_percent[match(
      data$patch_id,
      names(cover_percent)
    )]
    data$core_pixel <- data$n_pixel - data$bordering

    data$date <- rep(d, nrow(data))

    data$ID <- ifelse(
      is.nan(data$patch_id),
      NA,
      paste0(seq_len(nrow(data)), "_", d, "_", data$patch_id)
    )
    polygon <- terra::as.polygons(x, aggregate = T, round = F)
    perim_m <- data.table::data.table(perimeter = terra::perim(polygon))

    perim_m$patch_id <- unique(data$patch_id[!is.na(data$patch_id)]) # in meters

    data$perimeter <- perim_m$perimeter[match(data$patch_id, perim_m$patch_id)]

    data <- data |>
      dplyr::select(-cell_size)

    if (length(polygon) == 0) {
      data$centroid_x <- rep(NA_real_, nrow(data))
      data$centroid_y <- rep(NA_real_, nrow(data))
      data$perim_pixel_ratio <- rep(NA_real_, nrow(data))
      data$perim_area_ratio_m <- rep(NA_real_, nrow(data))
      data$shape_index <- rep(NA_real_, nrow(data))
      data$fractal_cor_index <- rep(NA_real_, nrow(data))
      data$perim_area_ratio_m <- rep(NA_real_, nrow(data))
      data$circle_ratio_index <- rep(NA_real_, nrow(data))
      data$contiguity_index <- rep(NA_real_, nrow(data))
    } else {
      centro <- terra::centroids(polygon, TRUE)
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
      Perim_area_ratios_m <- landscapemetrics::lsm_p_para(x)
      data$perim_area_ratio_m <- Perim_area_ratios_m$value[match(
        data$patch_id,
        Perim_area_ratios_m$id
      )]
      #Shape Index :
      x <- raster::raster(x)
      Shape_indexes <- landscapemetrics::lsm_p_shape(x)
      data$shape_index <- Shape_indexes$value[match(
        data$patch_id,
        Shape_indexes$id
      )]
      #Fractal correlation :
      Fractal_cor_indexes <- landscapemetrics::lsm_p_frac(x)
      data$fractal_cor_index <- Fractal_cor_indexes$value[match(
        data$patch_id,
        Fractal_cor_indexes$id
      )]

      #Ratio btw patch shape and circle
      circle_ratio <- data.frame(
        patch_id = as.numeric(names(sp)[names(sp) != "NaN"]),
        circle_ratio_index = as.numeric(patch_circle_ratio(polygon, crs = crs))
      )
      data <- merge(
        data,
        circle_ratio,
        by = "patch_id",
        all.x = TRUE,
        sort = FALSE
      )

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

  # Summarise by patch
  data_summarised <- lapply(dist_list, function(x) {
    x <- unique(x, by = "patch_id")
    x[c("x", "y")] <- NULL
    x
  })

  data_summarised$pixel_id <- NULL
  # Add a list of patch_ID per pixel (/!\ there are several patch_ID per EE
  # because we cannot do spatio-temporal analysis here)
  if (per_pix == TRUE) {
    pixels_nb <- seq(1, nrow(dist_list[[1]]))
    dist_tab <- dplyr::bind_rows(dist_list)
    data.table::setDT(dist_tab)
    dist_tab <- dist_tab |>
      dplyr::mutate(
        valid_patch_ID = ifelse(!is.nan(patch_id), as.Date(dist_tab$date), NA)
      ) # keeps only
    # the patch that represents an EE
    patch_list <- terra::rast(patch_list)
    names(patch_list) <- terra::time(rasters)
    output <- list(dist_tab, patch_list)
    return(output)
  }
  patch_list <- terra::rast(patch_list)
  names(patch_list) <- terra::time(rasters)
  data_summarised <- data.table::rbindlist(data_summarised)
  output <- list(data_summarised, patch_list)
  return(output)
}


#' Compute de smallest circle in which a patch fits and compare patch are to circle area
#'
#' @noRd
#'
patch_circle_ratio <- function(polygon, crs = crs) {
  # WORK one raster per one raster, not on spat raster
  patch_sf <- st_as_sf(polygon, crs = crs)
  patch_sf <- sf::st_transform(patch_sf, crs = crs)
  patch_area <- terra::expanse(polygon, unit = "m")
  circle <- sf::st_minimum_bounding_circle(patch_sf)
  circle_area <- as.numeric(sf::st_area(circle))
  patch_area / circle_area
}
