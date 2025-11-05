#'Compute metrics spatial features of the extreme event (only in surface, 2D)
#'
#'@description This function give the morphological metrics for an extreme 
#' event, Corrected_raster is the SpatRaster you want to analyze (it needs to be
#' binarise with the mhw function)
#'@param Corrected_raster the SpatRaster you want to analyse
#'@param start_date and @param end_date allows to perform the analysis on a 
#' specific  time_frame, this allows to save computation time.
#' @param per_pix use TRUE if you want a list with one dt per pixel as an output
#'@note Your SpatRaster need to be binarised with the mhw function before you
#' put it in this function 
#' BEE.calc.metrics_morpho is not designed to work on 4D data (time+spatial 3D).
#' 
#-------------------------------------------------------------------------------

#start_date \<- "2022-01-01" ; end_date \<- "2023-12-31"



BEE.calc.metrics_morpho <- function(Corrected_rasters,
                                    start_date = NULL,
                                    end_date = NULL,
                                    per_pix = FALSE) {
# Retrive the dataframe if required
  if (!is.null(start_date) | !is.null(end_date)) {
    ## Check that date format matches
    format_start_date <- lubridate::guess_formats(start_date,
                                       orders = c("dmy", "ymd", "mdy"))
    format_end_date <- lubridate::guess_formats(end_date, orders = c("dmy", "ymd", "mdy"))
    format_layers_names <- lubridate::guess_formats(names(Corrected_rasters[[1]]),
                                         orders = c("dmy", "ymd", "mdy"))
    same_format <- all(format_start_date %in% format_layers_names)
    same_format <- c(same_format,
                     all(format_end_date %in% format_layers_names))
    if (any(same_format == FALSE)) {
      warning(
        "The format of start_date or end_date does not match the format of the
        dates in Corrected_raster. Please run names(Corrected_raster[[1]]) to 
        check the required format."
      )
      stop()
    }
    else{
      ## Retrive
      start_date <- ifelse(is.null(start_date),
                           names(Corrected_rasters[[1]]),
                           start_date)
      end_date <- ifelse(is.null(end_date), 
                         names(Corrected_rasters[[terra::nlyr(Corrected_rasters)]]), 
                         end_date)
      rasters <- Corrected_rasters[[which(names(Corrected_rasters) >= start_date
                                    & names(Corrected_rasters) <= end_date)]]
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
    fun = function(x)
      sum(is.na(x))
  ))
  if (length(nb_NA) > 1) {
    message(
      "The number of NA pixels varies through the layers, this may leads to 
      unexpected results in the output of that funciton."
    )
  }
  nb_pixel_studied <- terra::ncell(rasters[[1]]) - nb_NA
  dist_list <- lapply(patch_list, function(x) {
    # x <- patch_list[[1]]
    # Create a df with on row per pixel
    vals <- terra::values(x)
    cell_size <- terra::values(terra::cellSize(x))[, 1]
    boundary <- terra::boundaries(x, directions = 8)
    coordonates <- xyFromCell(x, 1:terra::ncell(x)) #each pixel coordinates
    d <- terra::time(x) # PAS CERTAINE QUE ÇA SOIT LA FONCTION DE TERRA QUI SOIT
    #APPELLEE ICI
    
    data <- data.table::data.table(
      x = coordonates[, 1],
      y = coordonates[, 2],
      values = vals,
      # value of the patch
      cell_size = cell_size,
      boundary = terra::values(terra::boundaries(x, directions = 8))
    ) %>% dplyr::rename(patch_id = values.patches, boundary = boundary.patches)
    
    data[, length_boundary := ifelse(is.nan(patch_id),
                                     NA_real_, sum(boundary),
                                     by = patch_id)]
    data[, patch_area := ifelse(is.nan(patch_id),
                                NA_real_, sum(cell_size)),
         by = patch_id]
    data[, core_area := ifelse(is.nan(patch_id), NA_real_,
                               sum(cell_size * (boundary == 0))),
         by = patch_id]
    data[, core_area_index := core_area / as.numeric(nb_pixel_studied)]
    data[, n_pixel := ifelse(is.nan(patch_id), NA_real_, .N),
         by = patch_id]
    data[, pixel_EE_tot := ifelse(is.nan(patch_id), NA_real_,
                                  sum(unique(n_pixel), na.rm = T))]
    data[, cover_percent := n_pixel / as.numeric(nb_pixel_studied),
         by = patch_id]
    data[, core_pixel := n_pixel - length_boundary]
    date <- terra::time(x)
    data[, ID := ifelse(is.nan(patch_id), NA, paste0(date, "_", patch_id))]
    data[, pixel_id := .I]
    polygon <- terra::as.polygons(x, aggregate = T, round = F)
    perim <- data.table::data.table(perimeter = terra::perim(polygon))
    
    id <- data[!is.nan(patch_id), .(patch_id = unique(patch_id))]
    
    perim[, patch_id := id]
    
    data <- merge(data, perim, by = "patch_id", all.x = TRUE)
    
    data <- data %>%
      dplyr::select(-cell_size)
    
    if (length(polygon) == 0) {
      data[, centroid_x := rep(NA_real_, nrow(data))]
      data[, centroid_y := rep(NA_real_, nrow(data))]
      data[, perim_area_ratio := length_boundary / n_pixel, by = patch_id]
      data[, shape_index := (0.25 * length_boundary / (n_pixel)^0.5),
           by = patch_id]
      data[, fract_corel_dim := 
(2 * log(0.25 * length_boundary, base = exp(1)) / log(n_pixel, base = exp(1))),
by = patch_id]
      data[, date := d]
      
    } else {
      centro <- terra::centroids(polygon, TRUE)
      centro <- terra::extract(x, centro, xy = TRUE)
      centro <- centro %>%
        dplyr::select(patches, x, y) %>%
        dplyr::rename(
          centroid_x = x,
          centroid_y = y,
          patch_id = patches
        )
      data <- merge(data, centro, by = "patch_id", all.x = TRUE)
      
      data[, perim_area_ratio := length_boundary / n_pixel, by = patch_id]
      data[, shape_index := (0.25 * length_boundary / (n_pixel)^0.5),
           by = patch_id]
      data[, fract_corel_dim := (2 * log(0.25 * length_boundary,
                            base = exp(1)) / log(n_pixel, base = exp(1))),
           by = patch_id]
      data[, date := d]
    }
    return(data)
  }) # one dt per pixel
  
  names(dist_list) <- terra::time(rasters)
  
  # Summarise by patch
  data_summarised <- lapply(dist_list, function(x) {
    # x <- dist_list[[1]]
    x <- unique(x, by = "patch_id")
    x[, ':='(x = NULL, y = NULL)]
    return(x)
  })
  
  data_summarised$pixel_id <- NULL
  # Add a list of patch_ID per pixel (/!\ there are several patch_ID per EE 
  # because we cannot do spatio-temporal analysis here)
  if (per_pix == TRUE) {
    pixels_nb <- seq(1, nrow(dist_list[[1]]))
    dist_tab <- dplyr::bind_rows(dist_list)
    data.table::setDT(dist_tab)
    dist_tab <- dist_tab %>%
      dplyr::mutate(valid_patch_ID = ifelse(!is.nan(patch_id), date, NA)) # keeps only
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
