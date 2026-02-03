#' Merge the outputs of the metrics functions.
#'
#'@description To merge the daily outputs of at least two datasets of the
#' following fonctions from BEE package : bee.calc.metrics_point() ;
#' bee.calc.metrics_morpho() ; bee.calc.escape(). This function can also
#' summarise metrics over different time periods: extreme events, weak periods,
#' monthly periods and yearly periods.
#'
#'@param data_metrics_point the output of the BEE.calc.metrics_point computed
#' using the argument group_by_event = FALSE.
#'
#'@param data_metrics_morpho the output of the BEE.calc.metrics_morpho
#'
#'@param data_escape the output of the BEE.calc.escape computed
#' using the argument group_by_event = FALSE.
#'
#'@param crs a METRICS crs that suits the studdied area
#'
#'@return A dataframe with the metrics of all the datasets provided in column,
#' if you haven't choose a daily resolution, a mean, median, variance, minimum
#' and maximum will be computed for each variables over the choosen time step.
#'
#'@examples
#' # TO BE ADDED
#'
#'@export
#'
#-------------------------------------------------------------------------------

# data_metrics_point <- points_metrics ;
# data_metrics_morpho <- list_morpho_metrics[[1]] ;
# data_escape <- dist_to_escape ; crs = "EPSG:3035" ; summarize_by = "day"

BEE.data.merge <- function(
  data_metrics_point = NULL,
  data_metrics_morpho = NULL,
  data_escape = NULL,
  crs
) {
  ################## FORMAT TESTS AND WARNINGS #################################

  # Check that the provided crs is valid:
  crs_valid <- tryCatch(
    {
      terra::crs(crs, describe = TRUE)
      TRUE # if the crs is recognize -> T
    },
    error = function(e) {
      FALSE # if error from terra::crs -> F
    }
  )
  if (crs_valid == FALSE) {
    warnings(
      "The crs you have provided is not valid. Check this webpage to 
    select a suitable crs code : https://crs-explorer.proj.org/?ignoreWorld=false&allowDeprecated=false&authorities=EPSG&activeTypes=PROJECTED_CRS&map=osm"
    )
  }

  # Identify which datasets have been provided :
  ## Is there at least two datasets provided ?
  check_1 <- paste0(
    is.null(data_metrics_point),
    is.null(data_metrics_morpho),
    is.null(data_escape)
  )

  if (check_1 == "TRUETRUETRUE") {
    warnings(
      "You haven't provided any data to be merged and summarised. Please
    provide at least two datasets."
    )
  }
  if (check_1 == "FALSETRUETRUE") {
    message(
      "As you have only provided the metrics_point dataset, 
    summarising through time is all that will be done."
    )
  }
  if (check_1 == "TRUEFALSETRUE") {
    message(
      "As you have only provided the metrics_morpho dataset, 
    summarising through time is all that will be done."
    )
  }
  if (check_1 == "TRUETRUEFALSE") {
    message(
      "As you have only provided the escape dataset, 
    summarising through time is all that will be done."
    )
  }

  ## Is it a daily resolution ?
  if (!is.null(data_metrics_point)) {
    dates_point <- as.Date(data_metrics_point[[1]]$date)
    if (
      length(dates_point) !=
        dates_point[length(dates_point)] - dates_point[1] + 1
    ) {
      warnings(
        "Some days are missing in data_metrics_point or the time 
      resolution is not daily. Please, make sure to use 'group_by_event = FALSE'
      when running BEE.calc.metrics_point"
      )
    }
  }
  if (!is.null(data_metrics_morpho)) {
    dates_morpho <- unique(data_metrics_morpho[[1]]$date)
    dates_morpho <- as.Date(dates_morpho)
    if (
      length(dates_morpho) !=
        dates_morpho[length(dates_morpho)] - dates_morpho[1] + 1
    ) {
      warnings(
        "Some days are missing in data_metrics_point or the time 
      resolution is not daily. Please, make sure to use 'group_by_event = FALSE'
      when running BEE.calc.metrics_point"
      )
    }
  }
  if (!is.null(data_escape)) {
    dates_escape <- unique(data_escape[[1]]$date)
    dates_escape <- as.Date(dates_escape)
    if (
      length(dates_escape) !=
        dates_escape[length(dates_escape)] - dates_escape[1] + 1
    ) {
      warnings(
        "Some days are missing in data_metrics_point or the time 
      resolution is not daily. Please, make sure to use 'group_by_event = FALSE'
      when running BEE.calc.metrics_point"
      )
    }
  }
  ## Are the time periods overlapping ?
  if (check_1 == "FALSEFALSETRUE") {
    shared <- length(intersect(dates_point, dates_morpho))
    if (shared == 0) {
      warnings(
        "data_metrics_point and data_metrics_morpho do not share the same
      time period (not a single day in common). Please provide data sharing some
      dates in common and comming from the same rasters."
      )
    }
    if (shared < length(dates_point)) {
      message(
        "data_metrics_point covers more dates than data_metrics_morpho, 
      the merging will only be done on the days they share."
      )
    }
    if (shared < length(dates_morpho)) {
      message(
        "data_metrics_morpho covers more dates than data_metrics_point, 
      the merging will only be done on the days they share."
      )
    }
  }
  if (check_1 == "FALSETRUEFALSE") {
    shared <- length(intersect(dates_point, dates_escape))
    if (shared == 0) {
      warnings(
        "data_metrics_point and data_escape do not share the same
      time period (not a single day in common). Please provide data sharing some
      dates in common and comming from the same rasters."
      )
    }
    if (shared < length(dates_point)) {
      message(
        "data_metrics_point covers more dates than data_escape, 
      the merging will only be done on the days they share."
      )
    }
    if (shared < length(data_escape)) {
      message(
        "data_escape covers more dates than data_metrics_point, 
      the merging will only be done on the days they share."
      )
    }
  }
  if (check_1 == "TRUEFALSEFALSE") {
    shared <- length(intersect(dates_morpho, dates_escape))
    if (shared == 0) {
      warnings(
        "data_metrics_morpho and data_escape do not share the same
      time period (not a single day in common). Please provide data sharing some
      dates in common and comming from the same rasters."
      )
    }
    if (shared < length(dates_morpho)) {
      message(
        "data_metrics_morpho covers more dates than data_escape, 
      the merging will only be done on the days they share."
      )
    }
    if (shared < length(data_escape)) {
      message(
        "data_escape covers more dates than data_metrics_morpho, 
      the merging will only be done on the days they share."
      )
    }
  }
  if (check_1 == "FLASEFALSEFALSE") {
    shared <- length(intersect(
      intersect(dates_point, dates_morpho),
      dates_escape
    ))
    if (shared == 0) {
      warnings(
        "The three dataset do not share a same time period (not a single day
      in common). Please provide data sharing some dates in common and comming 
      from the same rasters."
      )
    }
    if (shared < length(dates_morpho)) {
      message(
        "data_metrics_morpho covers more dates than the others datasets, 
      the merging will only be done on the days they share."
      )
    }
    if (shared < length(data_escape)) {
      message(
        "data_escape covers more dates than the others datasets, 
      the merging will only be done on the days they share."
      )
    }
    if (shared < length(dates_point)) {
      message(
        "data_metrics_point covers more dates than the others datasets, 
      the merging will only be done on the days they share."
      )
    }
  }

  ## Are the dataset spatially overlapping ?
  ### data_metrics_point polygon
  if (!is.null(data_metrics_point)) {
    # Get informations positions
    data_metrics_point <- data.table::rbindlist(data_metrics_point)
    data_metrics_point_xy <- data.frame(
      lon = data_metrics_point$x,
      lat = data_metrics_point$y
    )
    data_metrics_point_xy <- stats::na.omit(unique(data_metrics_point_xy))
    # Create polygones for each datasets :
    ## Convert df to sf (using a metric crs)
    sf_data_metrics_point_xy <- sf::st_as_sf(
      data_metrics_point_xy,
      coords = c("lon", "lat"),
      crs = crs
    )
    ## Gather the point into one object (otherwise it will create one polygone
    # per point)
    sf_data_metrics_point_xy <- sf::st_combine(sf_data_metrics_point_xy)
    ## Compute polygones
    pol_data_metrics_point <- sf::st_concave_hull(
      sf_data_metrics_point_xy,
      ratio = 0.8,
      allow_holes = FALSE
    )
  }
  ### data_morpho_point polygon
  if (!is.null(data_metrics_morpho)) {
    # Get informations positions
    data_metrics_morpho <- data.table::rbindlist(data_metrics_morpho)
    data_metrics_morpho_xy <- data.frame(
      lon = data_metrics_morpho$centroid_x,
      lat = data_metrics_morpho$centroid_y
    )
    data_metrics_morpho_xy <- stats::na.omit(unique(data_metrics_morpho_xy))
    # Create polygones for each datasets :
    ## Convert df to sf (using a metric crs)
    sf_data_metrics_morpho_xy <- sf::st_as_sf(
      data_metrics_morpho_xy,
      coords = c("lon", "lat"),
      crs = crs
    )
    ## Gather the morpho into one object (otherwise it will create one polygone
    # per morpho)
    sf_data_metrics_morpho_xy <- sf::st_combine(sf_data_metrics_morpho_xy)
    ## Compute polygones
    pol_data_metrics_morpho <- sf::st_concave_hull(
      sf_data_metrics_morpho_xy,
      ratio = 0.8,
      allow_holes = FALSE
    )
  }
  ### data_escape polygon
  if (!is.null(data_escape)) {
    # Get informations positions
    data_escape <- data.table::rbindlist(data_escape)
    data_escape_xy <- data.frame(
      lon = data_escape$x,
      lat = data_escape$y
    )
    data_escape_xy <- stats::na.omit(unique(data_escape_xy))
    # Create polygones for each datasets :
    ## Convert df to sf (using a metric crs)
    sf_data_escape_xy <- sf::st_as_sf(
      data_escape_xy,
      coords = c("lon", "lat"),
      crs = crs
    )
    ## Gather the escape into one object (otherwise it will create one polygone
    # per morpho)
    sf_data_escape_xy <- sf::st_combine(sf_data_escape_xy)
    ## Compute polygones
    pol_data_escape <- sf::st_concave_hull(
      sf_data_escape_xy,
      ratio = 0.8,
      allow_holes = FALSE
    )
  }

  ### Tests overlapping :
  if (
    !is.null(data_metrics_point) &
      !is.null(data_metrics_morpho) &
      is.null(data_escape)
  ) {
    inter_point_morpho <- sf::st_intersection(
      pol_data_metrics_point,
      pol_data_metrics_morpho
    )
    if (is.null(inter_point_morpho)) {
      # no intersection
      warnings(
        "data_metrics_point and data_metrics_morpho are not overlapping 
      spatially, merging is not possible. Please check the following points : 
        - computation of metrics with BEE.calc.metrics_point() and 
      BEE.calc.metrics_morpho() were done on the same files.
        - The crs of on of the datasets has been modified in between operations.
        - The grid resolution has been modified in between operations.
        "
      )
    }
  }
  if (
    is.null(data_metrics_point) &
      !is.null(data_metrics_morpho) &
      !is.null(data_escape)
  ) {
    inter_morpho_escape <- sf::st_intersection(
      pol_data_metrics_morpho,
      pol_data_escape
    )
    if (is.null(inter_morpho_escape)) {
      warnings(
        "data_metrics_morpho and data_escape are not overlapping 
      spatially, merging is not possible. Please check the following points : 
        - computation of metrics with BEE.calc.metrics_morpho() and 
      BEE.calc.escape() were done on the same files.
        - The crs of on of the datasets has been modified in between operations.
        - The grid resolution has been modified in between operations.
        "
      )
    }
  }
  if (
    !is.null(data_metrics_point) &
      is.null(data_metrics_morpho) &
      !is.null(data_escape)
  ) {
    inter_point_escape <- sf::st_intersection(
      pol_data_metrics_point,
      pol_data_escape
    )
    if (is.null(inter_point_escape)) {
      warnings(
        "data_metrics_point and data_escape are not overlapping 
       spatially, merging is not possible. Please check the following points : 
        - computation of metrics with BEE.calc.metrics_point() and 
      BEE.calc.escape() were done on the same files.
        - The crs of on of the datasets has been modified in between operations.
        - The grid resolution has been modified in between operations.
        "
      )
    }
  }
  if (
    !is.null(data_metrics_point) &
      !is.null(data_metrics_morpho) &
      !is.null(data_escape)
  ) {
    inter_point_morpho_escape <- sf::st_intersection(
      pol_data_metrics_point,
      pol_data_metrics_morpho,
      pol_data_escape
    )
    if (is.null(inter_point_morpho_escape)) {
      warnings(
        "at least one dataset among data_metrics_point, 
      data_metrics_morpho and data_escape is not overlapping spatially with the
      others, merging is not possible. Please check the following points : 
        - computation of metrics with BEE.calc.metrics_point(),
      BEE.calc.metrics_morpho() and BEE.calc.escape() were done on the same files.
        - There were no modifications of the crs in between operations.
        - There were no modifications of the grid resolution in between operations.
        "
      )
    }
  }

  ########################## CODE ##############################################
  ### Add a sufix to the column_names that are shared btw sevral dataset:
  if (!is.null(data_metrics_point)) {
    to_rename <- which(names(data_metrics_point) %in% c("ID", "x", "y"))
    colnames(data_metrics_point)[to_rename] <- paste0(
      colnames(data_metrics_point)[to_rename],
      "_",
      "df_point"
    )
  }
  if (!is.null(data_metrics_morpho)) {
    to_rename <- which(names(data_metrics_morpho) %in% c("ID", "x", "y"))
    colnames(data_metrics_morpho)[to_rename] <- paste0(
      colnames(data_metrics_morpho)[to_rename],
      "_",
      "df_morpho"
    )
  }
  if (!is.null(data_escape)) {
    to_rename <- which(names(data_escape) %in% c("ID", "x", "y"))
    colnames(data_escape)[to_rename] <- paste0(
      colnames(data_escape)[to_rename],
      "_",
      "df_escape"
    )
  }
  ### merge the dataset by pixel_id and date:
  ## point and morpho
  if (check_1 == "FALSEFALSETRUE") {
    data_metrics_point$date <- as.Date(data_metrics_point$date)
    data_metrics_morpho$date <- as.Date(data_metrics_morpho$date)
    data_metrics_point$pixel_id <- as.numeric(data_metrics_point$pixel_id)
    data_metrics_morpho$pixel_id <- as.numeric(data_metrics_morpho$pixel_id)
    #crop so that only shared pixel are computed:
    pixel_to_keep <- unique(data_metrics_point$pixel_id)[which(
      unique(data_metrics_point$pixel_id) %in%
        unique(data_metrics_morpho$pixel_id)
    )]
    data_metrics_point <- data_metrics_point[
      data_metrics_point$pixel_id %in% pixel_to_keep,
    ]
    data_metrics_morpho <- data_metrics_morpho[
      data_metrics_morpho$pixel_id %in% pixel_to_keep,
    ]
    #crop so that only shared dates are computed:
    date_to_keep <- unique(data_metrics_point$date)[which(
      unique(data_metrics_point$date) %in%
        unique(data_metrics_morpho$date)
    )]
    data_metrics_point <- data_metrics_point[
      data_metrics_point$date %in% date_to_keep,
    ]
    data_metrics_morpho <- data_metrics_morpho[
      data_metrics_morpho$date %in% date_to_keep,
    ]
    merged_df <- merge(
      data_metrics_point,
      data_metrics_morpho,
      by = c("pixel_id", "date"),
      all = T
    )
  }
  ## morpho and escape
  if (check_1 == "TRUEFALSEFALSE") {
    data_metrics_morpho$date <- as.Date(data_metrics_morpho$date)
    data_escape$date <- as.Date(data_escape$date)
    data_metrics_morpho$pixel_id <- as.numeric(data_metrics_morpho$pixel_id)
    data_escape$pixel_id <- as.numeric(data_escape$pixel_id)
    #crop so that only shared pixel are computed:
    pixel_to_keep <- unique(data_escape$pixel_id)[which(
      unique(data_escape$pixel_id) %in%
        unique(data_metrics_morpho$pixel_id)
    )]
    data_metrics_morpho <- data_metrics_morpho[
      data_metrics_morpho$pixel_id %in% pixel_to_keep,
    ]
    data_escape <- data_escape[data_escape$pixel_id %in% pixel_to_keep, ]
    #crop so that only shared dates are computed:
    date_to_keep <- unique(data_escape$date)[which(
      unique(data_escape$date) %in%
        unique(data_metrics_morpho$date)
    )]
    data_metrics_morpho <- data_metrics_morpho[
      data_metrics_morpho$date %in% date_to_keep,
    ]
    data_escape <- data_escape[data_escape$date %in% date_to_keep, ]
    merged_df <- merge(
      data_metrics_morpho,
      data_escape,
      by = c("pixel_id", "date"),
      all = T
    )
  }
  if (check_1 == "FALSETRUEFALSE") {
    data_metrics_point$date <- as.Date(data_metrics_point$date)
    data_escape$date <- as.Date(data_escape$date)
    data_metrics_point$pixel_id <- as.numeric(data_metrics_point$pixel_id)
    data_escape$pixel_id <- as.numeric(data_escape$pixel_id)
    #crop so that only shared pixel are computed:
    pixel_to_keep <- unique(data_escape$pixel_id)[which(
      unique(data_escape$pixel_id) %in%
        unique(data_metrics_point$pixel_id)
    )]
    data_metrics_point <- data_metrics_point[
      data_metrics_point$pixel_id %in% pixel_to_keep,
    ]
    data_escape <- data_escape[data_escape$pixel_id %in% pixel_to_keep, ]
    #crop so that only shared dates are computed:
    date_to_keep <- unique(data_escape$date)[which(
      unique(data_escape$date) %in%
        unique(data_metrics_point$date)
    )]
    data_metrics_point <- data_metrics_point[
      data_metrics_point$date %in% date_to_keep,
    ]
    data_escape <- data_escape[data_escape$date %in% date_to_keep, ]
    merged_df <- merge(
      data_metrics_point,
      data_escape,
      by = c("pixel_id", "date"),
      all = T
    )
  }
  if (check_1 == "FALSEFALSEFALSE") {
    data_metrics_point$date <- as.Date(data_metrics_point$date)
    data_metrics_morpho$date <- as.Date(data_metrics_morpho$date)
    data_escape$date <- as.Date(data_escape$date)
    data_metrics_point$pixel_id <- as.numeric(data_metrics_point$pixel_id)
    data_metrics_point$pixel_id <- as.numeric(data_metrics_point$pixel_id)
    data_escape$pixel_id <- as.numeric(data_escape$pixel_id)
    #crop so that only shared pixel are computed:
    pixel_to_keep <- Reduce(
      intersect,
      list(
        data_metrics_morpho$pixel_id,
        data_metrics_point$pixel_id,
        data_escape$pixel_id
      )
    )
    data_metrics_point <- data_metrics_point[
      data_metrics_point$pixel_id %in% pixel_to_keep,
    ]
    data_metrics_morpho <- data_metrics_morpho[
      data_metrics_morpho$pixel_id %in% pixel_to_keep,
    ]
    data_escape <- data_escape[data_escape$pixel_id %in% pixel_to_keep, ]
    #crop so that only shared dates are computed:
    date_to_keep <- Reduce(
      intersect,
      list(
        data_metrics_morpho$date,
        data_metrics_point$date,
        data_escape$date
      )
    )
    data_metrics_point <- data_metrics_point[
      data_metrics_point$date %in% date_to_keep,
    ]
    data_metrics_morpho <- data_metrics_morpho[
      data_metrics_morpho$date %in% date_to_keep,
    ]
    data_escape <- data_escape[data_escape$date %in% date_to_keep, ]

    merged_df <- merge(
      data_metrics_point,
      data_metrics_morpho,
      by = c("pixel_id", "date"),
      all = T
    )
    merged_df <- merge(
      merged_df,
      data_escape,
      by = c("pixel_id", "date"),
      all = T
    )
  }
  merged_df <- data.table::setDF(merged_df)
  merged_df <- split(merged_df, merged_df$pixel_id)
  
  return(merged_df)
}
