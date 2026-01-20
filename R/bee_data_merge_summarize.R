#' Merge the outputs of the metrics functions and summarize them trought time.
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
#'@param summarize_by takes the followings options, into "" :
#' - "extreme_event" for each metrics, a mean, median, variance, min and max will
#' be computed for each extrem event
#' - "day" metrics are not summarized through time, the function keeps a daily
#' resolution and just merge the datasets.
#' - "week" for each metrics, a mean, median, variance, min and max will
#' be computed weekly. Note that this imply to provide daily output
#' (BEE.calc.escape(only_days_EE = FALSE))
#' - "two_weeks" for each metrics, a mean, median, variance, min and max will
#' be computed every two weeks. Note that this imply to provide daily output
#' (BEE.calc.escape(only_days_EE = FALSE))
#' - "month" for each metrics, a mean, median, variance, min and max will
#' be computed monthly. Note that this imply to provide daily output
#' (BEE.calc.escape(only_days_EE = FALSE))
#' - "year" for each metrics, a mean, median, variance, min and max will
#' be computed yearly. Note that this imply to provide daily output
#' (BEE.calc.escape(only_days_EE = FALSE))
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
# data_metrics_morpho <- list_morpho_metrics ;
# data_escape <- dist_to_escape ; crs = "EPSG:3035"

BEE.data.merge_summarize <- function(
  data_metrics_point = NULL,
  data_metrics_morpho = NULL,
  data_escape = NULL,
  summarize_by,
  crs
) {
  ################## FORMAT TESTS AND WARNINGS ###################################
  # Check that the summarize_by is a valid option
  # ("extreme_event"/"weak"/"month"/"year") :
  if (!(summarize_by %in% c("extreme_event", "day", "weak", "month", "year"))) {
    warnings(
      "The summarize_by argument is not a suitable option, please 
    provide one from the following list : extreme_event, day, weak, month, year"
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
    dates_point <- as.Date(data_metrics_point[[1]]$Date)
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
    data_metrics_morpho <- data_metrics_morpho[[1]]
    dates_morpho <- unique(data_metrics_morpho$date)
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
    dates_escape <- unique(data_escape$date)
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
  ### data_metrics_point coordinates
  data_metrics_point_xy <- data.table::rbindlist(data_metrics_point)
  data_metrics_point_xy <- data.frame(
    lon = data_metrics_point_xy$x,
    lat = data_metrics_point_xy$y
  )
  data_metrics_point_xy <- stats::na.omit(unique(data_metrics_point_xy))

  ### data_metrics_morpho coordinates
  data_metrics_morpho_xy <- data_metrics_morpho[[1]]
  data_metrics_morpho_xy <- data.frame(
    lon = data_metrics_morpho_xy$centroid_x,
    lat = data_metrics_morpho_xy$centroid_y
  )
  data_metrics_morpho_xy <- stats::na.omit(unique(data_metrics_morpho_xy))

  ### data_escape coordinates
  data_escape_xy <- data.frame(
    lon = data_escape$from_x,
    lat = data_escape$from_y
  )
  data_escape_xy <- stats::na.omit(unique(data_escape_xy))

  ### Create polygones for each datasets :
  #### Convert df to sf (using a metric crs)
  sf_data_metrics_point_xy <- sf::st_as_sf(
    data_metrics_point_xy,
    coords = c("lon", "lat"),
    crs = crs
  )
  sf_data_metrics_morpho_xy <- sf::st_as_sf(
    data_metrics_morpho_xy,
    coords = c("lon", "lat"),
    crs = crs
  )
  sf_data_escape_xy <- sf::st_as_sf(
    data_escape_xy,
    coords = c("lon", "lat"),
    crs = crs
  )
  #### COmpute polygones
  pol_data_metrics_point <- sf::st_concave_hull(
    sf_data_metrics_point_xy,
    ratio = 0.8,
    allow_holes = FALSE
  )
  #### Test overlapping :

  ########################## CODE ################################################
}
