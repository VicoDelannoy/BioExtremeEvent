#' Plot event categories
#'
#' @description
#'  Use *BEE.calc.metric_point()* to plot extreme event categorie over the all
#' time period or over a subset of that period.
#'
#' @param metric_point_df the output of *BEE.calc.metric_point()* for a given
#' location. Please note that the output of BEE.calc.metric_point() is a list.
#' You only need to provide a dataframe from that list, not the whole list.
#' @param color_theme A vector of four colour codes of your choice (the first
#' one is for category 1 and the last one is for category 4), or 'red' to use an
#' automatic palette of red shades, or 'blue' for a colour-shaded palette.
#' @param start_date = NULL by default or the date at which you want to start
#' the plot. It mus be in the same format than metric_point_df$date.
#' @param end_date = NULL by default or the date at which you want to stop
#' the plot. It mus be in the same format than metric_point_df$date.
#' @param ... Customising the graph is possible by adding any general ggplot2
#' argument.
#'
#' @return A 'ggplot2' with the date on the x-axis, the studied parameter value
#' on the y-axis, and the extreme events coloured according to their category.
#'
#'
#' @examples
#' # Load data:
#' file_name_1 <- system.file(file.path("extdata", "metrics_points_day.rds"),
#'                                   package = "BioExtremeEvent")
#' metrics_points_day <- readRDS(file_name_1)
#' # Get data set for one GPS point
#' metrics_points_day <- metrics_points_day[[1]] # first GPS point
#' BEE.plot.categories(metric_point_df = metrics_points_day,
#'                     color_theme = "red",
#'                     start_date = NULL,
#'                     end_date = NULL)
#'
#' @export
#'
#---------------------------------------------------------------------------

BEE.plot.categories <- function(
  metric_point_df,
  color_theme = "red",
  start_date = NULL,
  end_date = NULL,
  ...
) {
  if (!is.data.frame(metric_point_df)) {
    warning(
      "metric_point_df must be a dataframe from the list built by 
    *BEE.calc.metrics_poin()*, not the whole list. Please provide one 
    dataframe."
    )
  }
  # subset_date
  if (!is.null(start_date)) {
    metric_point_df <- metric_point_df[
      which(
        metric_point_df$date >= start_date
      ),
    ]
  }
  if (!is.null(end_date)) {
    metric_point_df <- metric_point_df[
      which(metric_point_df$date <= end_date),
    ]
  }
  # calulate daily categories limit
  one_place_col <- metric_point_df |>
    dplyr::mutate(
      lim1 = ifelse(
        value > baseline_qt,
        pmin(baseline_qt + 1 * anomaly_unit, value),
        NA
      ),
      lim2 = ifelse(
        value > baseline_qt,
        pmin(baseline_qt + 2 * anomaly_unit, value),
        NA
      ),
      lim3 = ifelse(
        value > baseline_qt,
        pmin(baseline_qt + 3 * anomaly_unit, value),
        NA
      ),
      lim4 = ifelse(
        value > baseline_qt,
        pmin(baseline_qt + 4 * anomaly_unit, value),
        NA
      )
    )

  if (1 < length(color_theme) & 4 != length(color_theme)) {
    warning(
      "*color_theme* argument should be of length 1 ('red' or 'blue) or a
  vector of length 4, with one color per category."
    )
  }
  if (color_theme == "blue") {
    cat1 <- "#00e6c09f"
    cat2 <- "#00b7ffff"
    cat3 <- "#0047ccff"
    cat4 <- "#09136bff"
  }
  if (color_theme == "red") {
    cat1 <- "#ffda0aff"
    cat2 <- "#ff800aff"
    cat3 <- "#da0000ff"
    cat4 <- "#400364ff"
  }
  if (length(color_theme) == 4) {
    cat1 <- color_theme[1]
    cat2 <- color_theme[2]
    cat3 <- color_theme[3]
    cat4 <- color_theme[4]
  }

  # Graph
  plot <- ggplot2::ggplot(one_place_col, ggplot2::aes(x = date)) +

    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = baseline_qt, ymax = lim1, fill = "Category 1")
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lim1, ymax = lim2, fill = "Category 2")
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lim2, ymax = lim3, fill = "Category 3")
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lim3, ymax = lim4, fill = "Category 4")
    ) +

    ggplot2::geom_line(
      ggplot2::aes(y = value, color = "Parameter value")
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = baseline_qt, color = "Baseline (quantile)"),
      linetype = "dashed",
      linewidth = 2
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = baseline_mean, color = "Baseline (mean)"),
      linetype = "dashed",
      linewidth = 2
    ) +
    ggplot2::scale_fill_manual(
      name = "Event categories",
      values = c(
        "Category 1" = cat1,
        "Category 2" = cat2,
        "Category 3" = cat3,
        "Category 4" = cat4
      )
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Parameter value" = "#41383aff",
        "Baseline (quantile)" = "darkred",
        "Baseline (mean)" = "darkgreen"
      )
    ) +
    ggplot2::scale_x_date(date_labels = "%m-%Y") +
    ggplot2::labs(x = "Date", y = "Parameter value") +
    ggplot2::theme_minimal()

  extras <- list(...)
  if (length(extras) > 0) {
    plot <- Reduce(`+`, c(list(plot), extras))
  }
  print(plot)
}
