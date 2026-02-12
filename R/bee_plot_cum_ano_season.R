#' Cumulative anomaly against months of the years
#'
#' @description allows to see at which time of the year a place is explosed to
#' the bigest anomalies.
#' #' @param metric_point_df the output of *BEE.calc.metric_point()* for a given
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
#' @return A 'ggplot2' with the months on the x-axis, the cumulative anomaly per
#' extreme event of the studied parameter value on the y-axis.
#'
#'
#' @examples
#' # Load data:
#' file_name_1 <- system.file(file.path("extdata", "metrics_points_day.rds"),
#'                                   package = "BioExtremeEvent")
#' metrics_points_day <- readRDS(file_name_1)
#' # Get data set for one GPS point
#' metrics_points_day <- metrics_points_day[[1]] # first GPS point
#' BEE.plot.cumulative_anomaly(metric_point_df = metrics_points_day,
#'                     start_date = NULL,
#'                     end_date = NULL)
#'
#' @export
#-----------------------------------------------------------------------------
BEE.plot.cumulative_anomaly <- function(
  metric_point_df,
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

  # Month in english
  Sys.setlocale("LC_TIME", "en_US.UTF-8")
  # withdraw events that are not extreme:
  one_place_filtered <- one_place %>%
    dplyr::filter(daily_category != "No extreme event") %>%
    dplyr::mutate(
      # date fictive pour aligner tous les mois sur une même année
      month_day = as.Date(paste0("2000-", format(date, "%m-%d"))),
      year_num = lubridate::year(date) # année numérique pour le gradient
    )

  # color scales labels:
  years <- sort(unique(one_place_filtered$year_num))
  if (length(years) > 10) {
    breaks_years <- round(seq(
      from = min(years),
      to = max(years),
      length.out = 10
    ))
  } else {
    breaks_years <- years
  }

  # Graph
  plot <- ggplot2::ggplot(
    one_place_filtered,
    ggplot2::aes(
      x = month_day,
      y = cumulative_anomaly_qt,
      group = ID,
      color = year_num,
      fill = year_num
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = 0, ymax = cumulative_anomaly_qt),
      alpha = 1 / length(unique(one_place_filtered$year_num)) + 0.1,
      color = NA
    ) +
    ggplot2::geom_line(
      size = 1,
      alpha = 1 / length(unique(one_place_filtered$year_num)) + 0.1
    ) +
    ggplot2::scale_color_gradientn(
      colors = c("darkblue", "turquoise", "green", "yellow", "orange", "red"),
      name = "Years",
      breaks = breaks_years,
      labels = breaks_years
    ) +
    ggplot2::scale_fill_gradientn(
      colors = c("darkblue", "turquoise", "green", "yellow", "orange", "red"),
      name = "Years",
      breaks = breaks_years,
      labels = breaks_years
    ) +
    ggplot2::scale_x_date(date_labels = "%b", date_breaks = "1 month") +
    ggplot2::labs(
      x = "Months",
      y = "Cumulative anomaly (quantile)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  extras <- list(...)
  if (length(extras) > 0) {
    plot <- Reduce(`+`, c(list(plot), extras))
  }
  print(plot)
  Sys.sleep(1)

  # Second plot to compare cumulative anomalies profiles:
  ## prepare data
  one_place_filtered <- one_place %>%
    dplyr::filter(daily_category != "No extreme event") %>%
    dplyr::arrange(ID, date) %>%
    dplyr::group_by(ID) %>%
    dplyr::mutate(event_day = row_number() - 1)

  # Grap
  plot2 <- ggplot2::ggplot(
    one_place_filtered,
    ggplot2::aes(
      x = event_day,
      y = cumulative_anomaly_qt,
      group = ID,
      color = factor(ID)
    )
  ) +
    ggplot2::geom_line(size = 1) +
    ggplot2::labs(
      x = "Duration of event (days)",
      y = "Cumulative anomaly (quantile)",
      color = "ID"
    ) +
    ggplot2::theme_minimal()

  extras <- list(...)
  if (length(extras) > 0) {
    plot2 <- Reduce(`+`, c(list(plot2), extras))
  }
  print(plot2)
}
