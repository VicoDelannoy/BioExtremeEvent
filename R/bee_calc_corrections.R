#' Compute the impact of the selected definition of extreme events on the number
#' of days of extreme event.
#'
#' @description It computes the number of modifications between the raw
#' binerised Spatraster and the Spartraster produced using the
#' *BEE.calc.true_event()*.
#'
#' @param true_event_raster The first element of the *BEE.calc.true_event()*
#' ouput, which is a spatraster binarized that inclued the corrections by
#' *BEE.calc.true_event()*.
#' @param true_event_df_list The second element of the *BEE.calc.true_event()* ouput,
#' which is a list of data.table containing information about the value of each
#' pixel before and after definition criteria are applied to distinguish
#' isolated extreme days from extreme events.
#' @param plot Accepted values are TRUE of FALSE. Set to “FALSE” to not display
#' the graphs.
#'
#' @returns In R consol it returns a message with the number of 1 corrected to 0
#' and the number of 0 corrected to 1. If plot=TURE it shows the cumulative
#' percentage of changes through time.
#' The object return if a dataframe with one line percombinaison of date and
#' pixel_id. The dataframe has the same structure than the ones in
#' true_event_df_list.
#'
#' @examples
#' # to be added
#'
#' @export
#'
#-------------------------------------------------------------------------------
# For tests : true_event_raster <- Corrected_rasters;
# true_event_df_list <- Events_corrected
BEE.calc.corrections <- function(
  true_event_raster,
  true_event_df_list,
  plot = TRUE
) {
  df <- data.table::rbindlist(true_event_df_list)
  # total changes in the spatraster
  delta <- df$cleanned_value - df$original_value
  n_positive <- sum(delta > 0, na.rm = TRUE) # 0 -> 1
  n_negative <- sum(delta < 0, na.rm = TRUE) # 1 -> 0
  # total changes in percentages (among non NA values) in the spatraster:
  n_tot <- sum(!is.na(delta))
  n_pos_per <- n_positive * 100 / n_tot
  n_neg_per <- n_negative * 100 / n_tot
  message(paste0(
    "Number of 0 corrected to 1: ",
    n_positive,
    "\n",
    "Number of 1 corrected to 0: ",
    n_negative,
    "\n",
    "Percentage of non NA values corrected from 0 to 1: ",
    n_pos_per,
    "\n",
    "Percentage of non NA values corrected from 1 to 0: ",
    n_neg_per,
    "\n",
    n_pos_per + n_neg_per,
    " % of the non NA values have been modified."
  ))
  # plot number of modifications trough time:
  df$delta <- delta
  df$one_to_zero <- ifelse(delta < 0, 1, 0)
  df$zero_to_one <- ifelse(delta > 0, 1, 0)
  if (plot == TRUE) {
    # Number of corrections per date
    sum_one_to_zero <- tapply(df$one_to_zero, df$date, sum, na.rm = TRUE)
    sum_zero_to_one <- tapply(df$zero_to_one, df$date, sum, na.rm = TRUE)
    # Cumulated sum
    cum_one_to_zero_count <- cumsum(sum_one_to_zero)
    cum_zero_to_one_count <- cumsum(sum_zero_to_one)

    dates <- names(sum_one_to_zero)

    # Percentage cumulated and df for plot
    df_plot <- data.frame(
      date = as.Date(names(sum_one_to_zero)),
      cum_one_to_zero = 100 * cum_one_to_zero_count / n_tot,
      cum_zero_to_one = 100 * cum_zero_to_one_count / n_tot
    )

    cumulative_changes <- ggplot2::ggplot(df_plot) +
      ggplot2::geom_line(
        ggplot2::aes(
          x = date,
          y = cum_one_to_zero,
          color = "1 -> 0",
          group = 1
        ),
        size = 1
      ) +
      ggplot2::geom_line(
        ggplot2::aes(
          x = date,
          y = cum_zero_to_one,
          color = "0 -> 1",
          group = 1
        ),
        size = 1
      ) +
      ggplot2::scale_color_manual(
        name = "Type of change",
        values = c("1 -> 0" = "blue", "0 -> 1" = "red")
      ) +
      ggplot2::scale_x_date(
        date_breaks = "1 year",
        date_labels = "%Y"
      ) +
      ggplot2::labs(
        x = "Year",
        y = "Cumulative percentage change",
        title = "Cumulative percentage of changes by date among non NA values."
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
    print(cumulative_changes)
    # Sum of changes per pixel trough all the time period:
    sum_delta <- tapply(df$delta, df$pixel_id, sum, na.rm = TRUE)
    shell <- true_event_raster[[1]]
    vals <- terra::values(shell)
    pos <- as.integer(names(sum_delta))
    vals[pos] <- sum_delta
    terra::values(shell) <- vals
    terra::plot(
      shell,
      main = "Number of extreme days added or withdrawn per pixel."
    )
  }
  return(df)
}
