#'
#' Create an Hexagonal Sticker for the BEE Package
#'

#library(hexSticker)
#library(magick)
#library(sysfonts)
#library(tidyverse)

# https://www.youtube.com/watch?v=O34vzdHOaEk

mhw_image <- magick::image_read("inst/sticker/mhw_time_2D.png")

hexSticker::sticker(
  subplot = mhw_image,
  package = "BioExtremeEvent",
  s_width = 1.7,
  s_height = 1.8,
  s_x = 0.99,
  s_y = 0.82,
  p_size = 16,
  p_y = 1.33,
  p_color = "#020002ff",
  h_fill = "#96e0f7a9",
  h_color = "black",
  h_size = 1,
  #url = "https://vicodelannoy.github.io/BioExtremeEvent/",
) |>
  print()
