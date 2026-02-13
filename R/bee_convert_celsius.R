#'Convert temperature data to Celsius
#'
#' @description
#'  Function that detects the unit of measurement of a pixel value in a 
#'  SpatRaster (Celsius, Kelvin or Fahrenheit) and converts it to the
#'  international unit 'Celsius' if necessary.
#' 
#' @param yourspatraster :
#'  It is the SpatRaster containing temperature data, for which you want to
#'  check the units and convert the values to Celsius.
#' 
#' @return 
#'  yourspatraster with values corrected to celsius and 'unit' metadata updated
#'  to 'Celsisus'.
#' 
#' @examples
#'
#'### Load the example dataset in R environement :
#' file_name <- system.file(file.path("extdata", "copernicus_example_data.nc"),
#'                                   package = "BioExtremeEvent")
#' copernicus_data <- terra::rast(file_name)
#'
#'### Values before applying functions :
#' head(na.omit(terra::values(x = copernicus_data[[1]],
#'                    dataframe=TRUE))) # first pixels are on land so their value is NA
#'
#'### Apply function :
#' copernicus_data_celsius <- BioExtremeEvent::BEE.convert.celsius(copernicus_data)
#'
#'### Values after applying functions :
#' head(na.omit(terra::values(x = copernicus_data_celsius[[1]],
#'                   dataframe=TRUE))) # first pixels are on land so their value is NA
#'
#'### Checking the unit of the spatraster :
#' # Before using BEE.convert.celsius :
#' terra::units(copernicus_data[[1]])
#' # After using BEE.convert.celsius :
#' terra::units(copernicus_data_celsius[[1]])
#'
#'@export
#'
#-------------------------------------------------------------------------------

BEE.convert.celsius <- function(yourspatraster) {
  # Get the unit of each layer.
  units_count <- terra::units(yourspatraster)

  # Check that no layer contains several unit.
  if (any(sapply(units_count, length) != 1)) {
    if (any(sapply(units_count, length) > 1)) {
      warning(
        "Multiple units were detected for some layers.
              Plot 'yourspatraster[[1]]' to verify if a unit is provided, and
              ensure that there are no multiple values assigned to a single
              pixel in your raster. You can use terra::units to do so."
      )
      return(NULL)
    }
    # Script goes here if there is one layer or more with no unit provided :
    choice <- utils::menu(
      c(
        "Fahrenheit is the current unit of measurement for
                     all provided layers.",
        "Kelvin is the current unit of all provided layers.",
        "Stop the function without modification."
      ),
      title = "At least one of the layers of yourspatraster has no defined unit.
    Please choose one of the following options :"
    )
    if (choice == 1) {
      unit <- "fahrenheit"
    }
    if (choice == 2) {
      unit <- "kelvin"
    }
    if (choice == 3) {
      return(NULL)
    }
  }

  # Check that all layers have the same unit.
  unit <- terra::unique(units_count)
  if (length(unit) > 1) {
    warning(
      "Units differ between layers. If you have merged layers from
            different datasets, ensure that all layers measure the same
            parameter and that the units are consistent across datasets
            before merging."
    )
    return(NULL)
  }

  # Check that the unit is not already Celsius
  if (
    unit %in%
      c(
        "Celsius",
        "celsius",
        "<c2><b0>C",
        "<c2><b0>c",
        "c",
        "C",
        "degrees_C",
        "degrees_c"
      )
  ) {
    #buble symboles have been replaced
    #by ASCII code to avoid bugs.
    print(
      "Your dataset is already in celsius, there is no need to use this function."
    )
    return(yourspatraster)
  }
  # Store original metadata for later
  original_metadata <- list(
    time = terra::time(yourspatraster),
    crs = terra::crs(yourspatraster),
    extent = terra::ext(yourspatraster)
  )

  # Convert from kelvin to celsius
  if (tolower(unit) %in% c("kelvin", "k", "degrees_K", "degrees_k")) {
    # conversion en Celsius :
    yourspatraster <- terra::app(yourspatraster, function(x) {
      x - 273.15
    })
    terra::units(yourspatraster) <- "Celsius" # unfortunately this deletes other
    # metadata
    # restore other metadata
    terra::time(yourspatraster) <- original_metadata$time
    terra::crs(yourspatraster) <- original_metadata$crs
    terra::ext(yourspatraster) <- original_metadata$extent
    warning(
      "Your data were in Kelvin and have been converted to Celsius using:
            former value - 273.15 = new value."
    )
    return(yourspatraster)
  }

  # Convert from Fahrenheit to celsius
  if (tolower(unit) %in% c("fahrenheit", "f", "degrees_F", "degrees_f")) {
    yourspatraster <- terra::app(yourspatraster, function(x) {
      terra::round((x - 32) * (5 / 9), digits = 3)
    }) # conversion en celsius
    terra::units(yourspatraster) <- "Celsius"
    # unfortunately this deletes other metadata
    # restore other metadata
    terra::time(yourspatraster) <- original_metadata$time
    terra::crs(yourspatraster) <- original_metadata$crs
    terra::ext(yourspatraster) <- original_metadata$extent
    warning(
      "Your data were in Fahrenheit and have been converted to Celsius
            using: round((former value - 32)*(5/9), digits = 3) = new value."
    )
    return(yourspatraster)
  }

  warning("The unit of your dataset is not recognized.")
  return(NULL)
}
