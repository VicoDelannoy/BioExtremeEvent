#' Download spataraster from Copernicus to your local memory.
#' 
#' @details
#' You can select your downloading conditions in 
#'  https://data.marine.copernicus.eu/products : Product : Data access : 
#'  Form and then copy and past the argument of the copernicusmarine.subset 
#'  function in https://data.marine.copernicus.eu/products > Product > 
#'  Data access : Form : Automate : Pyhton API and use them as this function 
#'  arguments.
#' 
#' @param username is your Copernicus account username, you can create an 
#'  account here : 
#' https://data.marine.copernicus.eu/register?redirect=%2Fproducts%3Fdisc%3Dnone
#' 
#' @param password is your Copernicus account password.
#' 
#' @param dataset_id is available on the corpernicus page of your product, in 
#'  "Data access" rubric, under the temporal resolution name. It is usually 
#'  start by cmems_ ...
#'
#' @param dataset_version it the version of the dataset, you can found it in 
#'  the "Product user manual" or as subtitle in the Marine Data Store, it is 
#'  usually in capital letters.
#' 
#' @param variables is the variable you are interested in, you can only call 
#' one variable at a time. Variables name are visible in Marine Data Store :
#' Product : Data access : From (bellow subset) : use the name in grey before
#' square bracket
#' 
#' @param minimum_longitude You can use the "Draw on map" Copernicus tool to 
#'  select an area of interest in Marine Data Store : Product : Data access
#'  : Form (bellow subset) and copy the values using "N" (north) as the maximum 
#'  latitude and "W" (west) as the minimum longitude if working in the northwest 
#'  hemisphere. 
#' 
#' @param maximum_longitude
#' 
#' @param minimum_latitude
#' 
#' @param maximum_latitude
#' 
#' @param start_datetime First day wanted in the downloaded spatraster. Please 
#'  make sure that the date is available on copernicus. 
#' 
#' @param end_datetime Last day wanted in the downloaded spatraster. Please 
#'  make sure that the date is available on copernicus.
#' 
#' @param description This function is a wrapper from the CopernicusMarine 
#' function from pyhton. A pyhton setup must be present in your computer, thus 
#' this function is testing if some pyhton is present and if not it offer to 
#' download it from R.
#' 
#' @param coordinates_selection_method = stric-inside" by default, extraction 
#'  method, to see more info, check the python API help to get moer detail 
#'  https://help.marine.copernicus.eu/en/articles/7949409-copernicus-marine-
#'  toolbox-introduction
#' 
#' @param disable_progress_bar allows to pot a progress bar of the downloading,
#' by default it is set on 'FALSE'.
#' 
#' @param output_directory where you want the data to be save, by default it 
#'   goes to the "Data" file next to your script if there is one.
#' 
#' @return A SpatRaster within the timeframe provided and the spatial frame 
#'  given. 
#'
#' @export
#' 
#-------------------------------------------------------------------------------


BEE.data.load_copernicus <- function(username,
                                     password,
                                     dataset_id,
                                     dataset_version,
                                     variables,
                                     minimum_longitude,
                                     maximum_longitude,
                                     minimum_latitude,
                                     maximum_latitude,
                                     start_datetime,
                                     end_datetime,
                                coordinates_selection_method = "strict-inside",
                                     disable_progress_bar = FALSE,
                                     output_directory = here::here("Data")) {
  # Create a virtual Python space and connect to Copernicus API
  reticulate::virtualenv_create(envname = "CopernicusMarine")
  reticulate::virtualenv_install("CopernicusMarine", 
                                 packages = c("copernicusmarine"))
  reticulate::use_virtualenv("CopernicusMarine", required = TRUE)
  reticulate::py_install("copernicusmarine", pip = TRUE)
  # Check that it work, usually, if it doesn't works it is because Pyhton is 
  # nowhere in the system.
  version <- system("python3 --version", intern = TRUE)
  
  if (length(version) == 0) {
    version <- system("python3 --version", intern = TRUE)
    user_response <- tolower(
      readline(prompt = "No Python version was found on your device. Python is 
               required to use the CopernicusMarine interface, which allows 
               downloading Marine Copernicus datasets from R. Would you like to
               install Python in a virtual environment within the 'reticulate' 
               package? This process will take approximately 5 minutes and needs
               to be done only once. (Yes/No):")
    )
    if (user_response == "YES" || user_response == "yes") {
      reticulate::install_python(version = "3.10")
      version <- reticulate::py_config()
      print(version)
      # Create a virtual Python space and connect to Copernicus API
      reticulate::virtualenv_create(envname = "CopernicusMarine")
      reticulate::virtualenv_install("CopernicusMarine",
                                     packages = c("copernicusmarine"))
      reticulate::use_virtualenv("CopernicusMarine", required = TRUE)
      reticulate::py_install("copernicusmarine", pip = TRUE)
    }
    if (user_response == "NO" || user_response == "No") {
      return(
        "At last update, Copernicus Marine data center has not provided an API
        that allows downloding directly from R, thus creating a virutal pyhton
        environment is necessary. If you want to download a small area or a 
        short time window, you may consider downloding rasters directly on 
        https://data.marine.copernicus.eu/ . Both methods work well with the 
        other functions."
      )
    }
  }
  if (length(version) > 0) {
    cm <- reticulate::import("copernicusmarine") #not called in R
    cm$login(username, password)
    
    # Get the data
    ds <- cm$subset(
      dataset_id = dataset_id,
      dataset_version = dataset_version,
      variables = list(gsub(
        pattern = "[\\[\\]]", replacement = "", variables
      )),
      minimum_longitude = minimum_longitude,
      maximum_longitude = maximum_longitude,
      minimum_latitude = minimum_latitude,
      maximum_latitude = maximum_latitude,
      start_datetime = start_datetime,
      end_datetime = end_datetime,
      coordinates_selection_method = coordinates_selection_method,
      disable_progress_bar = disable_progress_bar,
      output_directory = output_directory
    )
  }
}
