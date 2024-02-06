#' Show the code for an HBAM or FBAM model
#'
#' Show the Stan code for one of the models in the package.
#'
#' @export
#' @param model Character: Name of the model to show. Defaults to "HBAM". See the documentation for the `hbam()` function for a list of the available models.
#' @return The function prints Stan code.
#'
#' @examples
#' show_code("HBAM")

show_code <- function(model = "HBAM") {
  if(!model %in% names(stanmodels)) {
    stop(paste(model, "is not a valid model name."))
  }
  show(stanmodels[[model]])
}
