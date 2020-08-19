#' Combines multiple autoregressive models in a least sqares regression model
#'
#' Need to fill in description data
#'
#'
#' @param vec A vector of numeric data
#' @param x A data frame containing covariates with which to generate predictive model, if unspecified, defaults to vec
#' @param predict_variable A string indicating which outcome measure should be predicted. Must be one of Min, FirstQu, Median, Mean, ThirdQu, or Max
#' @param wsize Number of prior observations to use for averaging
#' @param method Type of weighting to use
#' @param pdays Number of days into the future to make predictions
#' @param nsim  Number of simulations
#' @param skip Number of input values to skip
#' @param seed Seed for random number generator
#' @param output_type Type of output to predict, currently max, mean, or all are valid
#' @param regression_weights An array specifying regression weights. Currently not implemented

#' @return A data frame containing the specified output statistics for each sim
#'
#' @import dplyr
#' @import purrr
#' @import reshape2
#'
#' @export
#'
#'
regressAR <- function(vec, x = NULL, predict_variable = "Max", wsize, method = c("unweighted", "equal", "triangle"), pdays, nsim, skip  = 0, seed = NULL, output_type = "all", regression_weights = NULL){

  output_type <- "all"

  if( !is.null(seed) ){
    set.seed(seed)
  }

  if( is.null(x) ) {
    x = data.frame(default = vec)
  }

  if( class(x) != "data.frame") {
    stop("if specified, x must be a data frame")
  }

  if( nrow(x) != length(vec) ) {
    stop("if specified, x must be the same length as vec")
  }

  if( !(predict_variable %in% c("Min", "FirstQu", "Median", "Mean", "ThirdQu", "Max") ) ) {
    stop("predict_var not valid, must be one of Min, FirstQu, Median, Mean, ThirdQu, or Max")
  }

  n_predictors <- ncol(x)
  n_vec <- length(vec)

  build_ar_object_list <- map(x, buildAR, vec = vec, wsize = wsize, method = method)

  model_data <- map2_dfr(build_ar_object_list, names(build_ar_object_list), function(.build_ar_object, .name, .vec){

    data.frame( variable = .name,
                t = 1:(n_vec - 1),
                vec = vec[-1],
                fits = .build_ar_object$fits)
  }, .vec = vec) %>%
    dcast(t + vec ~ variable, value.var = "fits")

  if( is.null( regression_weights ) ) {
    # model vec ~ fits
    variable_string <- paste0(names(x), collapse = " + ")

    model_formula <- as.formula( paste0( "vec ~ ", variable_string ) )

    lm_object <- lm(data = model_data, formula = model_formula)

    new_x <- predict(lm_object)

    # create Y-based model
    return_build_ar_object <- buildAR(vec = vec[-1], x = new_x, wsize = wsize, method = method)

    return_build_ar_object$x <- x
    return_build_ar_object$regression_weights <- "none"

  }

  ar_out <- predictAR(buildAR_obj = return_build_ar_object,
                      pdays = pdays,
                      nsim = nsim,
                      skip = skip,
                      output_type = "all")

  return(ar_out)
}
