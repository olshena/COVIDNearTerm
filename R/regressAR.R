#' Combines multiple autoregressive models into one prediction model using a least sqares regression model
#'
#' Need to fill in description data
#'
#'
#' @param vec A vector of numeric data
#' @param x A data frame containing covariates with which to generate predictive model, if unspecified, defaults to vec
#' @param output_type A string indicating which outcome measure should be predicted. Must be one of Min, FirstQu, Median, Mean, ThirdQu, or Max
#' @param wsize Number of prior observations to use for averaging
#' @param method Type of weighting to use in individual prediction models
#' @param pdays Number of days into the future to make predictions
#' @param nsim  Number of simulations
#' @param skip Number of input values to skip
#' @param seed Seed for random number generator
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
regressAR <- function(vec,
                      x = NULL,
                      output_type = "Max",
                      wsize,
                      method = c("unweighted", "equal", "triangle"),
                      pdays,
                      nsim,
                      skip = 0,
                      seed = NULL,
                      regression_weights = NULL){

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

  if( !(output_type %in% c("min", "max", "mean", "all") ) ) {
    stop("output_type not valid, must be one of min, max, mean, or all")
  }

  n_predictors <- ncol(x)
  n_vec <- length(vec)

  # build ar object for individual variables
  build_ar_object_list <- map(x,
                              buildAR,
                              vec = vec,
                              wsize = wsize,
                              method = method)

  # fit regression model on buildAR variables
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

    model_formula <- as.formula( paste0( "vec ~ ", variable_string, " - 1" ) )

    # remove intercept
    lm_object <- lm(data = model_data, formula = model_formula)

    errors_regression <- residuals(lm_object)
    lowess_fit_regression <- lowess(vec[-1], errors_regression^2)

  }

  # make predictions on individual variables
  predict_ar_object_list <- map(build_ar_object_list,
                                predictAR,
                                pdays = pdays,
                                nsim = nsim,
                                skip = skip,
                                output_type = "predictions")


  #combine predictions into a data frame
  predicted_data <- map2_dfr(predict_ar_object_list, names(predict_ar_object_list), function(.predict_ar_object, .name){

    .predict_ar_object %>%
      melt(id.var = "t") %>%
      mutate(pred_set = as.numeric( gsub("X", "", variable) ),
             variable = .name)

  }) %>%
    dcast(pred_set + t ~ variable) %>%
    arrange(pred_set, t)


  # generate regression predicted values
  predicted_values <- predicted_data %>%
    group_by(pred_set) %>%
    mutate( y_hat_mean = predict(object = lm_object, newdata = data.frame(x = x, w = w, z = z)),
            errors = addError(y_hat_mean, lowess_fit_regression, n_draws = pdays),
            predicted_value = y_hat_mean + errors) %>%
    as.data.frame

  output <- predicted_values %>%
    select(t, pred_set, predicted_value) %>%
    dcast(pred_set ~ t, value.var = "predicted_value") %>%
    select(-pred_set)
  # summarize predicted values

  if(output_type == "max") {
    return_stat <- data.frame(max = apply(output, 1, max) )
    return_stat$rep <- 1:nsim
  }
  if(output_type == "min") {
    return_stat <- data.frame(min = apply(output, 1, max) )
    return_stat$rep <- 1:nsim
  }
  if(output_type=="mean") {
    return_stat <- data.frame(mean = apply(output, 1, mean) )
    return_stat$rep <- 1:nsim
  }
  if(output_type == "all" ) {
    return_stat <- data.frame( t( apply(output, 1, summary) ) )
    names( return_stat ) <- gsub("\\.|X", "", names( return_stat ) )
    names( return_stat ) <- gsub("1st", "First", names( return_stat ) )
    names( return_stat ) <- gsub("3rd", "Third", names( return_stat ) )
    names( return_stat ) <- tolower(names( return_stat ))
    return_stat$rep <- 1:nsim

  }




  return_object <- list(wsize = wsize,
                        method = method,
                        vec = vec,
                        x = x,
                        predicted = predicted_values,
                        model_object = lm_object,
                        return_stats = return_stat)

  class(return_object) <- "regressAR"

  return(return_object)
}
