#' Combines multiple autoregressive models into one prediction model using a least sqares regression model
#'
#' Need to fill in description data
#'
#'
#' @param vec A vector of numeric data, will also be include as a default predictor
#' @param x A data frame containing covariates with which to generate predictive model, if unspecified, defaults to vec
#' @param output_type A string indicating which outcome measure should be predicted. Must be one of Min, FirstQu, Median, Mean, ThirdQu, or Max
#' @param wsize Number of prior observations to use for averaging, default is 14
#' @param method Type of weighting to use in individual prediction models, default is equal
#' @param pdays Number of days into the future to make predictions, default is 28
#' @param nsim  Number of simulations, default is 100
#' @param skip Number of input values to skip, default is 0
#' @param seed Seed for random number generator
#' @param regression_weights An array specifying regression weights. Currently not implemented
#' @param rhat_method Method for calculating rhat, if "none", rhat = 1 and has no effect
#' @param debug TRUE returns buildAR objects in addition to standard output

#' @return A list containing the specified output statistics for each sim
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
                      output_type = "max",
                      wsize = 14,
                      method = c("equal", "unweighted", "triangle"),
                      pdays =28,
                      nsim = 100,
                      skip = 0,
                      seed = NULL,
                      regression_weights = NULL,
                      rhat_method = c("none", "geometric", "arithmetic"),
                      debug = FALSE){

  if( !is.null(seed) ){
    set.seed(seed)
  }

  if( class(x) != "data.frame" & !is.null(x) ) {
    stop("if specified, x must be a data frame")
  }

  if( is.null(x) ) {
    x_names <- NULL
  } else {
    x_names <- names(x)
  }

  # if( nrow(x) != length(vec) ) {
  #   stop("if specified, x must be the same length as vec")
  # }

  if( !(output_type %in% c("min", "max", "mean", "all") ) ) {
    stop("output_type not valid, must be one of min, max, mean, or all")
  }

  if(length(rhat_method) > 1) {
    rhat_method <- rhat_method[1]
  }

  n_predictors <- length(x_names)
  n_vec <- length(vec)

  # build ar object for vec
  build_ar_object_list_y <- map(list(y_hat = vec),
                                buildAR,
                                vec = vec,
                                wsize = wsize,
                                method = method,
                                rhat_method = rhat_method)

  model_data <- map2_dfr(build_ar_object_list_y, names(build_ar_object_list_y), function(.build_ar_object, .name, .vec){

    data.frame( variable = .name,
                t = 1:(n_vec - 1),
                vec = vec[-1],
                fits = .build_ar_object$fits)
  }, .vec = vec) %>%
    dcast(t + vec ~ variable, value.var = "fits")

  # build data frame for fitting regression model on buildAR variables
  if( !is.null(x) ) {
    build_ar_object_list_x <- map(x,
                                  buildAR,
                                  # vec = vec,
                                  wsize = wsize,
                                  method = method,
                                  rhat_method = rhat_method)

    model_data_x <- map2_dfr(build_ar_object_list_x, names(build_ar_object_list_x), function(.build_ar_object, .name, .vec){

      data.frame( variable = .name,
                  t = 1:(n_vec - 1),
                  fits = .build_ar_object$x[-1])
    }, .vec = vec) %>%
      dcast(t  ~ variable, value.var = "fits")

    model_data <- model_data %>%
      merge(y = model_data_x,
            by = "t")
  }


  if( is.null( regression_weights ) ) {

    if( length( names(x) ) > 1 ) {
    cov_correlation <- model_data[ , names(x)] %>%
      # mutate( vec = vec) %>%
      cor
    } else {
      cov_correlation = "Not enough covariates for a correlation."
    }
    # model vec ~ fits
    variable_string <- paste0( c("y_hat", names(x) ), collapse = " + ")

    model_formula <- as.formula( paste0( "vec ~ ", variable_string, " - 1" ) )

    # remove intercept
    lm_object <- lm(data = model_data, formula = model_formula)

    errors_regression <- residuals(lm_object)
    lowess_fit_regression <- lowess(vec[-1], errors_regression^2)

  }


  if( !is.null(x) ) {
    build_ar_object_list <- c(build_ar_object_list_y, build_ar_object_list_x)
  } else {
    build_ar_object_list <- build_ar_object_list_y
  }

  # make predictions on individual variables
  predict_ar_object_list <- map(build_ar_object_list,
                                predictAR,
                                pdays = pdays,
                                nsim = nsim,
                                skip = skip,
                                output_type = "predictions",
                                debug = debug)

  # combined debug info
  if( debug ) {

    debug_predicted_data <- map2_dfr(predict_ar_object_list, names(predict_ar_object_list), function(.predict_ar_object, .name){

      .predict_ar_object$predict_debug_object %>%
        mutate(variable = .name)

    }) %>%
      select(variable, everything() ) %>%
      arrange(variable, sim, day, phi_index)

  }

  #combine predictions into a data frame
  predicted_data <- map2_dfr(predict_ar_object_list, names(predict_ar_object_list), function(.predict_ar_object, .name){

    .predict_ar_object$return_stat %>%
      melt(id.var = "t") %>%
      mutate(pred_set = as.numeric( gsub("X", "", variable) ),
             variable = .name)

  }) %>%
    dcast(pred_set + t ~ variable) %>%
    arrange(pred_set, t)

  estimated_effect_names <- lm_object %>% summary %>% coefficients %>% rownames
  missing_estimates <- setdiff(c("y_hat", names(x) ), estimated_effect_names)

  if( length( missing_estimates ) > 0 ) {

    cov_correlation
    warning(paste0("The following covariates were not estimable in the regression model: ", paste0(missing_estimates, collapse = ", "), "\n",
               "See debug data for more info.\n") )

  }
  # generate regression predicted values
  predicted_data$z_hat_mean = predict(object = lm_object, predicted_data)

  names(predicted_data)[which( names( predicted_data ) %in% x_names ) ] <- paste0(x_names, "_hat")

  predicted_values <- predicted_data %>%
    group_by(pred_set) %>%
    mutate( error = addError(z_hat_mean, lowess_fit_regression, n_draws = pdays),
            predicted_value = z_hat_mean + error) %>%
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
    return_stat <- data.frame(min = apply(output, 1, min) )
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
                        rhat_method = rhat_method,
                        predicted = predicted_values,
                        model_object = lm_object,
                        return_stats = return_stat,
                        buildAR_objects = NULL)

  if( debug ) {
    return_object[["debug_buildAR"]] <- build_ar_object_list
    return_object[["debug_lowess_object"]] <- lowess_fit_regression
    return_object[["debug_predictAR"]] <- debug_predicted_data
    return_object[["debug_regression_model"]] <- lm_object
    return_object[["debug_correlations"]] <- cov_correlation
  }

  class(return_object) <- "regressAR"

  return(return_object)
}
