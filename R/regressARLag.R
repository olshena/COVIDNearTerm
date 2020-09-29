#' Combines multiple autoregressive models into one prediction model using a least sqares regression model
#'
#' Need to fill in description data
#'
#'
#' @param vec A vector of numeric data, will also be include as a default predictor
#' @param x A data frame containing covariates with which to generate predictive model, if unspecified, defaults to vec
#' @param x_lag An array with an element for each non-date column in x, indicating the lag associated with that column
#' @param index_colname An character indicating the name of the column that temporally indexes data, e.g. "date". This column should be yyyy/mm/dd format
#' @param output_type A string indicating which outcome measure should be predicted. Must be one of Min, FirstQu, Median, Mean, ThirdQu, or Max
#' @param wsize Number of prior observations to use for averaging
#' @param method Type of weighting to use in individual prediction models
#' @param pdays Number of days into the future to make predictions
#' @param nsim  Number of simulations
#' @param skip Number of input values to skip
#' @param seed Seed for random number generator
#' @param regression_weights An array specifying regression weights. Currently not implemented
#' @param rhat_method Method for calculating rhat, if "none", rhat = 1 and has no effect
#' @param debug TRUE returns buildAR objects in addition to standard output

#' @return A list containing the specified output statistics for each sim
#'
#' @import dplyr
#' @import purrr
#' @import reshape2
#' @import lubridate
#'
#' @export
#'
#'
regressARLag <- function(vec,
                      x = NULL,
                      x_lag = NULL,
                      index_colname = "date",
                      output_type = "Max",
                      wsize,
                      method = c("unweighted", "equal", "triangle"),
                      pdays,
                      nsim,
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
    x_names <- x %>% select( -(!!!index_colname) ) %>% names

    if( is.null(x_lag) ) {
      x_lag <- rep(0, (ncol(x) - 1) )
    } else {
      if( length( x_lag )  != length(x_names) ) {
        stop("x_lag dimensionality doesn't match number of columns in x")
      }
    }
  }

  if( is.null(x) & !is.null(x_lag) ) {
    stop("if x is not specified, x_lag should not be specified")
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

  n_predictors <- x_names
  n_vec <- nrow(vec)

  y_data <- vec %>% select( -(!!!index_colname) )
  y_name <- names(y_data)

  predict_start_date_y <- vec %>% select( !!!index_colname ) %>% unlist %>% ymd %>%  max()
  # build ar object for vec
  build_ar_object_list_y <- map(y_data,
                                buildAR,
                                wsize = wsize,
                                method = method,
                                rhat_method = rhat_method)

  # .build_ar_object <- build_ar_object_list_y[[1]]
  # .name <- names(build_ar_object_list_y)[[1]]
  # .vec = vec
  model_data <- map2_dfr(build_ar_object_list_y, names(build_ar_object_list_y), function(.build_ar_object, .name, .vec){

    data.frame( variable = .name,
                t = ymd(.vec[[index_colname]][-1]),
                vec = .vec[-1, 2],
                fits = .build_ar_object$fits)
  }, .vec = vec) %>%
    dcast(t + vec ~ variable, value.var = "fits")

  # build data frame for fitting regression model on buildAR variables
  if( !is.null(x) ) {

    x_data <- x %>% select( -(!!!index_colname) )
    x_dates <- x %>% select( (!!!index_colname) )

    x_data_list <- map2(x_data, x_names, function(.x_data, .x_name, .dates) {

     .tor <- data.frame(date = .dates,
                        val = .x_data) %>%
       filter(!is.na(val)) %>%
       select(val)

     names(.tor)[1] <- .x_name

     .tor %>% unlist

    }, x_dates)

    x_date_list <- map2(x_data, x_names, function(.x_data, .x_name, .dates) {

      .tor <- data.frame(date = .dates,
                         val = .x_data) %>%
        filter(!is.na(val)) %>%
        select(date)


      .tor %>% unlist

    }, x_dates)

    x_model_data <- list(x_data, x_dates, x_lag)

    predict_start_date_x <- pmap_dfr(x_model_data, function(.data, .dates, .lag) {

      .dates[ !is.na(.data) ] %>% tail(1) %>% unlist %>% ymd %>%  max() + .lag

    } ) %>% as.data.frame %>% do.call(what = "c")

    build_ar_object_list_x <- map(x_data_list,
                                  buildAR,
                                  wsize = wsize,
                                  method = method,
                                  rhat_method = rhat_method)

    # .build_ar_object <- build_ar_object_list_x[[1]]
    # .name <- names(build_ar_object_list_x)[[1]]
    # .vec = x
    # .lag = x_lag[1]
    # .dates <- x_date_list[[1]]

    x_model_args <- list(build_ar_object_list_x,
                         names(build_ar_object_list_x),
                         x_lag,
                         x_date_list)


    model_data_x <- pmap_dfr(x_model_args, function(.build_ar_object, .name, .lag, .dates){
      data.frame( variable = .name,
                  t = ymd(.dates)[-1] + .lag,
                  fits = .build_ar_object$x[-1])
    }) %>%
      dcast(t  ~ variable, value.var = "fits")

    model_data <- model_data %>%
      merge(y = model_data_x,
            by = "t",
            all = TRUE)
  }


  if( is.null( regression_weights ) ) {

    if( length( x_names ) > 1 ) {
      cov_correlation <- model_data[ , x_names] %>%
        # mutate( vec = vec) %>%
        cor(use = "pairwise.complete")
    } else {
      cov_correlation = "Not enough covariates for a correlation."
    }
    # model vec ~ fits
    variable_string <- paste0( c(y_name, x_names ), collapse = " + ")

    model_formula <- as.formula( paste0( "vec ~ ", variable_string, " - 1" ) )

    # remove intercept
    lm_object <- lm(data = model_data, formula = model_formula)

    errors_regression <- residuals(lm_object)
    lowess_fit_regression <- lowess(lm_object$model$vec, errors_regression^2)

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

  # .predict_ar_object <- predict_ar_object_list[[1]]
  # .name <- names(predict_ar_object_list)[1]
  # .predict_start_date <- predict_start_dates[1]

  # names(predict_ar_object_list)

  predicted_data <- predict_ar_object_list[[y_name]]$return_stat %>%
    mutate(t = t + predict_start_date_y) %>%
    melt(id.var = "t") %>%
    mutate(pred_set = as.numeric( gsub("X", "", variable) ),
           variable = .name) %>%
    dcast(pred_set + t ~ variable) %>%
    arrange(pred_set, t)


  if( !is.null(x) ) {

    predict_ar_object_x_list <- predict_ar_object_list
    predict_ar_object_x_list[[y_name]] <- NULL

    predict_model_args_x <- list(predict_ar_object_x_list,
                                 names(predict_ar_object_x_list),
                                 predict_start_date_x,
                                 x_lag)

    predicted_data_x <- pmap_dfr(predict_model_args_x, function(.predict_ar_object, .name, .predict_start_date, .lag){

      .predict_ar_object$return_stat %>%
        mutate(t = t + .predict_start_date) %>%
        melt(id.var = c("t") ) %>%
        mutate(pred_set = as.numeric( gsub("X", "", variable) ),
               variable = .name) %>%
        rename(predicted_value = value)

    })

    expanded_observed <- map_dfr(1:nsim, function(.pred_set, .x_obs){

      .x_obs %>%
        melt(id.vars = c("t")) %>%
        mutate(pred_set = .pred_set) %>%
        select(pred_set, t, everything()) %>%
        rename(observed_value = value)

    }, model_data_x)

    x_all <- expanded_observed %>%
      merge(y = predicted_data_x,
            by = c("t", "pred_set", "variable"),
            all = TRUE) %>%
      mutate(value = ifelse( !is.na(observed_value), observed_value, predicted_value ) ) %>%
      arrange(variable, pred_set, t) %>%
      select(pred_set, t, variable, value) %>%
      dcast(pred_set + t ~ variable)

    predicted_data <- predicted_data %>%
      merge(x = x_all,
            by = c("pred_set", "t"),
            all = TRUE) %>%
      arrange(pred_set, t) %>%
      filter(!is.na( (!!rlang::sym(y_name) ) ))
  }

  #########################################################
  # check if lm was able to estimate for each variable
  #########################################################
  estimated_effect_names <- lm_object %>% summary %>% coefficients %>% rownames
  missing_estimates <- setdiff(c(y_name, x_names), estimated_effect_names)

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
