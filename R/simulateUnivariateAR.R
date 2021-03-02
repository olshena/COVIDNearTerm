#' Builds an autoregressive model and uses it to predict
#'
#' A wrapper function for buildAR and predictAR.
#'
#'
#' @param vec A vector of numeric data
#' @param wsize Number of prior observations to use for averaging
#' @param method Type of weighting to use
#' @param pdays Number of days into the future to make predictions
#' @param nsim  Number of simulations
#' @param skip Number of input values to skip
#' @param seed Seed for random number generator
#' @param output_type Type of output
#' @param lambda Shrinkage parameter, if not specified, default is 0 which produces no shrinkage. If an array is specified, all values in the array
#' are evaluated and the optimal lambda is chosen based on residual sum of squares. Values should be between 0 and 1 inclusive.
#' @param alpha Alpha parameter, if not specified, default is 0 which produces no shrinkage. If an array is specified, all values in the array
#' are evaluated and the optimal alpha is chosen based on residual sum of squares. Values should be >=0.
#' @param rolling_mean rolling mean window for x. Default is 1
#' @param debug TRUE returns buildAR objects in addition to standard output

#' @return A data frame containing the specified output statistics for each sim
#'
#' @import dplyr
#' @import purrr
#' @import reshape2
#' @import lubridate
#' @import zoo
#'
#' @export
#'
#'
simulateUnivariateAR <- function(vec,
                                 x = NULL,
                                 wsize,
                                 method = c("unweighted", "equal", "triangle"),
                                 pdays,
                                 nsim,
                                 skip  = 0,
                                 seed = NULL,
                                 output_type = "all",
                                 rhat_method = c("none", "geometric", "arithmetic"),
                                 lambda = 0,
                                 alpha = 0,
                                 rolling_mean = 1,
                                 debug = FALSE){

  if(length(rhat_method) > 1) {
    rhat_method <- rhat_method[1]
  }

  if(is.null(x)) {
    stop("x cannot be NULL")

  }

  if(rolling_mean < 1 | rolling_mean > length(x) ) {
    stop(paste0("rolling mean must be at least 1 and no more than the length of x. specified rolling mean is ", rolling_mean) )
  }

  if( round(rolling_mean, 0) != rolling_mean ) {
    stop(paste0("specified rolling mean, ", rolling_mean, ", is not an integer and should be") )
  }

  if( !is.numeric(lambda) ) {
    stop("lambda must be numeric")
  }

  if( any(lambda > 1 ) |  any(lambda < 0 ) ) {
    stop("lambda must be between 0 and 1 inclusive")
  }

  if( !is.numeric(alpha) ) {
    if( alpha != "c") {
      stop("if alpha is not 'c' it must be numeric")
    }
  } else {
    if( any(alpha < 0 ) ) {
      stop("alpha must be >= 0")
    }
  }

  #########################################################
  # transform x using rolling mean
  #########################################################
  ### generate rolling means for first k values

  if(rolling_mean != 1){
    first_means <- map(1:(rolling_mean - 1), function(.k, .x) {

      zoo::rollmean(x = .x[1:.k], k = .k, align = "right" )

    }, x) %>% unlist

    last_means <- zoo::rollmean(x = x, k = rolling_mean, align = "right" )

    x <- c(first_means, last_means)
  }
  # x_rolling <- c(first_means, last_means)

  #########################################################
  # build training models
  #########################################################
  build_ar_object_y <- buildAR(vec = vec, x = vec, wsize, method = method, seed = seed, rhat_method = rhat_method, lambda = lambda)
  build_ar_object_x <- buildAR(vec = x, x = x, wsize, method = method, seed = seed, rhat_method = rhat_method, lambda = lambda)
  # build_ar_object_x_rolling <- buildAR(vec = x_rolling, x = x_rolling, wsize, method = method, seed = seed, rhat_method = rhat_method, lambda = lambda)

  ###Fit lowess for relationship between data and error so error is proportional to magnitude
  errors <- build_ar_object_y$errors
  dat <- build_ar_object_y$vec[-1]
  ###Skip the first skip features in the vector
  if( skip > 0 ) {
    dat <- dat[-(1:skip)]
    errors <- errors[-(1:skip)]
    diff_phis <- diff_phis[-(1:(skip-1))]
  }

  ###Fit lowess for relationship between data and error so error is proportional to magnitude
  lowess_fit <- lowess(dat, errors^2)

  #########################################################
  # determine alpha to use
  #########################################################

  y_obs <-  build_ar_object_y$vec
  y_phis <- build_ar_object_y$phi_s
  x_phis <- build_ar_object_x$phi_s
  # x__rolling_phis <- build_ar_object_x_rolling$phi_s

  x_phis_standardized <- ( x_phis - mean( x_phis ) ) * ( sd( y_phis ) / sd( x_phis ) ) + mean( y_phis )

  grid_search <- FALSE

  if( length(alpha) > 1) {

    alpha_grid <- alpha
    grid_search <- TRUE

  } else {
    if(  alpha == "c" ) {

      # C <- 1 / mean( x_phis_standardized - y_phis )
      C <- 1/50
      alpha_grid = seq(0, C, length.out = 50)
      grid_search <- TRUE

    }
  }

  if( grid_search ){

    grid_alphas <- map_dfr(alpha_grid, function(.alpha, .y_phis, .x_phis_standardized, .y_obs) {

      .alpha_pred <- .y_phis * head(.y_obs, -1) +
        .alpha * ( .x_phis_standardized - y_phis ) * head(.y_obs, -1)

      .ssr <- sum( (.alpha_pred - .y_obs[-1] ) ^ 2 )

      data.frame(alpha = .alpha,
                 SSR_alpha = .ssr)

    }, y_phis, x_phis_standardized, y_obs) %>%
      arrange(SSR_alpha)

    alpha <- grid_alphas$alpha[1]
  }

  #########################################################
  # get predicted x and y values and corresponding phis
  #########################################################

  ar_out_y <- predictAR(buildAR_obj = build_ar_object_y,
                        pdays = pdays,
                        nsim = nsim,
                        skip = skip,
                        output_type = "predictions",
                        debug = TRUE)

  y_debug <- ar_out_y$predict_debug_object %>%
    filter(!is.na(predicted)) %>%
    select(sim, day, phi_s, predicted, current_vec)

  predict_y <- split(y_debug, y_debug$sim) %>%
    map2_dfr(unique(y_debug$sim),  function(.debug_sim, .sim) {

      data.frame(sim = .sim,
                 day = .debug_sim$day,
                 y_t1 = c(.debug_sim$current_vec[1], head(.debug_sim$predicted, -1)),
                 phi_y = .debug_sim$phi_s)


    })

  ar_out_x <- predictAR(buildAR_obj = build_ar_object_x,
                        pdays = pdays,
                        nsim = nsim,
                        skip = skip,
                        output_type = "predictions",
                        debug = TRUE)

  # ar_out_x_rolling <- predictAR(buildAR_obj = build_ar_object_x_rolling,
  #                       pdays = pdays,
  #                       nsim = nsim,
  #                       skip = skip,
  #                       output_type = "predictions",
  #                       debug = TRUE)

  x_debug <- ar_out_x$predict_debug_object %>%
    filter(!is.na(predicted)) %>%
    select(sim, day, phi_s, predicted, current_vec)

  # x_rolling_debug <- ar_out_x$predict_debug_object %>%
  #   filter(!is.na(predicted)) %>%
  #   select(sim, day, phi_s, predicted, current_vec)

  predict_x <- split(x_debug, x_debug$sim) %>%
    map2_dfr(unique(x_debug$sim), function(.debug_sim, .sim) {

      data.frame(sim = .sim,
                 day = .debug_sim$day,
                 phi_p = .debug_sim$phi_s)


    })

  # predict_x_rolling <- split(x_debug, x_debug$sim) %>%
  #   map2_dfr(unique(x_debug$sim), function(.debug_sim, .sim) {
  #
  #     data.frame(sim = .sim,
  #                day = .debug_sim$day,
  #                phi_p = .debug_sim$phi_s)
  #
  #
  #   })

  predict_data_sub <- predict_y %>%
    merge(y = predict_x,
          by = c("sim", "day") ) %>%
    arrange(sim, day)


  predict_data <- split(predict_data_sub, predict_data_sub$sim) %>%
    map2_dfr(unique(predict_data_sub$sim), function(.predict_sim, .sim, .alpha) {

      .predict_sim  %>%
        mutate(sim = .sim,
               phi_p_standardized = ( phi_p - mean( phi_p ) ) * ( sd( phi_y ) / sd( phi_p ) ) + mean( phi_y ),
               f_phi = phi_p_standardized - phi_y,
               value = phi_y * y_t1 + .alpha * f_phi * y_t1)
     }, alpha)


  predict_data <- pmap_dfr(predict_data, function(...) {
    current <- tibble(...)

    current$error <- addError(xpred = current$value, lowess_obj = lowess_fit)
    current
  }) %>%
    mutate(predicted = value + error)


  output <- predict_data %>%
    select(sim, day, predicted) %>%
    reshape2::dcast(sim ~ day, value.var = "predicted") %>%
    select(-sim)

  if(output_type == "max") {
    return_stat <- data.frame(max = apply(output, 1, max) )
    # return_stat$rep <- 1:nsim
  }
  if(output_type=="mean") {
    return_stat <- data.frame(mean = apply(output, 1, mean) )
    # return_stat$rep <- 1:nsim
  }
  if(output_type == "all" ) {
    return_stat <- data.frame( t( apply(output, 1, summary) ) )
    names( return_stat ) <- gsub("\\.|X", "", names( return_stat ) )
    names( return_stat ) <- gsub("1st", "First", names( return_stat ) )
    names( return_stat ) <- gsub("3rd", "Third", names( return_stat ) )
    # return_stat$rep <- 1:nsim

  }
  if(output_type == "predictions" ) {
    return_stat <- data.frame( t( output ) )
    return_stat$t <- 1:pdays
    return_stat <- return_stat %>% select(t, everything())
  }

  if( debug ) {

    return_object <- list("return_stat" = return_stat)

    return_object[["build_object_y"]] <- build_ar_object_y
    return_object[["build_object_x"]] <- build_ar_object_x
    return_object[["predict_debug_object_y"]] <- ar_out_y
    return_object[["predict_debug_object_x"]] <- ar_out_x
    return_object[["predict_info_matrix"]] <-  as.data.frame(predict_data)


    if(grid_search) {
      return_object[["alpha_grid"]] <- grid_alphas
    }

    return_object[["alpha_used"]] <- alpha

    class(return_object) <- "simulateAR"

  } else {
    return_object <- return_stat
  }
  return( return_object )
}
