#' Builds an autoregressive model and uses it to predict
#'
#' A wrapper function for buildAR and predictAR.
#'
#'
#' @param vec A vector of numeric data
#' @param x A vector of numeric data, used to predict vec. If not specified, defaults to vec
#' @param wsize Number of prior observations to use for averaging, default is 14
#' @param method Type of weighting to use, default is equal
#' @param pdays Number of days into the future to make predictions, default is 28
#' @param nsim  Number of simulations, default is 100
#' @param skip Number of input values to skip, default is 0
#' @param seed Seed for random number generator
#' @param output_type Type of output, default is all
#' @param rhat_method Method for calculating rhat, if "none", rhat = 1 and has no effect
#' @param lambda Shrinkage parameter, if not specified, default is a grid search from 0 to 1 by 0.05. A value of 0 produces no shrinkage. If an array is specified, all values in the array
#' are evaluated and the optimal lambda is chosen based on residual sum of squares. Values should be between 0 and 1 inclusive.
#' @param debug TRUE returns buildAR objects in addition to standard output

#' @return A data frame containing the specified output statistics for each sim
#'
#' @export
#'
#'
simulateAR <- function(vec, x = NULL,
                       wsize = 14,
                       method = c("equal", "unweighted", "triangle"),
                       pdays = 28,
                       nsim = 100,
                       skip  = 0,
                       seed = NULL,
                       output_type = "all",
                       rhat_method = c("none", "geometric", "arithmetic"),
                       lambda = seq(0, 1, 0.05),
                       debug = FALSE){

  if(length(rhat_method) > 1) {
    rhat_method <- rhat_method[1]
  }

  if( !is.numeric(lambda) ) {
    stop("lambda must be numeric")
  }

  if( any(abs(lambda) > 1 ) ) {
    stop("lambda must be between 0 and 1")
  }

  build_ar_object <- buildAR(vec, x, wsize, method = method, seed = seed, rhat_method = rhat_method, lambda = lambda)

  ar_out <- predictAR(buildAR_obj = build_ar_object,
                      pdays = pdays,
                      nsim = nsim,
                      skip = skip,
                      output_type = output_type,
                      debug = debug)

  if( debug ) {

    return_object <- list("return_stat" = ar_out$return_stat,
                          "predict_object" = ar_out,
                          "buildAR_object" = build_ar_object)

    class(return_object) <- "simulateAR"

  } else {
    return_object <- ar_out$return_stat
  }
  return( return_object )
}
