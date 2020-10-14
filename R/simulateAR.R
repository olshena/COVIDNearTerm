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
#' @param lambda Shrinkage parameter, if not specified, default is 0 which produces no shrinkage. If an array is specified, all values in teh array
#' are evaluated and the optimal lambda is chosen based on residual sum of squares. Values should be between 0 and 1 inclusive.
#' @param debug TRUE returns buildAR objects in addition to standard output

#' @return A data frame containing the specified output statistics for each sim
#'
#' @export
#'
#'
simulateAR <- function(vec, x = NULL,
                       wsize,
                       method = c("unweighted", "equal", "triangle"),
                       pdays,
                       nsim,
                       skip  = 0,
                       seed = NULL,
                       output_type = "all",
                       rhat_method = c("none", "geometric", "arithmetic"),
                       lambda = 0,
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
