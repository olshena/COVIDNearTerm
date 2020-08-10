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

#' @return A data frame containing the specified output statistics for each sim
#'
#' @export
#'
#'
simulateAR <- function(vec, wsize, method = c("unweighted", "equal", "triangle"), pdays, nsim, skip  = 0, seed = NULL, output_type = "all"){

  build_ar_object <- buildAR(vec, wsize, method = method, seed = seed)

  ar_out <- predictAR(buildAR_obj = build_ar_object,
                      pdays = pdays,
                      nsim = nsim,
                      skip = skip,
                      output_type = output_type)

  return( ar_out )
}
