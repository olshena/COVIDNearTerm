#' Builds an autoregressive model
#'
#' Description text
#'
#'
#'
#' @param vec A vector of numeric data
#' @param x A vector containing a covariate with which to generate AR model, if unspecified, defaults to vec
#' @param wsize Number of prior observations to use for averaging
#' @param method Type of weighting to use
#' @param seed Seed for random number generator

#' @return A list with the following components:
#' @return wsize is the number of prior observations used for averaging
#' @return vec is a numeric vector containing the input data
#' @return fit is a numeric vector containing the autoregressive model fitted data
#' @return errors is a numeric vector containing the difference between fit and vec
#' @return initphi is a numeric vector containing the initial phi
#' @return finalphi is a numeric vector containing the final phi
#'
#' @export
#'

buildAR <- function(vec, x = NULL, wsize, method = c("unweighted", "equal", "triangle"), seed = NULL ) {

  if( !is.null(seed) ){
    set.seed(seed)
  }

  if( is.null(x) ) {
    x <- vec
  }

  if( length(x) != length(vec) ) {
    stop("vec and x must be the same length")
  }

  if(length(method) > 1) {
    cat("Method has length > 1. Using first item:", method[1], "\n")
    method <- method[1]
  }

  if( !(method %in% c("unweighted", "equal", "triangle") ) ) {
    stop("Method must be either unweighted, equal or triangle")
  }

  if ( method == "equal" ) {
    weights_phi <- rep(1, wsize) / wsize
  }

  if( method == "triangle" ) {
    weights_phi <- ( 1:wsize ) / sum( 1:wsize )
  }

  if( method == "unweighted" ) {
    weights_phi <- rep(1, wsize) / wsize
  }

  if( wsize > n_vec ){
    stop("wsize can't be larger than the length of vec")
  }

  n_vec <- length(vec)

  ###Calculate initial phis
  initphi <- c(1, vec[-1] / x[-n_vec] )
  # initphi <- 1:n_vec

  ###Take weighted sum of initia phis to get final phi
  ###Fill in 1s as needed before estimating phis
  finalphi <- rep(NA,n_vec)

  # if( method == "unweighted" ) {
  #   for(i in 1:n_vec) {
  #     if(i == 1) {
  #       to_remove <- 1
  #     } else {
  #       #determine which random day to remove
  #       to_remove <- sample(x = wsize, size = 1)
  #     }
  #
  #     if(i < wsize) {
  #       to_remove <- 1
  #       current_data <- c( rep(1, wsize - i), initphi[1:i] )
  #     } else {
  #       to_remove <- sample(x = wsize, size = 1)
  #       current_data <- c(current_data, initphi[i])
  #     }
  #     # cat("phi calc data\n")
  #     # print(current_data)
  #     finalphi[i] <- sum(weights_phi * current_data)
  #
  #     # cat("Removing", current_data[ to_remove ], "\n")
  #     current_data <- current_data[ -to_remove ]
  #     # print(current_data)
  #
  #   }
  # } else {
  for(i in 1:n_vec) {
    if(i < wsize) {
      new_data <- c( rep(1, wsize - i), initphi[1:i] )
    } else {
      new_data <- initphi[ ( i - wsize + 1 ):i]
    }
    finalphi[i] <- sum(weights_phi * new_data)
  }
  # }
  ###Fits of the data are phi times data, subtract last point because phi only at t-1
  fits <- finalphi[-n_vec] * x[-n_vec]
  ###Errors are data minus fits, subtract first point because it cannot be predicted
  errors <- vec[-1] - fits

  return_object <- list(wsize = wsize,
                        method = method,
                        vec = vec,
                        x = x,
                        fits = fits,
                        errors = errors,
                        initphi = initphi,
                        finalphi = finalphi)

  class(return_object) <- "buildAR"

  return(return_object)
}
