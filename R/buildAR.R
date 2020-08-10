#' Builds an autoregressive model
#'
#' Description text
#'
#'
#'
#' @param vec A vector of numeric data
#' @param wsize Number of prior observations to use for averaging
#' @param method Type of weighting to use

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

buildAR <- function(vec, wsize, method = c("equal","triangle") ) {
  n_vec <- length(vec)

  ###Use equal or triangle weighting
  if ( method == "equal" ) {
    weights_phi <- rep(1, wsize) / wsize
  }
  if( method == "triangle" ) {
    weights_phi <- ( 1:wsize ) / sum( 1:wsize )
  }
  if( !(method %in% c("equal","triangle") ) ) {
    stop("Method must be either equal or triangle")
  }

  if( wsize > n_vec ){
    stop("wsize can't be larger than the length of vec")
  }

  ###Calculate initial phis
  initphi <- c(1, vec[-1] / vec[-n_vec] )
  ###Take weighted sum of initia phis to get final phi
  ###Fill in 1s as needed before estimating phis
  finalphi <- rep(NA,n_vec)

  for(i in 1:n_vec) {
    if(i < wsize) {
      new_data <- c( rep(1, wsize - i), initphi[1:i] )
    } else {
      new_data <- initphi[ ( i - wsize + 1 ):i]
    }
    finalphi[i] <- sum(weights_phi * new_data)

  }
  ###Fits of the data are phi times data, subtract last point because phi only at t-1
  fits <- finalphi[-n_vec] * vec[-n_vec]
  ###Errors are data minus fits, subtract first point because it cannot be predicted
  errors <- vec[-1] - fits

  return_object <- list(wsize = wsize,
                        method = method,
                        vec = vec,
                        fits = fits,
                        errors = errors,
                        initphi = initphi,
                        finalphi = finalphi)

  class(return_object) <- "buildAR"

  return(return_object)
}
