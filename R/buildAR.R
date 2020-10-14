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
#' @param rhat_method Method for calculating rhat, if "none", rhat = 1 and has no effect
#' @param lambda Shrinkage parameter, if not specified, default is 0 which produces no shrinkage. If an array is specified, all values in teh array
#' are evaluated and the optimal lambda is chosen based on residual sum of squares. Values should be between 0 and 1 inclusive.

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

buildAR <- function(vec,
                    x = NULL,
                    wsize,
                    method = c("unweighted", "equal", "triangle"),
                    seed = NULL,
                    rhat_method = c("none", "geometric", "arithmetic"),
                    lambda = 0) {

  if( !is.null(seed) ){
    set.seed(seed)
  }

  if( any(vec == 0)) {
    stop("vec can not have zero-values")
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

  if(length(rhat_method) > 1) {
    rhat_method <- rhat_method[1]
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

  if ( rhat_method == "none" ) {
    rhat <- 1
  }

  n_vec <- length(vec)

  if ( rhat_method == "geometric" ) {

    r_i <- vec[-1] / x[-n_vec]
    if( any(r_i <= 0 ) ){
      cat("geometric mean specified for rhat_method, but r_i value(s) <= 0, using arithmetic mean\n")
      rhat <- mean( r_i, na.rm = TRUE )
    } else {
      rhat <-  exp( sum( log( r_i ) ) / length( r_i ) )
    }
  }

  if ( rhat_method == "arithmetic" ) {
    r_i <- vec[-1] / x[-n_vec]
    rhat <- mean( r_i )
  }

  if( wsize > n_vec ){
    stop("wsize can't be larger than the length of vec")
  }


  if( !is.numeric(lambda) ) {
    stop("lambda must be numeric")
  }

  if( any(abs(lambda) > 1 ) ) {
    stop("lambda must be between 0 and 1")
  }

  ###Calculate initial phis
  initphi <-  c(1,  (vec[-1] /  x[-n_vec] ) / rhat )

  ###Take weighted sum of initia phis to get final phi
  ###Fill in 1s as needed before estimating phis
  finalphi <- rep(NA,n_vec)

  if( length(lambda) == 1 ) {

    for(i in 1:n_vec) {
      if(i < wsize) {
        new_data <- c( rep(1, wsize - i), initphi[1:i] )
      } else {
        new_data <- initphi[ ( i - wsize + 1 ):i]
      }

      finalphi[i] <- sum(weights_phi * new_data)
    }

    phi_s <- (lambda + (1 - lambda) * finalphi[-n_vec] )

    ###Fits of the data are phi times data, subtract last point because phi only at t-1
    fits <- phi_s *
      x[-n_vec] *
      rhat

    ###Errors are data minus fits, subtract first point because it cannot be predicted
    errors <- vec[-1] - fits

    return_object <- list(wsize = wsize,
                          method = method,
                          vec = vec,
                          x = x,
                          rhat_method = rhat_method,
                          rhat = rhat,
                          lambda = lambda,
                          fits = fits,
                          errors = errors,
                          initphi = initphi,
                          finalphi = finalphi,
                          phi_s = phi_s)

    class(return_object) <- "buildAR"

    return(return_object)

  }

  if( length( lambda ) > 1 ) {

    names(lambda) <- lambda

    buildAR_objects <- purrr:::map(lambda, function(.lambda, .vec, .x, .wsize, .method, .rhat_method) {

      # print(.lambda)
      buildAR(vec = .vec,
              x = .x,
              wsize = .wsize,
              method = .method,
              rhat_method = .rhat_method,
              lambda = .lambda)

    }, vec, x, wsize, method, rhat_method)

    SSRs <- purrr:::map_dfr(buildAR_objects, function(.buildAR_object) {

      data.frame(SSR = sum( .buildAR_object$error^2  ) ,
                 lambda = .buildAR_object$lambda)

    })

    lambda_optimal <- SSRs %>%
      arrange(SSR) %>%
      head(1) %>%
      select(lambda) %>%
      unlist %>%
      paste0

    return_object <- buildAR_objects[[lambda_optimal]]

    return_object$SSR <- SSRs

    return(return_object)

  }

}
