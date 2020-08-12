#' Makes predictions using an autogregressive model from the buildAR function
#'
#' Description text
#'
#' @param buildAR_obj An object output from the buildAR function
#' @param pdays Number of days into the future to make predictions
#' @param nsim  Number of simulations
#' @param skip Number of input values to skip
#' @param seed Seed for random number generator
#' @param output_type Type of output

#' @return A data frame containing the specified output for each sim
#'
#' @export
#'
predictAR <- function(buildAR_obj, pdays, nsim, skip=0, seed=NULL, output_type="max") {

  if( class(buildAR_obj) != "buildAR") {
    stop("buildAR_obj must be a buildAR object" )
  }
  if( !is.null(seed) ){
    set.seed(seed)
  }

  method <- buildAR_obj$method
  wsize <- buildAR_obj$wsize

  n_vec <- length(buildAR_obj$vec)
  initphi <- buildAR_obj$initphi
  errors <- buildAR_obj$errors
  diff_phis <- initphi[-1] - initphi[-n_vec]

  ###Adding variable dat that will be used in lowess of error, dropping first point since no prediction
  dat <- buildAR_obj$vec[-1]

  ###Skip the first skip features in the vector
  if( skip > 0 )
  {
    dat <- dat[-(1:skip)]
    errors <- errors[-(1:skip)]
    diff_phis <- diff_phis[-(1:(skip-1))]
  }
  ###Fit lowess for relationship between data and error so error is proportional to magnitude
  lowess_fit <- lowess(dat, errors^2)
  ###Use the median absolute deviation divided by the square root of two to estimate robustly the standard deviation of the phis
  sdphi <- mad(diff_phis) / sqrt(2)

  ###Use unweighted by default, if not, triangle

  if( !(method %in% c("unweighted", "equal", "triangle") ) ) {
    stop("Method must be unweighted, equal or triangle")
  }
  if ( method == "equal" ) {
    weights_phi <- rep(1, wsize) / wsize
  }
  if( method == "triangle" ) {
    weights_phi <- ( 1:wsize ) / sum( 1:wsize )
  }
  if ( method == "unweighted" ) {
    weights_phi <- rep(1, wsize) / wsize
  }


  ###The output is a matrix of nsim by pdays with each row a potential path
  ###Right now loop over nsim then pdays, should be done in a function
  ###Look to separate out errors from phis
  ###Set up output
  # i <- 1
  # j <- 1
  output <- matrix(NA, nsim, pdays)
  for(i in 1:nsim) {
    for(j in 1:pdays)
    {
      ###Need to treat the first day specially because already have phi estimate
      if(j==1) {
        indices <- ( n_vec - wsize + 1 ):n_vec
        current_phis <- initphi[indices]
        weighted_phi <- sum( weights_phi * current_phis )
        current_vec <- buildAR_obj$vec[n_vec]
        new_value <- weighted_phi * current_vec
        new_error <- addError( new_value, lowess_fit )
        output[i,1] <- round( new_value + new_error )
        old_value <- new_value
        old_error <- new_error
      }
      else {
        current_phis <- c( current_phis[ -to_remove ], rnorm(1, weighted_phi, sdphi) )
        weighted_phi <- sum( weights_phi * current_phis )
        new_value <- weighted_phi * old_value
        new_error <- old_error + addError( new_value, lowess_fit )
        output[i,j] <- round( new_value + new_error )
        old_value <- new_value
        old_error <- new_error
      }
      if( method == "unweighted" ){
        to_remove <- sample(x = wsize, size = 1)
      } else {
        to_remove <- 1
      }

    }
  }
  if(output_type == "max") {
    return_object <- data.frame(max = apply(output, 1, max) )
    return_object$rep <- 1:nsim
  }
  if(output_type=="mean") {
    return_object <- data.frame(mean = apply(output, 1, mean) )
    return_object$rep <- 1:nsim
  }
  if(output_type == "all" ) {
    return_object <- data.frame( t( apply(output, 1, summary) ) )
    # return_object$type <- rownames(return_object)
    names( return_object ) <- gsub("\\.|X", "", names( return_object ) )
    names( return_object ) <- gsub("1st", "First", names( return_object ) )
    names( return_object ) <- gsub("3rd", "Third", names( return_object ) )
    return_object$rep <- 1:nsim

  }

  return(return_object)
}
