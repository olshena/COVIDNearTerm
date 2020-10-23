#' Makes predictions using an autogregressive model from the buildAR function
#'
#' Description text
#'
#' @param buildAR_obj An object output from the buildAR function
#' @param pdays Number of days into the future to make predictions
#' @param nsim  Number of simulations
#' @param skip Number of input values to skip
#' @param seed Seed for random number generator
#' @param output_type Type of output to predict, currently max, mean, or all are valid

#' @return A data frame containing the specified output for each sim
#'
#' @export
#'
predictAR <- function(buildAR_obj,
                      pdays,
                      nsim,
                      skip=0,
                      seed=NULL,
                      output_type="max",
                      debug = FALSE) {

  if( class(buildAR_obj) != "buildAR") {
    stop("buildAR_obj must be a buildAR object" )
  }
  if( !is.null(seed) ){
    set.seed(seed)
  }

  method <- buildAR_obj$method
  wsize <- buildAR_obj$wsize

  rhat <- buildAR_obj$rhat
  n_vec <- length(buildAR_obj$vec)
  initphi <- buildAR_obj$initphi
  errors <- buildAR_obj$errors
  diff_phis <- initphi[-1] - initphi[-n_vec]

  lambda <- buildAR_obj$lambda

  ###Adding variable dat that will be used in lowess of error, dropping first point since no prediction
  dat <- buildAR_obj$vec[-1]

  ###Skip the first skip features in the vector
  if( skip > 0 ) {
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
  # j <- 2
  if( debug ) {
    debug_output <- expand.grid(sim = 1:nsim,
                                day = 1:pdays,
                                phi_index = 1:wsize,
                                current_phi = NA,
                                rhat = rhat,
                                weighted_phi = NA,
                                lambda = lambda,
                                phi_s = NA,
                                value = NA,
                                error = NA,
                                predicted = NA) %>%
      arrange(sim, day, phi_index)
  }
  output <- matrix(NA, nsim, pdays)
  for(i in 1:nsim) {
    for(j in 1:pdays) {
      ###Need to treat the first day specially because already have phi estimate
      if(j==1) {
        indices <- ( n_vec - wsize + 1 ):n_vec
        current_phis <- initphi[indices]
        weighted_phi <- sum( weights_phi * current_phis ) * rhat
        phi_s <- lambda + ( 1 - lambda ) * weighted_phi
        current_vec <- buildAR_obj$x[n_vec] #buildAR_obj$vec[n_vec]
        new_value <- phi_s * current_vec #* rhat# added rhat product to match how buildAR generates predictions
        new_error <- addError( new_value, lowess_fit )
        output[i,1] <- round( new_value + new_error )
        old_value <- new_value
        old_error <- new_error
      } else {
        current_phis <- c( current_phis[ -to_remove ], rnorm(1, weighted_phi, sdphi) )
        weighted_phi <- sum( weights_phi * current_phis ) * rhat
        phi_s <- lambda + ( 1 - lambda ) * weighted_phi
        new_value <- phi_s * old_value #* rhat
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

      if( debug ) {
        # debug_output %>% head

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "current_phi" ] <- current_phis

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "current_vec" ] <- current_vec

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "weighted_phi" ][1] <- weighted_phi

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "phi_s" ][1] <- phi_s

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "value" ][1] <- new_value

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "error" ][1] <- new_error

        debug_output[debug_output$sim == i &
                       debug_output$day == j, "predicted" ][1] <- new_value + new_error

      }
    }
  }
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

  return_object <- list(pdays = pdays,
                        nsim = nsim,
                        skip = skip,
                        phis = initphi,
                        return_stat = return_stat)

  class(return_object) <- "predictAR"

  if( debug ) {
    return_object[["predict_debug_object"]] <- debug_output
  }

  return(return_object)
}
