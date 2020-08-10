#' simulateHoltMA
#'
#' Need a better description of this function
#'
#'
#' @param dat One-dimensional vector of interest like test positivity
#' @param h How far in the future to forecast
#' @param npaths Number of simulations
#' @param level Prediction level, not used but need for existing code
#' @param bootstrap Whether to resample errors (TRUE) or use Gaussian errors (FALSE)
#' @param ma Width of moving average window in days
#' @param smooth_type Type of smoothing used by movavg, see movavg documentation for detials
#' @param model_type 	From the ets documentation: Usually a three-character string identifying method using the framework terminology of Hyndman et al. (2002) and Hyndman et al. (2008). The first letter denotes the error type ("A", "M" or "Z"); the second letter denotes the trend type ("N","A","M" or "Z"); and the third letter denotes the season type ("N","A","M" or "Z"). In all cases, "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. So, for example, "ANN" is simple exponential smoothing with additive errors, "MAM" is multiplicative Holt-Winters' method with multiplicative errors, and so on.
#' @param seed Random seed to use

#' @return A list with the following components:
#' @return paths is raw matrix of paths including training data
#' @return paths_ma is moving average matrix of paths
#' @return smoothed_training is moving average of training data
#' @return max_predictions is maximum for each row of paths_ma
#' @import forecast
#' @import pracma

#' @export
#'
simulateHoltMA <- function(dat, h, npaths, level = 0.95, bootstrap, ma, smooth_type = "s", model_type = "AAN", seed = NULL){

  if( !is.null(seed) ){
    set.seed(seed)
  }
  n.dat <- length(dat)
  n.dat.h <- n.dat+h
  paths <- matrix(NA,npaths,n.dat.h)
  paths[,1:n.dat] <- matrix(rep(dat,npaths),ncol=n.dat,byrow=TRUE)
  model.ets <- forecast::ets(dat, model = model_type)
  paths[,(n.dat+1):n.dat.h] <- modifiedPegelsfcastC(h = h,
                                                    obj = model.ets,
                                                    npaths = npaths,
                                                    level = level,
                                                    bootstrap = bootstrap)$y_paths
  which.0 <- which(paths<0)
  paths[which.0] <- 0
  paths_ma <- t(apply(paths,1, pracma::movavg, n=(ma-1), type = smooth_type) )
  smoothed_training <- paths_ma[1,1:n.dat]
  max_predictions <- apply(paths_ma[,(n.dat+1):n.dat.h],1,max)
  return(list(paths = paths,
              paths_ma = paths_ma,
              smoothed_training = smoothed_training,
              max_predictions = max_predictions))
}
