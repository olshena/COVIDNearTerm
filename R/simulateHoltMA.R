#' simulateHoltMA
#'
#' Need a better description of this function
#'
#'
#' @param dat One-dimensional vector of interest like test positivity.
#' @param h How far in the future to forecast.
#' @param npaths Number of simulations.
#' @param level Prediction level, not used but need for existing code
#' @param bootstrap Whether to resample errors (TRUE) or use Gaussian errors (FALSE)
#' @param ma Width of moving average window in days

#' @return A list with the following components:
#' @return paths is raw matrix of paths including training data
#' @return paths.ma is moving average matrix of paths
#' @return smoothed.training is moving average of training data
#' @return max.predictions is maximum for each row of paths.ma


simulateHoltMA <- function(dat, h, npaths, level, bootstrap, ma){
  n.dat <- length(dat)
  n.dat.h <- n.dat+h
  paths <- matrix(NA,npaths,n.dat.h)
  paths[,1:n.dat] <- matrix(rep(dat,npaths),ncol=n.dat,byrow=TRUE)
  model.ets <- ets(dat,model="AAN")
  paths[,(n.dat+1):n.dat.h] <- modifiedPegelsfcastC(h=h,
                                                    obj=model.ets,
                                                    npaths=npaths,
                                                    level=level,
                                                    bootstrap=bootstrap)$y.paths
  which.0 <- which(paths<0)
  paths[which.0] <- 0
  paths.ma <- t(apply(paths,1,movavg,n=(ma-1)))
  smoothed.training <- paths.ma[1,1:n.dat]
  max.predictions <- apply(paths.ma[,(n.dat+1):n.dat.h],1,max)
  return(list(paths=paths,
              paths.ma=paths.ma,
              smoothed.training=smoothed.training,
              max.predictions=max.predictions))
}
