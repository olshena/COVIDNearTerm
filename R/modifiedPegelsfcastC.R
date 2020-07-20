#' Modified version of pegelsfcast.C
#'
#' This function, originally called pegelsfcast.C from the forecast package,
#' implements the work in Forecasting: Principles and Practice by
#' Rob J Hyndman and George Athanasopoulos (https://otexts.com/fpp2/)
#'
#' @param h How far in the future to forecast.
#' @param obj A mode object hard-coded to the Holt model.
#' @param npaths Number of simulations.
#' @param level Quantile of distribution. Not utilized but maintained for consistency.
#' @param bootstrap Whether to resample errors (TRUE) or use Gaussian errors (FALSE).

#' @return A list with the following components:
#' @return mu is the mean path
#' @return lower is the lower quantile
#' @return upper is the upper quantile
#' @return y_paths is a matrix of paths
#'
#' @import forecast
#'
#' @export
#'
modifiedPegelsfcastC <- function(h, obj, npaths, level, bootstrap)
{
  y_paths <- matrix(NA, nrow = npaths, ncol = h)
  obj$lambda <- NULL
  for (i in 1:npaths) {
    y_paths[i, ] <- forecast::simulate.ets(obj, h, future = TRUE, bootstrap = bootstrap)
  }
  y_f <- .C("etsforecast",
            as.double(obj$states[length(obj$x) +  1, ]),
            as.integer(obj$m),
            as.integer(switch(obj$components[2], N = 0, A = 1, M = 2)),
            as.integer(switch(obj$components[3], N = 0, A = 1, M = 2)),
            as.double(ifelse(obj$components[4] == "FALSE", 1, obj$par["phi"])),
            as.integer(h),
            as.double(numeric(h)),
            PACKAGE = "forecast")[[7]]
  if (abs(y_f[1] + 99999) < 1e-07) {
    stop("Problem with multiplicative damped trend")
  }
  lower <- apply(y_paths, 2, quantile, 0.5 - level/200, type = 8,
                 na.rm = TRUE)
  upper <- apply(y_paths, 2, quantile, 0.5 + level/200, type = 8,
                 na.rm = TRUE)
  if (length(level) > 1) {
    lower <- t(lower)
    upper <- t(upper)
  }
  return(list(mu = y_f, lower = lower, upper = upper, y_paths=y_paths))
}
