#' Generates error using SD derived from lowess model
#'
#' Description text
#'
#' @param xpred a predicted x value
#' @param lowess_obj a list containing output from the lowess function

#' @return A random draw from a normal distribution with mean zero and sd derived from lowess object
#'
#' @export
#'


addError <- function(xpred, lowess_obj, n_draws = 1)
{
  # removing duplicates from lowess_obj to avoid
  # duplicate warning from approx()
  lowess_obj_reduced <- as.list( unique( data.frame(lowess_obj) ) )

  est_var <- approx(lowess_obj_reduced, xout = xpred, rule = 2)$y

  # if lowess est_var is negative, use smallest positive observed value
  if( est_var < 0 ) {
    est_var <-  min( lowess_obj_reduced$y[ lowess_obj_reduced$y > 0 ] )
  }

  sd_dat <- sqrt( est_var )

  rnorm(n_draws, 0, sd_dat)
}

