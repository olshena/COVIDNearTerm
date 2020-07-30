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


addError <- function(xpred, lowess_obj)
{
  # removing duplicates from lowess_obj to avoid
  # duplicate warning from approx()
  lowess_obj_reduced <- as.list( unique( data.frame(lowess_obj) ) )

  sd_dat <- sqrt( approx(lowess_obj_reduced, xout = xpred, rule = 2)$y)
  rnorm(1, 0, sd_dat)
}

