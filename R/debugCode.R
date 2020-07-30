#' Title text
#'
#' Description text
#'
#'
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



library(forcast)
library(forcats)

library(devtools)
install.packages("forecast")

setwd("C:/Users/kkapphahn/Box/Data-Driven Reopening Strategies/kikapp/package/COVIDNearTerm")

# create("COVIDNearTerm")
library(roxygen2)
roxygenise()

install_github("olshena/COVIDNearTerm")
library(COVIDNearTerm)
?modifiedPegelsfcastC


devtools::document()

install.packages("pracma")
library(pracma) # for movavg

dat <- read.table("C:/Users/kkapphahn/Box/Data-Driven Reopening Strategies/kikapp/covid/scc_testing_jun15_taken_jun17.txt",
                  sep="\t",header=FALSE)
dat <- dat[nrow(dat):1,]
num.pos <- dat[,2]
num.neg <- dat[,3]
rate <- (100*num.pos/(num.pos+num.neg))
rate.7 <- movavg(rate, n=7)
rate.28 <- rate[1:28]
rate.35 <- rate[1:35]
rate.42 <- rate[1:42]

set.seed(12345)
test.noboot.28.35 <- simulateHoltMA(rate.28, h=7, npaths=100, bootstrap=FALSE, ma=7, seed = NULL)


dat <- rate.28
h <- 7
npaths <- 10
bootstrap <- FALSE
level <- 0.95
ma <- 7

n.dat <- length(dat)
n.dat.h <- n.dat+h
paths <- matrix(NA, npaths, n.dat.h)
paths[, 1:n.dat] <- matrix(rep(dat, npaths), ncol = n.dat, byrow = TRUE)
model.ets <- forecast::ets(dat, model = "AAN")
paths[,(n.dat+1):n.dat.h] <- modifiedPegelsfcastC(h = h,
                                                  obj = model.ets,
                                                  npaths = npaths,
                                                  level = level,
                                                  bootstrap = bootstrap)$y_paths


which.0 <- which(paths<0)
paths[which.0] <- 0
paths.ma <- t( apply(paths, 1, movavg, n = (ma-1) ) )
smoothed.training <- paths.ma[1,1:n.dat]
max.predictions <- apply(paths.ma[,(n.dat+1):n.dat.h],1,max)
return(list(paths=paths,
            paths.ma=paths.ma,
            smoothed.training=smoothed.training,
            max.predictions=max.predictions))




h = h
obj = model.ets
npaths = npaths
level = level
bootstrap = bootstrap

y_paths <- matrix(NA, nrow = npaths, ncol = h)
obj$lambda <- NULL


for (i in 1:npaths) {
  y_paths[i, ] <- forecast:::simulate.ets(obj, h, future = TRUE, bootstrap = bootstrap)
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






dat <- read.delim("C:/Users/kkapphahn/Box/Data-Driven Reopening Strategies/kikapp/covid/SF_COVID-19_Hospitalizations.csv",header=TRUE,sep=",",as.is=TRUE)

###Look at total covid hopitalizations, variable is chosp

hosp <- data.frame(matrix(dat$PatientCount,ncol=4))
names(hosp) <- c("covid.acute","covid.icu","suspected.acute","suspected.icu")
chosp <- hosp$covid.acute+hosp$covid.icu

vec <- chosp
wsize <- 7
method <- "equal"



test.buildAR.error <- buildAR(chosp, 7, method = "triangle")

buildAR(chosp, 7)



# buildAR_obj, pdays, wsize, nsim, skip=0, seed=12345, method=c("equal","triangle"), output_type="max"

buildAR_obj = test.buildAR.error
pdays = 7
wsize = 7
nsim = 10000
skip = 15
method = "equal"
output_type = "max"


test.predictAR.error <- predictAR(buildAR_obj = test.buildAR.error,
                                  pdays = 7,
                                  wsize = 7,
                                  nsim = 100,
                                  skip = 15,
                                  method = "equal",
                                  output_type = "all")
