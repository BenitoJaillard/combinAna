#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Set of internal functions                                                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


EPSILON <- 1e-15


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Basic statistics                                                         ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean
#'
#' @description `amean` takes a logical/numeric vector and returns the arithmetic mean.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the arithmetic mean of a vector, as standard 'mean' function.
#' @details `amean` function is as standard `mean` function, but with option `na.rm = TRUE` by default.
#'
#' @examples
#' x <- c(1.2, 2.2, 2.1, NA, 2.6)
#' amean(x)
#'x <- c(NA, NA, NA)
#' amean(x)
#'
#' @export

amean <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sum(x) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetical standard deviation
#'
#' @description `asd` takes a logical/numeric vector and returns the arithmetic standard deviation.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the arithmetic standard deviation of a vector, as standard `stats::sd` function.
#' @details `asd` function is as standard `sd` function, but with option `na.rm = TRUE` by default.
#'
#' @export

asd <- function(x) {

  x  <- x[!is.na(x)]

  return( ifelse(length(x), sqrt( sum((x - amean(x)) ^ 2) / length(x) ), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetical mean
#'
#' @description `wamean` takes a logical/numeric vector and returns the weighted arithmetic mean.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of numeric weights, of same length as x.
#'
#' @return Return the weighted arithmetic mean of a vector.
#' @details None.
#'
#' @export

wamean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), sum(w * x) / sum(w), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetic standard deviation
#' @description `wasd` takes a logical/numeric vector and returns the weighted arithmetic standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of numeric weights, of same length as x.
#'
#' @return Return the weighted arithmetic standard deviation of a vector.
#' @details None.
#'
#' @export

wasd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return(sqrt(ifelse(length(x) > 0,
                     sum( w * ((x - sum(x*w)/sum(w)) ^ 2) ) / sum(w),
                     NA)))
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#' @description `gmean` takes a logical/numeric vector and returns its geometric mean.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the geometric mean of a vector.
#' @details None.
#'
#' @export

gmean <- function(x) {

  x <- x[!is.na(x)]

  #  return( ifelse(length(x), prod(x) ^ (1/length(x)), NA) )
  return( ifelse(length(x), exp(amean(log(x))), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric standard deviation
#' @description `gsd` takes a logical/numeric vector and returns its geometric  standard deviation.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the geometric standard deviation of a vector.
#' @details None.
#'
#' @export

gsd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x),
                 exp( sqrt( sum( log(x/gmean(x)) ^ 2 ) / length(x) ) ),
                 NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric mean
#' @description `wgmean` takes a logical/numeric vector and returns its weighted geometric mean.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of weights, of same length as x.
#'
#' @return Return the weighted geometric standard deviation of a vector.
#' @details None.
#'
#' @export

wgmean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric standard deviation
#' @description `wgsd` takes a logical/numeric vector and returns its weighted geometric standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of numeric weights, of same length as x.
#'
#' @return Return the weighted geometric standard deviation of a vector.
#' @details None.
#'
#' @export

wgsd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Switch for arithmetic or geometric mean according to option "opt.mean"
#' @description `mean_fct` takes a logical/numeric vector and returns its weighted geometric standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param opt.mean The option `opt.mean` equals `amean` or `gmean`. It has no default value.
#'
#' @return Return the arithmetic or geometric mean of a vector.
#' @details None.
#'
#' @export

mean_fct <- function(x, opt.mean) {

  return(switch(opt.mean,
                amean = amean(x),
                gmean = gmean(x)
  ))
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#  Root mean square functions                                              ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the Residual Sum of Square (RSS) of two numeric vectors
#' @description `rss` returns the Residual Sum of Square (RSS) of two numeric vectors.
#'
#' @param x A logical/numeric vector.
#' @param y A logical/numeric vector.
#'
#' @return Return the Residual Sum of Square (RSS).
#' @details None.
#'
#' @export

rss <- function(x, y) {

  res     <- sum((x - y) ^ 2)

  return( ifelse(res > EPSILON, res, 0) )
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the Mean Square Error (MSE) of two numeric vectors
#' @description `mse` returns the Mean Square Error (MSE) of two numeric vectors.
#'
#' @param x A logical/numeric vector.
#' @param y A logical/numeric vector.
#'
#' @return Return the Mean Square Error (MSE).
#' @details None.
#'
#' @export

mse <- function(x, y) {

  index <- !(is.na(x) | is.na(y))
  x     <- x[index]
  y     <- y[index]

  return( ifelse(length(x), sum((x - y) ^ 2) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the Root Mean Square Error RMSE) of two numeric vectors
#' @description `rmse` returns the Root Mean Square Error (RMSE) of two numeric vectors.
#'
#' @param x A logical/numeric vector.
#' @param y A logical/numeric vector.
#'
#' @return Return the Root Mean Square Error (RMSE).
#' @details None.
#'
#' @export

rmse <- function(x, y) {

  return( sqrt(mse(x, y)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the Coefficient of Determination (R2) of two numeric vectors
#' @description `R2mse` returns the Coefficient of Determination (R2) of two numeric vectors.
#'
#' @param prd A logical/numeric vector of predicted values.
#' @param obs A logical/numeric vector of reference values.
#'
#' @return Return the Coefficient of Determination (R2).
#' @details Be careful, the `R2mse` function is not symmetrical. The first argument is vector of estimated, modelled or predicted values, the second argument is the vector of reference.
#'
#' @export

R2mse <- function(prd, obs) {

  index <- !(is.na(prd) | is.na(obs))
  prd   <- prd[index]
  obs   <- obs[index]

  if (length(obs) > 1) {
    tss <- sum((obs - amean(obs)) ^ 2)
    rss <- sum((obs - prd) ^ 2)
    res <- (tss - rss) / tss
  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the Probability associated with the Coeficient of Determination (R2) of two numeric vectors
#' @description `pmse` returns the Probability associated with the Coefficient of Determination (R2) of two numeric vectors.
#'
#' @param prd A logical/numeric vector.
#' @param obs A logical/numeric vector.
#' @param nbK A scalar. `nbK`should be higher than 1. Returns an error if not.
#'
#' @return Return the Probability associated with the Coeficient of Determination (R2).
#' @details Be careful, the `pmse` function is not symmetrical. The first argument `prd` is vector of estimated, modelled or predicted values, the second argument `obs` is the vector of reference.
#'
#' @export

pmse <- function(prd, obs, nbK) {

  if (nbK == 1) stop("nbK should be higher than one")

  index <- !(is.na(prd) | is.na(obs))
  prd   <- prd[index]
  obs   <- obs[index]

  res <- 1
  if (length(obs) > 1) {
    tss   <- sum((obs - amean(obs)) ^ 2)
    rss   <- sum((obs - prd) ^ 2)
    ess   <- tss - rss

    nbfdt <- length(obs)
    nbfde <- nbK - 1
    nbfdr <- nbfdt - nbfde

    if ( (nbfdr > 0) && (nbfde > 0) ) {
      Fratio <- (ess / nbfde) / (rss / nbfdr)
      res    <- 1 - stats::pf(Fratio, nbfde, nbfdr)
    }
  }

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
#  AIC and related functions                                               ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the AIC from the RSS value.
#' @description `AIC` returns the Akaike Information Criterion (AIC) computed from the Residuals Sum of Square and the number of used parameters.
#'
#' @param x A logical/numeric vector.
#' @param y A logical/numeric vector.
#' @param nbK A positive scalar.
#'
#' @return Return the the Akaike Information Criterion (AIC).
#' @details None.
#'
#' @export

AIC <- function(x, y, nbK) {

  n   <- length(y)
  RSS <- rss(x, y)

  if (RSS > EPSILON) {

    res <- n * log(rss(x, y) / n) + 2 * nbK
    if (res == -Inf) res <- NA

  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Compute the AICc from the RSS value.
#' @description `AICc` returns the Akaike Information Criterion (AIC) computed from the Residuals Sum of Square and the number of used parameters, and corrected for .
#'
#' @param x A logical/numeric vector.
#' @param y A logical/numeric vector.
#' @param nbK A scalar. `nbK`should be higher than `length(x)`. Returns an error if not.
#'
#' @return Return the the Akaike Information Criterion (AIC).
#' @details Be careful, the `pmse` function is not symmetrical. The first argument is vector of estimated, modelled or predicted values, the second argument is the vector of reference.
#'
#' @export

AICc <- function(x, y, nbK) {

  n   <- length(y)
  RSS <- rss(x, y)

  if ( (RSS > EPSILON) & (n + 1 > nbK) ) {

    res <- n * log(RSS / n) + 2 * nbK * (n + 2) / (n + 1 - nbK)
    if (res == -Inf) res <- NA

  } else {
    res <- NA
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Returns the index of the minimum values of a vector.
#' @description `argmin` returns the index of the minimum values of a vector.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the index of the minimum values of a vector.
#' @details `argmin` function is as the standard `which.min` function.
#'
#' @export

argmin <- function(x) {

  x   <- x[!is.na(x)]

  return( which(x == min(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Returns the index of the first minimum value of a vector.
#' @description `first.argmin` returns the index of the first (the lowest index) minimum values of a vector.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the index of the first (the lowest index) minimum values of a vector.
#' @details `first.argmin` function is as the standard `min(which.min)` function. #'
#' @export

first.argmin <- function(x) {

  x   <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] < x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Returns the index of the maximum values of a vector.
#' @description `argmax` returns the index of the maximum values of a vector.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the index of the maximum values of a vector.
#' @details `argmax` function is as the standard `which.max` function.
#'
#' @export

argmax <- function(x) {

  x   <- x[!is.na(x)]

  return( which(x == max(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Returns the index of the first maximum value of a vector.
#' @description `first.argmax` returns the index of the first (the lowest index) maximum values of a vector.
#'
#' @param x A logical/numeric vector.
#'
#' @return Return the index of the first (the lowest index) maximum values of a vector.
#' @details `first.argmax` function is as the standard `min(which.max)` function. #'
#' @export

first.argmax <- function(x) {

  x   <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] > x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


