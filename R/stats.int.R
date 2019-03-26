#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Set of internal functions                                                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#loadNamespace("stats")
#loadNamespace("multcompView")

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
#' @description Take a numeric vector and return the arithmetic mean.
#'
#' @usage amean(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the arithmetic mean of a vector.
#'
#' @details The function \code{amean} works as the standard function
#'  \code{stats::mean}, but with option \code{na.rm = TRUE} by default.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


amean <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sum(x) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetical standard deviation
#'
#' @description Take a numeric vector and return the arithmetic standard
#'  deviation.
#'
#' @usage asd(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the arithmetic standard deviation of a vector,
#'  as the standard function \code{stats::sd}.
#'
#' @details The fucntion \code{asd} works as the standard function \code{sd},
#'  but with option \code{na.rm = TRUE} by default.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

asd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sqrt( sum((x - amean(x)) ^ 2) / length(x) ), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetical mean
#'
#' @description Take a numeric vector and return the weighted arithmetic mean.
#'
#' @usage wamean(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of numeric weights, of same length as \code{x}.
#'
#' @return Return the weighted arithmetic mean of a vector.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wamean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), sum(w * x) / sum(w), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetic standard deviation
#'
#' @description Take a numeric vector and return the weighted
#' arithmetic standard deviation.
#'
#' @usage wasd(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of numeric weights, of same length as \code{x}.
#'
#' @return Return the weighted arithmetic standard deviation of a vector.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wasd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( sqrt(ifelse(length(x) > 0,
                      sum( w * ((x - sum(x*w)/sum(w)) ^ 2) ) / sum(w),
                      NA)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean
#'
#' @description Take a numeric vector and return its geometric mean.
#'
#' @usage gmean(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the geometric mean of the vector \code{x}.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmean <- function(x) {

  x <- x[!is.na(x)]

  #  return( ifelse(length(x), prod(x) ^ (1/length(x)), NA) )
  return( ifelse(length(x), exp(amean(log(x))), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric standard deviation
#'
#' @description Take a numeric vector and return its geometric
#' standard deviation.
#'
#' @usage gsd(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the geometric standard deviation of the vector \code{x}.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gsd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x),
                 exp( sqrt( sum( log(x/gmean(x)) ^ 2 ) / length(x) ) ),
                 NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric mean
#'
#' @description Take a numeric vector and return its weighted geometric mean.
#'
#' @usage wgmean(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of weights, of same length as x.
#'
#' @return Return the weighted geometric standard deviation of the vector \code{x}.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wgmean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric standard deviation
#'
#' @description Take a numeric vector and return its weighted geometric
#'  standard deviation.
#'
#' @usage wgsd(x, w)
#'
#' @param x a numeric vector.
#'
#' @param w a vector of numeric weights, of same length as \code{x}.
#'
#' @return Return the weighted geometric standard deviation
#' of the vector \code{x}.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

wgsd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Switch for arithmetic or geometric mean according
#' to option "opt.mean"
#'
#' @description Take the numeric vector and return its weighted
#' geometric standard deviation.
#'
#' @usage mean_fct(x, opt.mean = c("amean", "gmean"))
#'
#' @param x a numeric vector.
#'
#' @param opt.mean equals \code{"amean"} or \code{"gmean"}
#' according to that mean value must be computed with arithmetic
#' or geometric formula. There is no default value.
#'
#' @return Return the arithmetic or geometric mean of the vector
#'  \code{x} according to \code{opt.mean}.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mean_fct <- function(x, opt.mean = c("amean", "gmean")) {

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
#
#  Return the Pearson' R2 of a linear regression
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Return the Pearson' R2 of a linear regression
#' @description ddd
#' @usage cor2(x, y, na.rm = TRUE)
#'
#' @param x ff
#' @param y fff
#' @param na.rm ff
#'
#' @return ff
#'
#' @importFrom stats var cov
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cor2 <- function(x, y, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(y))
    x     <- x[index]
    y     <- y[index]
  }

  res <- NA
  if (length(x) > 0)
  {
    z1 <- stats::var(x)
    z2 <- stats::var(y)
    z3 <- stats::cov(x,y)

    if ((z1 > 0) && (z2 > 0)) res <- (z3 ^ 2) / (z1 * z2)
  }

  return(res)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Residual Sum of Square (RSS) of two numeric vectors
#'
#' @description Take two numeric vectors and return their
#' Residual Sum of Square (RSS).
#'
#' @usage rss(x, y)
#'
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Residual Sum of Square (RSS).
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rss <- function(x, y) {

  res <- sum((x - y) ^ 2)

  return( ifelse(res > EPSILON, res, 0) )
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Mean Square Error (MSE) of two numeric vectors
#'
#' @description Take two numeric vectors and return their
#' Mean Square Error (MSE).
#'
#' @usage mse(x, y)
#'
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Mean Square Error (MSE).
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

mse <- function(x, y) {

  index <- !(is.na(x) | is.na(y))
  x     <- x[index]
  y     <- y[index]

  return( ifelse(length(x), rss(x, y) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Root Mean Square Error RMSE) of two numeric vectors
#'
#' @description Take two numeric vectors and return their
#' Root Mean Square Error (RMSE).
#'
#' @usage rmse(x, y)
#'
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Root Mean Square Error (RMSE).
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rmse <- function(x, y) {

  return( sqrt(mse(x, y)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Coefficient of Determination (R2) of two numeric vectors
#'
#' @description Take two numeric vectors and return their Coefficient
#'  of Determination (R2).
#'
#' @usage R2mse(prd, obs)
#'
#' @param prd numeric vector of predicted values.
#'
#' @param obs numeric vector of reference values of \code{length(prd)}.
#'
#' @return Return the Coefficient of Determination (R2).
#'
#' @details Be careful, the function \code{R2mse} is not symmetrical.
#' The first argument \code{prd} is vector of estimated, modelled or
#' predicted values, the second argument \code{obs} is the vector
#' of reference.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
#' @title Probability associated with the Coefficient of Determination (R2)
#'  of two numeric vectors
#'
#' @description Take two numeric vectors and return the
#'  Probability associated with their Coefficient of Determination (R2).
#'
#' @usage pmse(prd, obs, nbK)
#'
#' @param prd a numeric vector of predicted values.
#'
#' @param obs a numeric vector of reference values of \code{length(prd)}.
#'
#' @param nbK an integer. \code{nbK}should be higher than 1: return an
#'  error if not.
#'
#' @return Return the Probability associated with the Coefficient
#' of Determination (R2).
#'
#' @details Be careful, the function \code{pmse} is not symmetrical.
#' The first argument \code{prd} is vector of estimated, modelled
#' or predicted values, the second argument \code{obs} is
#' the vector of reference.
#'
#' @importFrom stats pf
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title AIC of two numeric vectors
#'
#' @description Take two numeric vectors and return their Akaike Information
#'  Criterion (AIC) computed from the Residuals Sum of Square and
#'  the number of used parameters.
#'
#' @usage AIC(x, y, nbK)
#'
#' @param x,y numeric vectors of same length.
#'
#' @param nbK an integer. \code{nbK}should be higher than 1: return
#' an error if not.
#'
#' @return Return the the Akaike Information Criterion (AIC).
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
#' @title AICc of two numeric vectors
#'
#' @description Take two numeric vectors and return the Akaike Information
#'  Criterion (AICc) computed from the Residuals Sum of Square and
#'   the number of used parameters, and corrected for .
#'
#' @usage AICc(x, y, nbK)
#'
#' @param x,y numeric vectors of same length.
#'
#' @param nbK an integer. \code{nbK}should be higher than \code{length(x)}: return an error if not.
#'
#' @return Return the the Akaike Information Criterion (AIC).
#'
#' @details Be careful, the function \code{pmse} is not symmetrical.
#' The first argument is vector of estimated, modelled or predicted values,
#'  the second argument is the vector of reference.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
#' @title Index of the minimum values of a vector
#'
#' @description Take a numeric vector and return the index of
#' its minimum values.
#'
#' @usage argmin(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the minimum values of a vector.
#'
#' @details The function \code{argmin} works as the standard
#' function \code{which.min}.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argmin <- function(x) {

  x <- x[!is.na(x)]

  return( which(x == min(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the first minimum value of a vector
#'
#' @description Take a numeric vector and return the index of the first
#' (the lowest index) minimum values.
#'
#' @usage first_argmin(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the first (the lowest index) minimum values of a vector.
#'
#' @details The function \code{first_argmin} works as the
#' standard function \code{min(which.min)}.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first_argmin <- function(x) {

  x <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] < x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the maximum values of a vector
#'
#' @description Take a numeric vector and return the index of the
#'  maximum values.
#'
#' @usage argmax(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the maximum values of a vector.
#'
#' @details The function \code{argmax} works as the standard function
#'  \code{which.max}.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

argmax <- function(x) {

  x <- x[!is.na(x)]

  return( which(x == max(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the first maximum value of a vector
#'
#' @description Take a numeric vector and return the index of
#' the first (the lowest index) maximum values.
#'
#' @usage first_argmax(x)
#'
#' @param x a numeric vector.
#'
#' @return Return the index of the first (the lowest index)
#'  maximum values of a vector.
#'
#' @details The function \code{first_argmax} works as the
#' standard function \code{min(which.max)}.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

first_argmax <- function(x) {

  x <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] > x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Statistical functions of high level                                      ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Number of Stirling of second kind
#'
#' @description Take an integer and return the number
#'  of Stirling of second kind, that is the number
#'  of all possible partitions into k clusters among n elements
#'
#' @usage stirling(n)
#'
#' @param n a positive integer.
#'
#' @return Return a vector of integer,
#'  of which the \code{names} are the number \code{k} of clusters,
#'   and the \code{values} are the number of partitions
#'   into \code{k} clusters possibly combinated with \code{n} elements.
#'
#' @details None.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

stirling <- function(n) {

  res        <- integer(n)
  names(res) <- seq_len(n)

  for (k in seq_len(n)) {
    tmp <- 0
    for (j in 0:k) tmp <- tmp + (-1) ^ (k - j) * choose(k, j) * j ^ n
    res[k] <- tmp / factorial(k)
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test posthoc of variance analysis
#'
#' @description Proceed to a variance analysis of the vector of data
#' \code{x} associated with the vector of factors \code{clusters},
#' and return a \code{data.frame} containing the name of the groups
#' \code{motif}, the group size \code{number}, the group mean
#' \code{mean} and standard deviation \code{sd}, and the two-to-two
#'  differences expressed in the form of letters \code{group}.
#'
#' @usage test_posthoc(x, clusters, pvalue = 0.05)
#'
#' @param x a numeric vector.
#'
#' @param clusters a vector of factors of same length as \code{x}.
#'
#' @param pvalue a marginal p-value.
#'
#' @return Return a \code{data.frame} containing
#'  the name of the groups \code{motif},
#'   the group size \code{number},
#'  the group mean \code{mean} and standard deviation \code{sd},
#'  and the two-to-two differences expressed
#'  in the form of letters \code{group}.
#'
#' @details \code{test_posthoc} uses Tukey method.
#' Different groups are sorted by decreasing means.
#' Letter rank increases with decreasing means.
#'
#' @importFrom stats aov TukeyHSD
#' @importFrom multcompView multcompLetters4

#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test_posthoc <- function(x, clusters, pvalue = 0.05) {

  # check the inputs
  if (sum(is.na(x)) > 0) {
    clusters <- clusters[!is.na(x)]
    x        <- x[!is.na(x)]
  }

  numbers <- table(clusters)
  levels  <- sort(unique(clusters))
  nblevel <- length(levels)

  index   <- which(table(clusters) > 1)
  if (length(index) > 1) {

    # compute the mean values of groups
    means <- sds <- numeric(nblevel)
    for (i in seq_len(nblevel)) means[i] <- amean(x[clusters == levels[i]])
    for (i in seq_len(nblevel)) sds[i]   <- asd(x[clusters == levels[i]])
    tmp   <- sort(means, decreasing = TRUE, index.return = TRUE)

    # analyse the variances
    if (!is.factor(clusters)) clusters <- as.factor(clusters)

    model <- stats::aov(x ~ clusters)
    test  <- stats::TukeyHSD(x = model, "clusters", conf.level = (1 - pvalue))
    vlet  <- multcompView::multcompLetters4(object = model, comp = test,
                                            threshold = pvalue,
                                            reversed = FALSE)

    tmp   <- sort(means, decreasing = TRUE, index.return = TRUE)

    res   <- data.frame(numbers[tmp$ix], tmp$x, sds[tmp$ix],
                        (vlet$clusters)$Letters)
    # Letters are always sorted out the function multcompLetters4

  } else {

    res   <- data.frame(levels, length(x), amean(x), asd(x), "a")
  }

  colnames(res) <- c("motif", "number", "mean", "sd", "group")
  res$motif     <- as.character(res$motif)
  rownames(res) <- NULL

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Confidence and prediction intervals for a linear model
#'
#' @description Compute the confidence and prediction intervals
#' for a linear regression model \code{model}.
#'
#' @usage confidence_intervals(model, pvalue = 0.05,
#'                             xlim = range(model$model[ , 2]), x.length = 20)
#'
#' @param model a linear model obtained using \code{lm} function.
#'
#' @param pvalue a p-value. \code{pvalue = 0.05} by default.
#'
#' @param xlim a vector of two numeric values.
#' \code{xlim = range(model$model[,2])} by default.
#'
#' @param x.length an integer, indicating the number of points to compute.
#'
#' @return Return a \code{list} containing a vector \code{x} of x-values
#' (ranging from \code{xlim[1]} to \code{xlim[2]}), the corresponding values
#'  \code{y.confidence} for confidence interval, the corresponding values
#'   \code{y.prediction} for prediction interval.
#'
#' @details None.
#'
#' @importFrom stats qt
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

confidence_intervals <- function(model, pvalue = 0.05,
                                 xlim = range(model$model[ , 2]),
                                 x.length = 20) {

  x   <- (model$model)[ , 2]

  n   <- length(x)
  xm  <- amean(x)
  xss <- sum(x ^ 2) - sum(x) ^ 2 / n
  tss <- stats::qt(1 - pvalue/2, (n - 2))

  xx  <- seq(xlim[1], xlim[2], length = 20)

  inter.confidence <-
    tss * sqrt(summary(model)$sigma ^ 2 * (1/n + (xx - xm) ^ 2 / xss))

  inter.prediction <-
    tss * sqrt(summary(model)$sigma ^ 2 * (1 + 1/n + (xx - xm) ^ 2/xss))

  res        <- list(xx, inter.confidence, inter.prediction)
  names(res) <- c("x", "y.confidence", "y.prediction")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two R2
#'
#' @description Tests for the dependence of two R2.
#'
#' @usage test_dependent_R2(v1, v2, v3, n = length(v1))
#'
#' @param v1,v2,v3 three numeric vectors of same length n,
#' or three coefficients of determination. If inputs are three coefficients
#'  of determination, length of vector n must be specified.
#'
#' @param n the length of vectors. \code{n} must be specified
#' if the inputs are three coefficients of determination.
#'
#' @return Return a p-value.
#'
#' @details Be careful, the three vectors are not symmetrical.
#' The function compare R2 obtained by the models \code{v1 ~ v3} and
#'  \code{v2 ~ v3}.
#'
#' @importFrom stats cor pt
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test_dependent_R2 <- function(v1, v2, v3, n = length(v1)) {

  if ((length(v1) > 1) || (length(v2) > 1) || (length(v3) > 1))
  {
    if  ((length(v1) != n) || (length(v2) != n) || (length(v3) != n))
      stop("The vectors v1, v2 and v3 should have a length of n")

    r13 <- abs(stats::cor(v1, v3, method = "pearson"))
    r23 <- abs(stats::cor(v2, v3, method = "pearson"))
    r12 <- abs(stats::cor(v1, v2, method = "pearson"))

  } else {

    r13 <- sqrt(v1)
    r23 <- sqrt(v2)
    r12 <- sqrt(v3)
  }

  R      <- diag(1,3)
  R[2,1] <- R[1,2] <- r13
  R[3,1] <- R[1,3] <- r23
  R[2,3] <- R[3,2] <- r12

  if ((r13 == r23) || (r12 == 1)) {
    p.value <- 1
  } else {
    tobs <- abs(r13 - r23) *
      sqrt((1 + r12) /
             (2*abs(det(R)) / (n - 3) + (r13 + r23) ^ 2 * (1 - r12) ^ 3 /
                (4*(n - 1))))
    p.value <- 2 * (1 - stats::pt(tobs, (n - 3)))
  }

  return(p.value)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two R2
#'
#' @description Tests for the dependence of two R2.
#'
#' @usage test_dependent_R2mse(prd2, prd1, obs)
#'
#' @param prd2,prd1,obs three numeric vectors of same length.
#'
#' @return Return a p-value.
#'
#' @details Be careful, the three vectors are not symmetrical.
#' The function compare R2 obtained by the models \code{prd1 ~ obs}
#' and \code{prd2 ~ obs}.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

test_dependent_R2mse <- function(prd2, prd1, obs) {

  R1 <- R2mse(prd2, obs)
  R2 <- R2mse(prd1, obs)
  R3 <- R2mse(prd2, prd1)

  test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
  test2 <- (is.na(R1) || is.na(R2) || is.na(R3))
  test3 <- ( !((test1 || test2) == TRUE) )

  return( ifelse(test3, test_dependent_R2(R1, R2, R3, length(obs)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two R2
#'
#' @description \code{pvalue_dependent_R2mse} tests for
#' the dependence of two R2.
#'
#' @usage pvalue_dependent_R2mse(mprd, obs)
#'
#' @param mprd a matrix of numeric vectors.
#'
#' @param obs the reference vector of length equals to \code{dim(mprd)[2]}.
#'
#' @return Return a vector of p-value of length equals
#'  to \code{dim(mprd)[1]}.
#'
#' @details The function compare R2 obtained by the models
#'  \code{mprd[i, ] ~ obs}.
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pvalue_dependent_R2mse <- function(mprd, obs) {

  nbline <- dim(mprd)[1]
  res    <- numeric(nbline)

  for (lin in 2:nbline) {

    index <- !is.na(mprd[lin, ])
    if (sum(index) > 0) {

      prd  <- mprd[lin,     index]
      pprd <- mprd[lin - 1, index]
      pobs <- obs[index]

    } else {

      prd  <- mprd[lin,     ]
      pprd <- mprd[lin - 1, ]
      pobs <- obs
    }

    R1 <- R2mse(prd,  pobs)
    R2 <- R2mse(pprd, pobs)
    R3 <- R2mse(pprd, prd)

    test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
    test2 <- (is.na(R1) || is.na(R2) || is.na(R3))

    res[lin] <- NA
    if (!((test1 || test2) == TRUE))
      res[lin] <- test_dependent_R2(R1, R2, R3, length(pobs))
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  END OF FILE
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


