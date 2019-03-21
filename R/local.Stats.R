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
#' @description Take a numeric vector and return the arithmetic mean.
#' @usage amean(x)
#' @param x a numeric vector.
#' @return Return the arithmetic mean of a vector, as standard 'mean' function.
#' @details `amean` function does as the standard `stats::mean` function, but with option `na.rm = TRUE` by default.
#'
#' @keywords internal

amean <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sum(x) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetical standard deviation
#' @description Take a numeric vector and return the arithmetic standard deviation.
#' @usage asd(x)
#' @param x a numeric vector.
#' @return Return the arithmetic standard deviation of a vector, as the standard `stats::sd` function.
#' @details `asd` function is as standard `sd` function, but with option `na.rm = TRUE` by default.
#'
#' @keywords internal

asd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x), sqrt( sum((x - amean(x)) ^ 2) / length(x) ), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetical mean
#' @description Take a numeric vector and return the weighted arithmetic mean.
#' @usage wamean(x, w)
#' @param x a numeric vector.
#' @param w a vector of numeric weights, of same length as x.
#'
#' @return Return the weighted arithmetic mean of a vector.
#' @details None.
#'
#' @keywords internal

wamean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), sum(w * x) / sum(w), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetic standard deviation
#' @description Take a numeric vector and return the weighted arithmetic standard deviation.
#' @usage wasd(x, w)
#' @param x a numeric vector.
#' @param w a vector of numeric weights, of same length as x.
#'
#' @return Return the weighted arithmetic standard deviation of a vector.
#' @details None.
#'
#' @keywords internal

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
#' @description Take a numeric vector and return its geometric mean.
#' @usage gmean(x)
#' @param x a numeric vector.
#'
#' @return Return the geometric mean of a vector.
#' @details None.
#'
#' @keywords internal

gmean <- function(x) {

  x <- x[!is.na(x)]

  #  return( ifelse(length(x), prod(x) ^ (1/length(x)), NA) )
  return( ifelse(length(x), exp(amean(log(x))), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric standard deviation
#' @description Take a numeric vector and return its geometricstandard deviation.
#' @usage gsd(x)
#' @param x a numeric vector.
#'
#' @return Return the geometric standard deviation of a vector.
#' @details None.
#'
#' @keywords internal

gsd <- function(x) {

  x <- x[!is.na(x)]

  return( ifelse(length(x),
                 exp( sqrt( sum( log(x/gmean(x)) ^ 2 ) / length(x) ) ),
                 NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric mean
#' @description Take a numeric vector and return its weighted geometric mean.
#' @usage wgmean(x, w)
#' @param x a numeric vector.
#' @param w a vector of weights, of same length as x.
#'
#' @return Return the weighted geometric standard deviation of a vector.
#' @details None.
#'
#' @keywords internal

wgmean <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric standard deviation
#' @description Take a numeric vector and return its weighted geometric standard deviation.
#' @usage wgsd(x, w)
#' @param x a numeric vector.
#' @param w a vector of numeric weights, of same length as x.
#'
#' @return Return the weighted geometric standard deviation of a vector.
#' @details None.
#'
#' @keywords internal

wgsd <- function(x, w) {

  index <- !(is.na(x) | is.na(w))
  x     <- x[index]
  w     <- w[index]

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Switch for arithmetic or geometric mean according to option "opt.mean"
#' @description Take the numeric vector and return its weighted geometric standard deviation.
#' @usage mean_fct(x, opt.mean = c("amean", "gmean"))
#' @param x a numeric vector.
#' @param opt.mean equals `amean` or `gmean` according to that mean value must be computed with arithmetic or geometric formula. There is no default value.
#'
#' @return Return the arithmetic or geometric mean of the vector `x` according to opt.mean.
#' @details None.
#'
#' @keywords internal

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
#'
#' @title Residual Sum of Square (RSS) of two numeric vectors
#' @description Take two numeric vectors and return their Residual Sum of Square (RSS).
#' @usage rss(x, y)
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Residual Sum of Square (RSS).
#' @details None.
#'
#' @export

rss <- function(x, y) {

  res <- sum((x - y) ^ 2)

  return( ifelse(res > EPSILON, res, 0) )
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Mean Square Error (MSE) of two numeric vectors
#' @description Take two numeric vectors and return their Mean Square Error (MSE).
#' @usage mse(x, y)
#' @param x,y numeric vectors of same length.
#'
#' @return Return the Mean Square Error (MSE).
#' @details None.
#'
#' @export

mse <- function(x, y) {

  index <- !(is.na(x) | is.na(y))
  x     <- x[index]
  y     <- y[index]

  return( ifelse(length(x), rss(x, y) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Root Mean Square Error RMSE) of two numeric vectors
#' @description Take two numeric vectors and return their Root Mean Square Error (RMSE).
#' @usage rmse(x, y)
#' @param x,y numeric vectors of same length.
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
#' @title Coefficient of Determination (R2) of two numeric vectors
#' @description Take two numeric vectors and return their Coefficient of Determination (R2).
#' @usage R2mse(prd, obs)
#' @param prd numeric vector of predicted values.
#' @param obs numeric vector of reference values of same `length` as `prd`.
#'
#' @return Return the Coefficient of Determination (R2).
#' @details Be careful, the `R2mse` function is not symmetrical. The first argument `prd` is vector of estimated, modelled or predicted values, the second argument `obs` is the vector of reference.
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
#' @title Probability associated with the Coefficient of Determination (R2) of two numeric vectors
#' @description Take two numeric vectors and return the Probability associated with their Coefficient of Determination (R2).
#' @usage pmse(prd, obs, nbK)
#' @param prd a numeric vector of predicted values.
#' @param obs a numeric vector of reference values of same `length` as `prd`.
#' @param nbK an integer. `nbK`should be higher than 1: return an error if not.
#'
#' @return Return the Probability associated with the Coefficient of Determination (R2).
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
#' @title AIC of two numeric vectors
#' @description Take two numeric vectors and return their Akaike Information Criterion (AIC) computed from the Residuals Sum of Square and the number of used parameters.
#' @usage .AIC(x, y, nbK)
#' @param x,y numeric vectors of same length.
#' @param nbK an integer. `nbK`should be higher than 1: return an error if not.
#'
#' @return Return the the Akaike Information Criterion (AIC).
#' @details None.
#'
#' @export

.AIC <- function(x, y, nbK) {

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
#' @description Take two numeric vectors and return the Akaike Information Criterion (AICc) computed from the Residuals Sum of Square and the number of used parameters, and corrected for .
#' @usage .AICc(x, y, nbK)
#' @param x,y numeric vectors of same length.
#' @param nbK an integer. `nbK`should be higher than `length(x)`: return an error if not.
#'
#' @return Return the the Akaike Information Criterion (AIC).
#' @details Be careful, the `pmse` function is not symmetrical. The first argument is vector of estimated, modelled or predicted values, the second argument is the vector of reference.
#'
#' @export

.AICc <- function(x, y, nbK) {

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
#' @description Take a numeric vector and return the index of its minimum values.
#' @usage argmin(x)
#' @param x a numeric vector.
#'
#' @return Return the index of the minimum values of a vector.
#' @details `argmin` function is as the standard `which.min` function.
#'
#' @export

argmin <- function(x) {

  x <- x[!is.na(x)]

  return( which(x == min(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the first minimum value of a vector
#' @description Take a numeric vector and return the index of the first (the lowest index) minimum values.
#' @usage first.argmin(x)
#' @param x a numeric vector.
#'
#' @return Return the index of the first (the lowest index) minimum values of a vector.
#' @details `first.argmin` function is as the standard `min(which.min)` function. #'
#' @export

first.argmin <- function(x) {

  x <- x[!is.na(x)]

  arg <- 1
  while ((x[arg + 1] < x[arg]) & (arg < length(x))) arg <- arg + 1

  return(arg)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the maximum values of a vector
#' @description Take a numeric vector and return the index of the maximum values.
#' @usage argmax(x)
#' @param x a numeric vector.
#'
#' @return Return the index of the maximum values of a vector.
#' @details `argmax` function is as the standard `which.max` function.
#'
#' @export

argmax <- function(x) {

  x <- x[!is.na(x)]

  return( which(x == max(x, na.rm = TRUE)) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Index of the first maximum value of a vector
#' @description Take a numeric vector and return the index of the first (the lowest index) maximum values.
#' @usage first.argmax(x)
#' @param x a numeric vector.
#'
#' @return Return the index of the first (the lowest index) maximum values of a vector.
#' @details `first.argmax` function is as the standard `min(which.max)` function. #'
#' @export

first.argmax <- function(x) {

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
#' @description Take an integer and return the number of Stirling of second kind, that is the number of all possible partitions into k clusters among n elements
#' @usage stirling(n)
#' @param n a positive integer.
#'
#' @return Return a vector of integer, of which the `names` are the number `k` of clusters, and the `values` are the number of partitions into `k` clusters possibly combinated with `n` elements.
#' @details None.
#'
#' @examples
#'  stirling(12)
#'
#' @export

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
#' @description Proceed to a variance analysis of the vector of data `x` associated with the vector of factors `clusters`, and return a `data.frame` containing the name of the groups `motif`, the group size `number`, the group mean `mean` and standard deviation `sd`, and the two-to-two differences expressed in the form of letters `group`.
#' @usage test.posthoc(x, clusters, pvalue = 0.05)
#' @param x a numeric vector.
#' @param clusters a vector of factors of same length as `x`.
#' @param pvalue a marginal p-value.
#'
#' @return Return a `data.frame` containing the name of the groups `motif`, the group size `number`, the group mean `mean` and standard deviation `sd`, and the two-to-two differences expressed in the form of letters `group`.
#' @details Posthoc test uses Tukey method. Different groups are sorted by decreasing means. Letter rank increases with decreasing means.
#'
#' @export

test.posthoc <- function(x, clusters, pvalue = 0.05) {

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
#' @description Compute the confidence and prediction intervals for a linear regression model `model`.
#' @usage confidence.intervals(model, pvalue = 0.05,
#'                             xlim = range(model$model[ , 2]), x.length = 20)
#' @param model a linear model obtained using `lm` function.
#' @param pvalue a p-value. `pvalue = 0.05` by default.
#' @param xlim a vector of two numeric values. `xlim = range(model$model[,2])` by default.
#' @param x.length an integer, indicating the number of points to compute.
#'
#' @return Return a `list` containing a vector `x` of x-values (ranging from `xlim[1]` to `xlim[2]`), the corresponding values `y.confidence` for confidence interval, the corresponding values `y.prediction` for prediction interval.
#'
#' @details None.
#'
#' @export

confidence.intervals <- function(model, pvalue = 0.05,
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
#' @description Tests for the dependence of two R2.
#' @usage test.dependent.R2(v1, v2, v3, n = length(v1))
#' @param v1,v2,v3 three numeric vectors of same length n, or three coefficients of determination. If inputs are three coefficients of determination, length of vector n must be specified.
#' @param n the length of vectors. n must be specified if the inputs are three coefficients of determination.
#'
#' @return Return a p-value.
#' @details None.
#'
#' @export

test.dependent.R2 <- function(v1, v2, v3, n = length(v1)) {

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
#' @description Tests for the dependence of two R2.
#' @usage test.dependent.R2mse(prd2, prd1, obs)
#' @param prd2,prd1 three numeric vectors of same length (opt = "x"), or three coefficients of determination if opt = "R2".
#' @param obs an option for indicating if the inputs are vectors (opt = "x") or coefficients of determination (opt = "R2")
#'
#' @return Return a p-value.
#' @details None.
#'
#' @export

#  Compute the differences between successive R2 are significant

test.dependent.R2mse <- function(prd2, prd1, obs) {

  R1 <- R2mse(prd2, obs)
  R2 <- R2mse(prd1, obs)
  R3 <- R2mse(prd2, prd1)

  test1 <- ((R1 < 0) || (R2 < 0) || (R3 < 0))
  test2 <- (is.na(R1) || is.na(R2) || is.na(R3))
  test3 <- ( !((test1 || test2) == TRUE) )

  return( ifelse(test3, test.dependent.R2(R1, R2, R3, length(obs)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Test for the dependence of two R2
#' @description `pvalue.dependent.R2mse` tests for the dependence of two R2.
#' @usage pvalue.dependent.R2mse(mprd, obs)
#' @param mprd v2 v3 reg A linear model obtained using `lm` function.
#' @param obs v2 v3 reg A linear model obtained using `lm` function.
#'
#' @return Return a `list` containing a vector `x` of x-values (ranging from `xlim[1]` to `xlim[2]`), the corresponding values `y.confidence` for confidence interval, the corresponding values `y.prediction` for prediction interval.
#'
#' @details None.
#'
#' @export
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute the differences between successive R2 are significant
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

pvalue.dependent.R2mse <- function(mprd, obs) {

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
      res[lin] <- test.dependent.R2(R1, R2, R3, length(pobs))
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


