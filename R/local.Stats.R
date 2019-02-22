#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Set of internal functions                                                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

# Basic statistics                                                         ####


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean
#'
#' @description `amean` takes a logical/numeric vector and returns the arithmetic mean.
#'
#' @param x A logical/numeric vector.
#' @param na.rm The option.
#'
#' @return Return the arithmetic mean of a vector, as standard 'mean' function.
#' @details None.
#'
#' @examples
#' x <- c(1.2, 2.2, 2.1, NA, 2.6)
#' amean(x, na.rm = FALSE)
#' amean(x, na.rm = TRUE)
#' amean(x)
#'
#' @export

amean <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  return( ifelse(length(x), sum(x) / length(x), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetical standard deviation
#'
#' @description `asd` takes a logical/numeric vector and returns the arithmetic standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the arithmetic standard deviation of a vector, as standard `stats::sd` function.
#' @details None.
#'
#' @export

asd <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  mu <- amean(x, na.rm)
  return( ifelse(length(x), sqrt( sum((x - mu) ^ 2) / length(x) ), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetical mean
#'
#' @description `wamean` takes a logical/numeric vector and returns the weighted arithmetic mean.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of numeric weights, of same length as x.
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the weighted arithmetic mean of a vector.
#' @details None.
#'
#' @export

wamean <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

  return( ifelse(length(x), sum(w * x) / sum(w), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted arithmetic standard deviation
#' @description `wasd` takes a logical/numeric vector and returns the weighted arithmetic standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of numeric weights, of same length as x.
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the weighted arithmetic standard deviation of a vector.
#' @details None.
#'
#' @export

wasd <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

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
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the geometric mean of a vector.
#' @details None.
#'
#' @export

gmean <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

  #  return( ifelse(length(x), prod(x) ^ (1/length(x)), NA) )
  return( ifelse(length(x), exp(amean(log(x))), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric standard deviation
#' @description `gsd` takes a logical/numeric vector and returns its geometric  standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the geometric standard deviation of a vector.
#' @details None.
#'
#' @export

gsd <- function(x, na.rm = TRUE) {

  if (na.rm == TRUE) x <- x[!is.na(x)]

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
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the weighted geometric standard deviation of a vector.
#' @details None.
#'
#' @export

wgmean <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

  return( ifelse(length(x), prod(x ^ w) ^ (1/sum(w)), NA) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Weighted geometric standard deviation
#' @description `wgsd` takes a logical/numeric vector and returns its weighted geometric standard deviation.
#'
#' @param x A logical/numeric vector.
#' @param w A vector of numeric weights, of same length as x.
#' @param na.rm The option `na.rm` indicates if `NA` values must be ignored (`na.rm = TRUE` by default) or not (`na.rm = FALSE`). Return `NA` if there are `NA` values and `na.rm = FALSE`.
#'
#' @return Return the weighted geometric standard deviation of a vector.
#' @details None.
#'
#' @export

wgsd <- function(x, w, na.rm = TRUE) {

  if (na.rm == TRUE) {
    index <- !(is.na(x) | is.na(w))
    x     <- x[index]
    w     <- w[index]
  }

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


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

