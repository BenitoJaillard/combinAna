% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicting.int.R
\name{predict_amean_byelt_jack}
\alias{predict_amean_byelt_jack}
\title{Arithmetic mean (amean) by motif (bymot) by jackknife (jack)}
\usage{
predict_amean_byelt_jack(fct, assMotifs, mOccur, jack)
}
\arguments{
\item{fct}{a vector of numeric values (assembly properties).}

\item{assMotifs}{a vector of labels of \code{length(fct)} (assembly motifs).}

\item{mOccur}{a matrix of occurrence (occurrence of elements).
Its first dimension equals to \code{length(fct)}. Its second
dimension equals to the number of elements.}

\item{jack}{a vector of two elements. The first one \code{jack[1]}
specifies the size of subset, the second one \code{jack[2]} specifies
 the number of subsets.}
}
\value{
Return a vector of \code{length(fct)}, of which values are
 computed as the arithmetic mean of all vector elements.
}
\description{
Take a numeric vector and return the predicted vector
 computed as the arithmetic mean of all elements belonging to the same motif.
}
\details{
Prediction is computed using arithmetic mean \code{amean}
by motif \code{bymot} in a whole (WITHOUT taking into account species
 contribution). The elements belonging to a same motif are divided
  into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
   is computed by excluding \code{jack[1]} elements, of which the
   element to predict. If the total number of elements belonging to
    the motif is lower than \code{jack[1]*jack[2]}, prediction is
    computed by Leave-One-Out (LOO).
}
\keyword{internal}
