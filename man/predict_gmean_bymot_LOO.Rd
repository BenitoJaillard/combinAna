% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicting.int.R
\name{predict_gmean_bymot_LOO}
\alias{predict_gmean_bymot_LOO}
\title{Geometric mean (gmean) by motif (bymot) by leave-one-out (LOO)}
\usage{
predict_gmean_bymot_LOO(fct, assMotifs)
}
\arguments{
\item{fct}{a vector of numeric values (assembly properties).}

\item{assMotifs}{a vector of labels of \code{length(fct)} (assembly motifs).}
}
\value{
Return a vector of \code{length(fctMot)}, of which values
 are computed as the geometric mean of all vector elements, excepted
  the value of element to predict that have been left out.
}
\description{
Take a numeric vector and return the predicted vector
 computed as the geometric mean of all elements belonging to a same motif.
}
\details{
Prediction computed using geometric mean \code{gmean}
 by motif \code{bymot} in a whole (WITHOUT taking into account
 species contribution) by excluding the element to predict.
}
\keyword{internal}
