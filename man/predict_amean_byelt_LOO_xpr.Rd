% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predicting.int.R
\name{predict_amean_byelt_LOO_xpr}
\alias{predict_amean_byelt_LOO_xpr}
\title{Arithmetic mean (amean) by elements occurring within
the assemblage and other assemblages belonging to a same motif
 (byelt), except the assemblage to predict (leave-one-out, LOO).
  Over several experiments.}
\usage{
predict_amean_byelt_LOO_xpr(fct, assMotifs, mOccur, xpr)
}
\arguments{
\item{fct}{a vector of numeric values (assembly properties).}

\item{assMotifs}{a vector of labels of \code{length(fct)} (assembly motifs).}

\item{mOccur}{a matrix of occurrence (occurrence of elements).
Its first dimension equals to \code{length(fct)}. Its second dimension
 equals to the number of elements.}

\item{xpr}{a vector of labels of same length as \code{fct}
(assembly experiments).}
}
\value{
Return a vector of same length as \code{fctMot}, of
which values are computed as an arithmetic mean.
}
\description{
The numeric vector \code{fct} get together the properties
of assemblages belonging to different assembly motifs. The properties
 \code{fct} of assemblages belonging to a given assembly motif and
 containing a given element are separately averaged. The property of
 each assemblage is computed as the average of mean values of
 assemblages containing the same elements as the considered assemblage,
  except the assemblage to predict (leave-one-out, LOO). The motif
  of each vector element is specified in the vector \code{assMotifs}.
   The elemental composition of each assemblage is specified in the
   binary matrix \code{mOccur}: \code{0} if the element does not occur,
    \code{1} if the element occurs.
}
\details{
Prediction is computed using arithmetic mean \code{amean}
 by element \code{byelt} occurring within the assemblage and other
  assemblages of a same motif, by including all assemblages belonging
   to a same motif, except the assemblage to predict. Over several
    experiments.
}
\keyword{internal}
