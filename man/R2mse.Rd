% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats.int.R
\name{R2mse}
\alias{R2mse}
\title{Coefficient of Determination (R2) of two numeric vectors}
\usage{
R2mse(prd, obs)
}
\arguments{
\item{prd}{numeric vector of predicted values.}

\item{obs}{numeric vector of reference values of \code{length(prd)}.}
}
\value{
Return the Coefficient of Determination (R2).
}
\description{
Take two numeric vectors and return their Coefficient
 of Determination (R2).
}
\details{
Be careful, the function \code{R2mse} is not symmetrical.
The first argument \code{prd} is vector of estimated, modelled or
predicted values, the second argument \code{obs} is the vector
of reference.
}
\keyword{internal}
