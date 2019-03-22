#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Set of internal functions                                                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include local.Stats.R
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by Motif                                                 ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot)
#'
#' @description Take the numeric vector \code{fct} and return a vector of same
#'  length, of which values are computed as the arithmetic mean of all vector
#'   elements belonging to a same motif. The motif of each vector element is
#'    specified in the vector \code{assMotifs}.
#'
#' @usage predictAmeanBymot(fct, assMotifs)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct} (assembly
#'  motifs).
#'
#' @return Return a vector of same length as \code{fct}, of which values are
#'  computed as the arithmetic mean of all vector elements belonging to a same
#'  motif.
#'
#' @details Prediction is computed using arithmetic mean \code{amean} by motif
#'  \code{bymot} in a whole (WITHOUT taking into account species contribution)
#'   by including all elements belonging to a same motif, even the one to
#'   predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanBymot <- function(fct, assMotifs) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- amean(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) over several experiments
#' (xpr)
#'
#' @description Take the numeric vector \code{fct} and return a vector of same
#' length, of which values are computed as the arithmetic mean of all vector
#' elements belonging to a same motif, over several experiments. The motif of
#'  each vector element is specified in the vector \code{assMotifs}. The
#'   experiment of each vector element is specified in the vector \code{xpr}.
#'
#' @usage predictAmeanBymotXpr(fct, assMotifs, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#' (assembly motifs).
#'
#' @param xpr a vector of labels of same length as \code{fct}
#' (assembly experiments).
#'
#' @return Return a vector of same length as \code{fct}, of which values
#' are computed as the arithmetic mean of all vector elements belonging to
#'  the same motif over several experiments.
#'
#' @details Prediction is computed using arithmetic mean \code{amean} by
#' motif \code{bymot} in a whole (WITHOUT taking into account species
#' contribution) by including all elements belonging to a same motif,
#'  even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanBymotXpr <- function(fct, assMotifs, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictAmeanBymot(fct[indXpr], assMotifs[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by leave-one-out (LOO)
#'
#' @description Take the numeric vector \code{fctMot} and return a vector of
#' same length, of which values are computed as the arithmetic mean of all
#' vector elements.
#'
#' @usage ameanBymotLOO(fctMot)
#'
#' @param fctMot a vector of numeric values of elements belonging to a same
#'  motif.
#'
#' @return Return a vector of same length as \code{fctMot}, of which values
#' are computed as the arithmetic mean of all vector elements, excepted the
#' value of element to predict that have been left out.
#'
#' @details Prediction computed using arithmetic mean \code{amean} by motif
#'  \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution) by excluding the element to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ameanBymotLOO <- function(fctMot) {

  nbass <- length(fctMot)

  if (nbass > 1) {
    fctPrd <- numeric(nbass)
    for (ind in seq_len(nbass)) fctPrd[ind] <- amean(fctMot[-ind])
  } else {
    fctPrd <- NA
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by leave-one-out (LOO)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predictAmeanBymotLOO(fct, assMotifs)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#'  (assembly motifs).
#'
#' @return Return a vector of same length as \code{fctMot}, of which
#' values are computed as the arithmetic mean of all vector elements,
#'  excepted the value of element to predict that have been left out.
#'
#' @details Prediction computed using arithmetic mean \code{amean} by
#'  motif \code{bymot} in a whole (WITHOUT taking into account species
#'   contribution) by excluding the element to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanBymotLOO <- function(fct, assMotifs) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- ameanBymotLOO(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by leave-one-out (LOO)
#'  over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predictAmeanBymotLOOXpr(fct, assMotifs, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#' (assembly motifs).
#'
#' @param xpr a vector of labels of same length as \code{fct}
#' (assembly experiments).
#'
#' @return Return a vector of same length as \code{fctMot}, of which
#' values are computed as the arithmetic mean of all vector elements,
#'  over  several experiments, excepted the value of element to predict
#'  that have been left out.
#'
#' @details Prediction computed using arithmetic mean \code{amean} by
#' motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution) by excluding the element to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanBymotLOOXpr <- function(fct, assMotifs, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictAmeanBymotLOO(fct[indXpr], assMotifs[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack)
#'  over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage ameanBymotJack(fctMot, jack)
#'
#' @param fctMot a vector of numeric values of elements belonging to a
#'  same motif.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'   the number of subsets.
#'
#' @return Return a vector of same length as \code{fctMot},
#'  of which values are computed as the arithmetic mean of all vector elements.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#' contribution). The elements belonging to a same motif are divided
#' into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'  is computed by excluding \code{jack[1]} elements, of which the element
#'   to predict. If the total number of elements belonging to the motif
#'   is lower than \code{jack[1]*jack[2]}, prediction is computed by
#'    Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ameanBymotJack <- function(fctMot, jack) {

  nbass <- length(fctMot)

  if (nbass > jack[1] * jack[2]) {

    fctPrd <- numeric(nbass)
    index  <- sample.int(nbass)
    size   <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack         <- index[(ind - 1) * size + (1:size)]
      fctPrd[indjack] <- amean(fctMot[-indjack])
    }

    indjack         <- index[(ind * size + 1):nbass]
    fctPrd[indjack] <- amean(fctMot[-indjack])

  } else {

    fctPrd <- ameanBymotLOO(fctMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage predictAmeanBymotJack(fct, assMotifs, jack)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#' (assembly motifs).
#'
#' @param jack a vector of two elements. The first one\code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]}
#'  specifies the number of subsets.
#'
#' @return Return a vector of same length as \code{fct},
#' of which values are computed as the arithmetic mean of all vector elements.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account
#'  species contribution). The elements belonging to a same motif are
#'  divided into \code{jack[2]} subsets of \code{jack[1]} elements.
#'  Prediction is computed by excluding \code{jack[1]} elements,
#'  of which the element to predict. If the total number of elements
#'   belonging to the motif is lower than \code{jack[1]*jack[2]},
#'    prediction is computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanBymotJack <- function(fct, assMotifs, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- ameanBymotJack(fct[indMot], jack)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack)
#' over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predictAmeanBymotJackXpr(fct, assMotifs, jack, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#'  (assembly motifs).
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'   the number of subsets.
#'
#' @param xpr a vector of labels of same length as \code{fct}
#' (assembly experiments).
#'
#' @return Return the arithmetic mean of a vector, as standard \code{mean}
#'  function.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'   is computed by excluding \code{jack[1]} elements, of which the
#'   element to predict. If the total number of elements belonging to
#'   the motif is lower than \code{jack[1]*jack[2]}, prediction is
#'   computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanBymotJackXpr <- function(fct, assMotifs, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictAmeanBymotJack(fct[indXpr],
                                            assMotifs[indXpr], jack)
  }

  return(fctPrd)
}


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Geometric mean by Motif                                                  ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot)
#'
#' @description Take the numeric vector \code{fct} and return a vector
#' of same length, of which values are computed as the geometric mean
#' of all vector elements belonging to a same motif. The motif of each
#' vector element is specified in the vector \code{assMotifs}.
#'
#' @usage predictGmeanBymot(fct, assMotifs)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)}
#' (assembly motifs).
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#' computed as the geometric mean of all vector elements belonging
#' to a same motif.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account
#' species contribution) by including all elements belonging
#' to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanBymot <- function(fct, assMotifs) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmean(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) over several experiments (xpr)
#'
#' @description Take the numeric vector \code{fct} and return a vector
#'  of same length, of which values are computed as the geometric mean
#'  of all vector elements belonging to a same motif, over several experiments.
#'   The motif of each vector element is specified in the vector
#'   \code{assMotifs}. The experiment of each vector element is
#'    specified in the vector \code{xpr}.
#'
#' @usage predictGmeanBymotXpr(fct, assMotifs, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#' (assembly motifs).
#'
#' @param xpr a vector of labels of same length as \code{fct}
#' (assembly experiments).
#'
#' @return Return a vector of same length as \code{fct}, of which values
#'  are computed as the geometric mean of all vector elements belonging
#'   to the same motif over several experiments.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#' contribution) by including all elements belonging to a same motif,
#' even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanBymotXpr <- function(fct, assMotifs, xpr) {

  fctPrd    <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictGmeanBymot(fct[indXpr], assMotifs[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) by leave-one-out (LOO)
#'
#' @description Take the numeric vector \code{fctMot} and return a vector
#'  of same length, of which values are computed as the geometric mean
#'  of all vector elements.
#'
#' @usage gmeanBymotLOO(fctMot)
#'
#' @param fctMot a vector of numeric values of elements belonging
#'  to a same motif.
#'
#' @return Return a vector of \code{length(fctMot)}, of which values
#' are computed as the geometric mean of all vector elements, excepted
#'  the value of element to predict that have been left out.
#'
#' @details Prediction computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account
#'  species contribution) by excluding the element to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmeanBymotLOO <- function(fctMot) {

  nbass <- length(fctMot)

  if (nbass > 1) {
    fctPrd <- numeric(nbass)
    for (ind in seq_len(nbass)) fctPrd[ind] <- gmean(fctMot[-ind])
  } else {
    fctPrd <- NA
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) by leave-one-out (LOO)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the geometric mean of all elements belonging to a same motif.
#'
#' @usage predictGmeanBymotLOO(fct, assMotifs)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @return Return a vector of \code{length(fctMot)}, of which values
#'  are computed as the geometric mean of all vector elements, excepted
#'   the value of element to predict that have been left out.
#'
#' @details Prediction computed using geometric mean \code{gmean}
#'  by motif \code{bymot} in a whole (WITHOUT taking into account
#'  species contribution) by excluding the element to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanBymotLOO <- function(fct, assMotifs) {

  fctPrd    <- numeric(length(assMotifs))

  setMot    <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmeanBymotLOO(fct[indMot])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) by leave-one-out (LOO)
#'  over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the gGeometric mean of all elements belonging to a same motif.
#'
#' @usage predictGmeanBymotLOOXpr(fct, assMotifs, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#' (assembly motifs).
#'
#' @param xpr a vector of labels of same length as \code{fct}
#' (assembly experiments).
#'
#' @return Return a vector of same length as \code{fct},
#' of which values are computed as the geometric mean of all
#' vector elements, over several experiments, excepted the value
#'  of element to predict that have been left out.
#'
#' @details Prediction computed using geometric mean \code{gmean}
#'  by motif \code{bymot} in a whole (WITHOUT taking into
#'  account species contribution) by excluding the element to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanBymotLOOXpr <- function(fct, assMotifs, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictGmeanBymotLOO(fct[indXpr], assMotifs[indXpr])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) by jackknife (jack)
#'  over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the geometric mean of all elements belonging to the same motif.
#'
#' @usage gmeanBymotJack(fctMot, jack)
#'
#' @param fctMot a vector of numeric values of elements belonging
#' to a same motif.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'  the number of subsets.
#'
#' @return Return a vector of same length as \code{fctMot}, of which
#'  values are computed as the geometric mean of all vector elements.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'   into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'   is computed by excluding \code{jack[1]} elements, of which
#'    the element to predict. If the total number of elements
#'    belonging to the motif is lower than \code{jack[1]*jack[2]},
#'     prediction is computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmeanBymotJack <- function(fctMot, jack) {

  nbass <- length(fctMot)

  if (nbass > jack[1] * jack[2]) {

    fctPrd <- numeric(nbass)
    index  <- sample.int(nbass)
    size   <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {
      indjack         <- index[(ind - 1) * size + (1:size)]
      fctPrd[indjack] <- gmean(fctMot[-indjack])
    }

    indjack         <- index[(ind * size + 1):nbass]
    fctPrd[indjack] <- gmean(fctMot[-indjack])

  } else {

    fctPrd <- gmeanBymotLOO(fctMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) by jackknife (jack)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the geometric mean of all elements belonging to the same motif.
#'
#' @usage predictGmeanBymotJack(fct, assMotifs, jack)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#' (assembly motifs).
#'
#' @param jack a vector of two elements. The first one\code{jack[1]}
#' specifies the size of subset, the second one \code{jack[2]} specifies
#'  the number of subsets.
#'
#' @return Return a vector of same length as \code{fct}, of which
#' values are computed as the geometric mean of all vector elements.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'   is computed by excluding \code{jack[1]} elements, of which the
#'   element to predict. If the total number of elements belonging
#'   to the motif is lower than \code{jack[1]*jack[2]}, prediction
#'    is computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanBymotJack <- function(fct, assMotifs, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmeanBymotJack(fct[indMot], jack)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by motif (bymot) by jackknife (jack)
#' over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the geometric mean of all elements belonging to a same motif.
#'
#' @usage predictAmeanBymotJackXpr(fct, assMotifs, jack, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#'  (assembly motifs).
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#' specifies the size of subset, the second one \code{jack[2]} specifies
#'  the number of subsets.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @return Return a vector of \code{length(fct)}, of which values
#' are computed as the geometric mean of all vector elements, over
#'  several experiments, excepted the value of element to predict
#'  that have been left out.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#' contribution). The elements belonging to a same motif are divided
#' into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#' is computed by excluding \code{jack[1]} elements, of which the
#' element to predict. If the total number of elements belonging to
#' the motif is lower than \code{jack[1]*jack[2]}, prediction is
#' computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanBymotJackXpr <- function(fct, assMotifs, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictGmeanBymotJack(assMotifs[indXpr],
                                               fct[indXpr], jack)
  }

  return(fctPrd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by Element within each Motif                             ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by elements occurring within the
#' assemblage and other assemblages, all assemblages belonging to
#'  a same motif (byelt)
#'
#' @description The numeric vector \code{fctMot} get together the
#' properties of assemblages belonging to a same assembly motif. The
#'  properties \code{fctMot} of assemblages containing a given
#'  element are separately averaged. The property of each assemblage
#'  is computed as the average of mean values of assemblages containing
#'   the same elements as the considered assemblage. The elemental
#'   composition of each assemblage is specified in the binary matrix
#'    \code{mOccurMot}: \code{0} if the element does not occur,
#'    \code{1} if the element occurs.
#'
#' @usage ameanByelt(fctMot, mOccurMot)
#'
#' @param fctMot a vector of numeric values (assembly properties).
#'
#' @param mOccurMot a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fctMot)}. Its second
#' dimension equals to the number of elements.
#'
#' @return Return a vector of \code{length(fctMot)}, of which values
#' are computed as an arithmetic mean.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#'  by element \code{byelt} occurring within the assemblage and other
#'  assemblages of a same motif, by including all assemblages belonging
#'   to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ameanByelt <- function(fctMot, mOccurMot) {

  fctPrd <- numeric(length(fctMot))

  setElt <- unique(which((mOccurMot[ , , drop = FALSE] == 1),
                         arr.ind = TRUE)[ , 2])

  for (elt in seq_along(setElt)) {
    indElt         <- which(mOccurMot[ , setElt[elt]] == 1)
    fctPrd[indElt] <- fctPrd[indElt] + amean(fctMot[indElt])
  }

  fctPrd <- fctPrd / apply(mOccurMot, MARGIN = 1, FUN = sum)

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by elements occurring within
#' the assemblage and other assemblages, all assemblages belonging
#' to a same motif (byelt)
#'
#' @description The numeric vector \code{fct} get together the
#' properties of assemblages belonging to different assembly motifs.
#'  The properties \code{fct} of assemblages belonging to a given
#'  assembly motif and containing a given element are separately averaged.
#'  The property of each assemblage is computed as the average of mean
#'  values of assemblages containing the same elements as the considered
#'  assemblage. The motif of each vector element is specified in
#'  the vector \code{assMotifs}. The elemental composition of
#'   each assemblage is specified in the binary matrix \code{mOccur}:
#'    \code{0} if the element does not occur, \code{1} if the element occurs.
#'
#' @usage predictAmeanByelt(fct, assMotifs, mOccur)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of same length as \code{fct}
#'  (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to length of \code{fct}. Its second
#'  dimension equals to the number of elements.
#'
#' @return Return a vector of same length as \code{fct}, of which
#'  values are computed as the arithmetic mean of all vector elements
#'   belonging to a same motif.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanByelt <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])

    if (length(indMot) > 1) {
      fctPrd[indMot] <- ameanByelt(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- fct[indMot]
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by elements occurring within the
#' assemblage and other assemblages, all assemblages belonging
#' to a same motif (byelt), over several experiments (xpr)
#'
#' @description The numeric vector \code{fct} get together the properties
#'  of assemblages belonging to different assembly motifs. The properties
#'  \code{fct} of assemblages belonging to a given assembly motif and
#'  containing a given element are separately averaged. The property
#'  of each assemblage is computed as the average of mean values
#'  of assemblages containing the same elements as the considered assemblage.
#'   The motif of each vector element is specified in the vector
#'    \code{assMotifs}. The elemental composition of each assemblage
#'     is specified in the binary matrix \code{mOccur}: \code{0}
#'      if the element does not occur, \code{1} if the element occurs.
#'
#' @usage predictAmeanByeltXpr(fct, assMotifs, mOccur, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)}
#'  (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fct)}. Its second
#'  dimension equals to the number of elements.
#'
#' @param xpr a vector of labels of same length as \code{fct}
#'  (assembly experiments).
#'
#' @return Return a vector of same length as \code{fct}, of which
#' values are computed as the arithmetic mean of all vector elements
#'  belonging to a same motif.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanByeltXpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictAmeanByelt(fct[indXpr],
                                          assMotifs[indXpr],
                                          mOccur[indXpr, ])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by elements occurring within
#' the assemblage and other assemblages belonging to a same motif (byelt),
#' except the assemblage to predict (leave-one-out, LOO)
#'
#' @description The numeric vector \code{fctMot} get together the
#'  properties of assemblages belonging to a same assembly motif.
#'  The properties \code{fctMot} of assemblages containing a given element
#'   are separately averaged, except the assemblage to predict. The property
#'    of each assemblage is computed as the average of mean values of
#'    assemblages containing the same elements as the considered assemblage.
#'     The elemental composition of each assemblage is specified in the
#'      binary matrix \code{mOccurMot}: \code{0} if the element does not
#'       occur, \code{1} if the element occurs.
#'
#' @usage ameanByeltLOO(fctMot, mOccurMot)
#'
#' @param fctMot a vector of numeric values (assembly properties).
#'
#' @param mOccurMot a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fctMot)}. Its second
#'  dimension equals to the number of elements.
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#'  computed as an arithmetic mean of all vector elements over several
#'   experiments, excepted the value of element to predict that have
#'   been left out.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by element \code{byelt} occurring within the assemblage and other
#'  assemblages of a same motif, by including all assemblages belonging
#'  to a same motif, except the assemblage to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ameanByeltLOO <- function(fctMot, mOccurMot) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)
  vfct   <- numeric(dim(mOccurMot)[2])

  for (ind in seq_len(nbass)) {

    vfct[] <- NA
    indOth <- seq_len(nbass)[-ind]
    setElt <- which(mOccurMot[ind, ] != 0)

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- amean(fctMot[indOth[indElt]])

      #      if (length(indElt) > 0) {
      #        vfct[setElt[elt]] <- amean(fctMot[indOth[indElt]])
      #    } else {
      #        vfct[setElt[elt]] <- amean(fctMot[indOth])
      #    }

    }

    fctPrd[ind] <- amean(vfct[setElt])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by elements occurring within the
#' assemblage and other assemblages belonging to a same motif (byelt),
#'  except the assemblage to predict (leave-one-out, LOO)
#'
#' @description The numeric vector \code{fct} get together the properties
#' of assemblages belonging to different assembly motifs. The properties
#' \code{fct} of assemblages belonging to a given assembly motif and
#' containing a given element are separately averaged. The property of
#'  each assemblage is computed as the average of mean values of assemblages
#'   containing the same elements as the considered assemblage, except the
#'    assemblage to predict (leave-one-out, LOO). The motif of each vector
#'    element is specified in the vector \code{assMotifs}. The elemental
#'     composition of each assemblage is specified in the binary matrix
#'     \code{mOccur}: \code{0} if the element does not occur, \code{1}
#'      if the element occurs.
#'
#' @usage predictAmeanByeltLOO(fct, assMotifs, mOccur)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fct)}. Its second
#'  dimension equals to the number of elements.
#'
#' @return Return a vector of same length as \code{fctMot}, of which
#'  values are computed as an arithmetic mean.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, except the assemblage to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanByeltLOO <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {
    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- ameanByeltLOO(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by elements occurring within
#' the assemblage and other assemblages belonging to a same motif
#'  (byelt), except the assemblage to predict (leave-one-out, LOO).
#'   Over several experiments.
#'
#' @description The numeric vector \code{fct} get together the properties
#' of assemblages belonging to different assembly motifs. The properties
#'  \code{fct} of assemblages belonging to a given assembly motif and
#'  containing a given element are separately averaged. The property of
#'  each assemblage is computed as the average of mean values of
#'  assemblages containing the same elements as the considered assemblage,
#'   except the assemblage to predict (leave-one-out, LOO). The motif
#'   of each vector element is specified in the vector \code{assMotifs}.
#'    The elemental composition of each assemblage is specified in the
#'    binary matrix \code{mOccur}: \code{0} if the element does not occur,
#'     \code{1} if the element occurs.
#'
#' @usage predictAmeanByeltLOOXpr(fct, assMotifs, mOccur, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of elements.
#'
#' @param xpr a vector of labels of same length as \code{fct}
#' (assembly experiments).
#'
#' @return Return a vector of same length as \code{fctMot}, of
#' which values are computed as an arithmetic mean.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, except the assemblage to predict. Over several
#'     experiments.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanByeltLOOXpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predictAmeanByeltLOO(fct[indXpr],
                                              assMotifs[indXpr],
                                              mOccur[indXpr, ] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack)
#'  over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage ameanByeltJack(fctMot, mOccurMot, jack)
#'
#' @param fctMot a vector of numeric values of elements belonging to
#' a same motif.
#'
#' @param mOccurMot a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fctMot)}. Its second
#' dimension equals to the number of elements.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#' specifies the size of subset, the second one \code{jack[2]} specifies
#' the number of subsets.
#'
#' @return Return a vector of \code{length(fctMot)}, of which values are
#' computed as the arithmetic mean of all vector elements.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#' contribution). The elements belonging to a same motif are divided
#' into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#' is computed by excluding \code{jack[1]} elements, of which the element
#'  to predict. If the total number of elements belonging to the motif
#'   is lower than \code{jack[1]*jack[2]}, prediction is computed by
#'   Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

ameanByeltJack <- function(fctMot, mOccurMot, jack) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {

      indjack <- index[(ind - 1) * size + (1:size)]
      indOth  <- seq_len(nbass)[-indjack]

      tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
      tmp[tmp == 0]   <- NA
      vfct            <- apply(tmp, MARGIN = 2, FUN = amean)

      tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
      tmp[tmp == 0]   <- NA
      fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = amean)
    }

    indjack <- index[(ind * size + 1):nbass]
    indOth  <- seq_len(nbass)[-indjack]

    tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
    tmp[tmp == 0]   <- NA
    vfct            <- apply(tmp, MARGIN = 2, FUN = amean)

    tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
    tmp[tmp == 0]   <- NA
    fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = amean)

  } else {

    fctPrd[ ] <- ameanByeltLOO(fctMot, mOccurMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage predictAmeanByeltJack(fct, assMotifs, mOccur, jack)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second
#' dimension equals to the number of elements.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'   the number of subsets.
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#'  computed as the arithmetic mean of all vector elements.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'   into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'    is computed by excluding \code{jack[1]} elements, of which the
#'    element to predict. If the total number of elements belonging to
#'     the motif is lower than \code{jack[1]*jack[2]}, prediction is
#'     computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanByeltJack <- function(fct, assMotifs, mOccur, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- ameanByeltJack(fct[indMot], mOccur[indMot, ], jack)
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack)
#' over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predictAmeanByeltJackXpr(fct, assMotifs, mOccur, jack, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of elements.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'   the number of subsets.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @return Return the arithmetic mean of a vector, as standard \code{mean}
#'  function.
#'
#' @details Prediction is computed using arithmetic mean \code{amean} by
#'  motif \code{bymot} in a whole (WITHOUT taking into account species
#'   contribution). The elements belonging to a same motif are divided
#'    into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'     is computed by excluding \code{jack[1]} elements, of which the
#'     element to predict. If the total number of elements belonging to
#'      the motif is lower than \code{jack[1]*jack[2]}, prediction is
#'       computed by Leave-One-Out (LOO).
#'
#" @keywords internal
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictAmeanByeltJackXpr <- function(fct, assMotifs, mOccur, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predictAmeanByeltJack(fct[index],
                                              assMotifs[index],
                                              mOccur[index, ],
                                              jack  )
  }

  return(fctPrd)
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Geometric mean by Element within each Motif                             ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by elements occurring within the assemblage
#' and other assemblages, all assemblages belonging to a same motif (byelt)
#'
#' @description The numeric vector \code{fctMot} get together the properties
#'  of assemblages belonging to a same assembly motif. The properties
#'  \code{fctMot} of assemblages containing a given element are separately
#'   averaged. The property of each assemblage is computed as the average
#'    of mean values of assemblages containing the same elements as the
#'     considered assemblage. The elemental composition of each assemblage
#'      is specified in the binary matrix \code{mOccurMot}: \code{0}
#'       if the element does not occur, \code{1} if the element occurs.
#'
#' @usage gmeanByelt(fctMot, mOccurMot)
#'
#' @param fctMot a vector of numeric values (assembly properties).
#'
#' @param mOccurMot a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fctMot)}. Its second
#'   dimension equals to the number of elements.
#'
#' @return Return a vector of \code{length(fctMot)}, of which values
#'  are computed as an geometric mean.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'   to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmeanByelt <- function(fctMot, mOccurMot) {

  fctPrd    <- numeric(length(fctMot))
  fctPrd[ ] <- 1

  setElt    <- unique(which((mOccurMot[ , , drop = FALSE] == 1),
                            arr.ind = TRUE)[ , 2])
  for (elt in seq_along(setElt)) {
    indElt <- which(mOccurMot[ , setElt[elt]] == 1)

    if (length(indElt) > 0)
      fctPrd[indElt] <- fctPrd[indElt] * gmean(fctMot[indElt])
  }

  fctPrd <- fctPrd ^ (1/apply(mOccurMot, MARGIN = 1, FUN = sum))

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by elements occurring within the assemblage
#'  and other assemblages, all assemblages belonging to a same motif (byelt)
#'
#' @description The numeric vector \code{fct} get together the properties
#' of assemblages belonging to different assembly motifs. The properties
#' \code{fct} of assemblages belonging to a given assembly motif and
#' containing a given element are separately averaged. The property of
#'  each assemblage is computed as the average of mean values of assemblages
#'   containing the same elements as the considered assemblage. The
#'   motif of each vector element is specified in the vector \code{assMotifs}.
#'    The elemental composition of each assemblage is specified in the
#'     binary matrix \code{mOccur}: \code{0} if the element does not occur,
#'      \code{1} if the element occurs.
#'
#' @usage predictGmeanByelt(fct, assMotifs, mOccur)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fct)}. Its second dimension
#'   equals to the number of elements.
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#'  computed as the geometric mean of all vector elements belonging to
#'   a same motif.
#'
#' @details Prediction is computed using geometric mean \code{gmean} by
#'  element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanByelt <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmeanByelt(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- fct[indMot]
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (gmean) by elements occurring within the assemblage
#'  and other assemblages, all assemblages belonging to a same motif
#'  (byelt), over several experiments (xpr)
#'
#' @description The numeric vector \code{fct} get together the properties
#' of assemblages belonging to different assembly motifs. The properties
#'  \code{fct} of assemblages belonging to a given assembly motif and
#'   containing a given element are separately averaged. The property of
#'    each assemblage is computed as the average of mean values of assemblages
#'     containing the same elements as the considered assemblage.
#'     The motif of each vector element is specified in the vector
#'     \code{assMotifs}. The elemental composition of each assemblage
#'      is specified in the binary matrix \code{mOccur}: \code{0} if
#'      the element does not occur, \code{1} if the element occurs.
#'
#' @usage predictGmeanByeltXpr(fct, assMotifs, mOccur, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second
#' dimension equals to the number of elements.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @return Return a vector of \code{length(fct)}, of which values
#' are computed as the geometric mean of all vector elements belonging
#' to a same motif.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, even the one to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanByeltXpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predictGmeanByelt( fct[index],
                                        assMotifs[index],
                                         mOccur[index, ] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (amean) by elements occurring within the assemblage
#'  and other assemblages belonging to a same motif (byelt),
#'  except the assemblage to predict (leave-one-out, LOO)
#'
#' @description The numeric vector \code{fctMot} get together the
#' properties of assemblages belonging to a same assembly motif.
#' The properties \code{fctMot} of assemblages containing a given element
#'  are separately averaged, except the assemblage to predict. The property
#'   of each assemblage is computed as the average of mean values of
#'   assemblages containing the same elements as the considered assemblage.
#'    The elemental composition of each assemblage is specified in the
#'     binary matrix \code{mOccurMot}: \code{0} if the element does not
#'      occur, \code{1} if the element occurs.
#'
#' @usage gmeanByeltLOO(fctMot, mOccurMot)
#'
#' @param fctMot a vector of numeric values (assembly properties).
#'
#' @param mOccurMot a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fctMot)}. Its second
#'  dimension equals to the number of elements.
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#'  computed as a geometric mean of all vector elements over several
#'  experiments, excepted the value of element to predict that have
#'  been left out.
#'
#' @details Prediction is computed using geomùetric mean \code{gmean}
#' by element \code{byelt} occurring within the assemblage and other
#' assemblages of a same motif, by including all assemblages belonging
#'  to a same motif, except the assemblage to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmeanByeltLOO <- function(fctMot, mOccurMot) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)
  vfct   <- numeric(dim(mOccurMot)[2])

  for (ind in seq_len(nbass)) {

    vfct[] <- NA
    indOth <- seq_len(nbass)[-ind]
    setElt <- which(mOccurMot[ind, ] != 0)

    for (elt in seq_along(setElt)) {
      indElt <- which(mOccurMot[indOth, setElt[elt]] == 1)
      if (length(indElt) > 0)
        vfct[setElt[elt]] <- gmean(fctMot[indOth[indElt]])
    }

    fctPrd[ind] <- gmean(vfct[setElt])
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (amean) by elements occurring within the assemblage
#' and other assemblages belonging to a same motif (byelt), except the
#'  assemblage to predict (leave-one-out, LOO)
#'
#' @description The numeric vector \code{fct} get together the properties
#'  of assemblages belonging to different assembly motifs. The properties
#'   \code{fct} of assemblages belonging to a given assembly motif and
#'   containing a given element are separately averaged. The property of
#'   each assemblage is computed as the average of mean values of
#'   assemblages containing the same elements as the considered assemblage,
#'    except the assemblage to predict (leave-one-out, LOO). The motif
#'     of each vector element is specified in the vector \code{assMotifs}.
#'      The elemental composition of each assemblage is specified in the
#'      binary matrix \code{mOccur}: \code{0} if the element does not occur,
#'       \code{1} if the element occurs.
#'
#' @usage predictGmeanByeltLOO(fct, assMotifs, mOccur)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements). Its first
#'  dimension equals to \code{length(fct)}. Its second dimension equals
#'  to the number of elements.
#'
#' @return Return a vector of \code{length(fctMot)}, of which values are
#' computed as a geometric mean.
#'
#' @details Prediction is computed using geometric mean \code{gmean} by
#'  element \code{byelt} occurring within the assemblage and other
#'  assemblages of a same motif, by including all assemblages belonging
#'  to a same motif, except the assemblage to predict.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanByeltLOO <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmeanByeltLOO(fct[indMot], mOccur[indMot, ])
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (amean) by elements occurring within the assemblage
#'  and other assemblages belonging to a same motif (byelt),
#'   except the assemblage to predict (leave-one-out, LOO).
#'   Over several experiments.
#'
#' @description The numeric vector \code{fct} get together the properties
#'  of assemblages belonging to different assembly motifs. The properties
#'   \code{fct} of assemblages belonging to a given assembly motif and
#'   containing a given element are separately averaged. The property of
#'    each assemblage is computed as the average of mean values of assemblages
#'    containing the same elements as the considered assemblage,
#'    except the assemblage to predict (leave-one-out, LOO). The motif
#'     of each vector element is specified in the vector \code{assMotifs}.
#'      The elemental composition of each assemblage is specified in the
#'       binary matrix \code{mOccur}: \code{0} if the element does not
#'        occur, \code{1} if the element occurs.
#'
#' @usage predictGmeanByeltLOOXpr(fct, assMotifs, mOccur, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#' equals to the number of elements.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @return Return a vector of \code{length(fctMot)}, of which
#' values are computed as a geometric mean.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#'  by element \code{byelt} occurring within the assemblage and other
#'   assemblages of a same motif, by including all assemblages belonging
#'    to a same motif, except the assemblage to predict. Over several
#'     experiments.
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanByeltLOOXpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd    <- numeric(length(assMotifs))

  setXpr    <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predictGmeanByeltLOO(fct[index],
                                          assMotifs[index],
                                          mOccur[index, ] )
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (amean) by motif (bymot) by jackknife (jack)
#'  over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage gmeanByeltJack(fctMot, mOccurMot, jack)
#'
#' @param fctMot a vector of numeric values of elements belonging
#'  to a same motif.
#'
#' @param mOccurMot a matrix of occurrence (occurrence of elements).
#'  Its first dimension equals to \code{length(fctMot)}. Its second
#'   dimension equals to the number of elements.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'  the number of subsets.
#'
#' @return Return a vector of same length as \code{fctMot}, of which values
#'  are computed as the geometric mean of all vector elements.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#'  by motif \code{bymot} in a whole (WITHOUT taking into account species
#'   contribution). The elements belonging to a same motif are divided
#'    into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'     is computed by excluding \code{jack[1]} elements, of which the
#'      element to predict. If the total number of elements belonging
#'      to the motif is lower than \code{jack[1]*jack[2]}, prediction
#'      is computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

gmeanByeltJack <- function(fctMot, mOccurMot, jack) {

  nbass  <- length(fctMot)
  fctPrd <- numeric(nbass)

  if (nbass > jack[1] * jack[2]) {

    index <- sample.int(nbass)
    size  <- floor(nbass / jack[2])

    for (ind in seq_len(jack[2] - 1)) {

      indjack <- index[(ind - 1) * size + (1:size)]
      indOth  <- seq_len(nbass)[-indjack]

      tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
      tmp[tmp == 0]   <- NA
      vfct            <- apply(tmp, MARGIN = 2, FUN = gmean)

      tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
      tmp[tmp == 0]   <- NA
      fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = gmean)
    }

    indjack <- index[(ind * size + 1):nbass]
    indOth  <- seq_len(nbass)[-indjack]

    tmp             <- mOccurMot[indOth, ] * fctMot[indOth]
    tmp[tmp == 0]   <- NA
    vfct            <- apply(tmp, MARGIN = 2, FUN = gmean)

    tmp             <- t(t(mOccurMot[indjack, ]) * vfct)
    tmp[tmp == 0]   <- NA
    fctPrd[indjack] <- apply(tmp, MARGIN = 1, FUN = gmean)

  } else {

    fctPrd[] <- gmeanByeltLOO(fctMot, mOccurMot)
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (amean) by motif (bymot) by jackknife (jack)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the arithmetic mean of all elements belonging to the same motif.
#'
#' @usage predictGmeanByeltJack(fct, assMotifs, mOccur, jack)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second
#' dimension equals to the number of elements.
#'
#' @param jack a vector of two elements. The first one\code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'   the number of subsets.
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#' computed as the geometric mean of all vector elements, excepted the
#'  subset left out.
#'
#' @details Prediction is computed using geometric mean \code{gmean} by
#'  motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'   is computed by excluding \code{jack[1]} elements, of which the
#'   element to predict. If the total number of elements belonging to
#'    the motif is lower than \code{jack[1]*jack[2]}, prediction is
#'    computed by Leave-One-Out (LOO).
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanByeltJack <- function(fct, assMotifs, mOccur, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmeanByeltJack(fct[indMot], mOccur[indMot, ], jack)
    } else {
      fctPrd[indMot] <- NA
    }
  }

  return(fctPrd)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Geometric mean (amean) by motif (bymot) by jackknife (jack)
#' over several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector
#'  computed as the geometric mean of all elements belonging to a same motif.
#'
#' @usage predictGmeanByeltJackXpr(fct, assMotifs, mOccur, jack, xpr)
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#' equals to the number of elements.
#'
#' @param jack a vector of two elements. The first one \code{jack[1]}
#'  specifies the size of subset, the second one \code{jack[2]} specifies
#'   the number of subsets.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @return Return a vector of \code{length(fct)}, of which values are
#'  computed as the geometric mean of all vector elements,
#'  over several experiments.
#'
#' @details Prediction is computed using geometric mean \code{gmean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'   into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'    is computed by excluding \code{jack[1]} elements, of which
#'    the element to predict. If the total number of elements belonging
#'     to the motif is lower than \code{jack[1]*jack[2]},
#'     prediction is computed by Leave-One-Out (LOO).
#'
#" @keywords internal
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictGmeanByeltJackXpr <- function(fct, assMotifs, mOccur, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predictGmeanByeltJack(fct[index],
                                           assMotifs[index],
                                           mOccur[index, ],
                                           jack  )
  }

  return(fctPrd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Public Functions                                                         ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack) over
#'  several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector computed
#'  as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predictCal(fct, assMotifs, mOccur, xpr,
#'                   opt.mean = "amean",
#'                   opt.mod  = "bymot"  )
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of elements.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @param opt.mean equal to \code{"amean"} (by default) or \code{"gmean"}.
#'
#' @param opt.mod  equal to \code{"bymot"} (by default) or \code{"byelt"}.
#'
#' @return Return the arithmetic mean of a vector, as standard \code{mean}
#'  function.
#'
#' @details Prediction is computed using arithmetic mean \code{amean} by motif
#'  \code{bymot} in a whole (WITHOUT taking into account species contribution).
#'   The elements belonging to a same motif are divided into \code{jack[2]}
#'    subsets of \code{jack[1]} elements. Prediction is computed by excluding
#'     \code{jack[1]} elements, of which the element to predict. If the total
#'      number of elements belonging to the motif is lower than
#'       \code{jack[1]*jack[2]}, prediction is computed by Leave-One-Out (LOO).
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictCal <- function(fct, assMotifs, mOccur, xpr,
                       opt.mean = "amean",
                       opt.mod  = "bymot") {

  optmean <- "amean"
  if (opt.mean == "gmean") optmean <- "gmean"

  optmod <- "bymot"
  if (opt.mod == "byelt") optmod <- "byelt"

  optxpr <- ""
  if (length(unique(xpr)) != 1) optxpr <- "xpr"

  option <- paste(optmean, optmod, optxpr, sep = ".")

  return(
    switch(option,
           amean.bymot. =
             predictAmeanBymot(fct, assMotifs),
           amean.bymot.xpr =
             predictAmeanBymotXpr(fct, assMotifs, xpr),
           gmean.bymot. =
             predictGmeanBymot(fct, assMotifs),
           gmean.bymot.xpr =
             predictGmeanBymotXpr(fct, assMotifs, xpr),

           amean.byelt. =
             predictAmeanByelt(fct, assMotifs, mOccur),
           amean.byelt.xpr =
             predictAmeanByeltXpr(fct, assMotifs, mOccur, xpr),
           gmean.byelt. =
             predictGmeanByelt(fct, assMotifs, mOccur),
           gmean.byelt.xpr =
             predictGmeanByeltXpr(fct, assMotifs, mOccur, xpr)  )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predition of assembly property by cross-validation
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predictPrd(fct, assMotifs, mOccur, xpr,
#'                   opt.mean = "amean",
#'                   opt.mod  = "bymot",
#'                   opt.jack = FALSE, jack = c(2, 5)  )
#'
#' @param fct a vector of numeric values (assembly properties).
#'
#' @param assMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param mOccur a matrix of occurrence (occurrence of elements). Its first
#' dimension equals to \code{length(fct)}. Its second dimension equals to
#' the number of elements.
#'
#' @param xpr a vector of labels of \code{length(fct)} (assembly experiments).
#'
#' @param opt.mean equal to \code{"amean"} (by default) or to \code{"gmean"}.
#'
#' @param opt.mod equal to \code{"bymot"} (by default) or to \code{"byelt"}.
#'
#' @param opt.jack a logical to block or switch to jackknife cross-validation.
#'
#' @param jack \code{jack = c(2, 5)} by default.

#' @return Return the arithmetic mean of a vector, as standard \code{mean}
#' function.
#'
#' @details Prediction is computed using arithmetic mean \code{amean}
#' by motif \code{bymot} in a whole (WITHOUT taking into account species
#'  contribution). The elements belonging to a same motif are divided
#'  into \code{jack[2]} subsets of \code{jack[1]} elements. Prediction
#'  is computed by excluding \code{jack[1]} elements, of which the
#'  element to predict. If the total number of elements belonging to
#'  the motif is lower than \code{jack[1]*jack[2]}, prediction is
#'  computed by Leave-One-Out (LOO).
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predictPrd <- function(fct, assMotifs, mOccur, xpr,
                       opt.mean = "amean",
                       opt.mod  = "bymot",
                       opt.jack = FALSE,   jack = c(2, 5)) {

  optmean <- "amean"
  if (opt.mean == "gmean") optmean <- "gmean"

  optmod <- "bymot"
  if (opt.mod == "byelt") optmod <- "byelt"

  optjack <- ""
  if (opt.jack == TRUE) optjack <- "jack"

  optxpr <- ""
  if (length(unique(xpr)) != 1) optxpr <- "xpr"

  option <- paste(optmean, optmod, optjack, optxpr, sep = ".")

  return(
    switch(option,
           amean.bymot.. =
             predictAmeanBymotLOO(fct, assMotifs),
           amean.bymot..xpr =
             predictAmeanBymotLOOXpr(fct, assMotifs, xpr),
           amean.bymot.jack. =
             predictAmeanBymotJack(fct, assMotifs, jack),
           amean.bymot.jack.xpr =
             predictAmeanBymotJackXpr(fct, assMotifs, jack, xpr),

           gmean.bymot.. =
             predictGmeanBymotLOO(fct, assMotifs),
           gmean.bymot..xpr =
             predictGmeanBymotLOOXpr(fct, assMotifs, xpr),
           gmean.bymot.jack. =
             predictGmeanBymotJack(fct, assMotifs, jack),
           gmean.bymot.jack.xpr =
             predictGmeanBymotJackXpr(fct, assMotifs, jack, xpr),

           amean.byelt.. =
             predictAmeanByeltLOO(fct, assMotifs, mOccur),
           amean.byelt..xpr =
             predictAmeanByeltLOOXpr(fct, assMotifs, mOccur, xpr),
           amean.byelt.jack. =
             predictAmeanByeltJack(fct, assMotifs, mOccur, jack),
           amean.byelt.jack.xpr =
             predictAmeanByeltJackXpr(fct, assMotifs, mOccur, jack, xpr),

           gmean.byelt.. =
             predictGmeanByeltLOO(fct, assMotifs, mOccur),
           gmean.byelt..xpr =
             predictGmeanByeltLOOXpr(fct, assMotifs, mOccur, xpr),
           gmean.byelt.jack. =
             predictGmeanByeltJack(fct, assMotifs, mOccur, jack),
           gmean.byelt.jack.xpr =
             predictGmeanByeltJackXpr(fct, assMotifs, mOccur, jack, xpr) )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


