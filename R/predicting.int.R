#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Set of internal functions                                                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

figures  <- c(21, 22, 24, 23, 25)
# 21 = "circle", 22 = "square", 23 = "diamond",
# 24 = triangle point-up", 25 = "triangle point-down"

couleurs <- c("red3", "blue2", "orange2", "turquoise3", "magenta", "green4",
              "pink", "violet", "salmon4", "skyblue2", "sienna3", "olivedrab3")

myLetters <- c(letters, LETTERS)


#' @include stats.int.R
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Arithmetic mean by Motif                                                 ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by assembly motif (bymot)
#'
#' @description Take the numeric vector \code{fct} and return a vector of same
#'  length, of which values are computed as the arithmetic mean of all vector
#'   elements belonging to a same motif. The motif of each vector element is
#'    specified in the vector \code{assMotifs}.
#'
#' @usage predict_amean_bymot(fct, assMotifs)
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

predict_amean_bymot <- function(fct, assMotifs) {

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
#' @usage predict_amean_bymot_xpr(fct, assMotifs, xpr)
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

predict_amean_bymot_xpr <- function(fct, assMotifs, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_amean_bymot(fct[indXpr], assMotifs[indXpr])
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
#' @usage amean_bymot_LOO(fctMot)
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

amean_bymot_LOO <- function(fctMot) {

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
#' @usage predict_amean_bymot_LOO(fct, assMotifs)
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

predict_amean_bymot_LOO <- function(fct, assMotifs) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- amean_bymot_LOO(fct[indMot])
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
#' @usage predict_amean_bymot_LOO_xpr(fct, assMotifs, xpr)
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

predict_amean_bymot_LOO_xpr <- function(fct, assMotifs, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_amean_bymot_LOO(fct[indXpr], assMotifs[indXpr])
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
#' @usage amean_bymot_jack(fctMot, jack)
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

amean_bymot_jack <- function(fctMot, jack) {

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

    fctPrd <- amean_bymot_LOO(fctMot)
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
#' @usage predict_amean_bymot_jack(fct, assMotifs, jack)
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

predict_amean_bymot_jack <- function(fct, assMotifs, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- amean_bymot_jack(fct[indMot], jack)
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
#' @usage predict_amean_bymot_jack_xpr(fct, assMotifs, jack, xpr)
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

predict_amean_bymot_jack_xpr <- function(fct, assMotifs, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_amean_bymot_jack(fct[indXpr],
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
#' @usage predict_gmean_bymot(fct, assMotifs)
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

predict_gmean_bymot <- function(fct, assMotifs) {

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
#' @usage predict_gmean_bymot_xpr(fct, assMotifs, xpr)
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

predict_gmean_bymot_xpr <- function(fct, assMotifs, xpr) {

  fctPrd    <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_gmean_bymot(fct[indXpr], assMotifs[indXpr])
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
#' @usage gmean_bymot_LOO(fctMot)
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

gmean_bymot_LOO <- function(fctMot) {

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
#' @usage predict_gmean_bymot_LOO(fct, assMotifs)
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

predict_gmean_bymot_LOO <- function(fct, assMotifs) {

  fctPrd    <- numeric(length(assMotifs))

  setMot    <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmean_bymot_LOO(fct[indMot])
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
#' @usage predict_gmean_bymot_LOO_xpr(fct, assMotifs, xpr)
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

predict_gmean_bymot_LOO_xpr <- function(fct, assMotifs, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_gmean_bymot_LOO(fct[indXpr], assMotifs[indXpr])
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
#' @usage gmean_bymot_jack(fctMot, jack)
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

gmean_bymot_jack <- function(fctMot, jack) {

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

    fctPrd <- gmean_bymot_LOO(fctMot)
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
#' @usage predict_gmean_bymot_jack(fct, assMotifs, jack)
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

predict_gmean_bymot_jack <- function(fct, assMotifs, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot         <- which(assMotifs == setMot[mot])
    fctPrd[indMot] <- gmean_bymot_jack(fct[indMot], jack)
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
#' @usage predict_amean_bymot_jack_xpr(fct, assMotifs, jack, xpr)
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

predict_gmean_bymot_jack_xpr <- function(fct, assMotifs, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_gmean_bymot_jack(assMotifs[indXpr],
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
#' @usage amean_byelt(fctMot, mOccurMot)
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

amean_byelt <- function(fctMot, mOccurMot) {

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
#' @name predict_amean_byelt
#'
#' @usage predict_amean_byelt(fct, assMotifs, mOccur)
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
# @keywords internal
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_amean_byelt <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])

    if (length(indMot) > 1) {
      fctPrd[indMot] <- amean_byelt(fct[indMot], mOccur[indMot, ])
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
#' @usage predict_amean_byelt_xpr(fct, assMotifs, mOccur, xpr)
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

predict_amean_byelt_xpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_amean_byelt(fct[indXpr],
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
#' @usage amean_byelt_LOO(fctMot, mOccurMot)
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

amean_byelt_LOO <- function(fctMot, mOccurMot) {

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
#' @usage predict_amean_byelt_LOO(fct, assMotifs, mOccur)
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

predict_amean_byelt_LOO <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {
    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- amean_byelt_LOO(fct[indMot], mOccur[indMot, ])
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
#' @usage predict_amean_byelt_LOO_xpr(fct, assMotifs, mOccur, xpr)
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

predict_amean_byelt_LOO_xpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    indXpr         <- which(xpr == setXpr[ix])
    fctPrd[indXpr] <- predict_amean_byelt_LOO(fct[indXpr],
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
#' @usage amean_byelt_jack(fctMot, mOccurMot, jack)
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

amean_byelt_jack <- function(fctMot, mOccurMot, jack) {

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

    fctPrd[ ] <- amean_byelt_LOO(fctMot, mOccurMot)
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
#' @usage predict_amean_byelt_jack(fct, assMotifs, mOccur, jack)
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

predict_amean_byelt_jack <- function(fct, assMotifs, mOccur, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- amean_byelt_jack(fct[indMot], mOccur[indMot, ], jack)
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
#' @usage predict_amean_byelt_jack_xpr(fct, assMotifs, mOccur, jack, xpr)
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
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_amean_byelt_jack_xpr <- function(fct, assMotifs, mOccur, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict_amean_byelt_jack(fct[index],
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
#' @usage gmean_byelt(fctMot, mOccurMot)
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

gmean_byelt <- function(fctMot, mOccurMot) {

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
#' @usage predict_gmean_byelt(fct, assMotifs, mOccur)
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

predict_gmean_byelt <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmean_byelt(fct[indMot], mOccur[indMot, ])
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
#' @usage predict_gmean_byelt_xpr(fct, assMotifs, mOccur, xpr)
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

predict_gmean_byelt_xpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict_gmean_byelt( fct[index],
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
#' @usage gmean_byelt_LOO(fctMot, mOccurMot)
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

gmean_byelt_LOO <- function(fctMot, mOccurMot) {

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
#' @usage predict_gmean_byelt_LOO(fct, assMotifs, mOccur)
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

predict_gmean_byelt_LOO <- function(fct, assMotifs, mOccur) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmean_byelt_LOO(fct[indMot], mOccur[indMot, ])
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
#' @usage predict_gmean_byelt_LOO_xpr(fct, assMotifs, mOccur, xpr)
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

predict_gmean_byelt_LOO_xpr <- function(fct, assMotifs, mOccur, xpr) {

  fctPrd    <- numeric(length(assMotifs))

  setXpr    <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict_gmean_byelt_LOO(fct[index],
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
#' @usage gmean_byelt_jack(fctMot, mOccurMot, jack)
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

gmean_byelt_jack <- function(fctMot, mOccurMot, jack) {

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

    fctPrd[] <- gmean_byelt_LOO(fctMot, mOccurMot)
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
#' @usage predict_gmean_byelt_jack(fct, assMotifs, mOccur, jack)
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

predict_gmean_byelt_jack <- function(fct, assMotifs, mOccur, jack) {

  fctPrd <- numeric(length(assMotifs))

  setMot <- unique(assMotifs)
  for (mot in seq_along(setMot)) {

    indMot <- which(assMotifs == setMot[mot])
    if (length(indMot) > 1) {
      fctPrd[indMot] <- gmean_byelt_jack(fct[indMot], mOccur[indMot, ], jack)
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
#' @usage predict_gmean_byelt_jack_xpr(fct, assMotifs, mOccur, jack, xpr)
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
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_gmean_byelt_jack_xpr <- function(fct, assMotifs, mOccur, jack, xpr) {

  fctPrd <- numeric(length(assMotifs))

  setXpr <- unique(xpr)
  for (ix in seq_along(setXpr)) {

    index         <- which(xpr == setXpr[ix])
    fctPrd[index] <- predict_gmean_byelt_jack(fct[index],
                                              assMotifs[index],
                                              mOccur[index, ],
                                              jack  )
  }

  return(fctPrd)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Supplementary assemblages to predict                                     ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      amean = by using arithmetic mean
#      bymot = by motif in a whole (WITHOUT taking into account
#                                                       species contribution)
#      by including all assemblies, even the one to predict
#'
#' @title Prediction of supplementary assemblages
#'
#' @description Take a numeric f.
#'
#' @usage predict_amean_bymot_supp(appFct, appMotifs, supMotifs)
#'
#' @param appFct a vector of numeric values (assembly properties).
#'
#' @param appMotifs a vector of labels of \code{length(appFct)}
#'  (assembly motifs).
#'
#' @param supMotifs a vector of labels of assembly motifs of which values must be predicted.
#'
#' @return Return a vector of \code{length(supMotifs)}. The values are
#'  computed using arithmetic mean of elements belonging to \code{appMotifs}
#'  and sharing a same motif.
#'
#' @details Prediction ...
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_amean_bymot_supp <- function(appFct, appMotifs, supMotifs) {

  supFct   <- numeric(length(supMotifs))
  supFct[] <- NA

  setMot <- unique(supMotifs)
  for (mot in seq_along(setMot)) {

    indSup <- which(supMotifs == setMot[mot])
    indApp <- which(appMotifs == setMot[mot])
    if ( (length(indSup) > 0) & (length(indApp) > 0) )
      supFct[indSup] <- amean(appFct[indApp])
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      gmean = by using geometric mean
#      bymot = by motif in a whole
#                            (WITHOUT taking into account species contribution)
#      by including all the assemblies, even the one to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      amean = by using geometric mean
#      bymot = by motif in a whole (WITHOUT taking into account
#                                                       species contribution)
#      by including all assemblies, even the one to predict
#'
#' @title Prediction of supplementary assemblages
#'
#' @description Take a numeric f.
#'
#' @usage predict_gmean_bymot_supp(appFct, appMotifs, supMotifs)
#'
#' @param appFct a vector of numeric values (assembly properties).
#'
#' @param appMotifs a vector of labels of \code{length(appFct)}
#'  (assembly motifs).
#'
#' @param supMotifs a vector of labels of assembly motifs of which values
#' must be predicted.
#'
#' @return Return a vector of \code{length(supMotifs)}. The values are
#'  computed using arithmetic mean of elements belonging to \code{appMotifs}
#'  and sharing a same motif.
#'
#' @details Prediction ...
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_gmean_bymot_supp <- function(appFct, appMotifs, supMotifs) {

  supFct   <- numeric(length(supMotifs))
  supFct[] <- NA

  setMot <- unique(supMotifs)
  for (mot in seq_along(setMot)) {

    indSup <- which(supMotifs == setMot[mot])
    indApp <- which(appMotifs == setMot[mot])
    if ( (length(indSup) > 0) & (length(indApp) > 0) )
      supFct[indSup] <- gmean(appFct[indApp])
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      amean = by using arithmetic mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#
#' @title Prediction of supplementary assemblages computed
#'
#' @description Take a numeric f.
#'
#' @usage predict_amean_byelt_supp(appFct, appMotifs, appOccur,
#'                              supMotifs, supOccur )
#'
#' @param appFct  cccc
#'
#' @param appMotifs  cccc
#'
#' @param appOccur  cccc
#'
#' @param supMotifs  cccc
#'
#' @param supOccur  cccc
#'
#' @details dd
#'
#' @return  cccc
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_amean_byelt_supp <- function(appFct, appMotifs, appOccur,
                                     supMotifs, supOccur ) {

  setAppMot <- unique(appMotifs)
  setSupMot <- unique(supMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2],
                      dimnames = list(setMot, colnames(appOccur)))

  for (mot in seq_along(setAppMot)) {

    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- amean(appFct[indApp[indElt]])
    }
  }


  supFct  <- numeric(length(supMotifs))
  sizeSup <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {
    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)

    if (length(indSupMot) > 0) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)

        if (length(indSupElt) > 0) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] + mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] / sizeSup[indSupMot]
    }
  }

  return(supFct)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction of supplementary assemblages computed
#      gmean = by using geometric mean
#      byelt = by motif WITH taking into account species contribution
#      by including all the assemblies, even the one to predict
#      for any Function (for instance Fobs)
#
#' @title Prediction of supplementary assemblages computed
#'
#' @description Take a numeric f.
#'
#' @usage predict_gmean_byelt_supp(appFct, appMotifs, appOccur,
#'                              supMotifs, supOccur )
#'
#' @param appFct   cccc
#'
#' @param appMotifs  cccc
#'
#' @param appOccur  cccc
#'
#' @param supMotifs  cccc
#'
#' @param supOccur  cccc
#'
#' @details dd
#'
#' @return  cccc
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_gmean_byelt_supp <- function(appFct, appMotifs, appOccur,
                                     supMotifs, supOccur  ) {

  setAppMot <- unique(appMotifs)
  setSupMot <- unique(supMotifs)
  setMot    <- sort(union(setAppMot, setSupMot))
  nbMot     <- length(setMot)

  mfct      <- matrix(NA, nrow = nbMot, ncol = dim(appOccur)[2],
                      dimnames = list(setMot, colnames(appOccur)))

  for (mot in seq_along(setAppMot)) {

    motif  <- setAppMot[mot]
    indApp <- which(appMotifs == motif)
    setElt <- unique(which((appOccur[indApp, , drop = FALSE] == 1),
                           arr.ind = TRUE)[ , 2])

    for (elt in seq_along(setElt)) {
      element <- setElt[elt]
      indElt  <- which(appOccur[indApp, element] == 1)
      if (length(indElt) > 0)
        mfct[motif, element] <- gmean(appFct[indApp[indElt]])
    }
  }


  supFct   <- numeric(length(supMotifs))
  supFct[] <- 1
  sizeSup  <- apply(supOccur, MARGIN = 1, FUN = sum)

  for (mot in seq_along(setSupMot)) {

    motif     <- setSupMot[mot]
    indSupMot <- which(supMotifs == motif)

    if (length(indSupMot) > 0) {
      setSupElt <- unique(which((supOccur[indSupMot, , drop = FALSE] == 1),
                                arr.ind = TRUE)[ , 2])

      for (elt in seq_along(setSupElt)) {
        element   <- setSupElt[elt]
        indSupElt <- which(supOccur[indSupMot, element] == 1)

        if (length(indSupElt) > 0) {
          index         <- indSupMot[indSupElt]
          supFct[index] <- supFct[index] * mfct[motif, element]
        }
      }

      supFct[indSupMot] <- supFct[indSupMot] ^ (1/sizeSup[indSupMot])
    }
  }

  return(supFct)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Functions for switch on different options                                ####
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
#' @usage predict_cal(fct, assMotifs, mOccur, xpr,
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
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_cal <- function(fct, assMotifs, mOccur, xpr,
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
             predict_amean_bymot(fct, assMotifs),
           amean.bymot.xpr =
             predict_amean_bymot_xpr(fct, assMotifs, xpr),
           gmean.bymot. =
             predict_gmean_bymot(fct, assMotifs),
           gmean.bymot.xpr =
             predict_gmean_bymot_xpr(fct, assMotifs, xpr),

           amean.byelt. =
             predict_amean_byelt(fct, assMotifs, mOccur),
           amean.byelt.xpr =
             predict_amean_byelt_xpr(fct, assMotifs, mOccur, xpr),
           gmean.byelt. =
             predict_gmean_byelt(fct, assMotifs, mOccur),
           gmean.byelt.xpr =
             predict_gmean_byelt_xpr(fct, assMotifs, mOccur, xpr)  )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Predition of assembly property by cross-validation
#'
#' @description Take a numeric vector and return the predicted vector
#' computed as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predict_prd(fct, assMotifs, mOccur, xpr,
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
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_prd <- function(fct, assMotifs, mOccur, xpr,
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
             predict_amean_bymot_LOO(fct, assMotifs),
           amean.bymot..xpr =
             predict_amean_bymot_LOO_xpr(fct, assMotifs, xpr),
           amean.bymot.jack. =
             predict_amean_bymot_jack(fct, assMotifs, jack),
           amean.bymot.jack.xpr =
             predict_amean_bymot_jack_xpr(fct, assMotifs, jack, xpr),

           gmean.bymot.. =
             predict_gmean_bymot_LOO(fct, assMotifs),
           gmean.bymot..xpr =
             predict_gmean_bymot_LOO_xpr(fct, assMotifs, xpr),
           gmean.bymot.jack. =
             predict_gmean_bymot_jack(fct, assMotifs, jack),
           gmean.bymot.jack.xpr =
             predict_gmean_bymot_jack_xpr(fct, assMotifs, jack, xpr),

           amean.byelt.. =
             predict_amean_byelt_LOO(fct, assMotifs, mOccur),
           amean.byelt..xpr =
             predict_amean_byelt_LOO_xpr(fct, assMotifs, mOccur, xpr),
           amean.byelt.jack. =
             predict_amean_byelt_jack(fct, assMotifs, mOccur, jack),
           amean.byelt.jack.xpr =
             predict_amean_byelt_jack_xpr(fct, assMotifs, mOccur, jack, xpr),

           gmean.byelt.. =
             predict_gmean_byelt_LOO(fct, assMotifs, mOccur),
           gmean.byelt..xpr =
             predict_gmean_byelt_LOO_xpr(fct, assMotifs, mOccur, xpr),
           gmean.byelt.jack. =
             predict_gmean_byelt_jack(fct, assMotifs, mOccur, jack),
           gmean.byelt.jack.xpr =
             predict_gmean_byelt_jack_xpr(fct, assMotifs, mOccur, jack, xpr) )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Prediction computed by excluding (LOO) the assembly to predict
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Arithmetic mean (amean) by motif (bymot) by jackknife (jack) over
#'  several experiments (xpr)
#'
#' @description Take a numeric vector and return the predicted vector computed
#'  as the arithmetic mean of all elements belonging to a same motif.
#'
#' @usage predict_supp(appFct, appMotifs, appOccur,
#'                            supMotifs, supOccur,
#'                            opt.mean = "amean",
#'                            opt.mod  = "bymot"  )
#'
#' @param appFct a vector of numeric values (assembly properties).
#'
#' @param appMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param appOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of elements.
#'
#' @param supMotifs a vector of labels of \code{length(fct)} (assembly motifs).
#'
#' @param supOccur a matrix of occurrence (occurrence of elements).
#' Its first dimension equals to \code{length(fct)}. Its second dimension
#'  equals to the number of elements.
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
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_supp <- function(appFct, appMotifs, appOccur,
                         supMotifs, supOccur,

                         opt.mean = "amean",
                         opt.mod  = "bymot"  ) {

  option <- paste(opt.mean, opt.mod, sep = ".")

  return(
    switch(option,
           amean.bymot =
             predict_amean_bymot_supp(appFct, appMotifs, supMotifs) ,
           gmean.bymot =
             predict_gmean_bymot_supp(appFct, appMotifs, supMotifs) ,

           amean.byelt =
             predict_amean_byelt_supp(appFct, appMotifs, appOccur,
                                      supMotifs, supOccur) ,
           gmean.byelt =
             predict_gmean_byelt_supp(appFct, appMotifs, appOccur,
                                      supMotifs, supOccur)
    )
  )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                   END of FILE
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


