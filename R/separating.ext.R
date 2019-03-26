#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                                 myCOMBINAT.R
#
#  set of functions for combinatorial analysis of community diversity effects
#
#                         Benoit JAILLARD, summer 2016
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include stats.int.R  clustering.int.R


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                               mySEPARATING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Remove species assemblies that were in several exemplar
#     and return the right mat (with only one assemblage)
#' @description cc
#' @usage rm_dual_assemblies(mat, fct,
#'                           xpr = rep(1, length(fct)),
#'                           opt.mean = "amean")
#'
#' @param mat bb
#' @param fct  vv
#' @param xpr vv
#' @param opt.mean vv
#'
#' @return ff
#' @export
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm_dual_assemblies <- function(mat, fct,
                               xpr      = rep(1, length(fct)),
                               opt.mean = "amean") {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setXpr <- unique(xpr)
  setAss <- which(xpr == setXpr[1])
  pas    <- length(setAss)
  mut    <- mat[setAss, ]


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  index <- which(is.na(fct) == FALSE)
  fct   <- fct[index]
  mut   <- mut[index, ]

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  index  <- which(apply(mut, MARGIN = 1, FUN = sum) != 0)
  fct   <- fct[index]
  mut   <- mut[index, ]

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  AssMotifs  <- affect_motifs(c(1:dim(mut)[2]), mut)
  index      <- table(AssMotifs)
  bool       <- !logical(length(AssMotifs))

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setDoublon <- which(index > 1)
  for (i in seq_along(setDoublon)) {
    doublon <- which(AssMotifs == setDoublon[i])

    for (ipr in seq_along(setXpr)) {
      off <- pas * (ipr - 1)
      fct[doublon[1] + off] <- mean_fct(fct[doublon + off], opt.mean)
      bool[doublon[2:length(doublon)] + off] <- FALSE
    }
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  res        <- list(mat[bool, ], fct[bool])
  names(res) <- c("mat", "fct")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title  Remove one or more species from data
#     and return the index.row and index.col to select
#' @description ddd
#' @usage rm_elements(mOccur, fct, elements)
#'
#' @param mOccur dd
#' @param fct dd
#' @param elements dd
#'
#' @return dd
#' @export
#'
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm_elements <- function(mOccur, fct, elements) {

  index.row <- c(1:dim(mOccur)[1])
  index.col <- c(1:dim(mOccur)[2])

  for (elt in seq_along(elements)) {
    index.row <- setdiff(index.row, which(mOccur[ ,elements[elt]] == 1))
    index.col <- setdiff(index.col, which(elements[elt] == colnames(mOccur)))
  }

  mOccur <- mOccur[index.row, index.col]
  fct    <- fct[index.row]

  res        <- list(mOccur, fct)
  names(res) <- c("mOccur", "fct")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Decompose the observed performance into a.inter and a.comp ratios
#' @description bbb
#' @usage multiplicative_decomposition(mOccur, Fobs, rm.mono = TRUE)

#' @param mOccur ff
#' @param Fobs dd
#' @param rm.mono dd
#'
#' @return dd
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


multiplicative_decomposition <- function(mOccur, Fobs,
                                         rm.mono = TRUE) {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check if all monocultures are observed
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  size      <- apply(mOccur, MARGIN = 1, FUN = sum)
  mono      <- which(size == 1)
  if (length(mono) == 0) stop("There is no monoculture in mOccur")

  elements  <- which(apply(mOccur[mono, ], MARGIN = 2, FUN = sum) != 0)

  if (length(elements) < dim(mOccur)[2]) {
    missing <- setdiff(colnames(mOccur), colnames(mOccur)[elements])
    res     <- rm_elements(mOccur, Fobs, missing)

    mOccur  <- res$mOccur
    Fobs    <- res$fct
    size    <- apply(mOccur, MARGIN = 1, FUN = sum)
    mono    <- which(size == 1)
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Look for and compute the monoculture performance
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  Fmono <-
    apply(mOccur[mono, ] * Fobs[mono], MARGIN = 2, FUN = sum) /
    apply(mOccur[mono, ], MARGIN = 2, FUN = sum)

  Fref    <- (mOccur %*% Fmono) / size
  Fscale  <- mean(Fmono)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Compute the interaction and composition effects
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  alpha     <- as.vector(Fobs/Fref)
  beta      <- as.vector(Fref/Fscale)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Remove the monocultures
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (rm.mono == TRUE) {
    index  <- which(size != 1)
    size   <- size[index]
    mOccur <- mOccur[index, ]
    Fobs   <- Fobs[index]
    alpha  <- alpha[index]
    beta   <- beta[index]
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Identify the transgressive performances
  # 1: transgressive under-yielding;  2: under-yielding
  # 3: over-yielding ;                4: transgressive over-yielding
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  trans <- numeric(dim(mOccur)[1])
  tmp   <- t( t(mOccur) * Fmono )
  for (i in 1:dim(mOccur)[1]) {
    ttp         <- tmp[i, (tmp[i,] != 0)]
    if (Fobs[i] <  mean(ttp)) trans[i] <- 2
    if (Fobs[i] <  min(ttp))  trans[i] <- 1
    if (Fobs[i] >= mean(ttp)) trans[i] <- 3
    if (Fobs[i] >  max(ttp))  trans[i] <- 4
  }


  res        <- list(size, trans, mOccur,
              Fobs, alpha, beta, Fmono, Fscale)
  names(res) <- c("size", "trans", "mOccur",
                  "Fobs", "alpha", "beta", "Fmono", "Fscale")

  return(res)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                         END of FILE mySEPARATING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
