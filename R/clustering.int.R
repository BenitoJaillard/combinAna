#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#               CLUSTERING OF ASSEMBLY ELEMENTS AND MOTIFS                 ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#' @include stats.int.R
NULL



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Return a vector continuously indexed from 1 to max(index)
#'
#' @description ddd
#'
#' @usage compact_index(v)
#'
#' @param v gg
#'
#' @details dd
#'
#' @return gg
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compact_index <- function(v) {

  res <- v
  set <- sort(unique(v))
  for (i in seq_along(set)) res[v == set[i]] <- i

  return(res)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title kkk
#' @description dddd
#' @usage cut_tree(res, nbcl)

#' @param res gg
#' @param nbcl gg
#'
#' @details dd
#'
#' @return gg
#' @keywords internal
#'
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cut_tree <- function(res, nbcl) {

  return( compact_index(res$aff[nbcl, ]) )
}


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Sort the resulting file of tree
#
#  Input:  X : matrix of leaves : affectation of species to classes
#          Y : vector of distance
#          Z : matrix of distance from the best element
#  Output: files of leaves sorted for plotting
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title  Sort the resulting file of tree
#' @description ddd
#' @usage sort_tree(tree, index.return = FALSE)
#' @param tree ff
#' @param index.return ff
#'
#' @details dd
#'
#' @return ff
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

sort_tree <- function(tree, index.return = FALSE) {

  X      <- shift_affectElt(tree$aff)

  nbline <- dim(X)[1]
  nbitem <- dim(X)[2]


  # Separate first the sorted from unsorted elements
  mask <- ordre <- seq_len(nbitem)
  for (lin in 2:nbline) {
    index       <- sort(x = X[lin, mask], index.return = TRUE)
    index       <- index$ix + nbitem - length(mask)
    X[ , mask]  <- X[ , index]
    ordre[mask] <- ordre[index]
    mask        <- mask[X[lin, mask] == max(X)]
  }

  # Sort the leaves of tree
  new <- old <- seq_len(nbitem)
  for (lin in seq_len(nbline)) {
    for (item in seq_len(nbitem))
      old[item] <- which(unique(X[lin,]) == X[lin, item])
    new             <- sort(old, index.return = TRUE)
    X[ , seq_len(nbitem)]  <- X[ , new$ix]
    ordre[seq_len(nbitem)] <- ordre[new$ix]
  }

  # Re-numeration of leaves
  old <- integer(nbline)
  for (lin in seq_len(nbline)) {
    v  <- setdiff(unique(X[lin, ]), unique(old))
    if (length(v) == 1) old[lin] <- v
  }

  newX <- matrix(0, nrow = nbline, ncol = nbitem)
  for (lin in seq_len(nbline))
    for (j in seq_len(lin)) newX[lin, which(old[j] == X[lin, ])] <- j

  colnames(newX) <- colnames(X)[ordre]

  if (index.return == TRUE) {
    res        <- list(newX, ordre, colnames(newX))
    names(res) <- c("x", "ix", "noms")
  } else {
    res <- newX
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#          Plot a Divisive or Agglomerative Hierarchical Clustering
#    obtained with the options ["divisive" and "tree"] or "agglomerative"
#
#  Inputs: X: matrix of classes of elements
#          Y: vector of distance
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Plot a Divisive or Agglomerative Hierarchical Clustering
#' @description ddd
#' @usage plot_tree(tree, signifElt = NULL, col = "black", titre = "")
#' @param tree  kk
#' @param signifElt  ll
#' @param col  ll
#' @param titre ll
#'
#' @details dd
#'
#' @return ll
#' @importFrom graphics plot points lines axis title
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_tree <- function(tree,
                      signifElt = NULL,
                      col       = "black", titre = "") {

  if (length(signifElt) == 0) signifElt <- rep(TRUE, length(tree$cor))

  vpch            <- numeric(dim(tree$aff)[2])
  vpch[ ]         <- 1
  vpch[signifElt] <- 19

  if (length(tree$cor) > 1) {
    tmpx <- sort_tree(tree, index.return = TRUE)
    X    <- tmpx$x
    vpch <- vpch[tmpx$ix]
    if (length(col) > 1) col <- col[tmpx$ix]
  } else {
    X    <- shift_affectElt(tree$aff)
    vpch <- vpch[1]
    if (length(col) > 1) col <- col[1]
  }

  Y   <- tree$cor
  Y[Y < 0] <- 0


  nbline <- dim(X)[1]
  nbitem <- dim(X)[2] + 1
  YY     <- c(0,Y)

  xx     <- matrix(0, nrow = nbline, ncol = nbline)
  xy     <- NULL

  # plot the framework
  graphics::plot(x = seq_len(nbitem), y = seq(1/nbitem, 1, 1/nbitem),
       xlab = "element",    ylab = "R2-value",
       xlim = c(1, nbitem), ylim = c(0, 1),
       type = "n", tck = 0.02, las = 1)

  # plot the vertical lines
  for (lin in 1:nbline)
    for (leave in 1:lin) {
      xx[lin, leave] <- 1 + mean(which(X[lin, ] == leave))
      graphics::lines(x = rep(xx[lin, leave], 2), y = c(YY[lin], YY[lin + 1]),
            lty = "solid")
    }

  # plot the horizontal jonction lines
  if (nbline > 1) for (lin in 2:nbline) {
    XX <- which(xx[lin, ] != xx[lin - 1, ])
    if (length(XX) == 1) XX <- rep(XX, 2)
    graphics::lines(x = xx[lin, XX], y = rep(YY[lin], 2), lty = "solid" )
  }

  # plot the final horizontal jonction lines
  setLeave <- sort(unique(X[nbline, ]))
  for (leave in seq_along(setLeave)) {
    tmp <- which(X[nbline, ] == setLeave[leave])
    graphics::lines(x = c(min(tmp), max(tmp)) + 1,
          y = rep(YY[nbline + 1], 2),
          lty = "solid" )
  }

  # plot the symbols of the most likely partition
  xpos <- 1 + c(1:(nbitem - 1))
  ypos <- YY[length(YY)] + 0.02

  if (length(col) > 1) {
    graphics::points(x = xpos, y = rep(ypos, (nbitem - 1)), pch = vpch,
           cex = 2, col = col, bg = "white")
  } else {
    graphics::points(x = xpos, y = rep(ypos, (nbitem - 1)), pch = vpch,
                     bg = "white")
  }

  #  plot the names of elements
  graphics::axis(side = 3, at = 2:nbitem, labels = colnames(X),
       tick = FALSE, las = 2, pos = ypos)

  graphics::title(titre)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Cluster the assemblages by assembly motifs
#
#  Inputs  : affectElt and mat Occurrence
#  Outputs : the vector of motifs to which belong the assemblages
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Cluster the assemblages by assembly motifs
#' @description fff
#' @usage affect_motifs(affectElt, mOccur)
#'
#' @param affectElt gg
#' @param mOccur gg
#'
#' @details dd
#'
#' @return ff
#' @importFrom plyr alply llply
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

affect_motifs <- function(affectElt, mOccur) {

  nbAss <- dim(mOccur)[1]

  EltClasses <- t( t(mOccur) * compact_index(affectElt) )
  rownames(EltClasses) <- seq_len(nbAss)

  myfct1  <- function(x) { x[x != 0] }
  tmp1    <- plyr::alply(EltClasses[ , , drop = FALSE],
                         .margins = 1, .fun = myfct1)

  myfct2  <- function(x) { sort(unique(x)) }
  tmp2    <- plyr::llply(tmp1, .fun = myfct2)

  size    <- unlist(plyr::llply(tmp2, .fun = length))
  setSize <- sort(unique(size))

  motif     <- 1
  AssMotifs <- integer(nbAss)
  for (siz in seq_along(setSize)) {

    index <- which(size == setSize[siz])

    while (length(index) > 0) {
      ind <- 1
      ref <- tmp2[[index[ind]]]
      for (i in 2:length(index))
        if (identical(ref, tmp2[[index[i]]])) ind <- c(ind, i)

      AssMotifs[index[ind]] <- motif

      motif <- motif + 1
      index <- index[-ind]
    }
  }

  return(AssMotifs)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Cluster the assemblages by assembly motifs
#
#  Inputs  : a tree and mat Occurrence
#  Outputs : the matrix of motifs to which belong the assemblages
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Cluster the assemblages by assembly motifs
#' @description fff
#' @usage maffect_motifs(tree, mOccur)
#'
#' @param tree ff
#' @param mOccur ff
#'
#' @details dd
#'
#' @return gg
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

maffect_motifs <- function(tree, mOccur) {

  nbclTot <- length(tree$cor)
  nbclMax <- which(tree$cor == max(tree$cor))[1]

  mAssMotifs <- matrix(as.integer(0), nrow = nbclTot, ncol = dim(mOccur)[1],
                       dimnames = list(c(1:nbclTot), rownames(mOccur)))

  for (nbcl in seq_len(nbclMax)) {
    affectElt          <- cut_tree(tree, nbcl)
    mAssMotifs[nbcl, ] <- affect_motifs(affectElt, mOccur)
  }

  if (nbclMax < nbclTot)
    for (nbcl in (nbclMax + 1):nbclTot)
      mAssMotifs[nbcl, ] <- mAssMotifs[nbclMax, ]

  return(mAssMotifs)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the names (in lowercase letters) of elements, clusters of elements,
#              and all possible assembly motifs associated to element clusters
#
#  Input:  mOccur    : matrix of occurrence
#          affectElt : affectation of species into classes
#  Output: names.elements         : names of elements
#          names.element.clusters : names of clusters of elements
#          names.assembly.motifs  : names of elemental assembly motifs
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title names (in minus letters) of elements,
#' @description bbb
#' @usage names_assembly(affectElt, mOccur)

#' @param affectElt tt
#' @param mOccur ttt
#'
#' @details dd
#'
#' @return tt
#' @importFrom plyr alply llply
#' @keywords internal
#'
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

names_assembly <- function(affectElt, mOccur) {

  nbAss <- dim(mOccur)[1]
  nbElt <- dim(mOccur)[2]

  # label the element clusters
  affectElt <- shift_affectElt(affectElt)

  nomAffect <- myLetters[affectElt]

  # label the assembly motifs
  EltClasses <- t( t(mOccur) * affectElt )
  rownames(EltClasses) <- seq_len(nbAss)

  myfct1  <- function(x) { x[x != 0] }
  tmp1    <- plyr::alply(EltClasses[ , , drop = FALSE],
                         .margins = 1, .fun = myfct1)

  myfct2  <- function(x) { sort(unique(x)) }
  tmp2    <- plyr::llply(tmp1, .fun = myfct2)

  nomMotif <- character(nbAss)
  for (ass in seq_len(nbAss))
    for (elt in seq_along(tmp2[[ass]]))
      nomMotif[ass] <- paste(nomMotif[ass],
                             myLetters[tmp2[[ass]][elt]], sep = "")

  noms        <- list(colnames(mOccur), nomAffect, rownames(mOccur), nomMotif)
  names(noms) <- c("elements", "clusters", "assemblies", "motifs")

  return(noms)
}


#system.time(for (i in 1:100000) affect_motif.name(x)    )



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title mmm
#' @description fff
#' @usage nbMax_tests(s)
#'
#' @param s ff
#'
#' @details dd
#'
#' @return gg
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

nbMax_tests <- function(s) {

  return( s*(s + 1)*(2*s + 1)/12 + s*(s + 1)/4 - s)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Indexes the species clusters by decreasing structuring effect
#
#   Input and output: a vector of affectation
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#'
#' @title Indexes the species clusters by decreasing structuring effect
#' @description ddd
#' @usage shift_affectElt(affectElt)
#'
#' @param affectElt gg
#'
#' @details dd
#'
#' @return ff
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

shift_affectElt <- function(affectElt) {

  affectElt[affectElt == 1] <- max(affectElt) + 1
  affectElt <- affectElt - 1

  return(affectElt)
}




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                   END of FILE myCLUSTERING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
