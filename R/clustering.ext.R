#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#               CLUSTERING OF ASSEMBLY ELEMENTS AND MOTIFS                 ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



#' @include stats.int.R predicting.int.R clustering.int.R
NULL



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Return the R2-value of the element clustering and Fpredicted values
#
#  The function "R2.cluster.elt" is autonomous.
#
#  An alternative way is to run the function such as :
#         R2.cluster.elt(mOccur, fct, affect, nbclElt) or
#
#  Input:  mOccur  : matrix of occurrence
#          fct    : vector of performance
#          affect  : affectation of elements into classes
#          nbclElt : number of element classes
#
#  Output: R2      : R2-value between fpred and fct
#          fpred   : vector of predicted performance
#
#' @title  ddd
#' @description lll
#' @usage rss_cluster_elt(affectElt, mOccur, fct, xpr,
#'                        opt.mean, opt.mod)
#'
#'
#' @param affectElt ff
#' @param mOccur ff
#' @param fct ff
#' @param xpr ff
#' @param opt.mean ff
#' @param opt.mod gg
#'
#' @return gg
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rss_cluster_elt <- function(affectElt, mOccur, fct, xpr,
                            opt.mean, opt.mod) {

  AssMotifs <- affect_motifs(affectElt, mOccur)
  fprd      <- predict_cal(AssMotifs, mOccur, fct, xpr,
                           opt.mean, opt.mod)

  return( rss(fprd, fct) )
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#               CLUSTERING of ELEMENTS of ASSEMBLIES
#       for their ability to shift an assembly function fct
#
# Input: mOccur: matrix of occurence of nbElt elements within the assemblies
#        fct   : observed performances of assemblies
#        nbclEltMax : number maximum of classes of elements
#        option = "divisive" or "agglomerative"
#        method = "sort", "cluster" or "tree"
#
# if option="divisive" and method="sort": we recommand to left nbclEltMax=nbElt
# if option="divisive" and method="cluster": nbclEltMax is important
# if option="divisive" and method="tree": we recommand to left nbclEltMax=nbElt
# if option="agglomerative": method="tree" and nbclEltMax=nbElt.
#
# Outputs : elt : successive affectations of elements,
#           cor : the corresponding correlations
# if option="divisive" and method="sort" or "cluster":
#
#' @title fffff
#' @description ddd
#' @usage cluster_elements(mOccur, fct,
#'                         xpr = rep(1, length(fct)),
#'                         opt.meth = "divisive",
#'                         opt.mean = "amean",
#'                         opt.mod  = "byelt",
#'                         opt.prn  = TRUE)
#'
#'
#' @param mOccur ff
#' @param fct ff
#' @param xpr ff
#' @param opt.meth gg
#' @param opt.mean gg
#' @param opt.mod gg
#' @param opt.prn gg
#'
#' @return gg
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


cluster_elements <- function(mOccur, fct,
                             xpr = rep(1, length(fct)),

# options for computing
                             opt.meth = "divisive",
                             opt.mean = "amean",
                             opt.mod  = "byelt",
# options for plotting
                             opt.prn  = TRUE)
{

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Checking the inputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  tmp <- names(table(mOccur))
  if ( (sum(is.na(mOccur)) != 0) || (length(tmp) != 2) ||
       ( (tmp[1] != "0") & (tmp[1] != "FALSE") ) ||
       ( (tmp[2] != "1") & (tmp[2] != "TRUE") ) )
        stop("The matrix of occurence should be a binary matrix")

  nbAss <- dim(mOccur)[1]
  nbElt <- dim(mOccur)[2]

  if ( (length(fct) != nbAss) | (sum(!is.na(fct)) != nbAss) )
    stop("dim(mOccur)[1] should be equal to length(fct)")

  if ( (opt.meth != "sort") & (opt.meth != "divisive") &
       (opt.meth != "agglomerative") )
    stop(" 'opt.meth' can be 'sort', 'divisive' or 'agglomerative' ")

  if ( (opt.mean != "amean") & (opt.mean != "gmean") )
    stop(" 'opt.mean' can be 'amean' or 'gmean' ")

  if ( (opt.mod != "bymot") & (opt.mod != "byelt") )
    stop(" 'opt.mod' can be 'bymot' or 'byelt' ")

  res           <- matrix(Inf, nrow = nbElt, ncol = nbElt)
  rownames(res) <- c(1:nbElt)
  colnames(res) <- colnames(mOccur)



  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Choice of methods
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  if (opt.meth == "sort") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #               One-by-one Divisive SORTING of SPECIES
    #
    # Initial: all elements are gathered into a single big cluster. R2 is null.
    # Iteration:
    #      each element belonging to the big cluster is successively removed,
    #      then put as a new cluster containing only a single species.
    # Each new clustering is evaluated by R2
    #                                    of resulting community classification.
    # The most likely clustering is kept.
    # A new element is isolated in a singleton.
    # The process stops when each element is isolated into a singleton.
    #
    # Input : matrix of occurence
    #         the corresponding vector of observations
    # Outputs : successive affectations,
    #           the corresponding correlations
    #           the whole used matrix of correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # First line : Trunk
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    nbclElt   <- 1
    affectElt <- as.integer(rep(1, nbElt))
    RES.aff   <- matrix(affectElt, nrow = nbElt, ncol = nbElt, byrow = TRUE,
                        dimnames = list(seq_len(nbElt), colnames(mOccur)))

    RES.rss          <- numeric(nbElt)
    RES.rss[nbclElt] <- oldRes <-
      rss_cluster_elt(affectElt, mOccur, fct, xpr,
                      opt.mean, opt.mod)

    res[,] <- Inf

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Other lines : Branches
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    if (nbElt > 1) for (nbclElt in 2:nbElt) {

#      nbclElt <- nbclElt + 1

      setElt  <- which(affectElt == 1)

      # Test the elements one-after-one into a new singleton
      for (elt in setElt) {
        oldAffect         <- affectElt
        affectElt[elt]    <- nbclElt
        res[nbclElt, elt] <-
          rss_cluster_elt(affectElt, mOccur, fct, xpr,
                          opt.mean, opt.mod)
        affectElt         <- oldAffect
      }

      # Decision: one keeps the best fit
      minRes <- min(res[nbclElt, ], na.rm = TRUE)

      #     with predict_cal(), oldRes is always < Inf,
      # but with predict.prd(), oldRes is == Inf when all Pred == NA
      if (minRes < Inf) {
        elt    <- which(res[nbclElt, ] == minRes, arr.ind = TRUE)[1]
        affectElt[elt] <- nbclElt

        RES.aff[nbclElt, ] <- affectElt
        RES.rss[nbclElt]   <-
          rss_cluster_elt(affectElt, mOccur, fct, xpr,
                          opt.mean, opt.mod)
      } else {
        nbclElt <- nbclElt - 1
      }

      # to avoid unuseful time-consuming computations...
      if ( (minRes == Inf) | (RES.rss[nbclElt] == 0) |
           (RES.rss[nbclElt] == RES.rss[nbclElt - 1]) ) {

        tmp   <- table(RES.aff[nbclElt, ])
        index <- as.integer(names(tmp))[tmp > 1]
        while (length(index) > 0) {
          nbclElt <- nbclElt + 1
          affectElt[ which(affectElt == index[1])[1] ] <- nbclElt

          RES.aff[nbclElt, ] <- affectElt
          RES.rss[nbclElt]   <- RES.rss[nbclElt - 1]

          tmp     <- table(RES.aff[nbclElt, ])
          index   <- as.integer(names(tmp))[tmp > 1]
        }
        break
      }

      # print to follow the computation
      if (opt.prn) {
        if (nbAss < 30) {
          print(cbind(RES.aff[1:nbclElt, , drop = FALSE ], RES.rss[1:nbclElt]))
        } else {
          print(RES.rss[1:nbclElt])
        }
      } else {
        cat(".", fill = FALSE)
      }

    }    #  End of the FOR

  }    #    END of IF on DIVISIVE SORTING of ELEMENTS
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx




  if (opt.meth == "divisive") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #    DIVISIVE HIERARCHICAL CLUSTERING of ELEMENTS by DICHOTOMY (ROUX, 2005)
    #
    # Initial: all elements are gathered into a single big cluster. R2 is null.
    # Iteration:
    #   each element belonging to the big cluster is successively removed,
    #            then put in other clusters while R2 increased.
    #
    # Input : matrix of occurence,
    #         the corresponding vector of observations,
    #         a maximum number of element clusters.
    # Outputs : successive affectations,
    #           corresponding correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # First line : Trunk
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    nbclElt   <- 1
    affectElt <- as.integer(rep(1, nbElt))
    RES.aff   <- matrix(affectElt, nrow = nbElt, ncol = nbElt, byrow = TRUE,
                        dimnames = list(seq_len(nbElt), colnames(mOccur)))

    RES.rss    <- numeric(nbElt)
    RES.rss[1] <- oldRes <- rss_cluster_elt(affectElt, mOccur, fct, xpr,
                                            opt.mean, opt.mod)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines : Leaves
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    oldRes.affectElt <- affectElt

    if (nbElt > 1) for (nbclElt in 2:nbElt) {

#      nbclElt <- nbclElt + 1

      oldRes  <- res[ , ] <- Inf

      # One tests split of each leaves of tree
      for (clElt in 1:(nbclElt - 1)) {
        affectElt <- RES.aff[nbclElt - 1, ]
        setElt    <- which(affectElt == clElt)

        test2     <- (length(setElt) > 1)
        while (test2) {
          # One tests elements one-after-one
          for (elt in setElt) {
            oldAffect       <- affectElt
            affectElt[elt]  <- nbclElt
            res[clElt, elt] <-
              rss_cluster_elt(affectElt, mOccur, fct, xpr,
                              opt.mean, opt.mod)
            affectElt       <- oldAffect
          }

          # Decision: one keeps the best local fit
          minRes <- min(res, na.rm = TRUE)
          if (minRes < oldRes) {
            coord                 <- which(res == minRes, arr.ind = TRUE)
            affectElt[coord[1,2]] <- nbclElt
            oldRes                <- minRes
            oldRes.affectElt      <- affectElt
          } else {
            test2 <- FALSE
          }
        }
        # END of WHILE on test
      }
      # END of LOOP on LEAVE

      #     with predict_cal(), oldRes is always < Inf,
      # but with predcit.prd(), oldRes is == Inf when all Pred == NA
      if (minRes < Inf) {
        affectElt          <- oldRes.affectElt
        RES.aff[nbclElt, ] <- affectElt
        RES.rss[nbclElt]   <-
          rss_cluster_elt(affectElt, mOccur, fct, xpr,
                          opt.mean, opt.mod)
      } else {
        nbclElt <- nbclElt - 1
      }

      # to avoid unuseful time-consuming computations...
      if ( (minRes == Inf) | (RES.rss[nbclElt] == 0) |
           (RES.rss[nbclElt] == RES.rss[nbclElt - 1]) ) {

        tmp   <- table(RES.aff[nbclElt, ])
        index <- as.integer(names(tmp))[tmp > 1]
        while (length(index) > 0) {
          nbclElt <- nbclElt + 1
          affectElt[ which(affectElt == index[1])[1] ] <- nbclElt

          RES.aff[nbclElt, ] <- affectElt
          RES.rss[nbclElt]   <- RES.rss[nbclElt - 1]

          tmp     <- table(RES.aff[nbclElt, ])
          index   <- as.integer(names(tmp))[tmp > 1]
        }
        break
      }

      # print to follow the computation
      if (opt.prn) {
        print(cbind(RES.aff[1:nbclElt, , drop = FALSE ], RES.rss[1:nbclElt]))
      } else {
        cat(".", fill = FALSE)
      }

    }    # END of LOOP

  }   # END of DIVISIVE HIERARCHICAL CLUSTERING of ELEMENTS by DICHOTOMY
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



  if (opt.meth == "agglomerative") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #     Single linkage AGGLOMERATIVE HIERARCHICAL CLUSTERING of ELEMENTS
    #      #              according to Legendre et Legendre, 1998    #
    # Initial: each element is  singleton.
    # Iteration: the elements are gathered two-by-two
    #            until  they are all togetger in a big cluster.
    #
    # Input : matrix of occurence,
    #         the corresponding vector of observations,
    #         a maximum number of element clusters.
    # Outputs : successive affectations,
    #           corresponding correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    nbclElt <- nbElt
    RES.aff <- matrix(as.integer(0), nrow = nbElt, ncol = nbElt,
                      dimnames = list(seq_len(nbElt), colnames(mOccur)))
    RES.aff[nbElt, ] <- affectElt <- as.integer(seq(1, nbElt))
    names(affectElt) <- colnames(mOccur)

    RES.rss <- numeric(nbElt)
    RES.rss[nbElt] <- oldRes <-
      rss_cluster_elt(affectElt, mOccur, fct, xpr,
                      opt.mean, opt.mod)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines : branchs and trunc
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Test adding on a element to each assembly class
    if (nbElt > 1) for (nbclElt in (nbElt - 1):1) {

#      nbclElt  <- nbclElt - 1

      res[ , ] <- 0
      setElt   <- sort(unique(affectElt))

      set1 <- setdiff(setElt, max(setElt))
      for (elt1 in set1) {
        set2 <- setElt[setElt > elt1]
        # Test the elements one-after-one
        for (elt2 in set2) {
          oldAffect       <- affectElt
          affectElt[affectElt == elt2] <- elt1
          res[elt1, elt2] <-
            rss_cluster_elt(affectElt, mOccur, fct, xpr,
                            opt.mean, opt.mod)
          affectElt       <- oldAffect
        }
      }

      # Decision: one keeps the best fit
      affectElt <- RES.aff[nbclElt + 1, ]

      maxRes <- min(res[res != 0])
      coord  <- which(res == maxRes, arr.ind = TRUE)

      affectElt[affectElt == coord[1, 2]] <- coord[1, 1]
      RES.aff[nbclElt, ] <- affectElt
      RES.rss[nbclElt]   <-
        rss_cluster_elt(affectElt, mOccur, fct, xpr,
                        opt.mean, opt.mod)

      # print to follow the computation
      if (opt.prn) {
        print(cbind(RES.aff[nbclElt:nbElt, , drop = FALSE ],
                    RES.rss[nbclElt:nbElt]))
      } else {
        cat(".", fill = FALSE)
      }

    }
    # END of LOOP
  }
  #   END of IF on AGGLOMERATIVE HIERARCHICAL CLUSTERING of ELEMENTS
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



  if (opt.meth == "cluster") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #      One-by-one Divisive CLUSTERING of SPECIES (McNaughton-Smith, 1964)
    #
    # Initial: all elements are gathered into a single big cluster. R2 is null.
    # Iteration:
    # each element belonging to the big cluster is successively removed,
    #            then put in other clusters or a new one.
    # Each new clustering is evaluated by R2
    #                                    of resulting community classification.
    # The most likely clustering is kept. A new element is clustered.
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # First line : Trunk
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    nbclElt   <- 1
    affectElt <- as.integer(rep(1, nbElt))
    RES.aff   <- matrix(affectElt, nrow = nbElt, ncol = nbElt, byrow = TRUE,
                        dimnames = list(seq_len(nbElt), colnames(mOccur)))
    RES.aff[nbclElt, ] <- affectElt

    RES.rss            <- numeric(nbElt)
    RES.rss[1] <- oldRes <- rss_cluster_elt(affectElt, mOccur, fct, xpr,
                                            opt.mean, opt.mod)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines : Leaves
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    oldRes.affectElt <- affectElt

    nbclEltMax <- 8

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    test    <- TRUE
    while (test) {
      res[ , ] <- Inf
      setElt   <- which(affectElt == 1)

      nbclElt  <- min(length(unique(affectElt)) + 1, nbclEltMax)

      # Put each element in a new cluster,
      #       then test if the R2 increase is significantly better
      #       than in yet existing cluster
      for (clElt in nbclElt:2) {
        for (elt in setElt) {
          oldAffect      <- affectElt
          affectElt[elt] <- clElt
          res[clElt, elt] <-
            rss_cluster_elt(affectElt, mOccur, fct, xpr,
                            opt.mean, opt.mod)
          affectElt      <- oldAffect
        }
      }

      # Decision: one keeps the best fit
      minRes <- min(res)
      if (minRes < oldRes) {
        coord                 <- which(res == minRes, arr.ind = TRUE)
        affectElt[coord[1,2]] <- coord[1,1]
        nbclElt               <- min(length(unique(affectElt)), nbclEltMax)
        oldRes                <- minRes

        RES.aff[nbclElt, ]    <- affectElt
        RES.rss[nbclElt]      <- minRes
      } else {
        test <- FALSE
      }

      # print to follow the computation
      if (opt.prn) {
        print(cbind(RES.aff[1:nbclElt, , drop = FALSE ], RES.rss[1:nbclElt]))
      } else {
        cat(".", fill = FALSE)
      }

    }
    #  End of the WHILE

  }
  #    END of IF on DIVISIVE CLUSTERING of ELEMENTS
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  RES.cor <- 1 - RES.rss / rss(fct, mean_fct(fct, opt.mean))

  rownames(RES.aff) <- names(RES.cor) <- c(1:dim(RES.aff)[1])

  res        <- list(RES.aff, RES.cor)
  names(res) <- c("aff", "cor")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   END of cluster_elements
#
#' @title rrr
#' @description ddd
#' @usage cluster_elements_LOO(mOccur, fct,
#'                             xpr      = rep(1, length(fct)),
#'                             opt.meth = "divisive",
#'                             opt.mean = "amean",
#'                             opt.mod  = "byelt",
#'                             opt.prn  = TRUE)
#'
#'
#'
#' @param mOccur gg
#' @param fct gf
#' @param xpr gg
#' @param opt.meth gg
#' @param opt.mean gg
#' @param opt.mod gg
#' @param opt.prn gg
#'
#' @return gg
#' @importFrom clusterCrit extCriteria
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

cluster_elements_LOO <- function(mOccur, fct,
                                 xpr      = rep(1, length(fct)),

                                 # options for computing
                                 opt.meth = "divisive",
                                 opt.mean = "amean",
                                 opt.mod  = "byelt",
                                 # options for plotting
                                 opt.prn  = TRUE)   {

  nbAss <- dim(mOccur)[1]
  nbElt <- dim(mOccur)[2]

  tree.ref <- cluster_elements(mOccur, fct, xpr,
                               opt.meth, opt.mean, opt.mod, opt.prn)

  mJaccard <- mCor <- matrix(0, nrow = nbAss, ncol = nbElt,
                             dimnames = list(rownames(mOccur), seq_len(nbElt)))

  ass <- 1
  for (ass in seq_len(nbAss)) {

    tree <- cluster_elements(mOccur[-ass, ], fct[-ass], xpr[-ass],
                             opt.meth, opt.mean, opt.mod, opt.prn)

    mCor[ass, ] <- tree$cor

    for (lev in 1:(nbElt - 1))
      mJaccard[ass, lev] <- as.vector(
        unlist(clusterCrit::extCriteria(
          as.integer(cut_tree(tree,     lev)),
          as.integer(cut_tree(tree.ref, lev)), crit = "Jaccard")))
  }

  apply(mCor, MARGIN = 2, FUN = mean) / tree.ref$cor

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Look for the elements of which the deletion
#             changes the clustering of other elements at the "nbcl"-level
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title mmm
#' @description fff
#' @usage significant_elements(tree.ref, mOccur, fct, nbOpt,
#'                             xpr      = rep(1, length(fct)),
#'                             opt.ind  = "Jaccard",
#'                             opt.meth = "divisive",
#'                             opt.mean = "amean",
#'                             opt.mod  = "byelt",
#'                             opt.plot = TRUE)
#'
#' @param tree.ref ff
#' @param mOccur gg
#' @param fct gg
#' @param nbOpt gg
#' @param xpr gg
#' @param opt.ind gg
#' @param opt.meth gg
#' @param opt.mean gg
#' @param opt.mod gg
#' @param opt.plot gg
#' @details "Jaccard",
#' opt.ind can be : "Czekanowski_Dice", "Folkes_Mallows", "Hubert", "Jaccard",
#' "Kulczynski", "McNemar", "Phi", "Precision", "Rand", "Recall",
#' "Rogers_Tanimoto", "Russel_Rao", "Sokal_Sneath1", "Sokal_Sneath2"
#' @importFrom clusterCrit getCriteriaNames extCriteria
#' @return gg
#' @export


significant_elements <- function(tree.ref, mOccur, fct, nbOpt,

                                 xpr      = rep(1, length(fct)),

                                 opt.ind  = "Jaccard",
# opt.ind can be : "Czekanowski_Dice", "Folkes_Mallows", "Hubert", "Jaccard",
# "Kulczynski", "McNemar", "Phi", "Precision", "Rand", "Recall",
# "Rogers_Tanimoto", "Russel_Rao", "Sokal_Sneath1", "Sokal_Sneath2"

                                 opt.meth = "divisive",
                                 opt.mean = "amean",
                                 opt.mod  = "byelt",
                                 opt.plot = TRUE) {


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  nbElt <- nbClu <- dim(tree.ref$aff)[2]
  setElt         <- seq_len(nbElt)
  names(setElt)  <- colnames(tree.ref$aff)

  lAffectElt <- vector(mode = "list", length = nbElt)
  for (nb in seq_len(nbClu)) lAffectElt[[nb]] <- cut_tree(tree.ref, nb)

  if (opt.plot == TRUE) plot_tree(tree.ref,
                                  col = couleurs[cut_tree(tree.ref, nbOpt)],
                                  titre = "tree.ref")

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  critNames <- clusterCrit::getCriteriaNames(FALSE)          #external criteria
  vCriteria <- numeric(length(critNames))
  names(vCriteria)  <- c(critNames)

  .aff      <- array(0, dim = c(nbClu, nbElt, length(opt.ind)))
  dimnames(.aff) <- list(seq_len(nbClu), colnames(mOccur), opt.ind)

  for (elt in seq_len(nbElt)) {

    element  <- setElt[elt]
    indElt   <- setdiff(setElt, element)
    tmp      <- rm_dual_assemblies(mOccur[ ,indElt], fct, xpr)

    tree     <- cluster_elements(tmp$mat, tmp$fct, xpr,
                                 opt.meth, opt.mean, opt.mod, opt.prn = FALSE)
    prd      <- predict_function(tree, tmp$mat, tmp$fct, opt.mean)
    tree$cor <- prd$tStats[ ,"R2cal"]

    for (nb in seq_len(nbClu - 1)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
                     as.integer(cut_tree(tree, nb) ),
                     as.integer(lAffectElt[[nb]][indElt]),
                     crit = "all")))
      .aff[nb, element, opt.ind] <- vCriteria[opt.ind]
    }

    if (opt.plot == TRUE) {
      plot_tree(tree, col = couleurs[cut_tree(tree, nbOpt)], titre = "")
      text(x = 3, y = 0, labels = names(element))
    }
  }

  coord <- which(.aff == "NaN", arr.ind = TRUE)
  if (amean(.aff[nbClu, , ]) < 0.1) {
    .aff[coord] <- 0
  } else {
    .aff[coord] <- 1
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  if (length(opt.ind) == 1) {
    res        <- list(.aff[ , , 1], tree.ref$cor)
    names(res) <- c("aff", "cor")
  } else {
    res <- vector(mode = "list", length = length(opt.ind))
    names(res) <- opt.ind
    for (ind in seq_along(opt.ind)) {
      res[[ind]]        <- list(.aff[ , , opt.ind[ind]], tree.ref$cor)
      names(res[[ind]]) <- c("aff", "cor")
    }
  }

  return(res)
}



#' @title dd
#' @description dd
#' @usage significant_assemblies(tree.ref, mOccur, fct, nbOpt,
#' xpr      = rep(1, length(fct)),
#' opt.ind  = "Jaccard",
#' opt.meth = "divisive",
#' opt.mean = "amean",
#' opt.mod  = "byelt",
#' opt.plot = TRUE)
#'
#' @param tree.ref  gghhh
#' @param mOccur gghhh
#' @param fct gghhh
#' @param nbOpt gghhh
#' @param xpr gghhh
#' @param opt.ind gghhh
#' @param opt.meth gghhh
#' @param opt.mean gghhh
#' @param opt.mod gghhh
#' @param opt.plot gghhh
#'
#' @return gg
#' @importFrom clusterCrit getCriteriaNames extCriteria
#' @export
#'

significant_assemblies <- function(tree.ref, mOccur, fct, nbOpt,

                                   xpr      = rep(1, length(fct)),

                                   opt.ind  = "Jaccard",
                                   opt.meth = "divisive",
                                   opt.mean = "amean",
                                   opt.mod  = "byelt",
                                   opt.plot = TRUE) {

  # opt.ind can be : "Czekanowski_Dice", "Folkes_Mallows", "Hubert", "Jaccard",
  # "Kulczynski", "McNemar", "Phi", "Precision", "Rand", "Recall",
  # "Rogers_Tanimoto", "Russel_Rao", "Sokal_Sneath1", "Sokal_Sneath2"


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  setXpr <- unique(xpr)
  setAss <- which(xpr == setXpr[1])
  names(setAss) <- rownames(mOccur[setAss, ])
  nbAss  <- length(setAss)

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  nbClu <- nbElt <- dim(tree.ref$aff)[1]
  lAffectElt <- vector(mode = "list", length = nbClu)
  for (nb in seq_len(nbClu)) lAffectElt[[nb]] <- cut_tree(tree.ref, nb)

  if (opt.plot == TRUE) plot_tree(tree.ref,
                                  col = couleurs[cut_tree(tree.ref, nbOpt)],
                                  titre = "tree.ref")

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  critNames <- clusterCrit::getCriteriaNames(FALSE)          #external criteria
  vCriteria <- numeric(length(critNames))
  names(vCriteria)  <- c(critNames)

  .aff      <- array(0, dim = c(nbClu, nbAss, length(opt.ind)))
  dimnames(.aff) <- list(seq_len(nbClu), names(setAss), opt.ind)

  for (ass in seq_along(setAss)) {

    indAss <- ass + nbAss * (seq_along(setXpr) - 1)

    tree   <- cluster_elements(mOccur[-indAss, ], fct[-indAss], xpr[-indAss],
                               opt.meth, opt.mean, opt.mod, opt.prn = FALSE)
    prd    <- predict_function(tree, mOccur[-indAss, ], fct[-indAss], opt.mean)
    tree$cor <- prd$tStats[ ,"R2cal"]

    for (nb in seq_len(nbClu)) {
      vCriteria[] <- as.vector(unlist(clusterCrit::extCriteria(
        as.integer(cut_tree(tree, nb) ),
        as.integer(lAffectElt[[nb]]),
        crit = "all")))
      .aff[nb, ass, opt.ind] <- vCriteria[opt.ind]
    }

    if (opt.plot == TRUE) {
      plot_tree(tree, col = couleurs[cut_tree(tree, nbOpt)], titre = "")
      text(x = 3, y = 0, labels = names(setAss)[ass])
    }
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  if (length(opt.ind) == 1) {
    res        <- list(.aff[ , , 1], tree.ref$cor)
    names(res) <- c("aff", "cor")
  } else {
    res <- vector(mode = "list", length = length(opt.ind))
    names(res) <- opt.ind
    for (ind in seq_along(opt.ind)) {
      res[[ind]]        <- list(.aff[ , , opt.ind[ind]], tree.ref$cor)
      names(res[[ind]]) <- c("aff", "cor")
    }
  }

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   External function for clustering of ASSEMBLY MOTIFS
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#               CLUSTERING of MOTIFS of ASSEMBLIES
#
# Input: mOccur: matrix of occurence of nbElt elements within the assemblies
#        fct : observed performances of assemblies
#        affectElt : affectation of elements into classes
#        nbclMotMax : maximum number of classes of assembly motifs
#        option = "divisive "or "agglomerative"
#        method = "sort", "cluster" or "tree"
#
# if option="divisive" and method="sort": we recommand to left nbclMotMax="NA"
# if option="divisive" and method="cluster": nbclMotMax is important
# if option="divisive" and method="tree": we recommand to left nbclMotMax="NA"
# if option="agglomerative": method="tree" and nbclMotMax="NA"
#          ( nbclMotMax="NA" induces that nbclMotMax=maximum of motif number)
#
# Outputs : elt : successive affectations of elements,
#           cor : the corresponding correlations
# if option="divisive" and method="sort" or "cluster":
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title kkk
#' @description fff
#' @usage cluster_motifs(affectElt, mOccur, fct,
#'                       opt.meth   = "divisive",
#'                       nbclMotMax = "NA",
#'                       opt.prn    = FALSE       )
#'
#'
#' @param affectElt gg
#' @param mOccur gg
#' @param fct gg
#' @param opt.meth ff
#' @param nbclMotMax gg
#' @param opt.prn gg
#'
#' @return gg
#' @export
#'


cluster_motifs <- function(affectElt, mOccur, fct,
                           opt.meth   = "divisive",
                           nbclMotMax = "NA",
                           opt.prn    = FALSE       ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Inputs checking
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  tmp <- names(table(mOccur))
  if ((length(tmp) != 2) | (tmp[1] != "0") | (tmp[2] != "1"))
              stop("The matrix of occurence should be a binary matrix")

  if (length(fct) != dim(mOccur)[1])
                stop("dim(mOccur)[1] should be equal to length(fct)")
  fpred <- fct

  if (length(affectElt) != dim(mOccur)[2])
                stop("dim(mOccur)[2] should be equal to length(affectElt)")

  tmp <- sort(unique(affectElt))
  for (i in 1:length(tmp)) affectElt[affectElt == tmp[i]] <- i

  if ( (opt.meth != "sort") && (opt.meth != "divisive") &&
       (opt.meth != "agglomerative") )
    stop(" 'opt.meth' can be 'sort', 'divisive' or 'agglomerative' ")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # First line : Trunk
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  AssMotifs <- affect_motifs(affectElt, mOccur)
  setMots   <- unique(AssMotifs)
  nomMots   <- unique(names_assembly(affectElt, mOccur)$motifs)
  nbMot     <- length(setMots)

  if (!is.numeric(nbclMotMax)) {
      nbclMotMax <- nbMot
  } else {
      if (nbclMotMax < 2)     nbclMotMax <- 2
      if (nbclMotMax > nbMot) nbclMotMax <- nbMot
  }

  res <- RES.aff <- matrix(0, nrow = nbMot, ncol = nbMot)
  rownames(res)  <- rownames(RES.aff) <- c(1:nbMot)
  colnames(res)  <- colnames(RES.aff) <- nomMots

  RES.rss        <- numeric(nbMot)
  fctPrd  <- numeric(length(fct))



  if (opt.prn)
        print(paste("Clustering of ASSEMBLY MOTIFS: option = '", opt.meth,
                    sep = ""))


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Choice of methods
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  myfct <- function(affectMot) {

    for (motif in unique(affectMot)) {
      msk  <- setMots[(affectMot == motif)]
      mask <- NULL
      for (i in seq_along(msk)) mask <- c(mask, which(AssMotifs == msk[i]))
      fctPrd[mask] <- amean(fct[mask])
    }

    return( rss(fctPrd, fct) )
  }

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # Choice of methods
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  if (opt.meth == "sort") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #               One-by-one Divisive SORTING of SPECIES
    #
    # Initial: all elements are gathered into a single big cluster. R2 is null.
    # Iteration:
    #      each element belonging to the big cluster is successively removed,
    #      then put as a new cluster containing only a single species.
    # Each new clustering is evaluated by R2
    #                                    of resulting community classification.
    # The most likely clustering is kept.
    # A new element is isolated in a singleton.
    # The process stops when each element is isolated into a singleton.
    #
    # Input : matrix of occurence
    #         the corresponding vector of observations
    # Outputs : successive affectations,
    #           the corresponding correlations
    #           the whole used matrix of correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    nbclMot            <- 1
    affectMot          <- rep(1, nbMot)
    RES.aff[nbclMot, ] <- affectMot
    RES.rss[nbclMot]   <- rss(amean(fct), fct)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    res[ , ] <- oldRes <- Inf

    if (nbMot > 1) for (nbclMot in 2:nbMot) {

#      nbclMot <- nbclMot + 1

      oldAffect <- affectMot
      setMot    <- which(affectMot == 1)

      # Test the elements one-after-one into a new singleton
      for (mot in setMot) {
        affectMot[mot]    <- nbclMot
        res[nbclMot, mot] <- myfct(affectMot)
        affectMot         <- oldAffect
      }

      # Decision: one keeps the best fit
      minRes         <- min(res[nbclMot, ])
      mot            <- which(res[nbclMot,] == minRes, arr.ind = TRUE)[1]
      affectMot[mot] <- nbclMot

      RES.aff[nbclMot, ] <- affectMot
      RES.rss[nbclMot]   <- minRes

      # print to follow the computation
      if (opt.prn) print(cbind(RES.aff[1:nbclMot, ], RES.rss[1:nbclMot]))

    }
    #  End of the FOR

  }
  #    END of IF on DIVISIVE SORTING of MOTIFS
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



  if (opt.meth == "divisive") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #    DIVISIVE HIERARCHICAL CLUSTERING of ELEMENTS by DICHOTOMY (ROUX, 2005)
    #
    # Initial: all elements are gathered into a single big cluster. R2 is null.
    # Iteration:
    #   each element belonging to the big cluster is successively removed,
    #            then put in other clusters while R2 increased.
    #
    # Input : matrix of occurence,
    #         the corresponding vector of observations,
    #         a maximum number of element clusters.
    # Outputs : successive affectations,
    #           corresponding correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    nbclMot            <- 1
    affectMot          <- rep(1, nbMot)
    RES.aff[nbclMot, ] <- affectMot
    RES.rss[nbclMot]   <- rss(amean(fct), fct)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines : Leaves
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    res[ , ] <- oldRes <- Inf
    oldRes.affect <- affectMot

#    nbclMot <- 0

    if (nbMot > 1) for (nbclMot in 1:(nbMot - 1)) {

#      nbclMot <- nbclMot + 1
      oldRes <- res[ , ] <- Inf

      # One tests split of each leaves of tree
      for (clMot in 1:nbclMot) {

        affectMot <- RES.aff[nbclMot, ]
        setMot    <- which(affectMot == clMot)

        test <- (length(setMot) > 1)
        while (test) {

          # One tests elements one-after-one
          oldAffect <- affectMot
          for (mot in setMot) {
            affectMot[mot]  <- nbclMot + 1
            res[clMot, mot] <- myfct(affectMot)
            affectMot       <- oldAffect
          }

          # Decision: one keeps the best local fit
          minRes <- min(res)
          if (minRes < oldRes) {
            coord                 <- which(res == minRes, arr.ind = TRUE)
            affectMot[coord[1,2]] <- nbclMot + 1
            oldRes                <- minRes
            oldRes.affect         <- affectMot
          } else {
            test <- FALSE
          }
        }
        # END of WHILE on test
      }
      # END of LOOP on LEAVE

      affectMot              <- oldRes.affect
      RES.aff[nbclMot + 1, ] <- affectMot
      RES.rss[nbclMot + 1]   <- oldRes

      # print to follow the computation
      if (opt.prn) print(cbind(RES.aff[1:(nbclMot + 1), , drop = FALSE],
                               RES.rss[1:(nbclMot + 1)]))
    }
    # END of MAIN WHILE

  }
  # END of DIVISIVE HIERARCHICAL CLUSTERING of MOTIFS by DICHOTOMY
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



  if (opt.meth == "agglomerative") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #   Single linkage AGGLOMERATIVE HIERARCHICAL CLUSTERING of ASSEMBLY MOTIFS
    #              according to Legendre et Legendre, 1998
    #
    # Initial: each assembly motif is a singleton.
    # Iteration: the assembly motifs are gathered two-by-two
    #            until  they are all togetger in a big cluster.
    #
    # Input : matrix of occurence,
    #         the corresponding vector of observations,
    #         a vector of element affectation.
    # Outputs : successive affectations,
    #           corresponding correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    nbclMot   <- nbMot
    affectMot <- c(1:nbMot)

    RES.aff[nbclMot, ] <- affectMot
    RES.rss[nbclMot]   <- oldRes <- myfct(affectMot)

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines : branchs and trunc
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Test adding on a element to each assembly class
    if (nbMot > 1) for (nbclMot in (nbMot - 1):2) {

#      nbclMot <- nbclMot - 1

      res[ , ] <- Inf
      setMot   <- sort(unique(affectMot))

      iter <- 0
      set1 <- setdiff(setMot, max(setMot))
      for (mot1 in set1) {

        iter <- iter + 1
        set2 <- setMot[setMot > mot1]
        oldAffect <- affectMot
        # Test the elements one-after-one
        for (mot2 in set2) {
          affectMot[affectMot == mot2] <- mot1
          res[mot1, mot2]              <- myfct(affectMot)
          affectMot                    <- oldAffect
        }
        # print to follow the computation
        if (opt.prn) print(res[1:iter, ])
      }

      # Decision: one keeps the best fit
      affectMot <- RES.aff[nbclMot + 1, ] # ou mot + 1 ou mot - 1

      minRes <- min(res)
      coord  <- which(res == minRes, arr.ind = TRUE)

      affectMot[affectMot == coord[1, 2]] <- coord[1, 1]

      RES.aff[nbclMot, ] <- affectMot
      RES.rss[nbclMot]   <- minRes

      # print to follow the computation
      if (opt.prn) print(cbind(RES.aff[nbclMot:nbMot, ],
                               RES.rss[nbclMot:nbMot]))
    }
    # END of LOOP

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Trunk
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    nbclMot   <- 1
    affectMot <- rep(1, nbMot)

    RES.aff[nbclMot, ] <- affectMot
    RES.rss[nbclMot]   <- myfct(affectMot)

    # print to follow the computation
    if (opt.prn) print(cbind(RES.aff[nbclMot:nbMot, ],
                             RES.rss[nbclMot:nbMot]))

    # re-organize the cluster notation
    tmp   <- RES.aff
    index <- numeric(nbMot)
    for (nbcl in 1:nbMot)
      index[nbcl] <- setdiff(unique(RES.aff[nbcl, ]), index[index != 0])
    for (nbcl in 1:nbMot) RES.aff[which(tmp == nbcl)] <- which(index == nbcl)

    # print to follow the computation
    if (opt.prn) print(cbind(RES.aff[nbclMot:nbMot, ],
                             RES.rss[nbclMot:nbMot]))
  }
  #   END of IF on AGGLOMERATIVE HIERARCHICAL CLUSTERING of ELEMENTS
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


  if (opt.meth == "cluster") {
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #
    #      One-by-one Divisive CLUSTERING of SPECIES (McNaughton-Smith, 1964)
    #
    # Initial: all elements are gathered into a single big cluster. R2 is null.
    # I teration:
    # each element belonging to the big cluster is successively removed,
    #            then put in other clusters or a new one.
    # Each new clustering is evaluated by R2
    #                                    of resulting community classification.
    # The most likely clustering is kept. A new element is clustered.
    # The process runs while R2-value increases.
    #
    # Input : matrix of occurence,
    #         the corresponding vector of observations,
    #         a maximum number of element clusters.
    # Outputs : successive affectations,
    #           corresponding correlations
    #           the whole used matrix of correlations.
    #
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    nbclMot   <- 1
    affectMot <- rep(1,nbMot)
    oldRes    <- 0

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Following lines
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    test <- TRUE
    while (test)
    {
      res[,]  <- 0
      setMot  <- which(affectMot == 1)
      nbclMot <- min(nbclMot + 1, nbclMotMax)

      # Test each element in each cluster
      if (nbMot > 1) for (clMot in 2:nbclMot)
      {
        for (mot in setMot)
        {
          oldAffect      <- affectMot
          affectMot[mot] <- clMot

          for (motif in sort(unique(affectMot)))
          {
            msk  <- setMots[(affectMot == motif)]
            mask <- NULL
            for (i in 1:length(msk))
              mask <- c(mask, which(AssMotifs == msk[i]))
            fpred[mask] <- mean(fct[mask])
          }
          res[clMot,mot] <- cor2(fpred, fct)

          affectMot      <- oldAffect
        }
      }

      # Decision: one keeps the best fit
      maxRes <- max(res)
      if (maxRes > oldRes)
      {
        coord                 <- which(res == maxRes, arr.ind = TRUE)
        affectMot[coord[1,2]] <- coord[1,1]
        oldRes                <- maxRes

        RES.aff[nbclMot, ] <- affectMot
        RES.rss[nbclMot]   <- maxRes
      } else test <- FALSE

      # print to follow the computation
     if (opt.prn) print(cbind(RES.aff, RES.rss))
    }
    #  End of the WHILE

    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # Sort the results
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #    index   <- sort.cluster(RES.aff, RES.wrk)
    #    RES.aff <- index$aff

  }
  #    END of IF on DIVISIVE CLUSTERING of ELEMENTS
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Compute pval1: the p.value that R2 are significantly different from 0.000
  #         pval2: the p.value that two close R2 are significantly different
  #         pval3: the p.value that R2 are significantly different from R2max
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  RES.cor <- 1 - RES.rss / rss(fct, amean(fct))

  res        <- list(RES.aff, RES.cor)
  names(res) <- c("aff", "cor")

  return(res)
}

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   END of CLUSTER.ASSEMBLY.MOTIFS
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                   END of FILE myCLUSTERING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
