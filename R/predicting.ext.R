#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Set of external functions                                                ####
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#' @include stats.int.R predicting.int.R clustering.int.R
NULL


#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Functions for computing statistics                                       ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Compute the matrix 3 x nbclusters of R2cal, R2prd, and missing values
#'
#' @description Take a numeric f.
#'
#' @usage compute_fit_stats(mCal, mPrd, fct, nbK)
#'
#' @param mCal gghhh
#'
#' @param mPrd gghhh
#'
#' @param fct gghhh
#'
#' @param nbK gghhh
#'
#' @details dd
#'
#' @return gghhh
#'
#' @keywords internal
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute_fit_stats <- function(mCal, mPrd, fct, nbK) {

  nbclMax <- dim(mCal)[1]

  tmp     <- c("missing", "R2cal", "R2prd", "pcal", "pprd", "AIC", "AICc")
  mStats  <- matrix(0, nrow = nbclMax, ncol = length(tmp),
                    dimnames = list(seq(1, nbclMax), tmp))

  for (nbcl in seq_len(nbclMax)) {

    mStats[nbcl, "missing"] <- sum(is.na(mPrd[nbcl, ])) / length(mCal[nbcl, ])
    mStats[nbcl, "R2cal"] <- R2mse(mCal[nbcl, ], fct)
    mStats[nbcl, "R2prd"] <- R2mse(mPrd[nbcl,],  fct)
    mStats[nbcl, "AIC"]   <- AIC(  mCal[nbcl, ], fct, nbK[nbcl])
    mStats[nbcl, "AICc"]  <- AICc( mCal[nbcl, ], fct, nbK[nbcl])

    if (nbK[nbcl] > 1) {
      mStats[nbcl, "pcal"] <- pmse( mCal[nbcl, ], fct, nbK[nbcl])
      mStats[nbcl, "pprd"] <- pmse( mPrd[nbcl,],  fct, nbK[nbcl])
    }
  }

  return(mStats)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title  Compute and plot the best Prediction
#'                        as a concatenation of the best predictions
#'
#' @description Take a numeric f.
#'
#' @usage compute_tree_stats(mCal, mPrd, mStats, fct, nbK)
#'
#' @param mCal gghhh
#'
#' @param mPrd gghhh
#'
#' @param mStats gghhh
#'
#' @param fct gghhh
#'
#' @param nbK gghhh
#'
#' @details dd
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute_tree_stats <- function(mCal, mPrd, mStats, fct, nbK) {

  nbclMax <- dim(mCal)[1]
  nbAss   <- dim(mCal)[2]

  tCal <- tPrd <- tNbcl <- matrix(NA, nrow = nbclMax, ncol = nbAss)

  for (nbcl in seq_len(nbclMax)) {
    for (ass in seq_len(nbAss)) {
      last  <- min(sum(!is.na(mPrd[1:nbcl, ass])),
                   sum(!is.na(mStats[1:nbcl, "R2prd"])),
                   nbcl)
      index <- first_argmax(mStats[1:last, "R2prd"])

      tNbcl[nbcl, ass] <- index
      tPrd[nbcl, ass]  <- mPrd[index, ass]
      tCal[nbcl, ass]  <- mCal[index, ass]
    }
  }
  tStats <- compute_fit_stats(tCal, tPrd, fct, nbK)


  if (dim(tStats)[1] > 1) for (nbcl in 2:nbclMax)
    if (tStats[nbcl, "R2prd"] <= tStats[nbcl - 1, "R2prd"]) {
      tStats[nbcl, "R2prd"] <- tStats[nbcl - 1, "R2prd"]
      tNbcl[nbcl, ]         <- tNbcl[nbcl - 1, ]
      tPrd[nbcl, ]          <- tPrd[nbcl - 1, ]
      tCal[nbcl, ]          <- tCal[nbcl - 1, ]
    }
  tStats <- compute_fit_stats(tCal, tPrd, fct, nbK)

  res        <- list(tCal, tPrd, tStats, tNbcl)
  names(res) <- c("tCal", "tPrd", "tStats", "tNbcl")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title  Compute the statistiques of each motifs
#'
#' @description Take a numeric f.
#'
#' @usage compute_motif_stats(tCal, tPrd, mMotifs)
#'
#' @param tCal gghhh
#'
#' @param tPrd gghhh
#'
#' @param mMotifs gghhh
#'
#' @details dd
#'
#' @return gghhh
#'
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

compute_motif_stats <- function(tCal, tPrd, mMotifs) {

  nbclMax <- dim(tCal)[1]

  # build the general structure of motif.Stats

  uTab <- list()
  for (nbcl in 1:nbclMax) uTab[[nbcl]] <- table(mMotifs[nbcl, ])


  # compute the values of motif.Stats

  uMean <- uSd <- uRmse <- uR2 <- uSlope <- uTab

  for (nbcl in 1:nbclMax)
    for (motif in names(uTab[[nbcl]])) {

      index <- which(mMotifs[nbcl, ] == motif)

      uMean[[nbcl]][motif] <- amean(tPrd[nbcl, index])

      if (length(index) > 1) {

        uSd[[nbcl]][motif]    <- asd(tPrd[nbcl, index])
        uRmse[[nbcl]][motif]  <- rmse(tPrd[nbcl, index], tCal[nbcl, index])
        uR2[[nbcl]][motif]    <- R2mse(tPrd[nbcl, index], tCal[nbcl, index])
        uSlope[[nbcl]][motif] <-
          stats::lm(tPrd[nbcl, index] ~ tCal[nbcl, index])$coef[2]

      } else {

        oldMotif <- mMotifs[nbcl - 1, index]

        uSd[[nbcl]][motif]    <- uSd[[nbcl - 1]][oldMotif]
        uRmse[[nbcl]][motif]  <- uRmse[[nbcl - 1]][oldMotif]
        uR2[[nbcl]][motif]    <- uR2[[nbcl - 1]][oldMotif]
        uSlope[[nbcl]][motif] <- uSlope[[nbcl - 1]][oldMotif]
      }
    }

  res        <- list(uTab, uMean, uSd, uRmse, uR2, uSlope)
  names(res) <- c("uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  return(res)
}




#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Public functions                                                         ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title Compute assembly functioning and associated statistiques
#'          using the Clustering model
#'
#' @description Take a numeric f.
#'
#' @usage predict_function(tree.cal, mOccur, fct,
#'                         xpr = rep(1, length(fct)),
#'                         opt.mean = "amean",
#'                         opt.mod  = "byelt",
#'                         opt.jack = FALSE,  jack = c(2,5))
#'
#' @param tree.cal gghhh
#'
#' @param mOccur gghhh
#'
#' @param fct gghhh
#'
#' @param xpr gghhh
#'
#' @param opt.mean gghhh
#'
#' @param opt.mod gghhh
#'
#' @param opt.jack gghhh
#'
#' @param jack gghhh
#'
#' @details dd
#'
#' @return gghhh
#'
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_function <- function(tree.cal, mOccur, fct,
                             xpr = rep(1, length(fct)),
                             # options for computing
                             opt.mean = "amean",
                             opt.mod  = "byelt",
                             opt.jack = FALSE,  jack = c(2,5)) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbElt      <- dim(mOccur)[2]
  nbAss      <- length(fct)

  mAssMotifs <- maffect_motifs(tree.cal, mOccur)

  # Compute the cross-validation predictions
  mCal <- mPrd <- matrix(NA, nrow = nbElt, ncol = nbAss,
                         dimnames = list(seq_len(nbElt), rownames(mOccur)))

  for (nbcl in seq_len(nbElt)) {
    assMotifs    <- mAssMotifs[nbcl, ]
    mCal[nbcl, ] <- predict_cal(fct, assMotifs, mOccur, xpr, opt.mean, opt.mod)
    mPrd[nbcl, ] <- predict_prd(fct, assMotifs, mOccur, xpr,
                                opt.mean, opt.mod, opt.jack, jack)
  }

  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mAssMotifs[nbcl, ]))

  mStats <- compute_fit_stats(mCal, mPrd, fct, nbK)

  # Compute the Global Predictions for all the number of clusters
  tree.prd <- compute_tree_stats(mCal, mPrd, mStats, fct, nbK)

  uStats   <- compute_motif_stats(tree.prd$tCal, tree.prd$tPrd, mAssMotifs)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  res <- list(rownames(mOccur), fct, xpr, opt.mean, opt.mod,
              mCal, mPrd, mStats, mAssMotifs,
              tree.prd$tCal, tree.prd$tPrd, tree.prd$tStats, tree.prd$tNbcl,
              uStats$uTab, uStats$uMean, uStats$uSd,
              uStats$uRmse, uStats$uR2, uStats$uSlope)


  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
                  "mCal", "mPrd", "mStats", "mMotifs",
                  "tCal", "tPrd", "tStats", "tNbcl",
                  "uTab", "uMean", "uSd", "uRmse", "uR2", "uSlope")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute assembly functioning and associated statistiques
#          using the Clustering model
#          for Interaction effet (alpha), Composition effect (beta),
#            and assembly Functioning as Alpha * Beta * Fscale
#'
#' @title Compute assembly functioning and associated statistiques
#'          using the Clustering model
#'          for Interaction effet (alpha), Composition effect (beta),
#'            and assembly Functioning as Alpha * Beta * Fscale
#'
#' @description Take a numeric f.
#'
#' @usage predict_twin_fct(tree.cal, mOccur, alpha, beta,
#'                         xpr       = rep(1, length(alpha)),
#'                         fscale    = 1,
#'                         titre     = "",
#'                         opt.alpha = "gmean",
#'                         opt.beta  = "amean",
#'                         opt.mod   = "byelt",
#'                         opt.jack  = FALSE,    jack = c(2,5) )
#'
#' @param tree.cal gghhh
#' @param mOccur gghhh
#' @param alpha gghhh
#' @param beta gghhh
#' @param xpr gghhh
#' @param fscale gghhh
#' @param titre gghhh
#' @param opt.alpha gghhh
#' @param opt.beta gghhh
#' @param opt.mod gghhh
#' @param opt.jack gghhh
#' @param jack gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_twin_fct <- function(tree.cal, mOccur, alpha, beta,
                             xpr       = rep(1, length(alpha)),
                             fscale    = 1,
                             titre     = "",
                             opt.alpha = "gmean",
                             opt.beta  = "amean",
                             opt.mod   = "byelt",
                             opt.jack  = FALSE,    jack = c(2,5) ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  fct     <- alpha * beta * fscale

  mAssMotifs <- maffect_motifs(tree.cal, mOccur)

  nbElt   <- dim(mOccur)[2]
  nbAss   <- length(fct)
  nbElt   <- dim(mAssMotifs)[1]

  # Compute the cross-validation predictions
  mCal <- mPrd <- matrix(NA, nrow = nbElt, ncol = nbAss,
                         dimnames = list(seq_len(nbElt), rownames(mOccur)))

  for (nbcl in seq_len(nbElt)) {
    assMotifs <- mAssMotifs[nbcl, ]

    mCal[nbcl, ] <- fscale *
      predict_cal(alpha, assMotifs, mOccur, xpr, opt.alpha, opt.mod) *
      predict_cal(beta,  assMotifs, mOccur, xpr, opt.beta,  opt.mod)

    mPrd[nbcl, ] <- fscale *
      predict_prd(alpha, assMotifs, mOccur, xpr,
                  opt.alpha, opt.mod, opt.jack, jack) *
      predict_prd(beta, assMotifs, mOccur, xpr,
                  opt.beta, opt.mod, opt.jack, jack)
  }

  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mAssMotifs[nbcl, ]))
  mStats <- compute_fit_stats(mCal, mPrd, fct, nbK)

  # Compute the Global Predictions for all the number of clusters
  tree.prd <- compute_tree_stats(mCal, mPrd, mStats, fct, nbK)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  opt.mean <- paste(opt.beta)

  res <- list(rownames(mOccur), fct, xpr, opt.mean, opt.mod,
              mCal, mPrd, mStats, mAssMotifs,
              tree.prd$tCal, tree.prd$tPrd, tree.prd$tStats, tree.prd$tNbcl)

  names(res) <- c("names", "fct", "xpr", "opt.mean", "opt.mod",
                  "mCal", "mPrd", "mStats", "mMotifs",
                  "tCal", "tPrd", "tStats", "tNbcl")

  return(res)
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Compute assembly functioning and associated statistiques
#          using the Clustering model
#          for Supplementary Assemblies by knowing their elemental Composition
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title #  Compute assembly functioning and associated statistiques
#'          using the Clustering model
#'          for Supplementary Assemblies by knowing their elemental Composition
#'
#'
#' @description Take a numeric f.
#'
#' @usage predict_function_supp(tree, res.prd,
#'                              supOccur, appOccur, appFct,
#'                              opt.mean = "amean",
#'                              opt.mod  = "byelt"    )
#'
#' @param tree gghhh
#' @param res.prd gghhh
#' @param supOccur gghhh
#' @param appOccur gghhh
#' @param appFct gghhh
#' @param opt.mean gghhh
#' @param opt.mod gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

predict_function_supp <- function(tree, res.prd,
                                  supOccur, appOccur, appFct,

                                  # options for computing
                                  opt.mean = "amean",
                                  opt.mod  = "byelt"    ) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Main calculations of Cal, Prd, R2cal and R2prd
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #  Cas 1 où on suppose le nombre d'éléments identiques
  #      dans apprentissage et supplémentaires'
  #  Prévoir le cas 2 où le nombre d'éléments est différent
  #      dans apprentissage et supplémentaire'

  nbElt      <- dim(appOccur)[2]    # ????

  nbAppAss   <- dim(appOccur)[1]
  nbSupAss   <- dim(supOccur)[1]

  mAssMotifs <- maffect_motifs(tree, rbind(appOccur, supOccur))
  mAppMotifs <- mAssMotifs[ , 1:nbAppAss]
  mSupMotifs <- mAssMotifs[ , (nbAppAss + 1):(nbAppAss + nbSupAss)]

  # Compute the raw predictions
  mSup <- mSd <- tNbcl <-
    matrix(NA, nrow = nbElt, ncol = nbSupAss,
           dimnames = list(seq_len(nbElt), rownames(supOccur)))

  for (nbcl in seq_len(nbElt)) {

    mSup[nbcl, ] <- predict_supp(appFct, mAppMotifs[nbcl, ], appOccur,
                                         mSupMotifs[nbcl, ], supOccur,
                                 opt.mean, opt.mod)

    mSd[nbcl, ] <- res.prd$uRmse[[nbcl]][mSupMotifs[nbcl, ]]

    index <- which(mSup[nbcl, ] > appFct)
    mSd[nbcl, index] <- -mSd[nbcl, index]
  }


  # Compute the associated statistiques
  tNbcl[1, ] <- 1
  for (nbcl in 2:nbElt) {
    tNbcl[nbcl, ] <- nbcl
    tNbcl[nbcl, is.na(mSup[nbcl, ])] <- tNbcl[(nbcl - 1), is.na(mSup[nbcl, ])]
  }

  tSup <- mSup
  for (nbcl in 2:nbElt)
    tSup[nbcl, is.na(mSup[nbcl, ])] <- tSup[(nbcl - 1), is.na(mSup[nbcl, ])]


  # Compute the associated statistiques
  nbK <- integer(nbElt)
  for (nbcl in seq_len(nbElt))
    nbK[nbcl] <- length(unique(mSupMotifs[nbcl, ]))

  mStats <- compute_fit_stats(mSup, mSup, appFct, nbK)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Outputs
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  #                  , "tError" A vérifier

  res        <- list(rownames(supOccur), appFct, opt.mean, opt.mod,
                     mSupMotifs, tSup, mSd, tNbcl, mStats)
  names(res) <- c("names", "fct", "opt.mean", "opt.mod",
                  "mMotifs", "tSup", "tSd", "tNbcl", "tStats")

  return(res)
}



#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
#
# Plot functions                                                           ####
#
#oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#' @title  Plot the matrix 3 x nbclusters of R2cal, R2prd, and missing values
#' @description Take a numeric f.
#'
#' @usage plot_stats(mStats, nbElt, titre = "")
#'
#' @param mStats gghhh
#' @param nbElt gghhh
#' @param titre gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @importFrom graphics plot points text abline
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_stats <- function(mStats, nbElt, titre = "") {

  nbclMax <- dim(mStats)[1]

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a first graph without any pvaluesx
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx-

  graphics::plot(x = seq_len(nbclMax), y = mStats[ ,"R2cal"],
       xlim = c(1, nbElt),   ylim = c(0,1),
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "R2 (in black), E (in red)",
       xlab = "number of clusters",
       main = titre)

  # plot R2 calibration
  graphics::points(x = seq_len(nbclMax), y = mStats[ ,"R2cal"],
         type = "b", pch = 1, cex = 2, bg = "white", col = "black")

  # plot R2 prediction with pvalues
  graphics::points(x = seq_len(nbclMax), y = mStats[ ,"R2prd"],
         type = "b", pch = 1, cex = 2, bg = "white", col = "red3")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a second graph without any pvalues
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  graphics::plot(x = seq_len(nbclMax), y = mStats[ ,"R2cal"],
       xlim = c(1, nbElt), ylim = c(0,1),
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "Predicting ratio",
       xlab = "number of clusters",
       main = titre)

  # plot predictiong ratio
  graphics::points(x = seq_len(nbclMax), y = 1 - mStats[ ,"missing"],
         type = "b", pch = 0, cex = 2, bg = "white", col = "blue")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # plot a third graph with AIC and AICc
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  ylim <- c(min(mStats[ , "AIC"], mStats[ , "AICc"], na.rm = TRUE),
            max(mStats[1, "AIC"], mStats[1, "AICc"], na.rm = TRUE))

  # plot AIC
  graphics::plot(x = seq_len(nbclMax), y = mStats[ ,"AIC"],
       xlim = c(1, nbElt), ylim = ylim,
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "AIC", xlab = "number of clusters",
       main = titre)

  fct <- mStats[ which(!is.na(mStats[ , "AIC"]) == TRUE), "AIC"]
  graphics::abline(v = first_argmin(fct)[1], col = "green3")
  graphics::points(x = seq_along(fct), y = fct,
         type = "b", pch = 1, cex = 2, bg = "white", col = "blue3")

  # plot AICc
  graphics::plot(x = seq_len(nbclMax), y = mStats[ ,"AIC"],
       xlim = c(1, nbElt), ylim = ylim,
       type = "n", pch = 1, cex = 2, tck = 0.02, las = 1,
       bg = "white", col = "black",
       ylab = "AICc", xlab = "number of clusters",
       main = titre)

  fct <- mStats[ which(!is.na(mStats[ , "AICc"]) == TRUE), "AICc"]
  graphics::abline(v = first_argmin(fct)[1], col = "green3")
  graphics::points(x = seq_along(fct), y = fct,
         type = "b", pch = 1, cex = 2, bg = "white", col = "red3")

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot Predictions vs Observations
#
#  Inputs : Fprd : vector of modelled or predicted values
#           Fobs : vector of observed values
#           assMotifs : vector of motifs of which belong the assemblies
#
#  Options : titre : main titre of Figure
#            opt.reg : technical information on the quality of the regression
#            opt.aov : variance analysis of assembly function by motif
#            pvalue : threshold for the aov analysis
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Plot Predictions vs Observations
#' @description gggg
#' @usage plot_prediction_simple(Fprd, Fobs, assMotifs, nbcl,
#'                               titre    = "",
#'                               xylim    = range(Fobs),
#'                               opt.mean = "amean",
#'                               opt.aov  = FALSE,  pvalue   = 0.025)
#'
#' @param Fprd gghhh
#' @param Fobs gghhh
#' @param assMotifs gghhh
#' @param nbcl gghhh
#' @param titre gghhh
#' @param xylim gghhh
#' @param opt.mean gghhh
#' @param opt.aov gghhh
#' @param pvalue gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @importFrom graphics plot points text abline lines
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_prediction_simple <- function(Fprd, Fobs, assMotifs, nbcl,
                                   titre    = "",
                                   xylim    = range(Fobs),
                                   opt.mean = "amean",
                                   opt.aov  = FALSE,
                                   pvalue   = 0.025) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the inputs and plot the figure
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbAss <- length(Fobs)
  tmp   <- mean_fct(Fobs, opt.mean)

  index <- which(is.na(Fprd) == TRUE)
  if (length(index) == nbAss) stop("the vector Fprd cannot be null")

  if (length(index) > 0) {
    Fprd      <- Fprd[-index]
    Fobs      <- Fobs[-index]
    assMotifs <- assMotifs[-index]
  }

  graphics::plot(x = Fobs, xlab = "Observations", xlim = xylim,
       y = Fprd, ylab = "Predictions",  ylim = xylim,
       main = titre,
       type = "n", tck = 0.02, las = 1)

  graphics::abline(v = tmp, h = tmp, lty = "dotted", col = "blue")

  graphics::lines(x = xylim, y = xylim, lty = "solid", col = "red")

  graphics::points(x = Fobs, y = Fprd,
         pch = figures[assMotifs], col = couleurs[assMotifs],
         bg = "white", cex = 2)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Adding of various useful informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  width <- xylim[2] - xylim[1]

  posX1 <- xylim[1] + 0.75*width
  posX2 <- xylim[1] + 1.02*width
  posX3 <- xylim[1] - 0.02*width
  posX4 <- xylim[1] + 0.02*width
  posX5 <- xylim[1] + 0.25*width

  posY1 <- xylim[1] + 0.05*width
  posY2 <- xylim[1]
  posY3 <- xylim[1] + 1.02*width
  posY4 <- xylim[1] + 0.98*width


  # predicting ratio
  tmp <- paste("predicting = ", sum(!is.na(Fprd)), "/", nbAss, sep = "")
  graphics::text(x = posX1, y = posY1, labels = tmp, col = "red")

  # R2 value
  if ((nbAss - length(index)) > 1) {
    tmp <- paste("R2 = ", signif(R2mse(Fprd, Fobs), digits = 3), sep = "")
    graphics::text(x = posX1, y = posY2, labels = tmp, col = "red")
  }

  # Number of clusters
  graphics::text(paste("Nb clusters = ", nbcl, sep = ""),
       x = posX5, y = posY4, col = "red")


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Optional informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  # option variance analysis

  if (opt.aov == TRUE) {
    setMot <- sort(unique(assMotifs))

    graphics::text(x = posX2, y = posY3, labels = "pred", col = "black")

    test <- test_posthoc(Fprd, assMotifs, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(test$motif == motif)
        graphics::text(x = posX2, y = test[index, "mean"],
             labels = as.character(test[index, "group"]),
             col =  couleurs[motif], font = 3)
      }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot Predictions vs Observations
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title mmm
#' @description  gghhh
#' @usage plot_prediction_LOO(Fcal, Fprd, Fobs, assMotifs, nbcl,
#'                            xylim    = range(Fobs),
#'                            titre    = "",
#'                            opt.mean = "amean",
#'                            opt.aov  = FALSE,  pvalue   = 0.05)
#'
#' @param Fcal gghhh
#' @param Fprd gghhh
#' @param Fobs gghhh
#' @param assMotifs gghhh
#' @param nbcl gghhh
#' @param xylim gghhh
#' @param titre gghhh
#' @param opt.mean gghhh
#' @param opt.aov gghhh
#' @param pvalue gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @importFrom graphics plot points text abline lines
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_prediction_LOO <- function(Fcal, Fprd, Fobs, assMotifs, nbcl,
                                xylim    = range(Fobs),
                                titre    = "",
                                opt.mean = "amean",
                                opt.aov  = FALSE,
                                pvalue   = 0.05) {

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Check the inputs and plot the figure
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbAss <- length(Fobs)
  #  xylim <- c(0.7, 2.5) # pour alpha
  #  xylim <- c(0.45, 1.73) # pour beta
  tmp   <- mean_fct(Fobs, opt.mean)

  index <- which(is.na(Fprd) == TRUE)
  if (length(index) == nbAss) return(FALSE)

  if (length(index) > 0) {
    Fprd      <- Fprd[-index]
    Fcal      <- Fcal[-index]
    Fobs      <- Fobs[-index]
    assMotifs <- assMotifs[-index]
  }

  graphics::plot(x = Fobs, xlab = "Observations", xlim = xylim,
       y = Fcal, ylab = "Predictions",  ylim = xylim,
       main = titre,
       las = 1, type = "n", tck = 0.02)

  for (elt in seq_along(Fobs))
    graphics::lines(x = c(Fobs[elt], Fobs[elt]), y = c(Fprd[elt], Fcal[elt]),
          col = couleurs[assMotifs[elt]], lty = "solid")

  graphics::abline(v = tmp, h = tmp, lty = "dotted", col = "blue")

  graphics::lines(x = xylim, y = xylim, lty = "solid", col = "red")

  graphics::points(x = Fobs, y = Fcal,
         pch = figures[assMotifs], col = couleurs[assMotifs],
         bg = "white", cex = 2)


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Adding of various useful informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  width <- xylim[2] - xylim[1]

  posX1 <- xylim[1] + 0.75*width
  posX2 <- xylim[1] + 1.02*width
  posX3 <- xylim[1] - 0.02*width
  posX4 <- xylim[1] + 0.02*width
  posX5 <- xylim[1] + 0.25*width

  posY1 <- xylim[1] + 0.05*width
  posY2 <- xylim[1]
  posY3 <- xylim[1] + 1.02*width
  posY4 <- xylim[1] + 0.98*width


  # Predicting assemblies

  tmp <- paste("predicting = ", sum(!is.na(Fprd)), "/", nbAss, sep = "")
  graphics::text(x = posX1, y = posY1, labels = tmp, col = "red")

  #  Number of clusters and R2 value

  if ((nbAss - length(index)) > 1) {
    tmp <- paste0("R2 = ",  signif(R2mse(Fcal, Fobs), digits = 3),
                  "  E = ", signif(R2mse(Fprd, Fobs), digits = 3))
    graphics::text(x = posX1, y = posY2, labels = tmp, col = "red")

    tmp <- paste0("Nb clusters = ", nbcl, "  E/R2 = ",
                  signif(R2mse(Fprd, Fobs)/R2mse(Fcal, Fobs), digits = 3))
    graphics::text(tmp, x = posX5, y = posY4, col = "red")
  }


  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  # Optional informations
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  # option variance analysis

  if (opt.aov == TRUE) {
    setMot <- sort(unique(assMotifs))

    graphics::text(x = posX3, y = posY3, labels = "cal", col = "black")
    test   <- test_posthoc(Fprd, assMotifs, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(test$motif == motif)
        graphics::text(x = posX3, y = test[index, "mean"],
             labels = as.character(test[index, "group"]),
             col =  couleurs[motif], font = 3)
      }

    graphics::text(x = posX2, y = posY3, labels = "prd", col = "black")
    test   <- test_posthoc(Fprd, assMotifs, pvalue)
    if (is.list(test))
      for (mot in seq_along(setMot)) {
        motif <- setMot[mot]
        index <- which(test$motif == motif)
        graphics::text(x = posX2, y = test[index, "mean"],
             labels = as.character(test[index, "group"]),
             col =  couleurs[motif], font = 3)
      }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Add the names of assemblies
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Add the names of assemblies
#' @description llll
#' @usage add_ass_names(Fcal, Fprd, Fobs, assMotifs, AssNames)
#'
#' @param Fcal gghhh
#' @param Fprd gghhh
#' @param Fobs gghhh
#' @param assMotifs gghhh
#' @param AssNames gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @importFrom graphics text
#' @keywords internal
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

add_ass_names <- function(Fcal, Fprd, Fobs,
                          assMotifs, AssNames) {

  index <- which(is.na(Fprd) == TRUE)

  if (length(index) != 0) {
    Fprd      <- Fprd[-index]
    Fcal      <- Fcal[-index]
    Fobs      <- Fobs[-index]
    assMotifs <- assMotifs[-index]
    AssNames  <- AssNames[-index]
  }

  thre   <- max(Fobs) - (max(Fobs) - min(Fobs)) / 5
  index1 <- which(Fobs <= thre)
  index2 <- which(Fobs >  thre)

  if (length(index1) != 0)
    graphics::text(x = Fobs[index1], y = Fcal[index1], labels = AssNames[index1],
         col = couleurs[assMotifs[index1]], pos = 4)          # to the right of

  if (length(index2) != 0)
    graphics::text(x = Fobs[index2], y = Fcal[index2], labels = AssNames[index2],
         col = couleurs[assMotifs[index2]], pos = 2)          # to the left of
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot all the necessary predicting figures
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Plot all the necessary predicting figures
#' @description ddd
#' @usage plot_prediction(res,
#'                         xpr.who  = 5,
#'                         nbOpt    = 0,
#'                         titre    = "",
#'                         clu.stat = FALSE,
#'                         clu.cal  = FALSE,
#'                         clu.prd  = FALSE,
#'                         tre.stat = FALSE,
#'                         tre.prd  = FALSE,
#'                         tre.pub  = FALSE,
#'                         tre.best = FALSE,
#'                         tre.opt  = TRUE,
#'                         tre.calvsprd = FALSE,
#'                         ass.loc  = FALSE,
#'                         opt.aov  = FALSE, pvalue = 0.025,
#'                         opt.all  = FALSE  )
#' @param res gg
#' @param xpr.who gg
#' @param nbOpt gg
#' @param titre gg
#' @param clu.stat gg
#' @param clu.cal gg
#' @param clu.prd gg
#' @param tre.stat gghhh
#' @param tre.prd gghhh
#' @param tre.pub gghhh
#' @param tre.best gghhh
#' @param tre.opt gghhh
#' @param tre.calvsprd gghhh
#' @param ass.loc gghhh
#' @param opt.aov gghhh
#' @param pvalue gghhh
#' @param opt.all gghhh
#'
#' @details dd
#'
#' @return gghhh
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_prediction <- function(res,
                            xpr.who  = 5,
                            nbOpt    = 0,

                            titre    = "",

                            clu.stat = FALSE,
                            clu.cal  = FALSE,
                            clu.prd  = FALSE,

                            tre.stat = FALSE,
                            tre.prd  = FALSE,
                            tre.pub  = FALSE,
                            tre.best = FALSE,
                            tre.opt  = TRUE,
                            tre.calvsprd = FALSE,

                            ass.loc  = FALSE,

                            opt.aov  = FALSE, pvalue = 0.05,

                            opt.all  = FALSE  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.all == TRUE) {
    clu.stat <- clu.cal <- clu.prd <-
      tre.stat <- tre.prd <- tre.pub <- tre.best <- tre.opt <-
      tre.calvsprd <- ass.loc <- TRUE
  }

  opt.xpr <- FALSE
  setXpr  <- unique(res$xpr)
  indxpr  <- seq_along(setXpr)
  if (length(indxpr) > 1) {
    opt.xpr <- TRUE
    setAss  <- res$names[which(res$xpr == setXpr[1])]
    pas     <- length(setAss)
  }

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Plot analysis cluster by cluster
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbElt <- dim(res$mStats)[1]
  nbAss <- length(res$fct)

  nbMax <- first_argmax(res$tStats[ ,"R2prd"])
  if (nbOpt == 0) { nbOpt <- first_argmin(res$tStats[ ,"AICc"])
  } else {
    if (nbOpt > nbMax) nbOpt <- nbMax
  }

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the optimum Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.opt == TRUE) {
    nbcl   <- nbOpt
    titre2 <- paste(titre, "Optimum Tree Prediction with", nbcl,
                    "clusters", sep = " ")
    plot_prediction_LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                        rep(1, nbAss), nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Stats
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (clu.stat == TRUE) {
    titre2 <- paste("R2 calibration and prediction:", titre,
                    paste(res$opt.mean, res$opt.mod, sep = "/"), sep = " ")
    plot_stats(res$mStats, nbElt, titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Calibrations for all the number of clusters
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (clu.cal == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Calibration with", nbcl, "clusters", sep = " ")
      plot_prediction_simple(res$mCal[nbcl, ], res$fct,
                             res$mMotifs[nbcl, ], nbcl,
                             xylim    = range(res$fct),
                             opt.mean = res$opt.mean,
                             opt.aov  = opt.aov,
                             pvalue   = pvalue,
                             titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Predictions for all the number of clusters
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (clu.prd == TRUE)
    for (nbcl in seq_len(nbMax))
      if (sum(!is.na(res$mPrd[nbcl, ])) > 0) {
        titre2 <- paste(titre, "Prediction with", nbcl, "clusters", sep = " ")
        plot_prediction_LOO(res$mCal[nbcl, ], res$mPrd[nbcl, ], res$fct,
                            res$mMotifs[nbcl, ], nbcl,
                            opt.mean = res$opt.mean,
                            opt.aov  = opt.aov,
                            pvalue   = pvalue,
                            titre    = titre2)
      }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Stats of Tree predictions
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.stat == TRUE) {
    titre2 <- paste("R2 tree calibration and tree prediction:", titre,
                    paste(res$opt.mean, res$opt.mod, sep = "/"), sep = " ")
    plot_stats(res$tStats, nbElt, titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #    using different colors for predictions with different cluster numbers
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.prd == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Tree Prediction with", nbcl,
                      "clusters", sep = " ")
      plot_prediction_LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                          res$tNbcl[nbcl, ], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      index <- which(res$tNbcl[nbcl, ] == (nbcl - 1))
      if (length(index) != 0)
        add_ass_names(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                      res$fct[index],
                      res$tNbcl[nbcl, index], res$names[index])
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #          with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.pub == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Tree Prediction with", nbcl,
                      "clusters", sep = " ")
      plot_prediction_LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                          rep(1, nbAss), nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the best Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.best == TRUE) {
    nbcl   <- nbMax
    titre2 <- paste(titre, "Best Tree Prediction with", nbcl,
                    "clusters", sep = " ")
    plot_prediction_LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                        rep(1, nbAss), nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the optimum (according to AICc) Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.opt == TRUE) {
    nbcl   <- nbOpt
    titre2 <- paste(titre, "Optimum Tree Prediction with", nbcl,
                    "clusters", sep = " ")
    plot_prediction_LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                        rep(1, nbAss), nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions versus Tree Calibrations for all numbers of clusters
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (tre.calvsprd == TRUE)
    for (nbcl in seq_len(nbMax)) {
      titre2 <- paste(titre, "Tree Prediction vs Calibration with", nbcl,
                      "clusters", sep = " ")
      plot_prediction_simple(res$tPrd[nbcl, ], res$tCal[nbcl, ],
                             res$mMotifs[nbcl, ], nbcl,
                             xylim    = range(res$fct),
                             opt.mean = res$opt.mean,
                             opt.aov  = opt.aov,
                             pvalue   = pvalue,
                             titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions by experiment with the names of Assemblages
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.xpr == TRUE) {
    nbcl   <- nbOpt

    for (ipr in seq_along(indxpr)) {
      titre2 <- paste0(titre,
                       " Tree Predictions for the Experiment ", indxpr[ipr])
      index  <- which(res$xpr == setXpr[ipr])
      plot_prediction_LOO(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                          res$fct[index], res$mMotifs[nbcl, index],
                          nbcl, xylim = range(res$fct),
                          opt.mean = res$opt.mean,
                          titre    = titre2)
      if (ass.loc == TRUE)
        add_ass_names(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                      res$fct[index], res$mMotifs[nbcl, index],
                      res$names[index])
    }
  } else {
    if (ass.loc == TRUE) {
      nbcl   <- nbOpt
      titre2 <- paste(titre, "Tree Prediction with", nbcl,
                      "clusters", sep = " ")
      plot_prediction_LOO(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                          res$mMotifs[nbcl, ], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      add_ass_names(res$tCal[nbcl, ], res$tPrd[nbcl, ], res$fct,
                    res$mMotifs[nbcl, ], res$names)
    }
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions by Assemblage over experiments
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.xpr == TRUE) {

    nbcl   <- nbOpt
    if (is.numeric(xpr.who) == TRUE) {
      setIndex <- sort(sample(1:pas, size = xpr.who[1]))
    } else {
      setIndex <- NULL
      for (elt in seq_along(xpr.who))
        setIndex <- c(setIndex, which(res$names == xpr.who[elt])[1])
    }

    for (ass in setIndex) {
      titre2 <-
        paste0(titre, " Tree Prediction of the Assemblage '",
               setAss[ass], "'")
      index <- ass + ((1:length(indxpr)) - 1) * pas
      plot_prediction_LOO(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                          res$fct[index], res$mMotifs[nbcl, index],
                          nbcl, xylim = range(res$fct),
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      if (ass.loc == TRUE)
        add_ass_names(res$tCal[nbcl, index], res$tPrd[nbcl, index],
                      res$fct[index], res$mMotifs[nbcl, index],
                      res$xpr[index])
    }
  }
}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#           FUNCTIONS OF PREDICTION of A SUPPLEMENTARY DATASET
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#  Plot all the necessary predicting figures
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#' @title Plot all the necessary predicting figures
#' @description ddd
#' @usage plot_prediction_supp(res, nbOpt,
#'                             titre    = "",
#'                             sup.prd  = FALSE,
#'                             sup.pub  = FALSE,
#'                             sup.opt  = TRUE,
#'                             ass.loc  = FALSE,
#'                             opt.all  = FALSE  )
#'
#' @param res gg
#' @param nbOpt gg
#' @param titre gg
#' @param sup.prd gg
#' @param sup.pub gg
#' @param sup.opt gg
#' @param ass.loc gg
#' @param opt.all gg
#'
#' @details dd
#'
#' @return gg
#' @export
#'
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

plot_prediction_supp <- function(res, nbOpt,

                                 titre    = "",

                                 sup.prd  = FALSE,
                                 sup.pub  = FALSE,
                                 sup.opt  = TRUE,
                                 ass.loc  = FALSE,
                                 opt.all  = FALSE  )  {

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Check the inputs
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (opt.all == TRUE)
    sup.prd <- sup.pub <- sup.opt <- ass.loc <- TRUE

  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  #
  #  Plot analysis cluster by cluster
  #
  #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  nbElt  <- dim(res$tStats)[1]
  nbAss  <- length(res$fct)

  indFct <- which(!is.na(res$fct))

  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot the optimum Prediction with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (sup.opt == TRUE) {
    nbcl   <- nbOpt
    titre2 <- paste(titre, "Optimum Prediction of Supplementary data with",
                    nbcl, "clusters", sep = " ")
    plot_prediction_LOO(res$tSup[nbcl, indFct],
                        res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                        res$fct[indFct],
                        rep(1, nbAss)[indFct], nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)
  }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #    using different colors for predictions with different cluster numbers
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (sup.prd == TRUE)
    for (nbcl in seq_len(nbOpt)) {

      titre2 <- paste(titre, "Prediction of Supplementary data with", nbcl,
                      "clusters", sep = " ")
      plot_prediction_LOO(res$tSup[nbcl, indFct],
                          res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                          res$fct[indFct],
                          res$tNbcl[nbcl, indFct], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)

      index <- which(res$tNbcl[nbcl, indFct] == (nbcl - 1))
      if (length(index) != 0) {
        index <- indFct[index]
        add_ass_names(res$tSup[nbcl, index],
                      res$tSup[nbcl, index] + res$tSd[nbcl, index],
                      res$fct[index],
                      res$tNbcl[nbcl, index], res$names[index])
      }
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions for all numbers of clusters
  #          with only one color*symbol
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (sup.pub == TRUE)
    for (nbcl in seq_len(nbOpt)) {

      titre2 <- paste(titre, "Prediction of Supplementary data with", nbcl,
                      "clusters", sep = " ")
      plot_prediction_LOO(res$tSup[nbcl, indFct],
                          res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                          res$fct[indFct],
                          rep(1, nbAss)[indFct], nbcl,
                          opt.mean = res$opt.mean,
                          titre    = titre2)
    }


  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  # Plot Tree Predictions by experiment with the names of Assemblages
  #oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  if (ass.loc == TRUE) {

    nbcl   <- nbOpt
    titre2 <- paste(titre, "Prediction of Supplementary data with", nbcl,
                    "clusters", sep = " ")

    plot_prediction_LOO(res$tSup[nbcl, indFct],
                        res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                        res$fct[indFct],
                        res$mMotifs[nbcl, indFct], nbcl,
                        opt.mean = res$opt.mean,
                        titre    = titre2)

    add_ass_names(res$tSup[nbcl, indFct],
                  res$tSup[nbcl, indFct] + res$tSd[nbcl, indFct],
                  res$fct[indFct],
                  res$mMotifs[nbcl, indFct], res$names[indFct])
  }

}



#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                   END of FILE myPREDICTING.R
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


