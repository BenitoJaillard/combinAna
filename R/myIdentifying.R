
################################################################################
################################################################################
#
#                         ASSEMBLY SAMPLING EFFECT
#
#     Computation of the best species affectation to functional classes
#               by knowing the best partition of regional pool
#               ==============================================
#
#                     by Beno?t JAILLARD and Camille RICHON
#                               Avril 2014
#
#
#
################################################################################
################################################################################


#-------------------------------------------------------------------------------
# Clear Memory, load files and global variables
#-------------------------------------------------------------------------------

rm(list=ls())

#setwd("H:/2.Modelisation.probabiliste/0.R.code/clean")
#setwd("E:/2.Modelisation.probabiliste/0.R.code/clean2")
setwd("C:/Users/jaillard/Documents/2.Modelisation.probabiliste/0.R.code/clean2")
#setwd("D:/2.Modelisation.probabiliste/0.R.code/clean2")

source("myIO.r")
source("myPlot.r")
source("myCombine.new.r")

library(combinat)

enregistre <- FALSE


#-------------------------------------------------------------------------------
# Look for the n!/(n1!n2!n3!...) affectations of n elements of regional pool
#       inside a given partition [n1,n2,n3,...] of this regional pool
#
#  Inputs : the partition [n1,n2,n3,...] of the regional pool
#           the identity of n elements (nomElt)
#  Output : the mAffect matrix of size n!/(n1!n2!n3!...) x n elements
#-------------------------------------------------------------------------------
affectation <- function(nbAffect, elt, v, fclset, nomfic)
{
  # fclset contains only the gf that are not yet to the quorum
  for (fcl in fclset)
    if (sum(v==fcl) == partition[fcl]) fclset <- setdiff(fclset, fcl)

  for (fcl in fclset)
  {
    v[elt]  <- fcl
    if (elt < nbElt)
    {
      nbAffect <- affectation(nbAffect, (elt+1), v, fclset, nomfic)
    }
    else
    {
      write(x=v, file=nomfic, ncolumns=nbElt, append=TRUE, sep=",")
      nbAffect <- nbAffect + 1
    }
  }

  return(nbAffect)
}


#-------------------------------------------------------------------------------
# Drawing of the linear regression "Predicted.values vs Measured.values"
#
#  Inputs : reg = regression model,
#           alpha = value of 1er risk of first species
#    Output : plots
#-------------------------------------------------------------------------------
confidence.intervals <- function(xlim, reg, alpha)
{
  x   <- (reg$model)[,2]

  n   <- length(x)
  xm  <- mean(x)
  ssx <- sum(x^2) - sum(x)^2/n
  sst <- qt(1-alpha/2, (n-2))

  xx  <- seq(xlim[1], xlim[2], length=20)

  inter.confidence <- sst * sqrt(summary(reg)$sigma^2 * (1/n+(xx-xm)^2/ssx))
  inter.prediction <- sst * sqrt(summary(reg)$sigma^2 * (1+1/n+(xx-xm)^2/ssx))

  res        <- list(xx, inter.confidence, inter.prediction)
  names(res) <- c("x", "confidence", "prediction")

  return(res)
}



#------------------------------------------------------------------------------
plot.regression.11 <- function(x, y, xlim, ylim, reg, inter,
                                           fig, col, titre="Regression")
{
  plot(x, y, main=titre, las=1, lwd=2,
       xlim=xlim, ylim=ylim,
       xlab="measured values", ylab="predicted values",
       type="p", pch=fig, col=col, asp=1)
  abline(a=0, b=1, col="blue", lty="solid")
  abline(h=0, v=0, col="blue", lty="dotted")
  abline(a=reg$coef[1], b=reg$coef[2], col="red", lty="solid")

  fitted.y <- reg$coef[1] + inter$x * reg$coef[2]
  matlines(x=inter$x,
           y=cbind(fitted.y+inter$confidence, fitted.y-inter$confidence,
                   fitted.y+inter$prediction, fitted.y-inter$prediction),
           type="l", col="red", lty=c(rep("dotted",2), rep("longdash",2)))

  text(0.95*xlim[2], 0.05*ylim[2],
       paste("R2=",signif(summary(reg)$"r.squared",3),sep=""))
}



#------------------------------------------------------------------------------
plot.regression.sd <- function(x, y, classe, nbcl, xlim, ylim, fig, col, titre)
{
  mx <- dx <- my <- dy <- rep(0,nbcl)
  for (fcl in 1:nbcl)
  {
    mx[fcl] <- mean(x[classe==fcl])
    dx[fcl] <- sd(x[classe==fcl])
    my[fcl] <- mean(y[classe==fcl])
    dy[fcl] <- sd(y[classe==fcl])
  }

  plot(x=mx, y=my, main=titre, las=1, lwd=2,
       xlim=xlim, ylim=ylim,
       xlab="measured values", ylab="predicted values",
       type="p", pch=fig, col=col, asp=1)
  arrows(x0=mx[dx!=0]-dx[dx!=0], y0=my[dx!=0],
         x1=mx[dx!=0]+dx[dx!=0], y1=my[dx!=0], col=col[dx!=0],
         length=0.1, angle=90, code=3, lwd=1, lty="solid")
  arrows(x0=mx[dy!=0], y0=my[dy!=0]-dy[dy!=0],
         x1=mx[dy!=0], y1=my[dy!=0]+dy[dy!=0], col=col[dy!=0],
         length=0.1, angle=90, code=3, lwd=1, lty="solid")
  abline(h=0, v=0, col="blue", lty="dotted")
  abline(a=0, b=1, col="blue", lty="solid")
}



#------------------------------------------------------------------------------
plot.normality <- function(vecVar, nomVar, titre="normality")
{
  hist(vecVar, main=titre,
       xlab=nomVar, prob=TRUE, las=1, ylab="Density")
  lines(density(vecVar), col="blue")
  rug(vecVar, col="blue")

  f <- function(t) { dnorm(t, mean=mean(vecVar), sd=sd(vecVar), log=FALSE) }
  curve(f, add=TRUE, col="red", lty="dotted")

  moy <- mean(vecVar)
  std <- sqrt(var(vecVar))
  abline(v=c(moy,moy-std,moy+std),lty=c("dashed",rep("longdash",2)),col="red")

  r <- ks.test(vecVar, "pnorm", mean=mean(vecVar), sd=sd(vecVar))
  text(0.9*min(vecVar),0.9*max(density(vecVar)$y),
                                  paste("p(Kolmogorov)=",signif(r$p.value,3)))

  r <- shapiro.test(vecVar)
  text(0.9*min(vecVar), 0.8*max(density(vecVar)$y),
                                  paste("p(Shapiro)=",signif(r$p.value,3)))
}



#===============================================================================
#===============================================================================
#   Main
#===============================================================================
#===============================================================================

#-------------------------------------------------------------------------------
# Reading of assembly data
#-------------------------------------------------------------------------------

# Reading of data-file of assembly performance
res    <- work.directory("Reading of data-file")
oldrad <- paste(res$radical,res$radical,sep="/")

nomfic <- paste(oldrad, "mFmes0.csv", sep=".")
mFmes0 <- read.table(nomfic, header=TRUE, sep=",", dec=".", na.string="NA")
mFmes0 <- as.matrix(mFmes0)

#-------------------------------------------------------------------------------
#  Options
#-------------------------------------------------------------------------------
nbclasses  <- 3
woody      <- FALSE
sorted     <- FALSE

ref       <- "int"  # "int": Ffcl from the data, "ext": Ffcl from the fit-model
enregistre <- TRUE
nbpoints   <- 250    # on the graphs


chrono     <- TRUE
year       <- 2005
rule       <- "diversity"
#rule       <- "presence"


#===============================================================================
#  LOOP on YEARS
#===============================================================================
yearset    <- c(2001:2007,2009:2011)
for (year in yearset)
{
#-------------------------------------------------------------------------------
#  Extraction of appropriated informations
#-------------------------------------------------------------------------------
if (chrono) { nomfic <- paste(oldrad, year, "mOccur.Fmes.csv", sep=".")
} else {      nomfic <- paste(oldrad,       "mOccur.Fmes.csv", sep=".")   }
data   <- read.table(nomfic, header=TRUE, sep=",", dec=".", na.string="NA")
nbAss  <- dim(data)[1]

# Extraction of Fmes-vector (nbAss) of measured productivity for each assembly
Fmes   <- data[,dim(data)[2]]

# Extraction de la matrice mOccur (nbAss x nbElt) de pr?sence dans les assemblages
mOccur <- data[,-dim(data)[2]]
mOccur <- as.matrix(mOccur, nrow=dim(mOccur)[1], ncom=dim(mOccur)[2])
nomElt <- colnames(mOccur)
nbElt  <- length(nomElt)

# Extraction du vecteur Fmes0 (nbElt) de fonctionnement en culture pure
if (chrono) { Fmes0   <- mFmes0[(year-2000),]
} else {      Fmes0   <- mFmes0  }

#-------------------------------------------------------------------------------
# Reading of best-partition
#-------------------------------------------------------------------------------
if (chrono) { nomfic <- paste(oldrad,year,rule,nbclasses,
                                      "no.moore.best.Comp.csv", sep=".")
} else {      nomfic <- paste(oldrad,     rule,nbclasses,
                                      "no.moore.best.Comp.csv", sep=".") }
fit    <- read.table(nomfic, header=TRUE, sep=",", dec=".", na.string="NA")
fit    <- fit[1,]    # One keeps only the first line

if (chrono) { rad <- paste(oldrad, year, rule, nbclasses, "id", ref, sep=".")
} else {      rad <- paste(oldrad,       rule, nbclasses, "id", ref, sep=".") }
if (woody==TRUE)  rad <- paste(rad, "wood", sep=".")
if (sorted==TRUE) rad <- paste(rad, "sorted", sep=".")

partition <- c(fit$s1, fit$s2, fit$s3, fit$s4)
#index     <- sort(partition, decreasing=TRUE, index.return=TRUE)
#partition <- index$x

nbElt     <- sum(partition)
nbfcl     <- length(partition)

Ffclext <- Ffclint <- Ffcl <- c(fit$F1, fit$F2, fit$F3, fit$F4)

if (nbElt != dim(mOccur)[2])
    stop("The number of elements are not equal to in partition and matOccur")

if (rule=="diversity")   # symmetrical rule
{
  # test if there can be doublons in mAffect
  symmetry <- FALSE
  sym      <- table(partition)
  for (i in 1:length(sym)) if (sym[i]!= 1)
  {
    # WARNING: we assume there is only ONE symmetry
    symmetry <- TRUE
    sizeCl   <- as.numeric(names(sym[i]))
    eltCl    <- which(partition==sizeCl)
  }
}

#-------------------------------------------------------------------------------
#  Compute the reference productivities (without interaction)
#                        of functional classes
#-------------------------------------------------------------------------------
AssSize <- apply(mOccur, MARGIN=1, FUN=sum)
Fref0   <- (mOccur %*% Fmes0) / AssSize                 # nbAss


#-------------------------------------------------------------------------------
#  Compute the n!/(n1!n2!n3!...) affectations to test
#-------------------------------------------------------------------------------
nomfic   <- paste(rad, "mAffect.csv", sep=".")
write(x=nomElt, file=nomfic, ncolumns=nbElt, append=FALSE, sep=",")
nbAffect <- affectation(nbAffect=0, elt=1, v=rep(0,nbElt), fclset=1:nbfcl, nomfic)

mAffect  <- read.table(file=nomfic, header=TRUE, row.names=NULL, sep=",")
mAffect  <- as.matrix(mAffect, nrow=nbAffect, ncol=nbElt)


#-------------------------------------------------------------------------------
#  Sort the doublons by symmetry of rule
#-------------------------------------------------------------------------------
if (rule=="diversity")   # symmetrical rule
  if (symmetry)
  {
    # compute the number of doublons and location to index of doublons
    nbCl      <- length(eltCl)
    nbDoublon <- factorial(nbCl)
    mperm     <- t(matrix(unlist(permn(nbCl)),nrow=nbCl))
    index     <- matrix(0, nrow=nbCl, ncol=sizeCl)

    # look for doublons and delete them
    for (aff in 1:(dim(mAffect)[1]/nbDoublon))
    {
      for (j in 1:nbCl) index[j,] <- which(mAffect[aff,]==eltCl[j])

      res <- NULL
      for (i in 1:nbDoublon)
      {
        test <- TRUE
        for (j in 1:nbCl)
           for (k in 1:sizeCl)
               test <- test & (mAffect[,index[mperm[i,j],][k]]==eltCl[j])
        res <- union(res, which(test==TRUE))
      }

     mAffect <- mAffect[-res[-1],]  # we keep the first doublon
    }
    nbAffect <- dim(mAffect)[1]
  }

#### WARNING: HERE nbEltclasses= nbFUNCTClasses

#===============================================================================
#  LOOP on the PERMUTATIONS
#===============================================================================
prediction <- function(affect, mOccur, rule, nbfcl, Ffcl, Fmes, Fref0, int)
{
   #  Tag the functional classe of elements in matOccur
   EltClasses <- t( t(mOccur) * affect)
   #  Tag the classe of functioning of assemblies according to prevalent rule
   AssClasses <- apply(as.matrix(EltClasses), MARGIN=1,
                    FUN=affect.functional.classe, rule=rule)
   #  Compute the productivity mean-value of assembly classes
   if (ref == "int")
       for (fcl in 1:nbfcl) Ffcl[fcl] <- mean(Fmes[AssClasses==fcl])
   #  Compute the expected productivity mean-value of assembly classes
   #                              (without interactions)
   Fref  <- Ffcl
   for (fcl in 1:nbfcl) Fref[fcl] <- mean(Fref0[AssClasses==fcl])
   #  Compute the predicted productivity of assemblies
   Fpred <- Fmes
   for (fcl in 1:nbfcl)
     Fpred[(AssClasses==fcl)] <- Fref0[(AssClasses==fcl)]*Ffcl[fcl]/Fref[fcl]

   return(Fpred)
}

#-------------------------------------------------------------------------------
titre <- c("R2","Fischer","pf","intercept","slope","logLike", "RMSE")
mRes  <- matrix(0, nrow=nbAffect, ncol=length(titre))
colnames(mRes) <- titre

for (aff in 1:nbAffect)
{
  print(aff)
  affect <- mAffect[aff,]
  Fpred  <- prediction(affect, mOccur, rule, nbfcl, Ffcl, Fmes, Fref0, int)

  res1   <- lm(Fpred ~ Fmes)
  res2   <- summary(res1)
  res3   <- anova(res1)

  mRes[aff,1] <- res2$"r.squared"[1]
  mRes[aff,2] <- res3$"F value"[1]
  mRes[aff,3] <- res3$"Pr(>F)"[1]
  mRes[aff,4] <- res1$coefficients[1]
  mRes[aff,5] <- res1$coefficients[2]
  mRes[aff,6] <- logLik(res1)
  mRes[aff,7] <- rmse(Fpred, Fmes)
}
remove(aff)


#-------------------------------------------------------------------------------
#  Sorting of results and Identifying of the best affectation
#-------------------------------------------------------------------------------
nomfic <- paste(rad, "mRes.csv", sep=".")

out    <- cbind(mRes, mAffect)
index  <- sort(out[,1], decreasing=TRUE, index.return=TRUE)
out    <- out[index$ix,]

write.table(x=out, file=nomfic, append=FALSE,
                                col.names=TRUE, row.names=FALSE, sep=",")


#===============================================================================
#
#  END OF COMPUTATIONS - BEGINNING OF OUT OF RESULTS
#
#===============================================================================
if (chrono) { rad <- paste(oldrad, year, rule, nbclasses, "id", ref, sep=".")
} else {      rad <- paste(oldrad,       rule, nbclasses, "id", ref, sep=".") }
if (woody==TRUE)  rad <- paste(rad, "wood", sep=".")
if (sorted==TRUE) rad <- paste(rad, "sorted", sep=".")
nomfic <- paste(rad, "mRes.csv", sep=".")
titre <- c("R2","Fischer","pf","intercept","slope","logLike", "RMSE")

out      <- read.table(nomfic, header=TRUE, sep=",", dec=".")
mRes     <- as.matrix(out[, c(1:length(titre))])
mAffect  <- as.matrix(out[,-c(1:length(titre))])
nbAffect <- dim(mAffect)[1]


#-------------------------------------------------------------------------------
#  Look for stabilizing the element affectation when there is a symmetry of rule
#-------------------------------------------------------------------------------
if (rule=="diversity")
     if (symmetry)   # symmetrical rule
{
  # test if there can be doublons in mAffect
  sym      <- table(partition)
  for (i in 1:length(sym)) if (sym[i]!= 1)
  {
    sizeCl   <- as.numeric(names(sym[i]))
    eltCl    <- which(partition==sizeCl)
  }

  nbCl      <- length(eltCl)
  mperm     <- t(matrix(unlist(permn(nbCl)),nrow=nbCl))
  indexprec <- index <- matrix(0, nrow=nbCl, ncol=sizeCl)

  aff <- 1
  for (j in 1:nbCl) indexprec[j,] <- which(mAffect[aff,]==eltCl[j])

  aff <- 2
  for (aff in 2:nbAffect) #dim(mAffect)[1])
    {
      for (j in 1:nbCl) index[j,] <- which(mAffect[aff,]==eltCl[j])
      meme <- intersect(index,indexprec)
      if (length(meme) > 0)
        if (mAffect[aff,meme] != mAffect[aff-1,meme])
          for (j in 1:nbCl) for (k in 1:sizeCl)
              mAffect[aff,index[mperm[2,j],][k]] <- eltCl[j]
      for (j in 1:nbCl) indexprec[j,] <- which(mAffect[aff,]==eltCl[j])
    }
}


#===============================================================================
# Compute THE BEST AFFECTATION
#===============================================================================
# Affectation en groupes fonctionnels
# affect <- c(1,3,2,4,1,3,3,2,1,2,1,4,2,3,4,4)
# avec 1: forb, 2: fabac?es, 3: C3 et 4:C4.

bestAffect <- mAffect[1,]
bestEltClasses <- t( t(mOccur) * bestAffect)
bestAssClasses <- apply(as.matrix(bestEltClasses), MARGIN=1,
                    FUN=affect.functional.classe, rule=rule)
Fref <- Ffcl
for (fcl in 1:nbfcl) Fref[fcl] <- mean(Fref0[bestAssClasses==fcl])
for (fcl in 1:nbfcl) Ffclint[fcl] <- mean(Fmes[bestAssClasses==fcl])
if (ref == "int") Ffcl <- Ffclint

bestFpred <- Fmes
for (fcl in 1:nbfcl) bestFpred[(bestAssClasses==fcl)] <-
                             Fref0[(bestAssClasses==fcl)]*Ffcl[fcl]/Fref[fcl]


# Save the results of the best model
nomfic <- paste(rad, "bestFpred.csv", sep=".")
write.table(x=cbind(Fmes, bestFpred, Fref0, mOccur), nomfic,
                    append=FALSE, col.names=TRUE, row.names=FALSE, sep=",")


#-------------------------------------------------------------------------------
if (nbfcl > length(couleurs))
    for (i in 1:floor(nbfcl/length(couleurs))) couleurs <- c(couleurs, couleurs)
if (nbfcl > length(figures))
    for (i in 1:floor(nbfcl/length(figures)))  figures  <- c(figures,  figures)
regcouleurs <- couleurs[bestAssClasses]
regfigures  <- figures[bestAssClasses]

xreglim <- c(0,max(Fmes))
yreglim <- c(0,max(bestFpred, Fref0))


#-------------------------------------------------------------------------------
# Regression of the best Affectation
#-------------------------------------------------------------------------------
bestReg   <- lm(bestFpred ~ Fmes)
bestInter <- confidence.intervals(xlim=xreglim, reg=bestReg, alpha=0.05)

#  Plot the regression of the best Affectation
open.window(enregistre,1, paste(rad,"bestRegression",sep="."))

plot.regression.11(x=Fmes, y=bestFpred,
                   xlim=xreglim, ylim=yreglim,
                   reg=bestReg, inter=bestInter,
                   fig=regfigures, col=regcouleurs,
                   paste(rad, "bestRegression", sep="."))

#  Plot the regression with element names
plot.regression.11(x=Fmes, y=bestFpred,
                   xlim=xreglim, ylim=yreglim,
                   reg=bestReg, inter=bestInter,
                   fig=".", col="black",
                   paste(rad, "bestElements", sep="."))
for (elt in 1:nbElt) text(mean(Fmes[(mOccur[,elt]==1)]),
                          mean(bestFpred[(mOccur[,elt]==1)]), nomElt[elt])

#  Analyse the residuals of the regression
layout(matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE))
plot(bestReg)
layout(1)

#  Test the normality of residuals
plot.normality(bestReg$residuals, "bestFpred", paste(rad,"bestNormality",sep="."))

#  Plot the regression for each element
for (elt in 1:nbElt)
{
  plot.regression.11(x=Fmes[(mOccur[,elt]!=1)], y=bestFpred[(mOccur[,elt]!=1)],
                     xlim=xreglim, ylim=yreglim,
                     reg=bestReg, inter=bestInter,
                     fig=".", col="black",
                     paste(rad, "best", nomElt[elt], sep="."))
  points(x=Fmes[(mOccur[,elt]==1)], y=bestFpred[(mOccur[,elt]==1)],
         xlim=xreglim, ylim=yreglim,
         pch=regfigures[(mOccur[,elt]==1)],
         col=regcouleurs[(mOccur[,elt]==1)], lwd=2)
}

#  Plot the mean-Values of classes on the best regression
plot.regression.sd(x=Fmes, y=bestFpred,
                   classe=bestAssClasses, nbcl=nbfcl,
                   xlim=xreglim, ylim=yreglim,
                   fig=figures, col=couleurs,
                   paste(rad,"bestMeanRegression",sep="."))

#  Plot the mean-values of classes for each element
for (elt in 1:nbElt)
{
  plot.regression.sd(x=Fmes[(mOccur[,elt]==1)], y=bestFpred[(mOccur[,elt]==1)],
                     classe=bestAssClasses[(mOccur[,elt]==1)], nbcl=nbfcl,
                     xlim=xreglim, ylim=yreglim,
                     fig=figures, col=couleurs,
                     paste(rad, "best", nomElt[elt], sep="."))
  points(x=Fmes[(mOccur[,elt]!=1)], y=bestFpred[(mOccur[,elt]!=1)],
         xlim=xreglim, ylim=yreglim,
         pch=".", col="black")
 }

close.window(enregistre)


#-------------------------------------------------------------------------------
# Regression of reference (without interaction)
#-------------------------------------------------------------------------------
refReg   <- lm(Fref0 ~ Fmes)
refInter <- confidence.intervals(xlim=xreglim, reg=refReg, alpha=0.05)

#  Plot reference regression (without interaction)
open.window(enregistre,1, paste(rad,"refRegression",sep="."))

plot.regression.11(x=Fmes, y=Fref0,
                   xlim=xreglim, ylim=yreglim,
                   reg=refReg, inter=refInter,
                   fig=regfigures, col=regcouleurs,
                   paste(rad, "ref.Regression", sep="."))

#  Plot the reference regression with element names
plot.regression.11(x=Fmes, y=Fref0,
                   xlim=xreglim, ylim=yreglim,
                   reg=refReg, inter=refInter,
                   fig=".", col="black",
                   paste(rad, "ref", nomElt[elt], sep="."))
for (elt in 1:nbElt) text(mean(Fmes[(mOccur[,elt]==1)]),
                          mean(Fref0[(mOccur[,elt]==1)]), nomElt[elt])

#  Analyse the residuals of the regression
layout(matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE))
plot(refReg)
layout(1)

#  Test the normality of residuals
plot.normality(refReg$residuals, "Fref0", paste(rad,"refNormality",sep="."))

#  Plot the reference regression for each element
for (elt in 1:nbElt)
{
  plot.regression.11(x=Fmes[(mOccur[,elt]!=1)], y=Fref0[(mOccur[,elt]!=1)],
                     xlim=xreglim, ylim=yreglim,
                     reg=refReg, inter=refInter,
                     fig=".", col="black",
                     paste(rad, "ref", nomElt[elt], sep="."))
  points(x=Fmes[(mOccur[,elt]==1)], y=Fref0[(mOccur[,elt]==1)],
         xlim=xreglim, ylim=yreglim,
         pch=regfigures[(mOccur[,elt]==1)],
         col=regcouleurs[(mOccur[,elt]==1)], lwd=2)
}

#  Plot the mean-Values of classes on the ref regression
plot.regression.sd(x=Fmes, y=Fref0,
                   classe=bestAssClasses, nbcl=nbfcl,
                   xlim=xreglim, ylim=yreglim,
                   fig=figures, col=couleurs,
                   paste(rad,"ref.Regression",sep="."))

#  Plot the mean-values of classes for each element
for (elt in 1:nbElt)
{
  plot.regression.sd(x=Fmes[(mOccur[,elt]==1)], y=Fref0[(mOccur[,elt]==1)],
                     classe=bestAssClasses[(mOccur[,elt]==1)], nbcl=nbfcl,
                     xlim=xreglim, ylim=yreglim,
                     fig=figures, col=couleurs,
                     paste(rad, "ref", nomElt[elt], sep="."))
  points(x=Fmes[(mOccur[,elt]!=1)], y=Fref0[(mOccur[,elt]!=1)],
         xlim=xreglim, ylim=yreglim,
         pch=".", col="black")
}

close.window(enregistre)


#-------------------------------------------------------------------------------
# Outputs of statistical properties of all affectations
#-------------------------------------------------------------------------------
nomfic <- paste(rad, "StatCurves", sep=".")
open.window(enregistre, 1, nomfic)

mask  <- seq(from=1, to=nbAffect, by=1+floor(nbAffect/nbpoints))

# Compute the p.value that (bestR2-followingR2) are signifiant
pvalue  <- rep(1, length(mask))
for (aff in mask[-1])
{
  affect  <- mAffect[aff,]
  Ffollow <- prediction(affect, mOccur, rule, nbfcl, Ffcl, Fmes, Fref0, int)
  res     <- cor(bestFpred, Ffollow)^2
  pvalue[aff] <- compare.dependent.R2(mRes[1,"R2"],mRes[aff,"R2"],res,nbAss)
}

titre <- paste(rad, "diff.R2.log10(pvalue)", "decreaseAllValues", sep=".")
alpha <- 0.05
plot(x=mask, y=log10(pvalue[mask]), xlab="Rank", ylab="diff.R2.log10(pvalue)")
abline(h=log10(alpha), col="red")
text(0.9*mask[length(mask)],log10(alpha)+1,paste("alpha",alpha,sep="="), col="red")
title(titre)

# Compute the p.value that (refR2-followingR2) are signifiant
refCor   <- summary(refReg)$"r.squared"[1]
refvalue <- rep(1, length(mask))
for (aff in mask)
{
  affect  <- mAffect[aff,]
  Ffollow <- prediction(affect, mOccur, rule, nbfcl, Ffcl, Fmes, Fref0, int)
  res     <- cor(Fref0, Ffollow)^2
  pvalue[aff] <- compare.dependent.R2(refCor,mRes[aff,"R2"],res,nbAss)
}

titre <- paste(rad, "diff.refR2.log10(pvalue)", "decreaseAllValues", sep=".")
alpha <- 0.05
plot(x=mask, y=log10(pvalue[mask]), xlab="Rank", ylab="diff.R2.log10(pvalue)")
abline(h=log10(alpha), col="red")
text(0.9*mask[length(mask)],log10(alpha)+1,paste("alpha",alpha,sep="="), col="red")
title(titre)

for (i in 1:dim(mRes)[2])
{
  variable <- colnames(mRes)[i]
  titre    <- paste(rad, variable, "decreaseAllValues", sep=".")
  plot(x=mask, y=mRes[mask,variable], xlab="Rank", ylab=variable)
  title(titre)
}

variable <- "pf"
titre    <- paste(rad, "log10Pf", "decreaseAllValues", sep=".")
plot(x=mask, y=log10(mRes[mask,variable]), xlab="Rank", ylab="log10(pf)")
title(titre)


#-------------------------------------------------------------------------------
# Outputs of statistical properties of (first) best affectations
#-------------------------------------------------------------------------------
mask  <- c(1:min(nbAffect,nbpoints))

# Compute the p.value that (bestR2-followingR2) are signifiant
pvalue  <- rep(1, length(mask))
for (aff in mask[-1])
{
  affect  <- mAffect[aff,]
  Ffollow <- prediction(affect, mOccur, rule, nbfcl, Ffcl, Fmes, Fref0, int)
  res     <- cor(bestFpred, Ffollow)^2
  pvalue[aff] <- compare.dependent.R2(mRes[1,"R2"],mRes[aff,"R2"],res,nbAss)
}

titre <- paste(rad, "diff.R2.log10(pvalue)", "decreaseFirstValues", sep=".")
alpha <- 0.05
plot(x=mask, y=log10(pvalue[mask]), xlab="Rank", ylab="diff.R2.log10(pvalue)")
abline(h=log10(alpha), col="red")
text(0.9*mask[length(mask)],log10(alpha)+1,paste("alpha",alpha,sep="="), col="red")
title(titre)

for (i in 1:dim(mRes)[2])
{
  variable <- colnames(mRes)[i]
  titre    <- paste(rad, variable, "decreaseFirstValues", sep=".")
  plot(x=mask, y=mRes[mask,variable], xlab="Rank", ylab=variable)
  title(titre)
}

variable <- "pf"
titre    <- paste(rad, "log10Pf", "decreaseFirstValues", sep=".")
plot(x=mask, y=log10(mRes[mask,variable]), xlab="Rank", ylab="log10(pf)")
title(titre)

close.window(enregistre)


#-------------------------------------------------------------------------------
# Write the best Affectation in ASCII format
#-------------------------------------------------------------------------------
# options, stats de ref, stats de best, p.value, Fext, Fint, affect

titre <- c("rad", "rule", "nbClasses", "nbAffect", "nbAss",
           "F0", paste("F",1:max(3,nbclasses),"int",sep=""),
           colnames(mAffect),
           paste("best", colnames(mRes),  sep=""),
           paste("ref",  colnames(mRes),  sep=""))

out     <- matrix(0, nrow=1, ncol=length(titre))
colnames(out) <- titre

out[1,] <-  unlist(c(rad, rule, nbclasses, nbAffect, nbAss,
            mean(Fref0), Ffclint,
            mAffect[1,],
            mRes[1,],
            (summary(refReg))$"r.squared"[1],(anova(refReg))$"F value"[1],
            (anova(refReg))$"Pr(>F)"[1], refReg$coefficients[1],
            refReg$coefficients[2], logLik(lm(Fref0~Fmes)), rmse(Fref0,Fmes)))

write.table(x=out, file=paste(rad, "bestAffectation.csv", sep="."),
            append=FALSE, col.names=TRUE, row.names=FALSE, sep=",")



#===============================================================================
#  Plot the relationship "Productivity vs Diversity"
#===============================================================================
#t.test(Fmes[(bestAssClasses==3)], Fmes[(bestAssClasses==2)])

offset  <- 0.15
sizeset <- as.numeric(names(table(AssSize)))
nbsize  <- length(sizeset)

Prop    <- array(0, dim=c(5, nbsize, nbfcl))
for (size in 1:nbsize) for (fcl in 1: nbfcl)
    {
    index  <- which((AssSize==sizeset[size]) & (bestAssClasses==fcl))
                              Prop[1,size,fcl] <- length(Fmes[index])
    if (Prop[1,size,fcl] > 0) Prop[2,size,fcl] <- mean(Fmes[index])
    if (Prop[1,size,fcl] > 1) Prop[3,size,fcl] <- sd(Fmes[index])
    if (Prop[1,size,fcl] > 0) Prop[4,size,fcl] <-
                                         Prop[3,size,fcl]/sqrt(Prop[1,size,fcl])
    }
Prop[5,,] <- Prop[1,,]/apply(as.matrix(Prop[1,,]), MARGIN=1, FUN=sum)

#for (fcl in 2:nbfcl)
#  (t.test(Fmes[bestAssClasses==(fcl-1)], Fmes[bestAssClasses==fcl]))
#t.test(Fmes[bestAssClasses==1], Fmes[bestAssClasses==2])
#t.test(Fmes[bestAssClasses==2], Fmes[bestAssClasses==3])


# Plot theoretical probabilities and experimental proportions
open.window(enregistre, 1, paste(rad, "bestCurves", sep="."))

titre <- paste(rad, "Proportions", sep=".")
pdom  <- exact.generic.proba(rule, partition)
plot.exact.solution(cbind(1:nbElt,pdom),"")

mat   <- cbind(sizeset,Prop[5,,])
matpoints(x=mat[,1], y=mat[,-1], ylim=c(0,1),
          xlab="Assembly size", ylab="Assembly proportions",
          type="p", pch=figures, las=1, lwd=2, col=couleurs)
mtext(titre, side=3, cex=0.8)


# Plot mean-values of classes of functioning
titre <- paste(rad, "Fclasses", sep=".")

x <- sizeset
y <- Prop[2,,]
z <- Prop[4,,]
xlim <- c(1, nbElt)
ylim <- c(0, max(Fmes))

fcl <- 1
Offset <- -offset
plot(x=x[(y[,fcl]!=0)]+Offset, y=y[(y[,fcl]!=0),fcl],
     xlab="Assembly size", ylab="Assembly productivity",
     xlim=xlim, ylim=ylim,
     type="p", las=1, lwd=2, lty="solid", pch=figures[fcl], col=couleurs[fcl])
arrows(x0=x[(z[,fcl]!=0)]+Offset, y0=y[(z[,fcl]!=0),fcl]-z[(z[,fcl]!=0),fcl],
       x1=x[(z[,fcl]!=0)]+Offset, y1=y[(z[,fcl]!=0),fcl]+z[(z[,fcl]!=0),fcl],
       length=0.1, angle=90, code=3,
       lwd=1, lty="solid", col=couleurs[fcl])

for (fcl in 2:nbfcl)
{
    Offset <- (fcl-2) * offset
    lines(x=x[(y[,fcl]!=0)]+Offset, y=y[(y[,fcl]!=0),fcl],
          xlim=xlim, ylim=ylim,
          type="p", las=1, lwd=2, lty="solid",
          pch=figures[fcl], col=couleurs[fcl])
    arrows(x0=x[(z[,fcl]!=0)]+Offset, y0=y[(z[,fcl]!=0),fcl]-z[(z[,fcl]!=0),fcl],
           x1=x[(z[,fcl]!=0)]+Offset, y1=y[(z[,fcl]!=0),fcl]+z[(z[,fcl]!=0),fcl],
           length=0.1, angle=90, code=3,
           lwd=1, lty="solid", col=couleurs[fcl])
}

for (fcl in 1:nbfcl) abline(h=Ffclint[fcl], col=couleurs[fcl], lty="dotted")

remove(x,y,z)


# Plot mean-values of F and the adjustements
titre <- paste(rad, "Fadjust", sep=".")
xlim <- c(1, nbElt)
ylim <- c(0, max(Fmes))

plot(x=1:nbElt, y=(pdom %*% Ffclext),
     xlim=xlim, ylim=ylim,
     xlab="Assembly size", ylab="Assembly productivity",
     type="l", lty="solid")

for (fcl in 1:nbfcl) abline(h=Ffclext[fcl], col=couleurs[fcl], lty="solid")

points(x=1:nbElt, y=(pdom %*% Ffclint), type="l", lty="dotted")
for (fcl in 1:nbfcl) abline(h=Ffclint[fcl], col=couleurs[fcl], lty="dotted")
title(titre)


# Plot raw-values of classes of functioning
titre <- paste(rad, "Ftags", sep=".")

AssPlot <- AssSize
for (fcl in 1:nbfcl)
{
    Offset <- (fcl-2) * offset
    AssPlot[(bestAssClasses==fcl)] <- AssPlot[(bestAssClasses==fcl)]+Offset
}
xlim <- c(1, nbElt)
ylim <- c(0, max(Fmes))

plot(x=AssPlot, y=Fmes,
     xlim=xlim, ylim=ylim,
     xlab="Assembly size", ylab="Assembly productivity",
     type="p", las=1, lwd=2, pch=regfigures, col=regcouleurs)



#  Plot raw-values of classes for each element
xlim <- c(1, nbElt)
ylim <- c(0, max(Fmes))

for (elt in 1:nbElt)
{
  plot(x=AssPlot[(mOccur[,elt]!=1)], y=Fmes[(mOccur[,elt]!=1)],
       xlim=xlim, ylim=ylim,
       xlab="Assembly size", ylab="Assembly productivity",
       type="p", las=1, lwd=2, pch=".", col="black")
  for (fcl in 1:nbfcl) abline(h=Ffcl[fcl], col=couleurs[fcl])
  points(x=AssPlot[(mOccur[,elt]==1)], y=Fmes[(mOccur[,elt]==1)],
       xlim=xlim, ylim=ylim,
       xlab="Assembly size", ylab="Assembly productivity",
       type="p", las=1, lwd=2, pch=regfigures[(mOccur[,elt]==1)],
       col=regcouleurs[(mOccur[,elt]==1)])
  title(paste(rad, nomElt[elt], sep="."))
  points(x=1:nbElt, y=(pdom %*% Ffcl),    type="l", lty="solid")

}

remove(AssPlot)
close.window(enregistre)


out           <- t(Prop[1,,])
colnames(out) <- sizeset
for (i in 2:5) out <- rbind(out, sizeset, t(Prop[i,,]))
write.table(x=out, file=paste(rad, "bestAssClasses.csv", sep="."),
              append=FALSE, col.names=TRUE, row.names=FALSE, sep=",")




#===============================================================================
# Decomposition in "Productivity effect" and "Interaction effect"
# Plot histograms of interaction factor
#===============================================================================
titre <- paste(rad, "Histo", sep=".")
open.window(enregistre,1,paste(rad,"Interaction",sep="."))

histogram <- function(vecVar, nomVar)
{
  hist(vecVar, prob=TRUE, las=1, ylab="Density", main=nomVar)
  for (fcl in 1:nbfcl) if (length(vecVar[bestAssClasses==fcl])>1)
  {
    lines(density(vecVar[bestAssClasses==fcl],
                weights=rep(1/length(vecVar),length(vecVar[bestAssClasses==fcl]))),
                col=couleurs[fcl])
    abline(v=mean(vecVar[bestAssClasses==fcl]), col=couleurs[fcl])
    rug(vecVar[bestAssClasses==fcl], col=couleurs[fcl])
  }
}

histogram(Fref0/mean(Fref0), "Productivity effect")
for (fcl in 2:nbfcl) if (length(Fref0[(bestAssClasses==fcl)])>0)
  print(t.test(Fref0[(bestAssClasses==(fcl-1))], Fref0[(bestAssClasses==fcl)]))


fInter <- Fmes/Fref0
histogram(fInter, "Interaction effect")
for (fcl in 2:nbfcl) if (length(Fref0[(bestAssClasses==fcl)])>0)
  print(t.test(fInter[(bestAssClasses==(fcl-1))], fInter[(bestAssClasses==fcl)]))

close.window(enregistre)


#===============================================================================
#  Plot the diagramm dominance-rank with the best affectation of elements
#===============================================================================
open.window(enregistre,1,paste(rad,"affectElements",sep="."))

index <- sort(Fmes0, decreasing=TRUE, index.return=TRUE)

plot(x=c(1:nbElt), y=index$x,
     xlim=c(1,nbElt), ylim=c(0,1.05*max(Fmes0)),
     xlab="element", ylab="productivity in pure culture",
     las=1, lwd=2, col=couleurs[bestAffect[index$ix]])
text(x=c(1:nbElt), y=index$x, nomElt[index$ix], pos=3)
title("best.Affectation")

#-------------------------------------------------------------------------------
#  Plot the element affectation frequencies
#                          in non-significantly different assembly
#-------------------------------------------------------------------------------
mask <- which(pvalue >= alpha)

sfreq <- freq <- matrix(0, nrow=nbfcl, ncol=nbElt)
for (fcl in 1:nbfcl)
  if (length(mask)>1) {
      freq[fcl,] <- apply((mAffect[mask,]==fcl), MARGIN=2, FUN=sum)/length(mask)
  } else {
      freq[fcl,] <- (mAffect[mask,]==fcl)
  }
freq           <- freq[,index$ix]
colnames(freq) <- nomElt[index$ix]
sfreq          <- apply(freq, MARGIN=2, FUN=cumsum)

fcl <- 1
plot(x=which(freq[fcl,]!=0), y=sfreq[fcl,(freq[fcl,]!=0)],
     xlim=c(1,nbElt), ylim=c(0,1),
     xlab="element", ylab="productivity in pure culture",
     las=1, lwd=2, col=couleurs[fcl])
for (fcl in 2:nbfcl)
  points(x=which(freq[fcl,]!=0), y=sfreq[fcl,(freq[fcl,]!=0)],
       xlim=c(1,nbElt), ylim=c(0,1),
       xlab="element", ylab="productivity in pure culture",
       las=1, lwd=2, col=couleurs[fcl])

text(x=c(1:nbElt), y=index$x/max(Fmes0), nomElt[index$ix], pos=1)
title(paste("mean.Affectation for alpha=", alpha,
            " and n=", length(mask), sep=""))

close.window(enregistre)

write.table(x=freq, file=paste(rad, "affectElements.csv", sep="."),
            append=FALSE, col.names=TRUE, row.names=FALSE, sep=",")

#===============================================================================
#  END of the LOOP on YEARS
#===============================================================================
}

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
