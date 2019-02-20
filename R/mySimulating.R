#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#                   SIMULATING ELEMENT ASSEMBLAGES
#
#                          Benoit JAILLARD
#                           Octobre 2015
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
# Clear Memory, load files and global variables
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rm(list = ls())

setwd("D:/1.BEF.Relation/0.R.code/Y1") # lab
setwd("C:/Users/jaillard/Documents/1.BEF.Relation/0.R.code/Y1") # home




#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#
#   Main
#
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

rad <- c("virtual")

partition <- c(2,2,3,3,6)
nbElt     <- sum(partition)
nbclElt   <- length(partition)

rule    <- "diversity"



# naming each element by class

name.elements <- function(partition) {
  v1 <- v2 <- NULL
  for (i in 1:nbclElt)
  {
    #        v1 <- c(v1, rep(i, partition[i]))
    #        v2 <- c(v2, letters[1:partition[i]])

    v1 <- c(v1, rep(letters[i], partition[i]))
    v2 <- c(v2, seq(1:partition[i]))
  }

  return(paste(v1, v2, sep = ""))
}


# building a random element partition

set.partition <- function(partition) {
  nbclElt   <- length(partition)

  affectElt <- NULL
  for (i in 1:nbclElt) affectElt <- c(affectElt, rep(i, partition[i]))
  names(affectElt) <- name.elements(partition)

  return(affectElt)
}

set.random.partition <- function(nbElt, nbclElt, rule) {
  nbclElt   <- length(partition)
  partition <- integer(nbclElt)
  for (i in 1:(nbclElt - 1))
      partition[i] <- sample(1:(nbElt - (nbclElt - i) - sum(partition)), 1)
  partition[nbclElt] <- nbElt - sum(partition)

  affectElt <- set.partition(sort(partition))
  names(affectElt) <- name.elements(partition)

  return(affectElt)
}

affectElt <- set.random.partition(nbElt, nbclElt, rule)
affectElt <- set.partition(partition)

nbElt     <- sum(partition)
nbclElt   <- length(partition)




# building mOccur

mOccur.M <- NULL
for (j in 1:nbElt) {
  nbComb <- choose(nbElt, j)
  b      <- t(combn(nbElt, j))
  a      <- matrix(0, nrow = nbComb, ncol = nbElt)
  for (i in 1:nbComb) a[i, b[i, 1:j]] <- 1
  mOccur.M <- rbind(mOccur.M, a)
}
colnames(mOccur.M) <- names(affectElt)
nbAss <- dim(mOccur.M)[1]


# calculating Fobs the assembly performance

AssMotifs <- affect.motifs(affectElt, mOccur.M)
setMotif  <- unique(AssMotifs)

perform <- runif(length(setMotif), min = 20, max = 120)
hist(perform)

Fobs.M    <- numeric(nbAss)
for (i in seq_along(setMotif)) Fobs.M[(AssMotifs == i)] <- perform[i]

error <- 10
noise <- runif(nbAss, min = -error, max = +error)
noise <- rnorm(nbAss, mean = 0, sd = error)
Fobs.M  <- Fobs.M + noise


# saving the dataset

#index  <- sample(1:nbSpe, nbSpe, replace = FALSE)
#mOccur <- mOccur[,index]
write.table(x = cbind(mOccur.M, Fobs.M),
            file = paste(rad, paste(experiment, "csv", sep = "."), sep = "/"),
            append = FALSE, col.names = TRUE, row.names = FALSE, sep = ",")


file   <- paste(rad, paste(experiment, "csv", sep = "."), sep = "/")
dat    <- read.table(file, header = TRUE, sep = ",")

nbElt  <- dim(dat)[2] - 1
nbAss  <- 2 ^ nbElt - 1
mOccur.M <- as.matrix(dat[,1:nbElt])
Fobs.M   <- as.vector(as.numeric(unlist(dat[,nbElt + 1])))




experiment <- c("V8")
enregistre <- TRUE
#enregistre <- FALSE

size   <- apply(mOccur.M, MARGIN = 1, FUN = sum )

nbAss  <- 60
index2  <- sample(which(size == 2), 50, replace = FALSE)
index4  <- sample(which(size == 4), 50, replace = FALSE)
index8  <- sample(which(size == 8), 50, replace = FALSE)
index16 <- which(size == 16)

index <- c(index2, index4, index8, index16)
nbAss <- length(index)

index <- c(11299, 36337, 55623, 29550, 18820, 19140, 46946, 55539, 33667,
           22137, 28268, 64769, 31921,  1562, 19595, 18963, 49214, 61281,
           42882, 24246, 39804, 63336, 61015, 48324, 55574, 45535, 43042,
           55306, 12071, 12549, 53237, 11281, 33790,  3249, 40950, 18473,
           48856, 14612, 17426,  3386, 30379, 41783, 58458, 37825,  7183,
           21540,  8408, 26846, 50874, 46124, 13812, 44966, 32695, 11103,
           37366, 31998, 32199, 18595, 52584, 39023, 34285, 44412, 42337,
           10832, 48786, 43081, 20822, 33990, 35949, 17040, 27894, 47444,
            4601, 47783, 10579, 26305, 57135, 25912, 41065, 15067, 18575,
             925, 16767, 62711, 57396, 10771, 57903, 22654, 57851, 49899,
           39456, 24827, 53931,  5143, 40119, 56824, 60287, 26220, 59708,
           31400, 22723, 26048, 32818,  3938, 42562, 49375, 10003, 60987,
           60351, 52930, 13829,  9363, 18923, 14005, 22696, 60450, 12219,
           51045, 27826, 11751, 21924, 27210, 51060, 54070, 14329,  1396,
            8658, 61828, 55178, 30686, 57441, 15307, 22473,  5346, 22167,
           31824, 43257, 36989, 26695, 30798, 28770, 59527, 54901, 13987,
           14823, 51067,  7726, 14487, 21503, 38532, 54215,  7126, 28878,
           23541, 50933, 38795,  6133,  9676, 24558, 52836, 20049,  3920,
            5621, 13424, 49484,  4679, 35299, 40736,  6683, 56239, 44206,
           57657, 28746, 22197, 54815, 45925, 62559, 56476, 20967, 16263,
           37736, 62195, 45548, 46149, 25402, 50565, 39513, 60744, 49250,
           58852, 12840, 15955, 50538, 31130, 14059, 50419,  9907, 18006,
           44731, 52490)

mOccur <- mOccur.M[index, ]
Fobs   <- Fobs.M[index]
figures  <- check.symbol(figures,  dim(mOccur)[1])
couleurs <- check.symbol(couleurs, dim(mOccur)[1])

apply(mOccur, MARGIN = 2, FUN = sum)


tmp    <- rm.dual.assemblies(mOccur, Fobs)
mOccur <- tmp$mat
Fobs   <- tmp$fct


dim(mOccur)

#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
fct <- Fobs
nom <- experiment <- c("V8")

nomfic <- paste(rad, nom, sep = "/")
if (!file.exists(nomfic)) dir.create(nomfic)

titre    <- paste(setYears[year], setMonths[year], "Field.size", sep = ".")
filename <- paste(rad, nom, titre, sep = "/")
open.pdf(enregistre, filename)
open.txt(enregistre, filename)

size   <- apply(mOccur, MARGIN = 1, FUN = sum )
plot(y = Fobs, x = jitter(size), xlim = c(1, 16),
     tck = 0.02, cex = 2, bg = "white")
plot(y = Fobs, x = jitter(AssMotifs[index]), xlim = c(1, 31),
     tck = 0.02, cex = 2, bg = "white")


res.cal <- cluster.elements(mOccur, fct,
                            opt.meth = "divisive",
                            opt.mean = "amean",
                            opt.mod  = "byelt",
                            opt.ssr  = "cal")

plot.tree(res.cal, col = "black", titre)
plot.tree(res.cal, col = couleurs[gf], titre)

write.txt(enregistre, filename,
          full.names[sort.tree(res.cal$aff, index.return = TRUE)$ix]   )


mAffectElt <- res.cal$aff
mAssMotifs <- mAffect.motifs(res.cal, mOccur)

res.Fobs <- plot.predict.fct(mAssMotifs, mOccur, Fobs,
                             opt.mean = "amean",
                             opt.mod  = "byelt",
                             opt.R2  = TRUE,
                             opt.cal = FALSE,
                             opt.prd = FALSE,
                             opt.glb = FALSE,
                             opt.pub = TRUE,
                             opt.aov = FALSE,
                             pvalue = pvalue,
                             titre = "Fobs")


res     <- res.cal
res$cor <- res.Fobs$tR2[ ,"R2cal"]

nbcl    <- first.optimum(res.Fobs$tR2[, "AICc"], opt = "min")
cutline <- mean(res$cor[(nbcl - 1):nbcl])

plot.tree(res, col = "black", titre)
lines(x = c(0, nbElt + 1), y = rep(cutline, 2), col = "blue")
lines(x = c(0, nbElt + 1),
      y = rep(last(res.Fobs$tR2[ , "R2prd"]), 2), col = "red")

plot.tree(res, col = couleurs[gf], titre)
lines(x = c(0, nbElt + 1), y = rep(cutline, 2), col = "blue")
lines(x = c(0, nbElt + 1),
      y = rep(last(res.Fobs$tR2[ , "R2prd"]), 2), col = "red")


R2toCVmse(res.Fobs$tR2[dim(res.Fobs$tR2)[1],"R2prd"], Fobs)

close.pdf(enregistre)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
graphics.off()


nbElt     <- 10
partition <- c(nbElt)
affectElt <- set.partition(partition)

nbElt     <- sum(partition)
nbclElt   <- length(partition)

mOccur <- diag(1, nrow = nbElt)
colnames(mOccur) <- names(affectElt)
rownames(mOccur) <- seq(1, nbElt)
nbAss <- dim(mOccur)[1]

Fobs   <- seq(1, 10, by = 1)
error  <- 1
noise  <- runif(nbAss, min = -error, max = +error)
#noise  <- rnorm(nbAss, mean = 0, sd = error)
Fobs   <- Fobs + noise

plot(y = Fobs, x = 1:10)
summary(lm(Fobs ~ c(1:10)))

res.cal <- cluster.elements(mOccur, Fobs,
                            opt.meth = "divisive",
                            opt.mean = "amean",
                            opt.mod  = "byelt",
                            opt.ssr  = "cal")

plot.tree(res.cal, col = "black", titre)
plot.tree(res.cal, col = couleurs[gf], titre)

write.txt(enregistre, filename,
          full.names[sort.tree(res.cal$aff, index.return = TRUE)$ix]   )


mAffectElt <- res.cal$aff
mAssMotifs <- mAffect.motifs(res.cal, mOccur)

res.Fobs <- plot.predict.fct(mAssMotifs, mOccur, Fobs,
                             opt.mean = "amean",
                             opt.mod  = "byelt",
                             opt.R2  = TRUE,
                             opt.cal = FALSE,
                             opt.prd = FALSE,
                             opt.glb = FALSE,
                             opt.pub = TRUE,
                             opt.aov = FALSE,
                             pvalue = pvalue,
                             titre = "Fobs")


res     <- res.cal
res$cor <- res.Fobs$tR2[ ,"R2cal"]

nbcl    <- first.optimum(res.Fobs$tR2[, "AICc"], opt = "min")
cutline <- mean(res$cor[(nbcl - 1):nbcl])

plot.tree(res, col = "black", titre)
lines(x = c(0, nbElt + 1), y = rep(cutline, 2), col = "blue")
lines(x = c(0, nbElt + 1),
      y = rep(last(res.Fobs$tR2[ , "R2prd"]), 2), col = "red")

plot.tree(res, col = couleurs[gf], titre)
lines(x = c(0, nbElt + 1), y = rep(cutline, 2), col = "blue")
lines(x = c(0, nbElt + 1),
      y = rep(last(res.Fobs$tR2[ , "R2prd"]), 2), col = "red")


R2toCVmse(res.Fobs$tR2[dim(res.Fobs$tR2)[1],"R2prd"], Fobs)











