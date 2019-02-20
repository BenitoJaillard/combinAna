


#==============================================================================
# Programme de calcul des volumes nécessaires pour obtenir une dilution donnée
#==============================================================================

# excipient
cexcip <- 0
names(cexcip) <- "miel"

# solutions-mère
compo <- c("permetrin", "pyperonyl", "alcool")
ci1 <- c(100, 100, 95)  # en g/L et %
ci2 <- c(100,   0, 95)  # en g/L et %
names(ci1) <- names(ci2) <- compo

# volume et concentrations finales
vt <- 50          # en mL
cf <- c(1, 1, 5)  # en g/L

# calcul
ci <- ci1
vi <- vt * cf/ci  # en mL
vexcip <- vt - sum(vi)

# résultats en volume
res <- c(vi, vexcip)
names(res) <- paste("V(", c(names(ci), names(cexcip)), ") (mL)", sep = "")
print(signif(c(res, sum(res)), 3))

# résultats en gouttes
res <- 20 * res
names(res) <- paste("V(", c(names(ci), names(cexcip)), ") (gouttes)", sep = "")
print(signif(c(res, sum(res)), 3))


#==============================================================================
# Programme de calcul des volumes nécessaires pour obtenir une dilution donnée
#   Manips du 26 Juillet 2016
#==============================================================================
biocide.mere.1 <- 0.2229/10  # en g par goutte
biocide.mere.2 <- 0.2703/10  # en g par goutte

# alcool
alcool.ci1 <- 0.95  # en %
alcool.cf1 <- 0.11  # en %


Vti <- c(13.03, 13.07, 13.03, 13.17, 13.39)           # g de miel
Vti <- c(12.96, 13.16, 13.36, 13.33, 13.20)           # g de miel

alcool.v1 <- alcool.cf1 / (alcool.ci1 - alcool.cf1) * Vti

alcool.v1 <- c(2.006, 1.704, 1.708, 1.658, 1.604)
alcool.v1 <- c(1.722, 1.752, 1.678, 1.685, 1.955)


# biocide
biocide.ci1 <- 100
biocide.cf1 <- c(1, 2, 4, 6, 8)

biocide.v1 <- biocide.cf1 / (biocide.ci1 - biocide.cf1) * (Vti + alcool.v1)

biocide.v1 <- c(0.158, 0.385, 0.622, 0.989, 1.274)
biocide.v1 <- c(0.129, 0.329, 0.6114, 0.9931, 1.254)


(alcool.v1 / (Vti + alcool.v1 + biocide.v1))
(biocide.v1 / (Vti + alcool.v1 + biocide.v1))
(alcool.v1 + biocide.v1) / (Vti + alcool.v1 + biocide.v1)


(biocide.v1 / biocide.mere.1) # en gouttes
(biocide.v1 / biocide.mere.2) # en gouttes


#==============================================================================
# Programme de calcul des volumes nécessaires pour obtenir une dilution donnée
#   Manips du 28 Juillet 2016
#==============================================================================
biocide.mere.1 <- 0.2229/10  # en g par goutte
biocide.mere.2 <- 0.2703/10  # en g par goutte

# alcool
alcool.ci1 <- 0.95  # en %
alcool.cf1 <- 0.55  # en %


Vti <- 22.9299 - 9.6640           # g de miel

dv <- (alcool.ci1 - alcool.cf1) / alcool.cf1 * Vti

dv

alcool.v1
alcool.v1 - Vti

alcool.v1 <- c(2.006, 1.704, 1.708, 1.658, 1.604)
alcool.v1 <- c(1.722, 1.752, 1.678, 1.685, 1.955)


# biocide
biocide.ci1 <- 100
biocide.cf1 <- c(0.1, 0.2, 0.4, 0.6, 0.8)

biocide.v1 <- biocide.cf1 / (biocide.ci1 - biocide.cf1) * (Vti + alcool.v1)
biocide.v1

biocide.v1 <- c(0.158, 0.385, 0.622, 0.989, 1.274)
biocide.v1 <- c(0.129, 0.329, 0.6114, 0.9931, 1.254)


(alcool.v1 / (Vti + alcool.v1 + biocide.v1))
(biocide.v1 / (Vti + alcool.v1 + biocide.v1))
(alcool.v1 + biocide.v1) / (Vti + alcool.v1 + biocide.v1)


(biocide.v1 / biocide.mere.1) # en gouttes
(biocide.v1 / biocide.mere.2) # en gouttes


#==============================================================================
# Programme de calcul des volumes nécessaires pour obtenir une dilution donnée
#   Manips du 28 Juillet 2016
#==============================================================================
Vti <- 38.2472 - 9.0638

alcool.ci1 <- 0.75
alcool.cf1 <- 0.08

alcool.v1 <- alcool.cf1 / (alcool.ci1 - alcool.cf1) * Vti

alcool.v1


# biocide
biocide.ci1 <- 100
biocide.cf1 <- 1.2

biocide.v1 <- biocide.cf1 / (biocide.ci1 - biocide.cf1) * (Vti + alcool.v1)
biocide.v1

#miel à 18.5 d'eau

revient <- 10:20
prix <- revient * 5
revendeur <- 1.25 * prix
tva <- revendeur * 1.195
plot (y = prix, x = revient)
benefice <- 4 * revient
remboursement <- 230000 / benefice
round(remboursement)


conc   <- 0.8   # g/L
litre  <- 1000  # g 
goutte <- 0.025 / litre # L
qgoutte <- conc * goutte # g/L * L
tuegoutte <- qgoutte / (100 * 10^-9)

qgoutte * 1000000
tuegoutte








