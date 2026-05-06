# Funkcie: funkcia_concordant_discordant
#          had1_taub
#          nespravne_phi_taub (s n)

library(DescTools)
DescTools::KendallTauB # chyba vo phi
getAnywhere(ConDisPairs)

r <- 3
c <- 4
n_jedinci <- 300

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

ppA <- matrix(
  c(
    0.30, 0.20, 0.05, 0.05,   # sum = 0.6
    0.05, 0.05, 0.10, 0.10,   # sum = 0.3
    0.01, 0.03, 0.03, 0.03    # sum = 0.1
  ),
  nrow = 3,
  byrow = TRUE
) # neplati H0

vector <- rmultinom(1, n_jedinci, ppA)
tbl <- matrix(vector, r, c)

# KendallTauB funkcia: sigma^2
tab <- tbl
x <- ConDisPairs(tab)
n <- sum(tab)
n0 <- n * (n - 1)/2
ti <- rowSums(tab)
uj <- colSums(tab)
n1 <- sum(ti * (ti - 1)/2)
n2 <- sum(uj * (uj - 1)/2)
taub <- (x$C - x$D)/sqrt((n0 - n1) * (n0 - n2))
pi <- tab/sum(tab)
pdiff <- (x$pi.c - x$pi.d)/sum(tab)
Pdiff <- 2 * (x$C - x$D)/sum(tab)^2
rowsum <- rowSums(pi)
colsum <- colSums(pi)
rowmat <- matrix(rep(rowsum, dim(tab)[2]), ncol = dim(tab)[2])
colmat <- matrix(rep(colsum, dim(tab)[1]), nrow = dim(tab)[1], 
                 byrow = TRUE)
delta1 <- sqrt(1 - sum(rowsum^2))
delta2 <- sqrt(1 - sum(colsum^2))
tauphi <- (2 * pdiff + Pdiff * colmat) * delta2 * delta1 + (Pdiff * rowmat * delta2)/delta1
sigma2_funkciaKendall <- ((sum(pi * tauphi^2) - sum(pi * tauphi)^2)/(delta1 * delta2)^4)/n

# Sigma taub 1 mult. vyber: SPRAVNE phi  =======================================
g <- had1_taub(tbl)
Nasa_sigma2_spravna <- g^2

# Sigma taub 1 mult. vyber: NESPRAVNE phi  =====================================
n <- sum(tbl)
p <- tbl/n
total_pairs <- choose(n, 2)

# C a D
cd <- funkcia_concordant_discordant(tbl)
C <- cd$C
D <- cd$D

Pi_c <- C/ total_pairs
Pi_d <- D/ total_pairs

row_marg <- rowSums(p)
col_marg <- colSums(p)
delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))

# tau b
taub <- (Pi_c - Pi_d) / delta

phi_wrong <- nespravne_phi_taub(tbl)

sigma2_wrong <- 1/n *(( sum(p * phi_wrong^2) - (sum(p * phi_wrong))^2)/ delta^4 )

sigma2_funkciaKendall
sigma2_wrong
Nasa_sigma2_spravna

sqrt(sigma2_funkciaKendall)
sqrt(sigma2_wrong)
sqrt(Nasa_sigma2_spravna)

relativna_chyb <- (sqrt(Nasa_sigma2_spravna)-sqrt(sigma2_wrong))/sqrt(Nasa_sigma2_spravna)
relativna_chyb
