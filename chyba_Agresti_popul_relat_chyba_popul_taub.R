
# Potrebne funkcie: spravne_phi_taub (bez n)
#                   nespravne_phi_taub (bez n)

# Populacna relativna chyba medzi sigmami ######################################
r <- 3
c <- 4

# phi
phi_ok <- spravne_phi_taub(ppA)
phi_wrong <- nespravne_phi_taub(ppA)
  
row_marg <- rowSums(ppA)
col_marg <- colSums(ppA)
delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
# sigma2
sigma2_ok <- (( sum(ppA * phi_ok^2) - (sum(ppA * phi_ok))^2)/ delta^4 )
sigma2_wrong <- (( sum(ppA * phi_wrong^2) - (sum(ppA * phi_wrong))^2)/ delta^4 )

# rozdiel medzi sigmami^2
diff_sigma <- sigma2_wrong - sigma2_ok
diff_sigma

# relativna chyba sigmy^2
diff_sigma/sigma2_ok

# Tau_b populacne
r <- 3
c <- 4
ppA <- matrix(ppA, nrow = r, ncol = c)
outC <- 0
for (i in 1:(r-1)) {
  for (k in (i+1):r) {
    for (j in 1:(c-1)) {
      for (l in (j+1):c) {
        outC <- outC + ppA[i,j] * ppA[k,l]
      }
    }
  }
}
Pi_c <- 2*outC
outD <- 0
for (i in 1:(r-1)) {
  for (k in (i+1):r) {
    for (j in 2:c) {
      for (l in 1:(j-1)) {
        outD <- outD + ppA[i,j] * ppA[k,l]
      }
    }
  }
}
Pi_d <- 2*outD
nu    <- Pi_c - Pi_d
row_marg <- rowSums(ppA)
col_marg <- colSums(ppA)
tau_b <- nu/sqrt((1-sum(row_marg^2))*(1-sum(col_marg^2)))
tau_b
