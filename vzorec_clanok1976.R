
# Potrebne funkcie: nespravne_phi_taub (bez n)
#                   spravne_phi_taub (bez n)

sigma2_agresti <- function(pp) {
  r <- nrow(pp)
  c <- ncol(pp)
  
  p_i <- rowSums(pp)
  p_j <- colSums(pp)
  
  A1 <- 1 - sum(p_i^2)
  A2 <- 1 - sum(p_j^2)
  
  outC <- 0
  for (i in 1:(r-1)) {
    for (k in (i+1):r) {
      for (j in 1:(c-1)) {
        for (l in (j+1):c) {
          outC <- outC + pp[i,j] * pp[k,l]
        }
      }
    }
  }
  Pc <- 2*outC
  
  outD <- 0
  for (i in 1:(r-1)) {
    for (k in (i+1):r) {
      for (j in 2:c) {
        for (l in 1:(j-1)) {
          outD <- outD + pp[i,j] * pp[k,l]
        }
      }
    }
  }
  Pd <- 2*outD

###############################################################################
# kvadranty rozbite na mensie  
  ### SE quadrant
  sumSE <- matrix(0, r, c)
  for (i in 1:(r-1)) {
    for (j in 1:(c-1)) {
      for (ip in (i+1):r) {
        for (jp in (j+1):c) {
          sumSE[i,j] <- sumSE[i,j] + pp[ip, jp]
        }
      }
    }
  }
  
  ### NW quadrant
  sumNW <- matrix(0, r, c)
  for (i in 2:r) {
    for (j in 2:c) {
      for (ip in 1:(i-1)) {
        for (jp in 1:(j-1)) {
          sumNW[i,j] <- sumNW[i,j] + pp[ip, jp]
        }
      }
    }
  }
  
  ### SW quadrant (i' > i, j' < j)
  sumSW <- matrix(0, r, c)
  for (i in 1:(r-1)) {
    for (j in 2:c) {
      for (ip in (i+1):r) {      # i' > i
        for (jp in 1:(j-1)) {    # j' < j
          sumSW[i,j] <- sumSW[i,j] + pp[ip, jp]
        }
      }
    }
  }
  
  ### NE quadrant (i' < i, j' > j)
  sumNE <- matrix(0, r, c)
  for (i in 2:r) {
    for (j in 1:(c-1)) {
      for (ip in 1:(i-1)) {      # i' < i
        for (jp in (j+1):c) {    # j' > j
          sumNE[i,j] <- sumNE[i,j] + pp[ip, jp]
        }
      }
    }
  }
  
################################################################
  
  tmp1  <- sumSE + sumNW
  Pcc  <- sum(pp * tmp1^2)
  
  tmp2  <- sumSW + sumNE
  Pdd  <- sum(pp * tmp2^2)
  
  Pcd <- sum(pp * (sumSE + sumNW) * (sumSW + sumNE))
  
  part1 <- (Pcc - 2 * Pcd + Pdd) 
  part2 <- (Pc - Pd)^2 * (1/A1 + 1/A2)^2 / 4

################################################################
  
  sum1 <- 0
  for (i in 1:r) {
    for (j in 1:c) {
      sum1 <- sum1 +
        pp[i, j] * (sumSE[i, j] + sumNW[i, j] - sumSW[i, j] - sumNE[i, j])*(p_i[i]/A1 + p_j[j]/A2)
    }
  }
  
  sum2 <- 0
  for (i in 1:r) {
    for (j in 1:c) {
      p_ij <- pp[i, j]
      sum2 <- sum2 + p_ij * (p_i[i]/A1 + p_j[j]/A2)^2
    }
  }
  
  sigma2 <- (4 / (A1*A2)) * ( part1 - part2 + (Pc-Pd)*sum1 + (Pc-Pd)^2 * sum2 / 4 )
  return(sigma2)
}

################################################################################
# Pravdepodobnostna tabulka, na ktorej si vyskusame, ci Agrestiho tvar z clanku 1976
# je rovnaky ako nas odvodeny tvar pre varianciu tau_b pre jeden multinomicky vyber 
r <- 3
c <- 4
pp <- matrix(c(
  0.000, 0.001, 0.004, 0.045,
  0.002, 0.010, 0.030, 0.058,
  0.008, 0.029, 0.066, 0.732
), nrow = 3, byrow = TRUE)

# Vypocet variancii ############################################################

phi_wrong <- nespravne_phi_taub(pp)

row_marg <- rowSums(pp)
col_marg <- colSums(pp)
delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))

# sigma2
sigma2_wrong <- (( sum(pp * phi_wrong^2) - (sum(pp * phi_wrong))^2)/ delta^4 )
sigma2_wrong # Agrestiho kniha (v ktorej sme nasli chybu)

phi_ok <- spravne_phi_taub(pp)
sigma2_ok <- (( sum(pp * phi_ok^2) - (sum(pp * phi_ok))^2)/ delta^4 )
sigma2_ok # nami odvodena variancia

sigma2_agresti(pp) # sigma2 z clanku 1976
