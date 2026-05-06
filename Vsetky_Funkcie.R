# Vsetky esencialne funkcie k diplomovej práci

# C a D ========================================================================
funkcia_concordant_discordant <- function(tab) {
  r <- nrow(tab)
  c <- ncol(tab)
  C <- 0
  D <- 0
  
  for (i in 1:(r - 1)) { # X a Y rastu
    for (j in 1:(c - 1)) {
      for (k in (i + 1):r) {
        for (l in (j + 1):c) {
          C <- C + tab[i, j] * tab[k, l]
        }
      }
    }
  }
  
  for (i in 1:(r - 1)) { # X rastie, Y klesa
    for (j in 2:c) {
      for (k in (i + 1):r) {
        for (l in 1:(j - 1)) {
          D <- D + tab[i, j] * tab[k, l]
        }
      }
    }
  }
  
  return(list(C = C, D = D))
}

# Sigma gamma 1 multinomicky vyber  ============================================
gamma_1hadik <- function(tbl) {
  n <- sum(tbl)
  r <- nrow(tbl)
  c <- ncol(tbl)
  p <- tbl / n
  total_pairs <- choose(n,2)
  
  medzivys <- funkcia_concordant_discordant(tbl)
  C <- medzivys$C
  D <- medzivys$D
  
  Pi_c <- C / total_pairs
  Pi_d <- D / total_pairs
  
  phi <- matrix(0, nrow = r, ncol = c)
  pi_c <- matrix(0, nrow = r, ncol = c)
  pi_d <- matrix(0, nrow = r, ncol = c)
  
  for (i in 1:r) for (j in 1:c) {
    up    <- if (i > 1) 1:(i-1) else integer(0)
    down  <- if (i < r) (i+1):r else integer(0)
    left  <- if (j > 1) 1:(j-1) else integer(0)
    right <- if (j < c) (j+1):c else integer(0)
    pi_c[i,j] <- sum(p[up,  left , drop=FALSE]) + sum(p[down, right, drop=FALSE])
    pi_d[i,j] <- sum(p[up,  right, drop=FALSE]) + sum(p[down, left , drop=FALSE])
    phi[i, j] <- 4 * (Pi_d * pi_c[i,j] - Pi_c * pi_d[i,j])
  }
  
  sigma2_gamma <- 1/n * ( sum(p * phi^2) / ( (Pi_c + Pi_d)^4 ) )
  sigma_gamma <- sqrt(sigma2_gamma)
  
  return(sigma_gamma = sigma_gamma)
}

# Sigma gamma viac multinomickych vyberov ======================================
gamma_3hadici <- function(tbl) {
  n <- sum(tbl)
  r <- nrow(tbl)
  c <- ncol(tbl)
  p <- tbl / n
  n_i <- rowSums(tbl)
  omega <- n_i / n
  
  # C, D, Pi_c, Pi_d, nu, delta
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  total_pairs <- choose(n, 2)
  Pi_c <- C / total_pairs
  Pi_d <- D / total_pairs
  nu    <- Pi_c - Pi_d
  delta <- Pi_c + Pi_d
  
  R_c <- matrix(0,r,c)
  R_d <- matrix(0,r,c)
  for (i in 1:r) {
    for (j in 1:c) {
      # KONKORDANTNÉ páry
      part1_c <- if (i > 1 && j > 1) sum(p[1:(i-1), 1:(j-1)]) else 0
      part2_c <- if (i < r && j < c) sum(p[(i+1):r, (j+1):c]) else 0
      R_c[i, j] <- part1_c + part2_c
      
      # DISKORDANTNÉ páry
      part1_d <- if (i > 1 && j < c) sum(p[1:(i-1), (j+1):c]) else 0
      part2_d <- if (i < r && j > 1) sum(p[(i+1):r, 1:(j-1)]) else 0
      R_d[i, j] <- part1_d + part2_d
    }
  }
  
  # F_ij a φ_ij
  F <- delta * (R_c - R_d) - nu * (R_c + R_d)
  phi <- 2 * omega * F
  
  p_podm <- p / omega 
  
  # sigma
  inner1 <- rowSums(p_podm * (F^2))
  inner2 <- rowSums(p_podm *  F)
  sigma2 <- (4/ (n * delta^4)) * sum( omega * (inner1 - inner2^2))
  sigma  <- sqrt(sigma2)
  
  return(sigma = sigma)
}

# Sigma gamma populacna 1 mult. vyber  =========================================
gamma_1hadik_pop <- function(P) {
  r <- nrow(P)
  c <- ncol(P)
  
  outC <- 0
  for (i in 1:(r-1)) {
    for (k in (i+1):r) {
      for (j in 1:(c-1)) {
        for (l in (j+1):c) {
          outC <- outC + P[i,j] * P[k,l]
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
          outD <- outD + P[i,j] * P[k,l]
        }
      }
    }
  }
  Pi_d <- 2*outD
  
  phi <- matrix(0, nrow = r, ncol = c)
  pi_c <- matrix(0, nrow = r, ncol = c)
  pi_d <- matrix(0, nrow = r, ncol = c)
  
  for (i in 1:r) for (j in 1:c) {
    up    <- if (i > 1) 1:(i-1) else integer(0)
    down  <- if (i < r) (i+1):r else integer(0)
    left  <- if (j > 1) 1:(j-1) else integer(0)
    right <- if (j < c) (j+1):c else integer(0)
    pi_c[i,j] <- sum(P[up,  left , drop=FALSE]) + sum(P[down, right, drop=FALSE])
    pi_d[i,j] <- sum(P[up,  right, drop=FALSE]) + sum(P[down, left , drop=FALSE])
    phi[i, j] <- 4 * (Pi_d * pi_c[i,j] - Pi_c * pi_d[i,j])
  }
  
  sigma2_gamma <-  sum(P * phi^2) / ( (Pi_c + Pi_d)^4 ) 
  sigma_gamma <- sqrt(sigma2_gamma)
  
  return(sigma_gamma = sigma_gamma)
}

# Sigma gamma populacna viac mult. vyberov =====================================
gamma_3hadici_pop <- function(P) {
  r <- nrow(P)
  c <- ncol(P)
  omega <- rowSums(P)
  
  outC <- 0
  for (i in 1:(r-1)) {
    for (k in (i+1):r) {
      for (j in 1:(c-1)) {
        for (l in (j+1):c) {
          outC <- outC + P[i,j] * P[k,l]
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
          outD <- outD + P[i,j] * P[k,l]
        }
      }
    }
  }
  Pi_d <- 2*outD
  
  nu    <- Pi_c - Pi_d
  delta <- Pi_c + Pi_d
  
  R_c <- matrix(0,r,c)
  R_d <- matrix(0,r,c)
  for (i in 1:r) {
    for (j in 1:c) {
      # KONKORDANTNÉ páry
      part1_c <- if (i > 1 && j > 1) sum(P[1:(i-1), 1:(j-1)]) else 0
      part2_c <- if (i < r && j < c) sum(P[(i+1):r, (j+1):c]) else 0
      R_c[i, j] <- part1_c + part2_c
      
      # DISKORDANTNÉ páry
      part1_d <- if (i > 1 && j < c) sum(P[1:(i-1), (j+1):c]) else 0
      part2_d <- if (i < r && j > 1) sum(P[(i+1):r, 1:(j-1)]) else 0
      R_d[i, j] <- part1_d + part2_d
    }
  }
  
  # F_ij a φ_ij
  F <- delta * (R_c - R_d) - nu * (R_c + R_d)
  phi <- 2 * omega * F
  
  p_podm <- P / omega
  
  # sigma
  inner1 <- rowSums(p_podm * (F^2))
  inner2 <- rowSums(p_podm *  F)
  sigma2 <- (4/ (delta^4)) * sum( omega * (inner1 - inner2^2))
  sigma  <- sqrt(sigma2)
  
  return(sigma = sigma)
}

# Sigma taub 1 mult. vyber  ====================================================
had1_taub <- function(tbl) {
  n <- sum(tbl)
  p <- tbl/n
  r <- nrow(p)
  c <- ncol(p)
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  total_pairs <- choose(n,2)
  
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  nu <- (C - D)/total_pairs
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
  # Inicializacia vyslednej matice
  phi <- matrix(0, r, c)
  
  sum_row_sq <- sum(row_marg^2)
  sum_col_sq <- sum(col_marg^2)
  R_c <- matrix(0,r,c)
  R_d <- matrix(0,r,c)
  
  for (i in 1:r) {
    for (j in 1:c) {
      # KONKORDANTNÉ páry
      part1_c <- if (i > 1 && j > 1) sum(p[1:(i-1), 1:(j-1)]) else 0
      part2_c <- if (i < r && j < c) sum(p[(i+1):r, (j+1):c]) else 0
      R_c[i, j] <- part1_c + part2_c
      
      # DISKORDANTNÉ páry
      part1_d <- if (i > 1 && j < c) sum(p[1:(i-1), (j+1):c]) else 0
      part2_d <- if (i < r && j > 1) sum(p[(i+1):r, 1:(j-1)]) else 0
      R_d[i, j] <- part1_d + part2_d
      
      # phi
      A <- 1 - sum(row_marg^2)
      B <- 1 - sum(col_marg^2)
      
      term1 <- -delta*nu*(row_marg[i]/A + col_marg[j]/B)
      term2 <- -delta*2*(R_c[i,j] - R_d[i,j])
      
      phi[i, j] <- term1 + term2
    }
  }
  
  variance <- 1/n *(( sum(p * phi^2) - (sum(p * phi))^2)/ delta^4 )
  sigma <- sqrt(variance)
  
  return(sigma)
}

# Sigma taub viac mult. vyberov  ===============================================
hady3_taub <- function(tbl) {
  n <- sum(tbl)
  n_i <-rowSums(tbl)
  p <- tbl/n
  r <- nrow(tbl)
  c <- ncol(tbl)
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  total_pairs <- choose(n,2)
  
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  Pi_c <- C / total_pairs
  Pi_d <- D / total_pairs
  nu <- (C - D)/total_pairs
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  omega <- n_i/n
  B <- 1 - sum(col_marg^2)
  
  # Inicializacia vyslednej matice
  phi <- matrix(0, r, c)
  R_c <- matrix(0,r,c)
  R_d <- matrix(0,r,c)
  
  for (i in 1:r) {
    for (j in 1:c) {
      # KONKORDANTNÉ páry
      part1_c <- if (i > 1 && j > 1) sum(p[1:(i-1), 1:(j-1)]) else 0
      part2_c <- if (i < r && j < c) sum(p[(i+1):r, (j+1):c]) else 0
      R_c[i, j] <- part1_c + part2_c
      
      # DISKORDANTNÉ páry
      part1_d <- if (i > 1 && j < c) sum(p[1:(i-1), (j+1):c]) else 0
      part2_d <- if (i < r && j > 1) sum(p[(i+1):r, 1:(j-1)]) else 0
      R_d[i, j] <- part1_d + part2_d
      
      # phi
      phi[i,j] <- delta*omega[i]*(2*(R_c[i,j] - R_d[i,j]) + nu*(col_marg[j]/B))
    }
  }
  
  # Sigma2
  theta <- matrix(0, r, c) # theta[i,j] = pi_{j|i}
  T <- matrix(0, r, c)
  for (i in 1:r) {
    for (j in 1:c) {
      theta[i, j] <- tbl[i, j] / n_i[i]
      
      T[i,j] <- 2*(R_c[i,j] - R_d[i,j]) + nu * col_marg[j]/B
    }
  }
  
  inner1 <- rowSums(theta * (T^2))
  inner2 <- rowSums(theta * T)
  sigma2 <- 1/n * (1 / delta^2) * sum(omega * (inner1 - inner2^2))
  sigma <- sqrt(sigma2)
  
  return(sigma)
}

# Nase spravne phi (bez n)  ====================================================
spravne_phi_taub <- function(ppAB) {
  r <- nrow(ppAB)
  c <- ncol(ppAB)
  row_marg <- rowSums(ppAB)
  col_marg <- colSums(ppAB)
  
  outC <- 0
  for (i in 1:(r-1)) {
    for (k in (i+1):r) {
      for (j in 1:(c-1)) {
        for (l in (j+1):c) {
          outC <- outC + ppAB[i,j] * ppAB[k,l]
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
          outD <- outD + ppAB[i,j] * ppAB[k,l]
        }
      }
    }
  }
  Pi_d <- 2*outD
  
  nu <- Pi_c - Pi_d
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  A <- 1 - sum(row_marg^2)
  B <- 1 - sum(col_marg^2)
  
  # Inicializacia vyslednej matice
  phi <- matrix(0, r, c)
  pi_c <- matrix(0, r, c)
  pi_d <- matrix(0, r, c)
  for (i in 1:r) {
    for (j in 1:c) {
      ## π_ij^(c)
      part1_c <- if (i > 1 && j > 1) sum(ppAB[1:(i-1), 1:(j-1)]) else 0   # a < i, b < j
      part2_c <- if (i < r && j < c) sum(ppAB[(i+1):r, (j+1):c]) else 0   # a > i, b > j
      pi_c[i, j] <- part1_c + part2_c
      
      ## π_ij^(d)
      part1_d <- if (i > 1 && j < c) sum(ppAB[1:(i-1), (j+1):c]) else 0   # a < i, b > j
      part2_d <- if (i < r && j > 1) sum(ppAB[(i+1):r, 1:(j-1)]) else 0   # a > i, b < j
      pi_d[i, j] <- part1_d + part2_d
      
      ## phi
      term1 <- (nu/delta) * (row_marg[i]*B + col_marg[j]*A)
      term2 <- delta * 2 * (pi_c[i, j] - pi_d[i, j])
      
      phi[i, j] <- term1 + term2
    }
  }
  
  return(phi)
}

# Agrestiho nespravne phi (bez n)  =============================================
nespravne_phi_taub <- function(ppAB) {
  r <- nrow(ppAB)
  c <- ncol(ppAB)
  row_marg <- rowSums(ppAB)
  col_marg <- colSums(ppAB)
  
  outC <- 0
  for (i in 1:(r-1)) {
    for (k in (i+1):r) {
      for (j in 1:(c-1)) {
        for (l in (j+1):c) {
          outC <- outC + ppAB[i,j] * ppAB[k,l]
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
          outD <- outD + ppAB[i,j] * ppAB[k,l]
        }
      }
    }
  }
  Pi_d <- 2*outD
  
  nu <- Pi_c - Pi_d
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  A <- 1 - sum(row_marg^2)
  B <- 1 - sum(col_marg^2)
  
  # Inicializacia vyslednej matice
  phi <- matrix(0, r, c)
  pi_c <- matrix(0,r,c)
  pi_d <- matrix(0,r,c)
  
  for (i in 1:r) {
    for (j in 1:c) {
      ## π_ij^(c)
      part1_c <- if (i > 1 && j > 1) sum(ppAB[1:(i-1), 1:(j-1)]) else 0   # a < i, b < j
      part2_c <- if (i < r && j < c) sum(ppAB[(i+1):r, (j+1):c]) else 0   # a > i, b > j
      pi_c[i, j] <- part1_c + part2_c
      
      ## π_ij^(d)
      part1_d <- if (i > 1 && j < c) sum(ppAB[1:(i-1), (j+1):c]) else 0   # a < i, b > j
      part2_d <- if (i < r && j > 1) sum(ppAB[(i+1):r, 1:(j-1)]) else 0   # a > i, b < j
      pi_d[i, j] <- part1_d + part2_d
      
      # phi
      term1 <- nu*(row_marg[i]*sqrt(B/A) + col_marg[j]*sqrt(A*B))
      term2 <- delta*2*(pi_c[i,j] - pi_d[i,j])
      
      phi[i, j] <- term1 + term2
    }
  }
  
  return(phi)
}

# Nase spravne phi (s n)  ======================================================
spravne_phi_taub <- function(tbl) {
  n <- sum(tbl)
  p <- tbl/n
  r <- nrow(p)
  c <- ncol(p)
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  total_pairs <- choose(n,2)
  
  medzivys <- funkcia_concordant_discordant(tbl)
  C <- medzivys$C
  D <- medzivys$D
  
  Pi_c <- C / total_pairs
  Pi_d <- D / total_pairs
  nu <- Pi_c - Pi_d
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
  # Inicializacia vyslednej matice
  phi <- matrix(0, r, c)
  
  sum_row_sq <- sum(row_marg^2)
  sum_col_sq <- sum(col_marg^2)
  pi_c <- matrix(0,r,c)
  pi_d <- matrix(0,r,c)
  
  for (i in 1:r) {
    for (j in 1:c) {
      pi_c[i,j] <- sum(p[0:(i-1), 0:(j-1)]) + tryCatch(sum(p[(i+1):r, (j+1):c]), error= function(e) 0)
      pi_d[i,j] <- tryCatch(sum(p[0:(i-1), (j+1):c]), error= function(e) 0) + tryCatch(sum(p[(i+1):r, 0:(j-1)]), error= function(e) 0)
      
      # phi
      A <- 1 - sum(row_marg^2)
      B <- 1 - sum(col_marg^2)
      
      term1 <- -delta*nu*(row_marg[i]/A + col_marg[j]/B)
      term2 <- -delta*2*(pi_c[i,j] - pi_d[i,j])
      
      phi[i, j] <- term1 + term2
    }
  }
  
  return(phi)
}

# Agrestiho nespravne phi (s n)  ===============================================
nespravne_phi_taub <- function(tbl) {
  n <- sum(tbl)
  p <- tbl/n
  r <- nrow(p)
  c <- ncol(p)
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  total_pairs <- choose(n,2)
  
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  Pi_c <- C / total_pairs
  Pi_d <- D / total_pairs
  nu <- Pi_c - Pi_d
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
  # Inicializacia vyslednej matice
  phi <- matrix(0, r, c)
  
  pi_c <- matrix(0,r,c)
  pi_d <- matrix(0,r,c)
  
  for (i in 1:r) {
    for (j in 1:c) {
      pi_c[i,j] <- sum(p[0:(i-1), 0:(j-1)]) + tryCatch(sum(p[(i+1):r, (j+1):c]), error= function(e) 0)
      pi_d[i,j] <- tryCatch(sum(p[0:(i-1), (j+1):c]), error= function(e) 0) + tryCatch(sum(p[(i+1):r, 0:(j-1)]), error= function(e) 0)
      
      # phi
      A <- 1 - sum(row_marg^2)
      B <- 1 - sum(col_marg^2)
      
      term1 <- -nu*(row_marg[i]*sqrt(B/A) + col_marg[j]*sqrt(A*B))
      term2 <- -delta*2*(pi_c[i,j] - pi_d[i,j])
      
      phi[i, j] <- term1 + term2
    }
  }
  
  return(phi)
}

# Rho funkcia ==================================================================
rho_hat <- function(tab, u_scores = NULL, v_scores = NULL) {
  n <- sum(tab)
  p <- tab / n
  r <- nrow(tab)
  c <- ncol(tab)
  
  if (is.null(u_scores)) u_scores <- 1:r
  if (is.null(v_scores)) v_scores <- 1:c
  
  p_i <- rowSums(p)
  p_j <- colSums(p)
  
  u_bar <- sum(u_scores * p_i)
  v_bar <- sum(v_scores * p_j)
  
  num <- 0
  for (i in 1:r) {
    for (j in 1:c) {
      num <- num + (u_scores[i] - u_bar) * (v_scores[j] - v_bar) * p[i, j]
    }
  }
  
  den_u <- sum((u_scores - u_bar)^2 * p_i)
  den_v <- sum((v_scores - v_bar)^2 * p_j)
  
  rho <- num / sqrt(den_u * den_v)
  return(rho)
}

# Rho_b funkcia ================================================================
rho_b_hat <- function(tab) {
  n <- sum(tab)
  p <- tab / n
  r <- nrow(p)
  c <- ncol(p)
  
  row_marg <- rowSums(p)  # p_{i+}
  col_marg <- colSums(p)  # p_{+j}
  
  # ridit skóre
  a_X <- numeric(r)
  for (i in 1:r) {
    sum_pred <- if (i == 1) 0 else sum(row_marg[1:(i - 1)])
    a_X[i] <- sum_pred + row_marg[i] / 2
  }
  
  a_Y <- numeric(c)
  for (j in 1:c) {
    sum_pred <- if (j == 1) 0 else sum(col_marg[1:(j - 1)])
    a_Y[j] <- sum_pred + col_marg[j] / 2
  }
  
  # čitateľ
  citatel <- 0
  for (i in 1:r) {
    for (j in 1:c) {
      citatel <- citatel + (a_X[i] - 0.5) * (a_Y[j] - 0.5) * p[i, j]
    }
  }
  
  # menovatele
  menovatel1 <- sum((a_X - 0.5)^2 * row_marg)
  menovatel2 <- sum((a_Y - 0.5)^2 * col_marg)
  
  # výsledok
  rho_b <- citatel / sqrt(menovatel1 * menovatel2)
  return(rho_b)
}
