# Potrebne funkcie: funkcia_concordant_discordant
#                   gamma_1hadik

# Porovnanie Gammy vs. Zeta transformacie pomocou pravdepodobnosti chyby 1. druhu - GRAF
# --- graf zavislosti chyby 1. druhu od n_jedincov ---

Gama_vs_Zeta <- function(nj, n_sim, r, c, pp, alpha = 0.05) {
  rej_g <- logical(n_sim)
  rej_z <- logical(n_sim)
  for (s in 1:n_sim) {
    vektor <- rmultinom(1, nj, pp)
    tbl <- matrix(vektor, r, c, byrow = TRUE)
    n <- sum(tbl)                 
    p <- tbl / n                   
    
    # Vypocet gammy a sigma_gamma
    vysledky <- funkcia_concordant_discordant(tbl)
    C <- vysledky$C 
    D <- vysledky$D
    
    gamma <- (C - D) / (C + D)
    
    sigma_gamma <- gamma_1hadik(tbl)
    
    # Z-test pre gamma
    z1 <- gamma / sigma_gamma
    p_value1 <- 2 * (1 - pnorm(abs(z1)))
    
    # Zeta transformacia
    zetaHAT <- 0.5 * log((1 + gamma) / (1 - gamma))
    zeta_var <- sigma_gamma^2 * (1 - gamma^2)^(-2)
    
    # Z-test pre zeta
    z2 <- zetaHAT / sqrt(zeta_var)
    p_value2 <- 2 * (1 - pnorm(abs(z2)))
    
    rej_g[s] <- p_value1 < alpha
    rej_z[s] <- p_value2 < alpha
  }
  c(alpha_gamma = mean(rej_g), alpha_zeta = mean(rej_z))
}

# pocty jedincov - grid
n_grid <- c(30, 40, 60, 80, 100, 150, 200)
n_sim <- 10000
r <- 3
c <- 4

pp <- rep(1/12, 12)

# pre kazdy bod na gride (nj) spocitame alpha (zamietnutie)
alpha_mat <- sapply(n_grid, function(nj) Gama_vs_Zeta(nj, n_sim, r, c, pp))

# graf
ylim_all <- range(c(alpha_mat, 0.05), na.rm = TRUE)
plot(n_grid, alpha_mat["alpha_gamma", ], type = "o", pch = 16,
     ylim = ylim_all, xlab = "n",
     ylab = "Chyba 1. druhu",
     main = "Gamma vs. Zeta")
lines(n_grid, alpha_mat["alpha_zeta", ], type = "o", pch = 1, lty = 2)
abline(h = 0.05, lty = 3)
legend("topright",
       legend = c("Gamma test", "Zeta transformácia", "α = 0.05"),
       lty = c(1, 2, 3), pch = c(16, 1, NA), bty = "n")

# tabulka ku grafu
tab <- data.frame(
  n_jedincov   = n_grid,
  alpha_gamma  = as.numeric(alpha_mat["alpha_gamma", ]),
  alpha_zeta   = as.numeric(alpha_mat["alpha_zeta", ])
)
tab_round <- transform(tab,
                       alpha_gamma = round(alpha_gamma, 3),
                       alpha_zeta  = round(alpha_zeta, 3)
)
print(tab_round)
