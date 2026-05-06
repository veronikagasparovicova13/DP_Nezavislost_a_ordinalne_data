# Funkcie: rho_hat
#          rho_b_hat

# Porovnanie Rho vs. Rho_b pomocou chyby 1. druhu - GRAF
# --- graf zavislosti pravdepodobnosti chyby 1. druhu od n_jedincov ---

Rho_vs_Rhob <- function(nj, n_sim, r, c, pp, alpha = 0.05) {
  rej_rho  <- logical(n_sim)
  rej_rhob <- logical(n_sim)
  
  for (s in 1:n_sim) {
    vektor <- rmultinom(1, nj, pp)
    tbl <- matrix(vektor, r, c, byrow = TRUE)
    n <- sum(tbl)
    
    # odhady korelácie
    rhoHAT   <- rho_hat(tbl)
    rhoBHAT  <- rho_b_hat(tbl)
    
    # χ² štatistika a p-hodnoty
    X2_rho  <- (n - 1) * rhoHAT^2
    p_rho   <- 1 - pchisq(X2_rho, df = 1)
    
    X2_rhob <- (n - 1) * rhoBHAT^2
    p_rhob  <- 1 - pchisq(X2_rhob, df = 1)
    
    rej_rho[s]  <- p_rho  < alpha
    rej_rhob[s] <- p_rhob < alpha
  }
  c(alpha_rho = mean(rej_rho), alpha_rho_b = mean(rej_rhob))
}

n_grid <- seq(15, 100, by = 5)

pp <- rep(1/12, 12)
n_sim <- 10000
r <- 3
c <- 4

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# pp <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)
# pp <- as.vector(t(pp))

# výpočet pre všetky n_jedincov
alpha_mat_rho <- sapply(n_grid, function(nj) Rho_vs_Rhob(nj, n_sim, r, c, pp))

# --- GRAFICKE POROVNANIE ---
ylim_rho <- range(c(alpha_mat_rho, 0.05), na.rm = TRUE)

plot(n_grid, alpha_mat_rho["alpha_rho", ], type = "o", pch = 16,
     ylim = ylim_rho, xlab = "n",
     ylab = expression("Odhad" ~ alpha),
     main = expression(hat(rho) ~ " vs. " ~ hat(rho)[b]))

lines(n_grid, alpha_mat_rho["alpha_rho_b", ], type = "o", pch = 1, lty = 2)

abline(h = 0.05, lty = 3)

legend("bottomright",
       legend = c(expression(hat(rho)),
                  expression(hat(rho)[b]),
                  expression(alpha == 0.05)),
       lty = c(1, 2, 3), pch = c(16, 1, NA), bty = "n")

# --- TABULKA S HODNOTAMI ---
tab_rho <- data.frame(
  n_jedincov = n_grid,
  alpha_rho = round(alpha_mat_rho["alpha_rho", ], 3),
  alpha_rho_b = round(alpha_mat_rho["alpha_rho_b", ], 3)
)

print(tab_rho)
