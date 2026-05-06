# Potrebne funkcie: funkcia_concordant_discordant
#                   spravne_phi_taub (s n)
#                   nespravne_phi_taub (s n)

library(tibble)
library(scales)
library(ggplot2)
library(dplyr)

n_sim <- 10000
n_jedinci <- 300
r <- 3
c <- 4
n_kategorie <- r * c

# Vplyv chyby na pravdepodobnost chyby prveho druhu - H0 neplati
# ppA <- matrix(
#   c(
#     0.14, 0.18, 0.16, 0.12,
#     0.05, 0.09, 0.08, 0.08,
#     0.01, 0.01, 0.06, 0.02
#   ),
#   nrow = 3,
#   byrow = TRUE
# ) # tau_b = 0.1224932; A7

ppA <- matrix(c(
  0.18, 0.12, 0.06, 0.04,
  0.10, 0.14, 0.04, 0.02,
  0.06, 0.08, 0.10, 0.06
), nrow = 3, byrow = TRUE) # tau_b = 0.2260211; A8

z_A <- numeric(n_sim) # wrong
rej_A <- logical(n_sim)

z_AA <- numeric(n_sim) # ok
rej_AA <- logical(n_sim)

for (s in 1:n_sim) {
  # Generovanie tabulky
  vector <- rmultinom(1, n_jedinci, prob = as.vector(t(ppA)))  
  tbl <- matrix(vector, r, c, byrow = TRUE)
  
  n <- sum(tbl)
  p <- tbl/n
  total_pairs <- choose(n, 2)
  
  # C a D
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  Pi_c <- C / total_pairs
  Pi_d <- D / total_pairs
  
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
  # tau b
  taub <- (Pi_c - Pi_d) / delta
  
  # phi
  phi_ok <- spravne_phi_taub(tbl)
  phi_wrong <- nespravne_phi_taub(tbl)
  
  # Sigma
  sigma2_ok <- 1/n *(( sum(p * phi_ok^2) - (sum(p * phi_ok))^2)/ delta^4 )
  sigma2_wrong <- 1/n *(( sum(p * phi_wrong^2) - (sum(p * phi_wrong))^2)/ delta^4 )
  
  sigma_ok <- sqrt(sigma2_ok)
  sigma_wrong <- sqrt(sigma2_wrong)
  
  # Z-test pre taub = 0
  z_A[s] <- taub / sigma_wrong
  p_value_A <- 2 * (1 - pnorm(abs(z_A[s])))
  rej_A[s] <- p_value_A < 0.05
  
  z_AA[s] <- taub / sigma_ok
  p_value_AA <- 2 * (1 - pnorm(abs(z_AA[s])))
  rej_AA[s] <- p_value_AA < 0.05
}

k_A <- sum(rej_A, na.rm = TRUE) # pocet zamietnuti
n_A <- sum(!is.na(rej_A))
k_AA <- sum(rej_AA, na.rm = TRUE)
n_AA <- sum(!is.na(rej_AA))

alpha_A <- k_A / n_A
alpha_AA <- k_AA / n_AA

df_alpha <- tibble(
  scenár = c("nespravne phi", "spravne phi"),
  alpha  = c(alpha_A, alpha_AA)
)

ggplot(df_alpha, aes(x = scenár, y = alpha)) +
  geom_col(width = 0.6, alpha = 0.7, fill = c("red", "green")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Vplzv chyby vo phi na silu testu",
       x = NULL, y = "Odhad α") +
  theme_minimal(base_size = 13) + geom_text(aes(label = percent(alpha, accuracy = 0.01)), 
                                            vjust = -0.5, size = 4)

############################################# Relativna chyba medzi sigmami graf
n_sim <- 10000
n_jedinci <- 300
r <- 3
c <- 4

# ppA <- matrix(1/12, nrow = 3, ncol = 4) # rovnomerna pp - H0 plati

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

# ppA <- matrix(c(
#   0,   0,      0,      0.1,
#   0,   0.1,    0.2,    0,
#   0.6, 0,      0,      0
# ), nrow=3, byrow=TRUE) # H0 neplati; A3; tau_b = -0.9649013; relat. chyba popul. = 57,63%

# ppA <- matrix(c(
#   0.6, 0,      0,      0,
#   0,   0,      0,      0.1,
#   0,   0.1,    0.2,    0
# ), nrow=3, byrow=TRUE) # A4; tau_b = 0.7504788; relat. chyba popul. = 34,69%

# ppA <- matrix(c(
#   0.025,   0.025,      0.025,  0.025,
#   0.05,   0.05,    0.1,    0.1,
#   0.6, 0,      0,      0
# ), nrow=3, byrow=TRUE) # relat. chyba popul. = 0.05511434 tau_b = -0.7258865; A6

ppA <- matrix(c(
  1,  2,  2,  0,
  1, 11,  1,  0,
  0,  3,  2,  7
), nrow = 3, byrow = TRUE) / 30 # tau_b = 0.5, popul relat. chyba = 11% - dobre aj hadiky

# ulozime si sigmy + rel chyby
sigma_ok_vec    <- numeric(n_sim)
sigma_wrong_vec <- numeric(n_sim)
rel_err_vec     <- numeric(n_sim)

for (s in 1:n_sim) {
  vector <- rmultinom(1, n_jedinci, prob = as.vector(t(ppA)))  
  tbl <- matrix(vector, r, c, byrow = TRUE)
  
  n <- sum(tbl)
  p <- tbl / n
  
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
  # phi
  phi_ok <- spravne_phi_taub(tbl)
  phi_wrong <- nespravne_phi_taub(tbl)
  
  # sigma^2 (ak by delta bolo 0, tak NA)
  if (is.na(delta) || delta == 0) {
    sigma_ok <- NA_real_
    sigma_wrong <- NA_real_
  } else {
    sigma2_ok <- (1/n) * ((sum(p * phi_ok^2) - (sum(p * phi_ok))^2) / delta^4)
    sigma2_wrong <- (1/n) * ((sum(p * phi_wrong^2) - (sum(p * phi_wrong))^2) / delta^4)
    
    sigma_ok <- sqrt(sigma2_ok)
    sigma_wrong <- sqrt(sigma2_wrong)
  }
  
  sigma_ok_vec[s] <- sigma_ok
  sigma_wrong_vec[s] <- sigma_wrong
  
  # relatívna chyba (ak sigma_ok = 0 alebo NA -> NA)
  rel_err_vec[s] <- ifelse(is.na(sigma_ok) || sigma_ok == 0, NA_real_,
                           (sigma_wrong - sigma_ok) / sigma_ok)
}

df_rel <- tibble(
  sigma_ok = sigma_ok_vec,
  sigma_wrong = sigma_wrong_vec,
  rel_err = rel_err_vec,
  abs_rel_err = abs(rel_err_vec)
)

ggplot(df_rel, aes(x = rel_err)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 0.01)) +
  labs(
    title = "Relatívna chyba",
    x = "(σ_wrong − σ_ok) / σ_ok",
    y = "Počet simulácií"
  ) +
  theme_minimal(base_size = 13)
