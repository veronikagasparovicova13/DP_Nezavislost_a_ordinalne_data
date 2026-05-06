library(dplyr)
library(tibble)
library(ggplot2)
library(scales)

# Funkcie: funkcia_concordant_discordant
#          had1_taub
#          gamma_1hadik
#          rho_hat
#          rho_b_hat

# Parametre
n_sim <- 10000
n_jedinci <- 30
r <- 3
c <- 4

# Pravdepodobnostne tabulky
# ppA <- rep(1/(r*c), r*c) # rovnomerne pp. rozdelenie

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

ppA <- matrix(c(
  0.04, 0.028, 0.018, 0.014,   # sum = 0.10
  0.075, 0.095, 0.075, 0.055,  # sum = 0.30
  0.12, 0.147, 0.177, 0.156    # sum = 0.60
), nrow = 3, byrow = TRUE) # tau b = 0.1430873

z_tau_b <- numeric(n_sim)
rej_tau_b <- logical(n_sim)

z_gamma <- numeric(n_sim)
rej_gamma <- logical(n_sim)

X2_rho <- numeric(n_sim)
rej_rho <- logical(n_sim)

X2_rho_b <- numeric(n_sim)
rej_rho_b <- logical(n_sim)

for (s in 1:n_sim) {
  # Generovanie tabulky
  vector <- rmultinom(1, n_jedinci, ppA)
  tbl <- matrix(vector, r, c)
  
  n <- sum(tbl)
  total_pairs <- choose(n, 2)
  p <- tbl/n
  
  # C a D
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  Pi_c <- C/total_pairs
  Pi_d <- D/total_pairs
  
  row_sucty <- rowSums(p)
  col_sucty <- colSums(p)
  
  # tau_b, gamma, rho, rho_b
  taub <- (Pi_c - Pi_d) / sqrt((1-sum(row_sucty^2)) * (1-sum(col_sucty^2)))
  gamma <- (C - D) / (C + D)
  rho <- rho_hat(tbl)
  rho_b <- rho_b_hat(tbl)
  
  # Sigma
  sigma_tau_b <- had1_taub(tbl)
  sigma_gamma <- gamma_1hadik(tbl)
  
  # Z-test
  z_tau_b[s] <- taub / sigma_tau_b
  p_value_tau_b <- 2 * (1 - pnorm(abs(z_tau_b[s])))
  rej_tau_b[s] <- p_value_tau_b < 0.05
  
  z_gamma[s] <- gamma / sigma_gamma
  p_value_gamma <- 2 * (1 - pnorm(abs(z_gamma[s])))
  rej_gamma[s] <- p_value_gamma < 0.05
  
  X2_rho[s] <- (n - 1) * rho^2
  p_value_rho <- 1 - pchisq(X2_rho[s], df = 1)
  rej_rho[s] <- p_value_rho < 0.05
  
  X2_rho_b[s] <- (n - 1) * rho_b^2
  p_value_rho_b <- 1 - pchisq(X2_rho_b[s], df = 1)
  rej_rho_b[s] <- p_value_rho_b < 0.05
}

df_alpha <- tibble(
  miera = c("tau_b", "gamma", "rho", "rho_b"),
  zamietnutia = c(
    sum(rej_tau_b, na.rm = TRUE),
    sum(rej_gamma, na.rm = TRUE),
    sum(rej_rho, na.rm = TRUE),
    sum(rej_rho_b, na.rm = TRUE)
  ),
  n_platne = c(
    sum(!is.na(rej_tau_b)),
    sum(!is.na(rej_gamma)),
    sum(!is.na(rej_rho)),
    sum(!is.na(rej_rho_b))
  )
) %>%
  mutate(
    alpha = zamietnutia / n_platne,
    miera = factor(miera, levels = c("tau_b", "gamma", "rho", "rho_b"))
  )

ggplot(df_alpha, aes(x = miera, y = alpha)) +
  #geom_hline(yintercept = 0.05,
  #           linetype = "dashed", linewidth = 0.8, color = "grey35") +
  geom_col(width = 0.62, fill = "#6baed6") +
  geom_text(aes(label = percent(alpha, accuracy = 0.1)),
            vjust = -0.35, size = 5) +
  scale_x_discrete(labels = c(
    "tau_b" = expression(hat(tau)[b]),
    "gamma" = expression(hat(gamma)),
    "rho"   = expression(hat(rho)),
    "rho_b" = expression(hat(rho)[b])
  )) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(df_alpha$alpha) + 0.01),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    title = "Porovnanie sily testov pre jednotlivé miery asociácie",
    x = "Miera asociácie",
    y = "Odhad alphy"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    plot.title = element_text(size = 18, hjust = 0.5)
  )
