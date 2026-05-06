library(tibble)
library(ggplot2)
library(scales)

# Potrebne funkcie: funkcia_concordant_discordant
#                   had1_taub
#                   hady3_taub

# Parametre ####################################################################
n_sim <- 10000
n_jedinci <- 300
r <- 3
c <- 4
n_kategorie <- r * c

# Simulácia: 1 mult. vyber #####################################################
# ppA <- rep(1/n_kategorie, n_kategorie) # rovnomerna pp

row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

z_A <- numeric(n_sim) # nespravne pouzitie sigmy
rej_A <- logical(n_sim)

z_AA <- numeric(n_sim) # spravne pouzitie sigmy
rej_AA <- logical(n_sim)

for (s in 1:n_sim) {
  # Generovanie tabulky
  vector <- rmultinom(1, n_jedinci, as.vector(t(ppA)))
  tbl <- matrix(vector, r, c, byrow = TRUE)
  
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
  
  # tau b
  taub <- (Pi_c - Pi_d) / sqrt(
    (1 - sum(row_sucty^2)) * (1 - sum(col_sucty^2))
  )
  
  # Sigma
  sigma_spravna <- had1_taub(tbl)
  sigma_nespravna <- hady3_taub(tbl)
  
  # Z-test pre taub = 0
  z_A[s] <- taub / sigma_nespravna
  p_value_A <- 2 * (1 - pnorm(abs(z_A[s])))
  rej_A[s] <- p_value_A < 0.05
  
  z_AA[s] <- taub / sigma_spravna
  p_value_AA <- 2 * (1 - pnorm(abs(z_AA[s])))
  rej_AA[s] <- p_value_AA < 0.05
}

# Simulácia: 3 hady ############################################################

# P <- matrix(ppA, nrow = r, ncol = c, byrow = TRUE)
# ppB <- P / rowSums(P)

row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
pp <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)
ppB <- pp / row_sums

z_B <- numeric(n_sim) # nespravne pouzitie sigmy
rej_B <- logical(n_sim)

z_BB <- numeric(n_sim) # spravne pouzitie sigmy
rej_BB <- logical(n_sim)

riadok_pocet_jedincov <- n_jedinci * rowSums(matrix(ppA, nrow = 3, ncol = 4))

for (s in 1:n_sim) {
  # Generovanie tabulky
  tbl <- matrix(0L, r, c)
  for (i in 1:r) {
    tbl[i, ] <- as.vector(rmultinom(1, riadok_pocet_jedincov[i], ppB[i, ]))
  }
  
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
  
  # tau b
  taub <- (Pi_c - Pi_d) / sqrt(
    (1 - sum(row_sucty^2)) * (1 - sum(col_sucty^2))
  )
  
  # Sigma
  sigma_spravna <- hady3_taub(tbl)
  sigma_nespravna <- had1_taub(tbl)
  
  # Z-test pre taub = 0
  z_B[s] <- taub / sigma_nespravna
  p_value_B <- 2 * (1 - pnorm(abs(z_B[s])))
  rej_B[s] <- p_value_B < 0.05
  
  z_BB[s] <- taub / sigma_spravna
  p_value_BB <- 2 * (1 - pnorm(abs(z_BB[s])))
  rej_BB[s] <- p_value_BB < 0.05
}

# Vykreslenie chyby ############################################################
k_A <- sum(rej_A, na.rm = TRUE) # pocet zamietnuti pre tabulku: 1 mult. vyber
n_A <- sum(!is.na(rej_A))
k_AA <- sum(rej_AA, na.rm = TRUE)
n_AA <- sum(!is.na(rej_AA))

k_B <- sum(rej_B, na.rm = TRUE)
n_B <- sum(!is.na(rej_B))
k_BB <- sum(rej_BB, na.rm = TRUE)
n_BB <- sum(!is.na(rej_BB))

# odhad empiricky pre alphu = podiel zamietnuti / platne zamietnutia
alpha_A <- k_A / n_A
alpha_AA <- k_AA / n_AA
alpha_B <- k_B / n_B
alpha_BB <- k_BB / n_BB

df_alpha <- tibble(
  generovanie = c("1 mult.", "1 mult.", "3 mult.", "3 mult."),
  variancia   = c("SE(1M)", "SE(3M)", "SE(1M)", "SE(3M)"),
  typ         = c("správna", "nesprávna", "nesprávna", "správna"),
  alpha       = c(alpha_AA, alpha_A, alpha_BB, alpha_B)
)

ggplot(df_alpha, aes(x = generovanie, y = alpha, fill = typ)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
  
  geom_text(aes(label = percent(alpha, accuracy = 0.1)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  
  facet_wrap(~variancia) +
  
  scale_fill_manual(values = c("správna" = "#2ecc71",
                               "nesprávna" = "#e74c3c")) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(df_alpha$alpha) + 0.05)) +
  
  labs(title = "Pravdepodobnosť chyby I. druhu pre Tau-b",
       x = "Generovanie dát",
       y = "Odhad alphy",
       fill = "") +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
