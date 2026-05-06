# Potrebne funkcie: funkcia_concordant_discordant
#                   gamma_1hadik
#                   gamma_3hadici

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tibble)

# Parametre
n_sim <- 10000
n_jedinci <- 300
r <- 3
c <- 4

### Gamma: generovanie jedneho hadika

# ppA <- rep(1/(r*c), r*c) # rovnomerne pp. rozdelenie

row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
ppA <- t(outer(row_sums, col_sums)) # nerovnomerna pp. (platnost H0)
ppA <- as.vector(ppA)

z_A <- numeric(n_sim) # nespravne pouzitie sigmy
rej_A <- logical(n_sim)

z_AA <- numeric(n_sim) # spravne pouzitie sigmy
rej_AA <- logical(n_sim)

for (s in 1:n_sim) {
  hadik <- rmultinom(1, n_jedinci, ppA)
  tbl <- matrix(hadik, r, c, byrow = TRUE)
  
  # Vypocet gamma a sigma_gamma
  sigma_nespravna <- gamma_3hadici(tbl)
  sigma_spravna <- gamma_1hadik(tbl)
  
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  gamma <- (C - D) / (C + D)
  
  # Z-test pre gamma = 0
  z_A[s] <- gamma / sigma_nespravna
  p_value_A <- 2 * (1 - pnorm(abs(z_A[s])))
  rej_A[s] <- p_value_A < 0.05
  
  z_AA[s] <- gamma / sigma_spravna
  p_value_AA <- 2 * (1 - pnorm(abs(z_AA[s])))
  rej_AA[s] <- p_value_AA < 0.05
}

### Gamma: generovanie troch hadikov
# ppA <- matrix(1/12, nrow = 3, ncol = 4)
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

riadok_pocet_jedincov <- n_jedinci * rowSums(matrix(ppA, nrow = r, ncol = c, byrow = TRUE))

for (s in 1:n_sim) {
  tbl <- matrix(0L, r, c)
  for (i in 1:r) {
    tbl[i, ] <- as.vector(rmultinom(1, riadok_pocet_jedincov[i], ppB[i, ]))
  }
  
  # C a D
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  # gamma
  gamma <- (C-D)/(C+D)
  
  # sigma 
  sigma_nespravna <- gamma_1hadik(tbl)
  sigma_spravna <- gamma_3hadici(tbl)
  
  # Z-test pre gamma = 0
  z_B[s] <- gamma / sigma_nespravna
  p_value_B <- 2 * (1 - pnorm(abs(z_B[s])))
  rej_B[s] <- p_value_B < 0.05
  
  z_BB[s] <- gamma / sigma_spravna
  p_value_BB <- 2 * (1 - pnorm(abs(z_BB[s])))
  rej_BB[s] <- p_value_BB < 0.05
}

################################################################################

# ======================== 5% chyba porovnanie =================================
# pravdepodobnost zamietnutia H0, ked H0 plati

k_A <- sum(rej_A, na.rm = TRUE) # pocet zamietnuti pre tabulku: 1 mult. vyber
n_A <- sum(!is.na(rej_A)) # pocet platnych zamietnuti pre tabulku: 1 mult. vyber (bez NA)
k_AA <- sum(rej_AA, na.rm = TRUE)
n_AA <- sum(!is.na(rej_AA))

k_B <- sum(rej_B, na.rm = TRUE)
n_B <- sum(!is.na(rej_B))
k_BB <- sum(rej_BB, na.rm = TRUE)
n_BB <- sum(!is.na(rej_BB))

# odhad empiricky pre alphu = podiel zamietnuti spomedzi platnych zamietnuti
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
  geom_hline(yintercept = 0.05,
             linetype = "dashed", linewidth = 0.7, color = "grey35") +
  
  geom_col(width = 0.58) +
  
  geom_text(aes(label = percent(alpha, accuracy = 0.1)),
            vjust = -0.35, size = 4.2) +
  
  facet_wrap(~variancia) +
  
  scale_fill_manual(values = c("správna" = "#39c56b",
                               "nesprávna" = "#e74c3c")) +
  
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, max(df_alpha$alpha) + 0.02),
    expand = expansion(mult = c(0, 0.03))
  ) +
  
  labs(
    title = "Pravdepodobnosť chyby I. druhu pre Gamma",
    x = "Generovanie dát",
    y = "Odhad alphy",
    fill = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold")
  )

