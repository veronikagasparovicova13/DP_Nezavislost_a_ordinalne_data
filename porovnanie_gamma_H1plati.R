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
n_jedinci <- 30
r <- 3
c <- 4

### Gamma: generovanie jedneho hadika
# ppA <- matrix(c(
#   0.04, 0.03, 0.02, 0.01,   # sum = 0.10
#   0.08, 0.10, 0.07, 0.05,   # sum = 0.30
#   0.14, 0.17, 0.16, 0.13    # sum = 0.60
# ), nrow = 3, byrow = TRUE) # gamma = 0.1752988

ppA <- matrix(c(
  0.03, 0.35, 0.20, 0.02,   # sum = 0.60
  0.02, 0.06, 0.13, 0.09,   # sum = 0.30
  0.00, 0.00, 0.01, 0.09    # sum = 0.10
), nrow = 3, byrow = TRUE) # tau_b = 0.5

z_A <- numeric(n_sim) # nespravne pouzitie sigmy
rej_A <- logical(n_sim)

z_AA <- numeric(n_sim) # spravne pouzitie sigmy
rej_AA <- logical(n_sim)

for (s in 1:n_sim) {
  hadik <- rmultinom(1, n_jedinci, as.vector(t(ppA)))
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
ppB <- sweep(ppA, 1, rowSums(ppA), "/")

z_B <- numeric(n_sim) # nespravne pouzitie sigmy
rej_B <- logical(n_sim)

z_BB <- numeric(n_sim) # spravne pouzitie sigmy
rej_BB <- logical(n_sim)

riadok_pocet_jedincov <- n_jedinci * rowSums(ppA)

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

k_A <- sum(rej_A, na.rm = TRUE) # pocet zamietnuti pre tabulku: 1 had
n_A <- sum(!is.na(rej_A)) # pocet platnych zamietnuti pre tabulku: 1 had (bez NA)
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
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +
  
  geom_text(aes(label = percent(alpha, accuracy = 0.1)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 4) +
  
  facet_wrap(~variancia) +
  
  scale_fill_manual(values = c("správna" = "#2ecc71",
                               "nesprávna" = "#e74c3c")) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(df_alpha$alpha) + 0.05)) +
  
  labs(title = "Sila testu pre Gamma",
       x = "Generovanie dát",
       y = "Odhad sily testu",
       fill = "") +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold"),  # to je presne ten efekt čo chceš
    panel.grid.minor = element_blank()
  )

