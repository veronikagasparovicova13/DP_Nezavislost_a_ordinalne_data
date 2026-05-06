# Funkcie: gamma_1hadik_pop
#          gamma_3hadici_pop
#          gamma_1hadik
#          gamma_3hadici
#          funkcia_concordant_discordant

library(ggplot2)
library(scales)

r <- 3
c <- 4

# Pp tabulky ###################################################################
# ppA <- rep(1/12, 12) # rovnomerna pp
# ppA <- matrix(rep(1/12, 12), nrow = 3, ncol = 4, byrow = TRUE)

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

# ppA <- matrix(c(
#   0.04, 0.028, 0.018, 0.014,   # sum = 0.10
#   0.075, 0.095, 0.075, 0.055,  # sum = 0.30
#   0.12, 0.147, 0.177, 0.156    # sum = 0.60
# ), nrow = 3, byrow = TRUE) # tau b = 0.1430873 gamma = 0.2239137

# ppA <- matrix(c(
#   3, 2, 4, 1,
#   4, 2, 1, 1,
#   0, 4, 2, 6
# ), nrow = 3, byrow = TRUE) / 30 # gamma = 0.4008439 ; relat. chyba = 8%

ppA <- matrix(c(
  0.03, 0.35, 0.20, 0.02,   # sum = 0.60
  0.02, 0.06, 0.13, 0.09,   # sum = 0.30
  0.00, 0.00, 0.01, 0.09    # sum = 0.10
), nrow = 3, byrow = TRUE) # tau_b = 0.5

# Poppulacne sigmy #############################################################
gamma_1had <- gamma_1hadik_pop(ppA)
gamma_1had

gamma_3hady <- gamma_3hadici_pop(ppA)
gamma_3hady

pop_rel_err <- (gamma_3hady - gamma_1had)/gamma_1had
pop_rel_err

pop_rel_err2 <- (gamma_1had - gamma_3hady)/gamma_3hady
pop_rel_err2

# Vyberova relativna chyba (sigma_3 - sigma_1)/sigma_1 #########################
n_sim <- 10000
n_jedinci <- 300

# ulozime si sigmy + rel chyby
sigma_1_vec    <- numeric(n_sim) # spravna
sigma_3_vec <- numeric(n_sim) # nespravna
rel_err_vec     <- numeric(n_sim)

for (s in 1:n_sim) {
  vector <- rmultinom(1, n_jedinci, prob = as.vector(t(ppA)))  
  tbl <- matrix(vector, r, c, byrow = TRUE)
  
  sigma_ok <- gamma_1hadik(tbl)
  sigma_wrong <- gamma_3hadici(tbl)
  
  sigma_1_vec[s] <- sigma_ok
  sigma_3_vec[s] <- sigma_wrong
  
  # relatívna chyba (ak sigma_ok = 0 alebo NA -> NA)
  rel_err_vec[s] <- ifelse(is.na(sigma_ok) || sigma_ok == 0, NA_real_,
                           (sigma_wrong - sigma_ok) / sigma_ok)
}

df_rel <- tibble(
  sigma_ok = sigma_1_vec,
  sigma_wrong = sigma_3_vec,
  rel_err = rel_err_vec,
  abs_rel_err = abs(rel_err_vec)
)

# dynamický posun (2 % rozsahu dát)
x_offset <- 0.02 * diff(range(df_rel$rel_err))

ggplot(df_rel, aes(x = rel_err)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 0.01)) +
  labs(
    title = expression("Relatívna chyba, " ~ gamma == 0.73),
    x = "(σ_3 − σ_1) / σ_1",
    y = "Počet simulácií"
  ) +
  theme_minimal(base_size = 13) +
  geom_vline(
    xintercept = pop_rel_err,
    color = "red",
    linetype = "dashed",
    size = 1.2
  ) +
  annotate(
    "text",
    x = pop_rel_err - x_offset,
    y = Inf,
    label = paste("Popul. chyba:", round(pop_rel_err * 100, 2), "%"),
    vjust = 2,
    hjust = 1,
    color = "red",
    fontface = "bold"
  ) 

# Relativna chyba (sigma_1 - sigma_3)/sigma_3 ##################################

# ulozime si sigmy + rel chyby
sigma_1_vec    <- numeric(n_sim) # nespravna
sigma_3_vec <- numeric(n_sim) # spravna
rel_err_vec     <- numeric(n_sim)

riadok_pocet_jedincov <- n_jedinci * rowSums(matrix(ppA, nrow = 3, ncol = 4))
ppB <- sweep(ppA, 1, rowSums(ppA), "/")

for (s in 1:n_sim) {
  tbl <- matrix(0L, r, c)
  for (i in 1:r) {
    tbl[i, ] <- as.vector(rmultinom(1, riadok_pocet_jedincov[i], ppB[i, ]))
  }
  
  sigma_ok <- gamma_3hadici(tbl)
  sigma_wrong <- gamma_1hadik(tbl)
  
  sigma_3_vec[s] <- sigma_ok
  sigma_1_vec[s] <- sigma_wrong
  
  # relatívna chyba (ak sigma_ok = 0 alebo NA -> NA)
  rel_err_vec[s] <- ifelse(is.na(sigma_ok) || sigma_ok == 0, NA_real_,
                           (sigma_wrong - sigma_ok) / sigma_ok)
}

df_rel <- tibble(
  sigma_ok = sigma_3_vec,
  sigma_wrong = sigma_1_vec,
  rel_err = rel_err_vec,
  abs_rel_err = abs(rel_err_vec)
)

rel_finite <- df_rel$rel_err[is.finite(df_rel$rel_err)]

x_range <- range(rel_finite, na.rm = TRUE)
x_offset <- 0.02 * diff(x_range)

if (is.finite(pop_rel_err2) && pop_rel_err2 > mean(x_range)) {
  x_text <- pop_rel_err2 - x_offset
  hjust_val <- 1
} else {
  x_text <- pop_rel_err2 + x_offset
  hjust_val <- 0
}

ggplot(df_rel, aes(x = rel_err)) +
  geom_histogram(bins = 60, alpha = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 0.01)) +
  labs(
    title = expression("Relatívna chyba, " ~ gamma == 0.73),
    x = "(σ_1 − σ_3) / σ_3",
    y = "Počet simulácií"
  ) +
  theme_minimal(base_size = 13) +
  geom_vline(
    xintercept = pop_rel_err2,
    color = "red",
    linetype = "dashed",
    size = 1.2
  ) +
  annotate(
    "text",
    x = x_text,
    y = Inf,
    label = paste("Popul. chyba:", round(pop_rel_err2 * 100, 2), "%"),
    vjust = 2,
    hjust = hjust_val,
    color = "red",
    fontface = "bold"
  )

# Populacna gamma ##############################################################
ppA <- matrix(ppA, nrow = 3, ncol = 4)
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
delta <- Pi_c + Pi_d
gamma <- nu/delta
gamma
