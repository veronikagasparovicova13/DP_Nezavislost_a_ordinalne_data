# Funkcie: spravne_phi_taub (bez n)
#          had1_taub
#          hady3_taub
#          funkcia_concordant_discordant

library(ggplot2)
library(scales)

r <- 3
c <- 4

# ppA <- rep(1/12, 12) # rovnomerna pp
# ppA <- matrix(rep(1/12, 12), nrow = 3, ncol = 4, byrow = TRUE)

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

# ppA <- matrix(c(
#   0.04, 0.028, 0.018, 0.014,   # sum = 0.10
#   0.075, 0.095, 0.075, 0.055,  # sum = 0.30
#   0.12, 0.147, 0.177, 0.156    # sum = 0.60
# ), nrow = 3, byrow = TRUE) # tau b = 0.1430873; popul. relat. chyba = 0.7

# ppA <- matrix(c(
#   0.10, 0.15, 0.20, 0.15,
#   0.10, 0.10, 0.05, 0.05,
#   0.05, 0.03, 0.01, 0.01
# ), nrow = 3, byrow = TRUE) # tau_b = -0.2534267, popul. relat. chyba = 1.07%

# ppA <- matrix(c(
#   1, 10, 1, 1,
#   1,  2, 4, 3,
#   1,  1, 1, 4
# ), nrow = 3, byrow = TRUE) / 30 # taub = 0.3988913; 
# 
ppA <- matrix(c(
  0.03, 0.35, 0.20, 0.02,   # sum = 0.60
  0.02, 0.06, 0.13, 0.09,   # sum = 0.30
  0.00, 0.00, 0.01, 0.09    # sum = 0.10
), nrow = 3, byrow = TRUE) # tau_b = 0.5
rowSums(ppA)

# vypocet sigma_1had (populacna) ###############################################
phi <- spravne_phi_taub(ppA)
row_marg <- rowSums(ppA)
col_marg <- colSums(ppA)
delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
sigma2_1had <- ( sum(ppA * phi^2) - (sum(ppA * phi))^2)/ delta^4
sigma_1had <- sqrt(sigma2_1had)
sigma_1had
# vypocet sigma_3hady (populacna) ##############################################

row_marg <- rowSums(ppA)
col_marg <- colSums(ppA)

# Pravdepodobnosti konkordantnych a diskordantnych parov
outC <- 0
for (i in 1:(r - 1)) {
  for (k in (i + 1):r) {
    for (j in 1:(c - 1)) {
      for (l in (j + 1):c) {
        outC <- outC + ppA[i, j] * ppA[k, l]
      }
    }
  }
}
Pi_c <- 2 * outC

outD <- 0
for (i in 1:(r - 1)) {
  for (k in (i + 1):r) {
    for (j in 2:c) {
      for (l in 1:(j - 1)) {
        outD <- outD + ppA[i, j] * ppA[k, l]
      }
    }
  }
}
Pi_d <- 2 * outD

nu <- Pi_c - Pi_d
delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
omega <- row_marg
B <- 1 - sum(col_marg^2)

# Kvadrantove sucty
R_c <- matrix(0, r, c)
R_d <- matrix(0, r, c)

for (i in 1:r) {
  for (j in 1:c) {
    part1_c <- if (i > 1 && j > 1) sum(ppA[1:(i - 1), 1:(j - 1)]) else 0
    part2_c <- if (i < r && j < c) sum(ppA[(i + 1):r, (j + 1):c]) else 0
    R_c[i, j] <- part1_c + part2_c
    
    part1_d <- if (i > 1 && j < c) sum(ppA[1:(i - 1), (j + 1):c]) else 0
    part2_d <- if (i < r && j > 1) sum(ppA[(i + 1):r, 1:(j - 1)]) else 0
    R_d[i, j] <- part1_d + part2_d
  }
}

# Podmienene pravdepodobnosti theta_{j|i}
theta <- matrix(0, r, c)
T <- matrix(0, r, c)

for (i in 1:r) {
  for (j in 1:c) {
    theta[i, j] <- ppA[i, j] / omega[i]
    T[i, j] <- 2 * (R_c[i, j] - R_d[i, j]) + nu * col_marg[j] / B
  }
}

inner1 <- rowSums(theta * (T^2))
inner2 <- rowSums(theta * T)

sigma2_3hady <- (1 / delta^2) * sum(omega * (inner1 - inner2^2))
sigma_3hady <- sqrt(sigma2_3hady)
sigma_3hady
pop_rel_err <- (sigma_3hady - sigma_1had)/sigma_1had
pop_rel_err

pop_rel_err2 <- (sigma_1had - sigma_3hady)/sigma_3hady
pop_rel_err2

# Relativna chyba (sigma_3 - sigma_1)/sigma_1 ##################################
n_sim <- 10000
n_jedinci <- 300

# ulozime si sigmy + rel chyby
sigma_1_vec    <- numeric(n_sim) # spravna
sigma_3_vec <- numeric(n_sim) # nespravna
rel_err_vec     <- numeric(n_sim)

for (s in 1:n_sim) {
  vector <- rmultinom(1, n_jedinci, prob = as.vector(t(ppA)))  
  tbl <- matrix(vector, r, c, byrow = TRUE)
  
  sigma_ok <- had1_taub(tbl)
  sigma_wrong <- hady3_taub(tbl)
  
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
    title = expression("Relatívna chyba, " ~ tau[b] == 0.5),
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
  #+
  # scale_x_continuous(
  #   limits = c(-0.025, 0),
  #   labels = percent_format(accuracy = 0.1)
  # )

# Relativna chyba (sigma_1 - sigma_3)/sigma_3 ##################################
n_sim <- 10000
n_jedinci <- 300

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
  
  sigma_ok <- hady3_taub(tbl)
  sigma_wrong <- had1_taub(tbl)
  
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

# rozsah dát
x_range <- range(df_rel$rel_err)
x_offset <- 0.02 * diff(x_range)

# automatické rozhodnutie strany
if (pop_rel_err2 > mean(x_range)) {
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
    title = expression("Relatívna chyba, " ~ tau[b] == 0.5),
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

