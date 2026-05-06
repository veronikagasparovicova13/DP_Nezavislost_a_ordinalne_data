# Funkcie: had1_taub
#          funkcia_concordant_discordant

library(tibble)
library(ggplot2)
library(scales)

collapse_3x4_to_2x3_r12_c12 <- function(tbl) {
  stopifnot(all(dim(tbl) == c(3,4)))
  out <- matrix(0, 2, 3)
  out[1,1] <- sum(tbl[1:2, 1:2])
  out[1,2] <- sum(tbl[1:2, 3])
  out[1,3] <- sum(tbl[1:2, 4])
  out[2,1] <- sum(tbl[3, 1:2])
  out[2,2] <- tbl[3,3]
  out[2,3] <- tbl[3,4]
  out
}
ppA <- collapse_3x4_to_2x3_r12_c12(ppA)

# collapse_3x4_to_2x3_r23_c12 <- function(tbl) {
#   stopifnot(all(dim(tbl) == c(3,4)))
#   out <- matrix(0, 2, 3)
#   # R1 ostáva, C1+C2
#   out[1,1] <- sum(tbl[1, 1:2])
#   out[1,2] <- tbl[1,3]
#   out[1,3] <- tbl[1,4]
#   # R2+R3
#   out[2,1] <- sum(tbl[2:3, 1:2])
#   out[2,2] <- sum(tbl[2:3, 3])
#   out[2,3] <- sum(tbl[2:3, 4])
#   out
# }

# Parametre
n_sim <- 10000
n_jedinci <- 300
r <- 3
c <- 4

# ppA <- rep(1/(r*c), r*c) # rovnomerna pp

row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
ppA <- outer(row_sums, col_sums) # nerovnomerna pp. (platnost H0)

# ppA <- matrix(
#   c(
#     0.30, 0.20, 0.05, 0.05,   # sum = 0.6
#     0.05, 0.05, 0.10, 0.10,   # sum = 0.3
#     0.01, 0.03, 0.03, 0.03    # sum = 0.1
#   ),
#   nrow = 3,
#   byrow = TRUE
# ) # H0 neplati

# ppA <- matrix(c(
#   4,  3,  1,  6,
#   4, 20,  6,  1,
#   7,  3, 16,  3
# ), nrow = 3, byrow = TRUE) / 74 # po zluceni na mensiu tab je sila testu vacsia

z_3x4 <- numeric(n_sim) # tab 3x4
rej_3x4 <- logical(n_sim)

z_2x3 <- numeric(n_sim) # tab 2x3
rej_2x3 <- logical(n_sim)

for (s in 1:n_sim) {
  # Generovanie tabulky
  vector <- rmultinom(1, n_jedinci, as.vector(t(ppA)))
  tbl1 <- matrix(vector, r, c, byrow = TRUE)
  tbl2 <- collapse_3x4_to_2x3_r12_c12(tbl1)
  
  n <- sum(tbl1)
  total_pairs <- choose(n, 2)
  p1 <- tbl1/n
  p2 <- tbl2/n
  
  # C a D
  cd1 <- funkcia_concordant_discordant(tbl1)
  C1 <- cd1$C
  D1 <- cd1$D
  Pi_c1 <- C1 / total_pairs
  Pi_d1 <- D1 / total_pairs
  
  cd2 <- funkcia_concordant_discordant(tbl2)
  C2 <- cd2$C
  D2 <- cd2$D
  Pi_c2 <- C2 / total_pairs
  Pi_d2 <- D2 / total_pairs
  
  row_sucty1 <- rowSums(p1)
  col_sucty1 <- colSums(p1)
  
  row_sucty2 <- rowSums(p2)
  col_sucty2 <- colSums(p2)
  
  # tau b
  taub_3x4 <- (Pi_c1 - Pi_d1) / sqrt((1-sum(row_sucty1^2)) * (1-sum(col_sucty1^2)))
  taub_2x3 <- (Pi_c2 - Pi_d2) / sqrt((1-sum(row_sucty2^2)) * (1-sum(col_sucty2^2)))
  
  # Sigma
  sigma_3x4 <- had1_taub(tbl1)
  sigma_2x3 <- had1_taub(tbl2)
  
  # Z-test pre taub = 0
  z_3x4[s] <- taub_3x4 / sigma_3x4
  p_value_3x4 <- 2 * (1 - pnorm(abs(z_3x4[s])))
  rej_3x4[s] <- p_value_3x4 < 0.05
  
  z_2x3[s] <- taub_2x3 / sigma_2x3
  p_value_2x3 <- 2 * (1 - pnorm(abs(z_2x3[s])))
  rej_2x3[s] <- p_value_2x3 < 0.05
}

# Vykreslenie
k_3x4 <- sum(rej_3x4, na.rm = TRUE) # pocet zamietnuti pre tabulku: 1 mult. vyber
n_3x4 <- sum(!is.na(rej_3x4)) # pocet platnych zamietnuti pre tabulku: 1 mult. vyber (bez NA)
k_2x3 <- sum(rej_2x3, na.rm = TRUE)
n_2x3 <- sum(!is.na(rej_2x3))

# odhad empiricky pre alphu = podiel zamietnuti spomedzi platnych zamietnuti
alpha_3x4 <- k_3x4 / n_3x4
alpha_2x3 <- k_2x3 / n_2x3

df_alpha <- tibble(
  scenár = c("3x4", "2x3"),
  alpha  = c(alpha_3x4, alpha_2x3)
)

ggplot(df_alpha, aes(x = scenár, y = alpha)) +
  geom_col(width = 0.6, alpha = 0.7, fill = c("#7f7f7f", "#4C72B0")) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Pravdepodobnosť chyby I. druhu pre obe tabuľky pre TauB",
       x = "Veľkosť tabuľky", y = "Odhad alphy") +
  theme_minimal(base_size = 13) + geom_text(aes(label = percent(alpha, accuracy = 0.1)), 
                                            vjust = -0.5, size = 4) +
  geom_hline(yintercept = 0.05, linetype = "dashed")

