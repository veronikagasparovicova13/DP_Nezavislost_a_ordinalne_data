library(ggplot2)
library(tibble)
library(scales)

# Potrebne funkcie: funkcia_concordant_discordant
#                   spravne_phi_taub (s n)
#                   nespravne_phi_taub (s n)

# Prejav chyby v simulaciach ###################################################

n_sim <- 10000
n_jedinci <- 300
r <- 3
c <- 4
n_kategorie <- r * c

# ppA <- rep(1/n_kategorie, n_kategorie)

# row_sums <- c(0.6, 0.3, 0.1) # riadkove sucty
# col_sums <- c(0.16, 0.28, 0.3, 0.26) # stlpcove sucty
# ppA <- t(outer(row_sums, col_sums)) # nerovnomerna pp. (platnost H0)
# ppA <- as.vector(ppA)

z_A <- numeric(n_sim) # wrong
rej_A <- logical(n_sim)

z_AA <- numeric(n_sim) # ok
rej_AA <- logical(n_sim)

for (s in 1:n_sim) {
  # Generovanie tabulky
  vector <- rmultinom(1, n_jedinci, ppA)
  tbl <- matrix(vector, r, c, byrow = TRUE)
  
  n <- sum(tbl)
  p <- tbl/n
  total_pairs <- choose(n, 2)
  
  # C a D
  cd <- funkcia_concordant_discordant(tbl)
  C <- cd$C
  D <- cd$D
  
  row_sucty <- rowSums(tbl)
  col_sucty <- colSums(tbl)
  TY <- sum(choose(col_sucty, 2))
  TX <- sum(choose(row_sucty, 2))
  row_marg <- rowSums(p)
  col_marg <- colSums(p)
  delta <- sqrt((1 - sum(row_marg^2)) * (1 - sum(col_marg^2)))
  
  # tau b
  taub <- (C - D) / sqrt((total_pairs - TX) * (total_pairs - TY))
  
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
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  labs(title = "Vplyv chyby vo phi na pravdepodobnosť chyby I. druhu",
       x = NULL, y = "Odhad α") +
  theme_minimal(base_size = 13) + geom_text(aes(label = percent(alpha, accuracy = 0.01)), 
                                            vjust = -0.5, size = 4)

# Funkcie, kde je Agrestiho chyba
DescTools::KendallTauB # chyba vo phi
ConDisPairs(tbl)
