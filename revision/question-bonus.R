################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Question Bonus - 9 février 2025 ###########################################
################################################################################

# CONTEXTE: Soit X une variable aléatoire discrète définit sur {0, ..., 100}
# telle que sa fmp est contenue dans le vecteur fx ci-dessous.
x <- 0:100
fx <- ((5/(5 + x))^2 - (5/(6 + x))^2)/(1 - (5/106)^2)

# QUESTION: Soit la variable aléatoire Z = max(min(X, 75) - 20, 0). Calculez
# l'espérance et la variance de Z. Calculez également la VaR et la TVaR de Z.

# SOLUTION: On observe que Z est définie sur le support{0, ..., 55} avec une
# fmp que nous pouvons obtenir avec les deux méthodes qui suivent.
z <- 0:55
fz <- tapply(fx, pmax(pmin(x, 75) - 20, 0), sum)        # Directement
fz <- c(sum(fx[1:21]), fx[22:75], sum(fx[76:101]))      # Alternative

# Espérance et variance
EspZ <- sum(z * fz)
VarZ <- sum((z - EspZ)^2 * fz) 

# Fonction stop-loss
SL <- function(d) sum(pmax(z - d, 0) * fz)

# Value-at-Risk
VaR <- function(k) z[min(which(cumsum(fz) >= k))]
VaR(0.99)

# Tail Value-at-Risk
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
TVaR(0.99)
