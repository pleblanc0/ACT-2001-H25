################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
################################################################################

################################################################################
### Partie I - Lois discrètes ##################################################
################################################################################

# Fonction de masse de probabilité
x <- 0:10
fx <- c(0.1, 0.05, 0.11, 0.15, 0.17, 0.08, 0.1, 0.05, 0.09, 0.05, 0.05)

# Espérance
EspX <- sum(x * fx)
EspX

# Variance
VarX <- sum((x - EspX)^2 * fx)
VarX

# Fonction de répartition
Fx <- cumsum(fx)
Fx

# Espérance limitée
EspLim <- function(d) sum(pmin(x, d) * fx)
EspLim(5)

# Espérance tronquée
EspTr <- function(d) sum(x * fx * I(x > d))
EspTr(5)

# Fonction stop-loss
SL <- function(d) sum(pmax(x - d, 0) * fx)
SL(5)

# Fonction génératrice des moments
FGM <- function(t) sum(exp(t * x) * fx)
FGM(0.1)

# Fonction génératrice des probabilités
FGP <- function(t) sum(t^x * fx)
FGP(0.1)

# Value at Risk
VaR <- function(k) x[min(which(Fx >= k))]
VaR(0.95)

# Tail Value at Risk
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
TVaR(0.95)

# Left Tail Value at Risk
LTVaR <- function(k) (EspX - (1 - k) * TVaR(k))/k
LTVaR(0.95)

# Range Value at Risk
RVaR <- function(k1, k2) ((1 - k1) * TVaR(k1) - (1 - k2) * TVaR(k2))/(k2 - k1)
RVaR(0.6, 0.95)


################################################################################
### Partie II - Généralisations des lois de fréquence ##########################
################################################################################

# Dans les deux cas, nous utiliserons des lois de Poisson avec une moyenne de 5
lam <- 5

# Une valeur maximale de 100 pour le support est suffisante
x <- 0:100

# Mélange de lois de Poisson ---------------------------------------------------
poids <- c(0.8, 0.2)
theta <- c(0.8, 1.8)

f_mix <- sapply(x, function(x) sum(poids * dpois(x, theta * lam)))
sum(f_mix)

# Modification de la masse à 0 -------------------------------------------------
q <- 0.3
p0 <- 1 - q + q * dpois(0, lam)

f_mod <- dzmpois(x, lam, p0)
sum(f_mod)


################################################################################
### Partie III - Convolution ###################################################
################################################################################

# Nous allons réutiliser les densités de la partie II
f1 <- f_mix
f2 <- f_mod

# Convolution pour obtenir la fmp de S
fs <- convolve(f1, rev(f2), type = 'open')
sum(fs)

# Afin d'obtenir l'espérance de S, nous devons obtenir le nouveau support
s <- seq_along(fs) - 1
sum(s * fs)


################################################################################
### Partie IV - Lois continues #################################################
################################################################################

### Loi Gamma ##################################################################

al <- 2
be <- 0.01

# Fonction de répartition
Fx <- function(x) pgamma(x, al, be)
Fx(100)

# Fonction stop-loss
SL <- function(d) al/be * pgamma(d, al + 1, be, lower = FALSE) - d * (1 - Fx(d))
SL(100)

# Value at Risk
VaR <- function(k) qgamma(k, al, be)
VaR(0.99)

# Tail Value at Risk
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
TVaR <- function(k) al/be * pgamma(VaR(k), al + 1, be, lower = FALSE)/(1 - k)
TVaR(0.99)


### Loi Erlang Généralisée #####################################################

be <- c(1/20, 1/50, 1/100)
ci <- sapply(1:3, function(i) prod(be[-i]/(be[-i] - be[i])))

# Fonction de répartition
Fx <- function(x) sum(ci * pexp(x, be))
Fx(100)

# Fonction stop-loss
SL <- function(d) sum(ci * exp(-be * d)/be)
SL(100)

# Value at Risk
VaR <- function(k) optimize(function(x) abs(Fx(x) - k), c(0, 1000))$min
VaR(0.99)

# Tail Value at Risk
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
TVaR(0.99)


### Loi mélange d'exponentielles ###############################################

p <- c(0.2, 0.3, 0.5)
be <- c(1/20, 1/50, 1/100)

# Fonction de répartition
Fx <- function(x) sum(p * pexp(x, be))
Fx(100)

# Fonction stop-loss
SL <- function(d) sum(p/be * exp(-be * d))
SL(100)

# Value at Risk
VaR <- function(k) optimize(function(x) abs(Fx(x) - k), c(0, 1000))$min
VaR(0.99)

# Tail Value at Risk
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
TVaR(0.99)


### Loi Poisson Composée avec sinistres Gamma ##################################

al <- 2
be <- 0.01
lam <- 2
k0 <- 1:1000

# Fonction de masse de probabilité
p0 <- dpois(0, lam)
pk <- dpois(k0, lam)

# Fonction de répartition
Fx <- function(x) p0 + sum(pk * pgamma(x, k0 * al, be))
Fx(100)

# Value at Risk
VaR <- function(k) optimize(function(x) abs(Fx(x) - k), c(0, 2000))$min
VaR(0.99)

# Tail Value at Risk
TVaR <- function(k) sum(pk * k0 * al/be * pgamma(VaR(k), k0 * al + 1,
                                          be, lower = FALSE))/(1 - k)
TVaR(0.99)


### Somme de lois Gamma avec beta différents ###################################

al <- c(0.5, 1, 1.5)
be <- c(0.1, 0.2, 0.3)

# On crée les densités des deux binomiales négatives nécessaires
fj1 <- dnbinom(0:100, al[1], be[1]/be[3])
fj2 <- dnbinom(0:100, al[2], be[2]/be[3])

# On obtient les poids du mélange avec la convolution des deux densités
pk <- convolve(fj1, rev(fj2), type = 'open')
sum(pk)

# Fonction de répartition en utilisant les poids et les paramètres ajustés
Fs <- function(x) sum(pk * pgamma(x, sum(al) + 0:200, max(be)))
Fs(10)
