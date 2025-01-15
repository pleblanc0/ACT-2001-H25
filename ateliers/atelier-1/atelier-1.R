################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Atelier 1 - 15 janvier 2025 ###############################################
################################################################################

###
### Exemple 1
###

rm(list = ls())

# Nous observons que la variable aléatoire X obéit à un mélange de deux
# lois binomiales de paramètres (n = 100, q = 0.1) et (n = 100, q = 0.6)
# avec des poids respectifs de 0.8 et 0.2.

n <- 100
w <- c(0.8, 0.2)
q <- c(0.1, 0.6)

# Option 1
fx <- w[1] * dbinom(0:n, n, q[1]) + w[2] * dbinom(0:n, n, q[2])

# Option 2
fx <- sapply(0:n, function(x) sum(w * dbinom(x, n, q)))

## 1)
sum(fx)

## 2)
EspX <- sum(0:n * fx)
EspX

################################################################################

###
### Exemple 1.5
###

rm(list = ls())

# Il est important d'être en mesure de déduire le support ainsi que
# les probabilités associées à une variable aléatoire discrète à partir
# de sa fonction génératrice des probabilités.

## 1)
x <- c(0, 1, 2, 3)
fx <- c(0.2, 0.4, 0.3, 0.1)

## 2)

# i)
Fx <- cumsum(fx)
Fx[1.5 + 1]             # Ajouter le + 1 afin d'indexer correctement

# ii)
EspTr <- function(d) sum(x * fx * I(x <= d))
EspTr(1.5)

# iii)
SL <- function(d) sum(pmax(x - d, 0) * fx)
SL(1.5)

# iv)
ExMoy <- function(d) SL(d)/(1 - Fx[d + 1])
ExMoy(1.5)

# v)
FGM <- function(t) sum(exp(t * x) * fx)
FGM(2)


################################################################################

###
### Exemple 2
###

rm(list = ls())

n <- 100
q <- 0.05

x <- 0:n
fx <- dbinom(0:n, n, q)

SL <- function(d) sum(pmax(x - d, 0) * fx)
barplot(sapply(0:10, SL), names.arg = 0:10, col = 'lightcoral')


################################################################################

###
### Exemple 3
###

rm(list = ls())

lam <- 5
xmax <- 100

x <- 0:xmax
fx <- dpois(x, lam)

SL <- function(d) sum(pmax(x - d, 0) * fx)
barplot(sapply(0:10, SL), names.arg = 0:10, col = 'lightcoral')


################################################################################

###
### Exemples 4, 5 et 6
###

rm(list = ls())

r <- c(1, 5, 0.5)
q <- c(1/6, 1/2, 1/11)
xmax <- 100

x <- 0:xmax
f1 <- dnbinom(x, r[1], q[1])
f2 <- dnbinom(x, r[2], q[2])
f3 <- dnbinom(x, r[3], q[3])

SL <- function(d, fx) sum(pmax(x - d, 0) * fx)

par(mfrow = c(3, 1))
barplot(sapply(0:20, SL, fx = f1), names.arg = 0:20, col = 'lightblue')
barplot(sapply(0:20, SL, fx = f2), names.arg = 0:20, col = 'lightcoral')
barplot(sapply(0:20, SL, fx = f3), names.arg = 0:20, col = 'lightgreen')
par(mfrow = c(1, 1))


################################################################################

###
### Exemple 7
###

rm(list = ls())

n <- 100
x <- 0:n
w <- c(0.8, 0.2)
q <- c(0.0125, 0.2)

# Option 1
fx <- w[1] * dbinom(x, n, q[1]) + w[2] * dbinom(x, n, q[2])

# Option 2
fx <- sapply(0:n, function(x) sum(w * dbinom(x, n, q)))

SL <- function(d) sum(pmax(x - d, 0) * fx)

barplot(sapply(0:30, SL), names.arg = 0:30, col = 'lightcoral')


################################################################################

###
### Exemple 8
###

rm(list = ls())

be <- 1/5

SL <- function(d) exp(-be * d)/be
curve(SL, from = 0, to = 30, col = 'lightcoral', lwd = 3)


################################################################################

###
### Exemples 9 et 10
###

rm(list = ls())

al <- c(5, 0.5)
be <- c(1, 1/10)

SL <- function(d, i) al[i]/be[i] * pgamma(d, al[i] + 1, be[i], lower = FALSE) -
                     d * pgamma(d, al[i], be[i], lower = FALSE)

curve(SL(x, 1), from = 0, to = 30, col = 'lightcoral', lwd = 3)
curve(SL(x, 2), from = 0, to = 30, col = 'lightblue', lwd = 3, add = TRUE)


################################################################################

###
### Exemples 11 et 12
###

rm(list = ls())

mu <- c(log(5) - 0.125, log(5) - 0.5)
sig <- c(0.5, 1)

SL <- function(d, i)
{
    exp(mu[i] + sig[i]^2/2) * pnorm(
        (log(d) - mu[i] - sig[i] ^2)/sig[i], lower = FALSE
    ) - d * pnorm((log(d) - mu[i])/sig[i], lower = FALSE)
}

curve(SL(x, 1), from = 0, to = 30, col = 'lightcoral', lwd = 3)
curve(SL(x, 2), from = 0, to = 30, col = 'lightblue', lwd = 3, add = TRUE)


################################################################################

###
### Exemple 13
###

rm(list = ls())

x <- c(6, 9, 10, 13, 16, 17, 20)
fx <- c(0.12, 0.4, 0.3, 0.06, 0.01, 0.07, 0.04)
Fx <- cumsum(fx)

FxInv <- function(u) x[min(which(Fx >= u))]
FxInv(0.95)


################################################################################

###
### Exemple 15
###

rm(list = ls())

u <- c(0.9, 0.99)

## 1)
lam_exp <- 10

FxInvExp <- function(u) -log(1 - u)/lam_exp     # Option 1
FxInvExp <- function(u) qexp(u, lam_exp)        # Option 2
FxInvExp(u)

## 2)
al_par <- 3
lam_par <- 0.6

FxInvPar <- function(u) lam_par * ((1 - u)^(-1/al) - 1)
FxInvPar(u)

## 3)
al_gam <- 4
be <- 0.3

FxInvGam <- function(u) qgamma(u, al_gam, be)
FxInvGam(u)


################################################################################

###
### Exemple 16
###

rm(list = ls())

# Nous observons que la variable aléatoire X obéit à un mélange
# de deux lois exponentielles de paramètres 1/10 et 1/60 avec
# des poids respectifs de 0.8 et 0.2.

be <- c(1/10, 1/60)
w <- c(0.8, 0.2)
u <- 1:99/100

## 1)
fx <- function(x) sum(w * dexp(x, be))
plot(0:200, sapply(0:200, fx), type = 'l', col = 'lightcoral', lwd = 3)

Fx <- function(x) sum(w * pexp(x, be))
plot(0:200, sapply(0:200, Fx), type = 'l', col = 'lightcoral', lwd = 3)

## 2)
FxInv <- function(u) optimize(function(x) abs(Fx(x) - u), c(0, 200))$min
plot(u, sapply(u, FxInv), type = 'l', col = 'lightcoral', lwd = 3)


################################################################################

###
### Exemple 17
###

rm(list = ls())

# On reprend une partie du code de l'exemple 7.

n <- 100
x <- 0:n
w <- c(0.8, 0.2)
q <- c(0.0125, 0.2)
u <- 1:99/100

# Option 1
fx <- w[1] * dbinom(x, n, q[1]) + w[2] * dbinom(x, n, q[2])
Fx <- function(k) w[1] * pbinom(k, n, q[1]) + w[2] * pbinom(k, n, q[2])

# Option 2
fx <- sapply(0:n, function(x) sum(w * dbinom(x, n, q)))
Fx <- function(k) sum(w * pbinom(k, n, q))

# Option 3
Fx <- function(k) sum(fx * I(x <= k))

FxInv <- function(u) x[min(which(Fx >= u))]
plot(u, sapply(u, FxInv), type = 'l', col = 'lightcoral', lwd = 3)


################################################################################

###
### Exemple 18
###

rm(list = ls())

x <- 0:100
fx <- ((5/(5 + x))^2 - (5/(6 + x))^2)/(1 - (5/(106))^2)
sum(fx)     # Vérification

## 1)
EspX <- sum(x * fx)
EspX

## 2)
VarX <- sum((x - EspX)^2 * fx)
sqrt(VarX)

## 3)
Fx <- cumsum(fx)
Fx[30 + 1]          # Ajouter le + 1 afin d'indexer correctement

## 4)

# Je ferai un retour sur cette question lorsque vous aurez vu la TVaR!
SL <- function(d) sum(pmax(x - d, 0) * fx)
VaR <- function(k) x[min(which(Fx >= k))]       # Identique à FxInv!
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)

u <- c(0.6, 0.99)
theta <- (1 - u[1]) * TVaR(u[1]) - (1 - u[2]) * TVaR(u[2])
theta

## 5)
rho <- c(0.001, 0.01, 0.1)
FGM <- function(t) sum(exp(t * x) * fx)

ENTR <- function(rho) log(FGM(rho))/rho
sapply(rho, ENTR)


################################################################################

###
### Exemple 19
###

rm(list = ls())

x <- c(0, 10, 50, 100, 200, 500, 800, 1500, 4000, 10000)
fc <- c(100, 4, 10, 20, 30, 50, 25, 12, 6, 3)

## 1)
fx <- fc/sum(fc)
sum(fx)   # Vérification

## 2)
EspX <- sum(x * fx)
EspX

## 3)
VarX <- sum((x - EspX)^2 * fx)
sqrt(VarX)

## 4)
k <- c(0, 300, 800, 1000, 12000)
Fx <- function(k) sum(fx * I(x <= k))
sapply(k, Fx)

## 5)
u <- c(0.1, 0.5, 0.7, 0.9, 0.99, 0.999)
FxInv <- function(u) x[min(which(cumsum(fx) >= u))]     # Identique à FxInv!
sapply(u, FxInv)

## 6)

# Je ferai un retour sur cette question lorsque vous aurez vu la TVaR!
SL <- function(d) sum(pmax(x - d, 0) * fx)
VaR <- function(k) x[min(which(cumsum(fx) >= k))]
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)

u <- c(0.6, 0.99)
theta <- (1 - u[1]) * TVaR(u[1]) - (1 - u[2]) * TVaR(u[2])
theta

## 7)
rho <- c(0.0001, 0.001, 0.01)
FGM <- function(t) sum(exp(t * x) * fx)

ENTR <- function(rho) log(FGM(rho))/rho
sapply(rho, ENTR)
