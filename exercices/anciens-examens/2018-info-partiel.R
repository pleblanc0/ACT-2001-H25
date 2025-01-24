################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Examen partiel informatique H18 ###########################################
################################################################################

###
### Question 1
###

rm(list = ls())

al <- 0.1
be <- 0.001
n <- c(4, 100)

# i)
alW <- n * al
beW <- n * be

# ii)
EspW <- alW/beW
EspW

VarW <- alW/beW^2
sqrt(VarW)

## iii)
VaRW <- function(k) qgamma(k, alW, beW)
VaRW(0.9999)

# iv)
TVaRW <- function(k) alW/beW * pgamma(VaRW(k), alW + 1, beW, low = F)/(1 - k)
TVaRW(0.9999)


###
### Question 2
###

rm(list = ls())

al <- c(0.5, 1.5, 2.5)
be <- c(0.05, 0.15, 0.25)

m <- 1e5
set.seed(2018)

U <- matrix(runif(3 * m), m, byrow = TRUE)
X1 <- qgamma(U[, 1], al[1], be[1])
X2 <- qgamma(U[, 2], al[2], be[2])
X3 <- qgamma(U[, 3], al[3], be[3])

## b)
k <- 0.99
VaR <- function(x) sort(x)[k * m]

TVaR <- function(x) mean(x[x > VaR(x)])
sapply(list(X1, X2, X3), TVaR)

## c)
S <- X1 + X2 + X3
S[c(3, 4)]

## d)
TVaR(S)


###
### Question 3
###

rm(list = ls())

lam <- 1
al <- 1.5
be <- 0.15
k0 <- 1:1000

## a)
EspX <- lam * mgamma(1, al, be)
VarX <- lam * mgamma(2, al, be)

## b)
p0 <- dpois(0, lam)
pk <- dpois(k0, lam)

Fx <- function(x) p0 + sum(pk * pgamma(x, k0 * al, be))
sapply(c(0, 60, 100), Fx)

## d)
VaR <- function(k) ifelse(
    k < p0, 0, optimize(function(x) abs(Fx(x) - k), c(0, 180))$min
)
sapply(c(0.01, 0.99), VaR)

## e)
TVaR <- function(k) 1/(1 - k) * sum(
    pk * k0 * al/be * pgamma(VaR(k), k0 * al + 1, be, lower = FALSE)
)
sapply(c(0.01, 0.99), TVaR)


###
### Question 4
###

rm(list = ls())

be <- c(1/2, 1/6, 1/12)
ci <- sapply(1:3, function(i) prod(be[-i]/(be[-i] - be[i])))

# i)
VaRX <- function(k) qexp(k, be)
VaRX(0.995)

# ii)
TVaRX <- function(k) mexp(1, be) + VaRX(k)
TVaRX(0.995)

# iii)
Fs <- function(x) sum(ci * pexp(x, be))
sapply(c(50, 80), Fs)

# iv)
VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 200))$min
VaRS(0.995)

# v)
TVaRS <- function(k) sum(ci * exp(-be * VaRS(k)) * (VaRS(k) + 1/be))/(1 - k)
TVaRS(0.995)

# vi)
BM <- function(k) sum(TVaRX(k)) - TVaRS(k)
BM(0.995)


###
### Question 5
###

rm(list = ls())

lam <- 0.28
mu <- 5.9
sig <- 1.22

m <- 1e5
set.seed(2018)
M <- X <- numeric(m)

for (i in 1:m)
{
    M[i] <- qpois(runif(1), lam)
    X[i] <- sum(qlnorm(runif(M[i]), mu, sig))
}

## b)
mean(X > 1000)

## c)
mean(pmax(X - 1000, 0))

## d)
VaR <- function(k) sort(X)[m * k]
VaR(0.99)

TVaR <- function(k) mean(X[X > VaR(k)])
TVaR(0.99)

## e)
n <- 100
EspN <- n * dpois(0:2, lam)
