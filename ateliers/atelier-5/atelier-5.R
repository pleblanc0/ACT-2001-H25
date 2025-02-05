################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Atelier 5 - 5 février 2025 ################################################
################################################################################

###
### Partiel Informatique 2018 Question 1
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
### Partiel Informatique 2018 Question 4
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
### Partiel Informatique 2016 Question 1
###

rm(list = ls())

al <- 0.5
be <- 0.5

## b)
Fs <- function(x, n) pgamma(x, n * al, be)
Fs(200, 200)

## c)
EspX <- al/be
EspX

## d)
ruine <- function(n, prime, u) 1 - Fs(n * prime + u, n)
n <- c(10, 400, 1000, 2000)

lapply(c(0.9, 1.1), function(prime) ruine(n, prime, 20))
optimize(function(prime) abs(ruine(1000, prime, 20) - 0.01), c(0, 2))$min
