################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Examen partiel informatique H16 ###########################################
################################################################################

###
### Question 1
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


###
### Question 2
###

rm(list = ls())

lam <- 0.006
al <- 1.5
be <- 1/1500
k0 <- 1:1000

p0 <- dpois(0, lam)
pk <- dpois(k0, lam)

Fx <- function(x) p0 + sum(pk * pgamma(x, k0 * al, be))

VaRX <- function(k) ifelse(k < p0, 0,
    optimize(function(x) abs(Fx(x) - k), c(0, 1e5))$min
)
VaRX(0.99)

TVaRX <- function(k) 1/(1 - k) * sum(
    pk * k0 * al/be * pgamma(VaRX(k), k0 * al + 1, be, lower = FALSE)
)
TVaRX(0.99)

## b)
n <- 1000
beW <- n * be

pw0 <- dpois(0, n * lam)
pkw <- dpois(k0, n * lam)

Fw <- function(x) pw0 + sum(pkw * pgamma(x, k0 * al, beW))

VaRW <- function(k) ifelse(
    k < pw0, 0, optimize(function(x) abs(Fw(x) - k), c(0, 50))$min
)
VaRW(0.99)

TVaRW <- function(k) 1/(1 - k) * sum(
    pkw * k0 * al/beW * pgamma(VaRW(k), k0 * al + 1, beW, lower = FALSE)
)
TVaRW(0.99)

## c)
BMVaR <- function(k) VaRX(k) - VaRW(k)
BMVaR(0.99)

BMTVaR <- function(k) TVaRX(k) - TVaRW(k)
BMTVaR(0.99)


###
### Question 3
###

rm(list = ls())

al <- c(4, 2)
be <- c(1/100, 1/200)
q <- 0.25

Fs <- function(x)
{
    kmax <- 100
    alpha <- sum(al)
    beta <- max(be)

    sig <- prod((be/beta)^al)
    zeta <- sapply(1:kmax, function(k) sum(al/k * (1 - be/beta)^k))

    xi <- 1
    for(k in 1:kmax) xi <- c(xi, sum(1:k * zeta[1:k] * rev(xi))/k)

    pk <- sig * xi
    sum(pk * pgamma(x, alpha + 0:kmax, beta))
}

Ft <- function(x) (1 - q) * pgamma(x, al[1], be[1]) + q * Fs(x)
sapply(c(100, 500, 800, 1000), Ft)
