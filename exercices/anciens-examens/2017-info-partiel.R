################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Examen partiel informatique H17 ###########################################
################################################################################

###
### Question 1
###

rm(list = ls())

r <- 0.4
q <- 10/11
al <- 1.5
be <- al/1000
k0 <- 1:1000

## a)
p0 <- dnbinom(0, r, q)
pk <- dnbinom(k0, r, q)

Fx <- function(x) p0 + sum(pk * pgamma(x, k0 * al, be))
Fx(30)

VaRX <- function(k) ifelse(
    k < p0, 0, optimize(function(x) abs(Fx(x) - k), c(0, 1e4))$min
)
VaRX(0.95)

TVaRX <- function(k) 1/(1 - k) * sum(
    pk * k0 * al/be * pgamma(VaRX(k), k0 * al + 1, be, lower = FALSE)
)
TVaRX(0.95)

## c)
n <- 200
beW <- n * be

pw0 <- dnbinom(0, n * r, q)
pwk <- dnbinom(k0, n * r, q)

Fw <- function(x) pw0 + sum(pwk * pgamma(x, k0 * al, beW))
Fw(30)

VaRW <- function(k) ifelse(
    k < pw0, 0, optimize(function(x) abs(Fw(x) - k), c(0, 500))$min
)
VaRW(0.95)

TVaRW <- function(k) 1/(1 - k) * sum(
    pwk * k0 * al/beW * pgamma(VaRW(k), k0 * al + 1, beW, lower = FALSE)
)
TVaRW(0.95)


###
### Question 2
###

rm(list = ls())

r <- 1.5
q <- 1/3
al <- 1.5
lam <- 50

## b)
m <- 1e5
set.seed(20170222)
M <- X <- numeric(m)

FbInv <- function(u) lam * ((1 - u)^(-1/al) - 1)

for (i in 1:m)
{
    M[i] <- qnbinom(runif(1), r, q)
    X[i] <- sum(FbInv(runif(M[i])))
}

## c)
mean(X > 1200)

## d)
mean(pmax(X - 1200, 0))

## e)
mean(exp(-0.001 * X))


###
### Question 3
###

rm(list = ls())

al <- 2.5
lam <- c(1/50, 20, 100)
tau <- c(1/2, 2.5, 2)

set.seed(20160419)
m <- 100000

U <- matrix(runif(3 * m), m, byrow = TRUE)

# il faut coder les fonctions de quantiles, le reste c'est comme d'habitude
