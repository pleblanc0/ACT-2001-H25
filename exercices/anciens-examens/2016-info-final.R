################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Examen final informatique H16 #############################################
################################################################################

###
### Question 1
###

rm(list = ls())

## a)
EspX <- 164.8721
VarX <- 46707.74

al1 <- 2/(1 - EspX^2/VarX)
lam <- (al1 - 1) * EspX

be <- EspX/VarX
al2 <- be * EspX

sig <- sqrt(log(VarX/EspX^2 + 1))
mu <- log(EspX) - sig^2/2

## b)
set.seed(20160419)
m <- 1000000

U <- matrix(runif(3 * m), nrow = m, byrow = TRUE)

X1 <- lam * ((1 - U[, 1])^(-1/al1) - 1)
X2 <- qgamma(U[, 2], al2, be)
X3 <- qlnorm(U[, 3], mu, sig)

## c)
S <- X1 + X2 + X3
S[c(3, 4)]

## d)
phi <- mean(S > 1000)
errs <- sqrt(phi * (1 - phi)/m)

IC <- phi + c(-1, 1) * qnorm(0.975) * errs
IC

## e)
VaR <- sort(S)[0.999 * m]
VaR

## f)
mean(S[S > VaR])

## g)
sapply(list(X1, X2, X3), function(x) x[which(S == VaR)])

## h)
sapply(list(X1, X2, X3), function(x) mean(x[S > VaR]))


###
### Question 4
###

rm(list = ls())

q <- 0.1
al <- 0.5
be <- 0.05

## a)
Fx <- function(x) 1 - q + q * pgamma(x, al, be)
Fx(40)

VaRX <- function(k) ifelse(k < 1 - q, 0, qgamma((k - (1 - q))/q, al, be))
sapply(c(0.8, 0.99), VaRX)
VaRX(0.99)

TVaRX <- function(k) q/(1 - k) * al/be * pgamma(VaRX(k), al + 1, be, lower = F)
sapply(c(0.8, 0.99), TVaRX)

## b)
n <- 200
k0 <- 1:n
beW <- n * be

p0 <- dbinom(0, n, q)
pk <- dbinom(k0, n, q)

Fw <- function(x) p0 + sum(pk * pgamma(x, k0 * al, beW))
Fw(1)

VaRW <- function(k) ifelse(
    k < p0, 0, optimize(function(x) abs(Fw(x) - k), c(0, 3))$min
)
sapply(c(0.8, 0.99), VaRW)

TVaRW <- function(k) 1/(1 - k) * sum(
    pk * k0 * al/beW * pgamma(VaRW(k), k0 * al + 1, beW, lower = FALSE)
)
sapply(c(0.8, 0.99), TVaRW)

## c)
BMVaR <- function(k) VaRX(k) - VaRW(k)
sapply(c(0.8, 0.99), BMVaR)

BMTVaR <- function(k) TVaRX(k) - TVaRW(k)
sapply(c(0.8, 0.99), BMTVaR)
