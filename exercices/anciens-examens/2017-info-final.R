################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Examen final informatique H17 #############################################
################################################################################

###
### Question 1
###

rm(list = ls())

## b)
set.seed(20160419)
m <- 100000
U <- matrix(runif(3 * m), m, byrow = TRUE)

X1 <- (-log(U[, 1]))^(-0.5) * 1000
X2 <- qllogis(U[, 2], shape = 2, scale = 1000)
X3 <- -2000 * log(2 - 2^(U[, 3]))

S <- X1 + X2 + X3
c(X1[3], X2[3], X3[3])

## c)
phi <- mean(S > 20000)
phi

errs <- sqrt(phi * (1 - phi)/m)
IC <- phi + c(-1, 1) * qnorm(0.975) * errs
IC

## d)
VaR <- sort(S)[0.9999 * m]
mean(S[S > VaR])

## e)
sapply(list(X1, X2, X3), function(x) mean(x[S > VaR]))

## f)
mean(X2 * (X1 + X3 <= 1500 | X1 + X3 > 150000))


###
### Question 3
###

rm(list = ls())

n <- 5
q <- 0.3

## a)
EspM <- 0.1 * (n * q + 1)
EspX <- 1000 * EspM

## À TERMINER!

###
### Question 5
###

rm(list = ls())

al <- 1.5
be <- 0.01

## a)
Fs <- function(x) (pgamma(x, al, be) + pgamma(x, 3 * al, be))/2
Fs(100)

## b)
VaRS <- function(k) optimize(function(x) abs(Fs(x) - k), c(0, 1e4))$min
VaRS(0.9)

## c)
TVaRS <- function(k) 1/(2 * (1 - k)) * (
    al/be * pgamma(VaRS(k), al + 1, be, lower = FALSE) +
    3 * al/be * pgamma(VaRS(k), 3 * al + 1, be, lower = FALSE)
)
TVaRS(0.9)
