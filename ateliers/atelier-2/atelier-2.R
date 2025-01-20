################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Atelier 2 - 16 janvier 2025 ###############################################
################################################################################

###
### Exemple 4
###

rm(list = ls())

q <- 0.5
lam <- 2
k <- 1:100

p0 <- 1 - q + q * dpois(0, lam)
pk <- q * dpois(k, lam)

fx <- c(p0, pk)
head(fx)


###
### Exemple 5
###

rm(list = ls())

n <- 10
q <- c(0.2, 0.3)
k <- 0:10

fx1 <- dbinom(k, n, q[1])
fx2 <- dbinom(k, n, q[2])

## 4)
s <- 0:20
fs <- convolve(fx1, rev(fx2), type = 'open')

## 5)
EspS <- sum(s * fs)
EspS


###
### Exemple 6
###

rm(list = ls())

lam <- c(2, 3)
k <- 0:100

fx1 <- dpois(k, lam[1])
fx2 <- dpois(k, lam[2])

## 4)
s <- 0:200
fs <- convolve(fx1, rev(fx2), type = 'open')

## 5)
EspS <- sum(s * fs)
EspS


###
### Exemple 7
###

rm(list = ls())

fx1 <- c(rep(0, 12), 0.3, rep(0, 1345 - 13), 0.3, rep(0, 19394 - 1346), 0.4)
fx2 <- c(0, 1/6, 0, 0, 1/6, 1/6, rep(0, 112 - 6), 1/6,
         rep(0, 4321 - 113), 1/6, rep(0, 55555 - 4322), 1/6)

fs <- convolve(fx1, rev(fx2), type = 'open')
s <- seq_along(fs) - 1

## 1)
EspS <- sum(s * fs)
EspS

## 3)
SL <- function(d) sum(pmax(s - d, 0) * fs)
SL(20000)


###
### Exemple 9
###

rm(list = ls())

w <- c(0.45, 0.55)
theta <- c(0.4, 5)
lam <- 8
k <- 0:50

fx1 <- sapply(k, function(x) sum(w * dpois(x, theta)))
fx2 <- dpois(k, lam)

## 1)
s <- 0:100
fs <- convolve(fx1, rev(fx2), type = 'open')

## 2)
Fs <- cumsum(fs)
FsInv <- function(u) s[min(which(Fs >= u))] 
FsInv(0.8)

## 3)
fs <- zapsmall(fs)     # éviter les erreurs numériques (valeur négative sinon)
FGM <- function(t) sum(exp(t * s) * fs)
FGM(0.67)
