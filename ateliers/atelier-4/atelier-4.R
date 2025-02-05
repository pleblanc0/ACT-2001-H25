################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Atelier 4 - 29 janvier 2025 ###############################################
################################################################################

####
#### Mutualisation de risques indépendants avec distributions identiques
####

rm(list = ls())	

al <- 0.5
be <- al/1000
EspW <- al/be
theta <- c(0.2, 0.1, 0.01)
n <- c(10, 100, 500)

# Value-at-Risk de W
VaRW <- function(k) qgamma(k, n * al, n * be)
VaRW(0.99)

# Tail-Value-at-Risk de W
TVaRW <- function(k) al/be * pgamma(VaRW(k), n * al + 1, n * be,
                             lower = FALSE)/(1 - k)
TVaRW(0.99)

# Probabilité d'être l'intervalle de confiance
intervalW <- function(n, theta)
{
    pgamma((1 + theta) * EspW, n * al, n * be, lower = FALSE) +
    pgamma((1 - theta) * EspW, n * al, n * be)
}
lapply(n, function(n) intervalW(n, theta))

# Obtenir la valeur de theta pour une probabilité donnée
lapply(n, function(n) sapply(c(0.1, 0.01), function(eta)
    optimize(function(theta) abs(intervalW(n, theta) - eta), c(0, 1.6))$min))

# Obtenir la plus petite valeur de n pour une probabilité donnée
sapply(c(0.1, 0.05), function(theta)
    optimize(function(n) abs(intervalW(n, theta) - 0.01), c(0, 1e4))$min)


####
#### Mutualisation de risques indépendants avec distributions non identiques
####

rm(list = ls())

al <- 0.5 * 1:4
be <- 1/1000
theta <- c(0.1, 0.01)

## Calculs pour X
EspX <- al/be
VarX <- al/be^2

# Value-at-Risk de X
VaRX <- function(k) qgamma(k, al, be)
VaRX(0.99)

# Tail-Value-at-Risk de X
TVaRX <- function(k) al/be * pgamma(VaRX(k), al + 1, be, lower = FALSE)/(1 - k)
TVaRX(0.99)

# Probabilité d'être l'intervalle de confiance
intervalX <- function(theta)
{
    pgamma((1 + theta) * EspX, al, be, lower = FALSE) +
    pgamma((1 - theta) * EspX, al, be)
}
lapply(theta, intervalX)

## Calculs pour Cj
alC <- sum(al)
beC <- sum(al)/al * be

EspC <- alC/beC
VarC <- alC/beC^2

# Value-at-Risk de C
VaRC <- function(k) qgamma(k, alC, beC)
VaRC(0.99)

# Tail-Value-at-Risk de C
TVaRC <- function(k) alC/beC * pgamma(VaRC(k), alC + 1, beC,
                               lower = FALSE)/(1 - k) 
TVaRC(0.99)

# Probabilité d'être l'intervalle de confiance
intervalC <- function(theta)
{
    pgamma((1 + theta) * EspC, alC, beC, lower = FALSE) +
    pgamma((1 - theta) * EspC, alC, beC)
}
lapply(theta, intervalC)

## Alternative avec S
EspS <- sum(EspX)
alS <- sum(al)

# Value-at-Risk de S
VaRS <- function(k) qgamma(k, alS, be)
EspX/EspS * VaRS(0.99)

# Tail-Value-at-Risk de S
TVaRS <- function(k) alS/be * pgamma(VaRS(k), alS + 1, be, low = FALSE)/(1 - k)
EspX/EspS * TVaRS(0.99)

# Probabilité d'être l'intervalle de confiance
intervalS <- function(theta)
{
    pgamma((1 + theta) * EspS, alS, be, lower = FALSE) +
    pgamma((1 - theta) * EspS, alS, be)
}
sapply(theta, intervalS)


####
#### Portefeuille avec 4 contrats et péril inondation
####

rm(list = ls())

x <- expand.grid(i1 = 0:1, i2 = 0:1, i3 = 0:1, i4 = 0:1)
probs <- c(0.4815, 0.1605, 0.0705, 0.0235, 0.0360, 0.0120, 0.0120, 0.0040,
           0.0360, 0.0120, 0.0120, 0.0040, 0.0090, 0.0030, 0.0930, 0.0310)

q <- colSums(x * probs)

## Densités bivariées
duos <- combn(names(x), 2)
f2 <- apply(duos, 2, function(v) sum(probs * x[[v[1]]] * x[[v[2]]]))
names(f2) <- apply(duos, 2, paste, collapse = '_')

## Densités trivariées
trios <- combn(names(x), 3)
f3 <- apply(trios, 2, function(v) sum(probs * apply(x[, v], 1, prod)))
names(f3) <- apply(trios, 2, paste, collapse = '_')

## Calculs pour I
EspI <- q
VarI <- q * (1 - q)
CovI <- zapsmall(f2 - q[duos[1, ]] * q[duos[2, ]])

VarCovI <- matrix(0, 4, 4)
diag(VarCovI) <- VarI
VarCovI[lower.tri(VarCovI)] <- CovI
VarCovI[upper.tri(VarCovI)] <- t(VarCovI)[upper.tri(VarCovI)]

## Calculs pour N
n <- 0:4
fn <- tapply(probs, rowSums(x), sum)

EspN <- sum(n * fn)                # sum(EspI)
VarN <- sum((n - EspN)^2 * fn)     # sum(VarCovI) = sum(VarI) + 2 * sum(CovI)

## Calculs pour X
b <- 4:1
bj <- expand.grid(i1 = c(0, 4), i2 = c(0, 3), i3 = c(0, 2), i4 = c(0, 1))

EspX <- b * q
VarX <- b^2 * q * (1 - q)
VarCovX <- VarCovI * outer(b, b)

## Calculs pour S
s <- 0:10
fs <- tapply(probs, rowSums(bj), sum)
Fs <- cumsum(fs)

EspS <- sum(s * fs)                # sum(EspX)
VarS <- sum((s - EspS)^2 * fs)     # sum(VarCovX)

VaRS <- function(k) s[min(which(Fs >= k))]
VaRS(0.95)

SLS <- function(d) sum(pmax(s - d, 0) * fs)
TVaRS <- function(k) SLS(VaRS(k))/(1 - k) + VaRS(k)
TVaRS(0.95)

## Calculs pour Cj
Cj <- outer(s, EspX/EspS)
EspC <- EspX
VarC <- (EspX/EspS)^2 * VarS
