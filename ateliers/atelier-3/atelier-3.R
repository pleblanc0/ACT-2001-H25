################################################################################
#### Introduction à l'actuariat II - Hiver 2025 ################################
#### École d'Actuariat - Université Laval ######################################
#### Atelier 3 - 22 janvier 2025 ###############################################
################################################################################

###
### Exemple pages 7-10
###

rm(list = ls())

x <- c(0, 100, 400, 1000, 2000)
fx <- c(0.3, 0.15, 0.4, 0.1, 0.05)
sum(fx)     # Vérification

kap <- c(0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999)
k1k2 <- list(c(0.8, 0.9), c(0.8, 0.995), c(0.9, 0.995))

## Value at Risk
VaR <- function(k) x[min(which(cumsum(fx) >= k))]
sapply(kap, VaR)

## Tail Value at Risk
SL <- function(d) sum(pmax(x - d, 0) * fx)
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
sapply(kap, TVaR)

## Range Value at Risk
RVaR <- function(k1, k2) ((1 - k1) * TVaR(k1) - (1 - k2) * TVaR(k2))/(k2 - k1)
sapply(k1k2, function(k) RVaR(k[1], k[2]))


###
### Exemple pages 11-14
###

rm(list = ls())

a1 <- 0.2; a2 <- 3
x <- 0:2000; k <- 1:2000
fx <- c(1 - a1, a1 * ((k/2000)^a2 - ((k - 1)/2000)^a2))
sum(fx)     # Vérification

kap <- c(0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.95, 0.99, 0.999, 0.9999)
k1k2 <- list(c(0.8, 0.9), c(0.8, 0.995), c(0.9, 0.995))

## Value at Risk
Fx <- cumsum(fx)
VaR <- function(k) x[min(which(Fx >= k))]
sapply(kap, VaR)

## Tail Value at Risk
SL <- function(d) sum(pmax(x - d, 0) * fx)
TVaR <- function(k) SL(VaR(k))/(1 - k) + VaR(k)
sapply(kap, TVaR)

## Range Value at Risk
RVaR <- function(k1, k2) ((1 - k1) * TVaR(k1) - (1 - k2) * TVaR(k2))/(k2 - k1)
sapply(k1k2, function(k) RVaR(k[1], k[2]))


###
### Exemple pages 16-17
###

rm(list = ls())

# Un support jusqu'à 300, par exemple, est largement suffisant
x <- 0:300

# Créer les vecteurs pour les fonctions de masse de probabilité
f1 <- dpois(x, lambda = 4)
f2 <- dnbinom(x, size = 0.4, prob = 1/11)
f3 <- dnbinom(x, size = 4, prob = 1/2)
f4 <- 0.8 * dpois(x, lambda = 2) + 0.2 * dpois(x, lambda = 12)

##
## Méthode 1: Utiliser des fonctions
##

# Créer une liste contenant les fonctions de masse de probabilité
f_list <- list(f1, f2, f3, f4)

## Espérances
sapply(f_list, function(f) sum(x * f))

## Variances
sapply(f_list, function(f) sum((x - sum(x * f))^2 * f))

## Stop-loss
SL <- function(d, f) sum(pmax(x - d, 0) * f)
sapply(f_list, function(f) SL(5, f))

## Espérance tronquée
EspTr <- function(d, f) sum(x * f * I(x > d))
sapply(f_list, function(f) EspTr(5, f))

## Value at Risk
VaR <- function(k, f) x[min(which(cumsum(f) >= k))]
sapply(f_list, function(f) VaR(0.99, f))

## Tail Value at Risk
TVaR <- function(k, f) SL(VaR(k, f), f)/(1 - k) + VaR(k, f)
sapply(f_list, function(f) TVaR(0.99, f))

## Range Value at Risk
RVaR <- function(k1, k2, f)
    ((1 - k1) * TVaR(k1, f) - (1 - k2) * TVaR(k2, f))/(k2 - k1)
sapply(f_list, function(f) RVaR(0.99, 0.995, f))

##
## Méthode 2: Manuellement remplacer f, d et k par leurs valeurs
##

## Espérances
sum(x * f)

## Variances
sum((x - sum(x * f))^2 * f)

## Stop-loss
SL <- sum(pmax(x - d, 0) * f)

## Espérance tronquée
sum(x * f * I(x > d))

## Value at Risk
Fx <- cumsum(f)
x[min(which(Fx >= k))]

## Tail Value at Risk
# On ne peut pas calculer avec les fonctions, il faut remplacer les valeurs.
SL(VaR(k))/(1 - k) + VaR(k)

## Range Value at Risk
# On ne peut pas calculer avec les fonctions, il faut remplacer les valeurs.
((1 - k1) * TVaR(k1) - (1 - k2) * TVaR(k2))/(k2 - k1)


###
### Exemple page 19
###

rm(list = ls())

al <- c(0.5, 1, 2.5, 5)
be <- al/100

## Espérance
EspX <- al/be
EspX

## Variance
VarX <- al/be^2
VarX

## Fonction stop-loss
SL <- function(d) al/be * pgamma(d, al + 1, be, lower = FALSE) -
                  d * pgamma(d, al, be, lower = FALSE)

## Espérance tronquée
EspTr <- function(d) al/be * pgamma(d, al + 1, be, lower = FALSE)

## Value at Risk
VaR <- function(k) qgamma(k, al, be)
VaR(0.99)

## Tail Value at Risk
TVaR <- function(k) EspTr(VaR(k))/(1 - k)
TVaR(0.99)

## Range Value at Risk
RVaR <- function(k1, k2) ((1 - k1) * TVaR(k1) - (1 - k2) * TVaR(k2))/(k2 - k1)
RVaR(0.9, 0.99)


###
### Exemple page 20
###

rm(list = ls())

sig <- c(0.2, 0.6, 0.8, 1.2)
mu <- log(100) - sig^2/2

## Espérance
EspX <- exp(mu + sig^2/2)
EspX

## Variance
VarX <- exp(2 * mu + sig^2) * (exp(sig^2) - 1)
VarX

## Fonction stop-loss
SL <- function(d) exp(mu + sig^2/2) * (1 - pnorm((log(d) - mu - sig ^2)/sig)) -
                  d * pnorm((log(d) - mu)/sig, lower = FALSE)

## Espérance tronquée
EspTr <- function(d) exp(mu + sig^2/2) * (1 - pnorm((log(d) - mu - sig^2)/sig))

## Value at Risk
VaR <- function(k) qlnorm(k, mu, sig)
VaR(0.99)

## Tail Value at Risk
TVaR <- function(k) EspTr(VaR(k))/(1 - k)
TVaR(0.99)

## Range Value at Risk
RVaR <- function(k1, k2) ((1 - k1) * TVaR(k1) - (1 - k2) * TVaR(k2))/(k2 - k1)
RVaR(0.9, 0.99)


###
### Loi Weibull
###

rm(list = ls())

## Paramétrisation en R
# tau --> shape
# beta --> 1 / scale
?pweibull

be <- 2
tau <- c(0.5, 2)

## Fonctions de densité
fx <- function(x, tau, be) be * tau * (be * x)^(tau - 1) * exp(-(be * x)^tau)
fx(5, tau[1], be)
dweibull(5, tau[1], 1/be)

# Courbes de densité
curve(dweibull(x, tau[1], 1/be), lwd = 3, col = 'blue', xlim = c(0, 10),
      main = 'Densités de la loi Weibull', ylab = 'f(x)', lty = 1)
curve(dweibull(x, tau[2], 1/be), lwd = 3, col = 'red', add = TRUE, lty = 1)
legend(x = 'topright', legend = c('tau = 2', 'tau = 0.5'),
       col = c('blue', 'red'), cex = 1.5, lwd = 3, lty = 1)

## Fonctions de répartition
Fx <- function(x, be, tau) 1 - exp(-(be * x)^tau)
Fx(5, be, tau[1])
pweibull(5, tau[1], 1/be)

# Courbes de répartition
curve(pweibull(x, tau[1], 1/be), lwd = 3, col = 'blue', xlim = c(0, 10),
      main = 'Répartitions de la loi Weibull', ylab = 'f(x)')
curve(pweibull(x, tau[2], 1/be), lwd = 3, col = 'red', add = TRUE, lty = 1)
legend(x = 'bottomright', legend = c('tau = 2', 'tau = 0.5'),
       col = c('blue', 'red'), cex = 1.5, lwd = 3, lty = 1)

## Fonction quantile
quantile <- function(u, tau, be) ((-log(1 - u))^(1/tau))/be
qweibull(0.5, tau[1], 1/be)
quantile(0.5, tau[1], be)

# Courbes de quantile
curve(qweibull(x, tau[1], 1/be), lwd = 3, col = 'blue', xlim = c(0, 1),
      main = 'Fonctions quantiles', ylab = 'VaR(X)', lty = 1, xlab = 'kappa')
curve(qweibull(x, tau[2], 1/be), lwd = 3, col = 'red', add = TRUE, lty = 1)
legend(x = 'bottomright', legend = c('tau = 2', 'tau = 0.5'),
       col = c('blue', 'red'), cex = 1.5, lwd = 3, lty = 1)

## Tail Value at Risk
VaR <- function(u, tau, be) ((-log(1 - u))^(1/tau))/be
TVaR <- function(u, tau, be) 1/(be * (1 - u)) * gamma(1 + 1/tau) * pgamma(
      -log(1 - u), 1 + 1/tau, 1, lower = FALSE
)
TVaR(0.95, tau[1], be)

# Courbes de TVaR
curve(sapply(x, TVaR, tau[1], be), lwd = 3, col = 'blue', xlim = c(0, 1),
      main = 'TVaR de la Weibull', ylab = 'TVaR(X)', lty = 1, xlab = 'kappa')
curve(sapply(x, TVaR, tau[2], be), lwd = 3, col = 'red', add = TRUE, lty = 1)
legend(x = 'topleft', legend = c('tau = 2', 'tau = 0.5'),
       col = c('blue', 'red'), cex = 1.5, lwd = 3, lty = 1)
