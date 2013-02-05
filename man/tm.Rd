\name{tm}
\alias{tm}
\alias{density.tm}
\alias{mean.tm}
\alias{print.tm}
\alias{quantile.tm}
\alias{summary.tm}
\title{Log-likelihood Estimation using the Time Machine}
\description{Estimates the log-likelihood for the given parameter values using
             the Time Machine.}
\usage{
    tm(transitions, pi=NULL, population, n, mu, samples, threads=NULL)

    \method{density}{tm}(x, ...)

    \method{mean}{tm}(x, ...)

    \method{print}{tm}(x, ...)

    \method{quantile}{tm}(x, ...)

    \method{summary}{tm}(object, ..., digits=max(3, getOption("digits")-3))
}
\arguments{
    \item{transitions}{Transition matrix between types}
    \item{pi}{Stationary distribution associated with the transition matrix, or
              \code{NULL} to compute it}
    \item{population}{Vector representing the initial population}
    \item{n}{Target population size}
    \item{mu}{Mutation rate}
    \item{samples}{Number of samples to simulate}
    \item{threads}{Number of threads used for simulations, or \code{NULL} to use
                   the default OpenMP setting}
    \item{x,object}{An object of class \code{tm} as returned by a call to
                  \code{tm}}
    \item{...}{Further arguments passed to the corresponding generic function}
    \item{digits}{Minimum number of significant digits for number formatting}
}
\value{
    An object of class \code{tm} representing the result of the simulations.
    \item{population}{Vector representing the initial population}
    \item{n}{Target population size}
    \item{mu}{Mutation rate}
    \item{logliks}{Vector of simulated log-likelihoods}
    \item{corrections}{Vector of correction terms}
    \item{coalescent.events}{Matrix of column vectors (one per simulation)
                             containing SENs when coalescent events were
                             simulated}
    \item{final.populations}{Matrix of columns vectors (one per simulation)
                             containing final population sizes for each type}
    \item{simulation.times}{Computation time for each simulation (in seconds)}
    \item{total.time}{Total computation time (in seconds)}
}
\examples{
# Load example dataset
data(pdm)

transitions <- full.transitions(pdm$unitary.transitions, pdm$loci)
pi <- stationary.dist(transitions)
n <- 10
samples <- 10

# Estimate log-likelihood for a fixed value of the mutation rate
mu <- 1
est.res <- tm(transitions, pi, pdm$population, n, mu, samples)
print(est.res)

# Estimate log-likelihood for different values of the mutation rate
mus <- seq(0.1, 10, 0.1)
estimates <- sapply(mus, function(mu) {
    tm(transitions, pi, pdm$population, n, mu, samples)
}, simplify=FALSE)

# Compute mean log-likelihood for each value of the mutation rate mu
mean.logliks <- do.call(rbind, lapply(estimates, function(x) {
    c(x$mu, mean(x))
}))

# Plot mean log-likelihoods
plot(mean.logliks, pch=20, xlab=expression(mu), ylab="Mean log-likelihood",
     main=paste("Mean log-likelihoods - Sample size:", samples))

# Estimate log-likelihood for different target population sizes and fixed
# mutation rate
ns <- 2:50
mu <- 1
estimates <- sapply(ns, function(n) {
    tm(transitions, pi, pdm$population, n, mu, samples)
}, simplify=FALSE)

# Compute mean correction term/log-likelihood ratios for each target population size
mean.corrections <- do.call(rbind, lapply(estimates, function(x) {
    c(x$n, mean(x$correction / x$logliks))
}))

# Plot mean correction term/log-likelihood ratios
plot(mean.corrections[,1], rev(mean.corrections[,2]), pch=20,
     xlab="Target population size", ylab="Mean correction term/log-likelihood",
     axes=FALSE)
axis(1, at=mean.corrections[,1], labels=rev(mean.corrections[,1]))
axis(2)
}

