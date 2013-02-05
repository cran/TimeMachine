\name{mu.mle}
\alias{mu.mle}
\title{Maximum-Likelihood Estimation of the Mutation Rate}
\description{Estimates the mutation rate using a maximum-likelihood approach.}
\usage{mu.mle(transitions, pi=NULL, population, n, mu.int, samples,
              threads=NULL, ...)}
\arguments{
    \item{transitions}{Transition matrix between types}
    \item{pi}{Stationary distribution associated with the transition matrix, or
              \code{NULL} to compute it}
    \item{population}{Vector representing the initial population}
    \item{n}{Target population size}
    \item{mu.int}{Vector specifying the end-points of the interval to be
                  searched for the optimal mutation rate}
    \item{samples}{Number of samples to simulate}
    \item{threads}{Number of threads used for simulations, or \code{NULL} to use
                   the default OpenMP setting}
    \item{...}{Further arguments passed to \code{optimize}}
}
\value{
    \item{mu.hat}{The maximum-likelihood estimate of the mutation rate}
    \item{loglik}{The log-likelihood value associated with the estimated
                  mutation rate}
}
\examples{
# Load example dataset
data(pdm)
transitions <- full.transitions(pdm$unitary.transitions, pdm$loci)
pi <- stationary.dist(transitions)

mu.int <- c(0.1, 10)
samples <- 10

# MLE of the mutation rate for a single target population size
n <- 10
mle.res <- mu.mle(transitions, pi, pdm$population, n, mu.int, samples)
print(mle.res)

# MLE of the mutation rate for 10 different target population sizes, including
# up to the MRCA (n = 1)
ns <- 1:10
mle.res <- sapply(ns, mu.mle, transitions=transitions, pi=pi,
                  population=pdm$population, mu.int=mu.int, samples=samples)
print(mle.res)
}

