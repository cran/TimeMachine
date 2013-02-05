\name{time-dist}
\alias{time-dist}
\alias{ttc}
\alias{ttm}
\alias{ttc.tm}
\alias{ttm.tm}
\title{Time to Coalescent Event or Mutation}
\description{Estimates the distribution of the times (expressed in SEN) between
             two coalescent events or two mutations.}
\usage{
    \method{ttc}{tm}(x)

    \method{ttm}{tm}(x)
}
\arguments{
    \item{x}{An object of class \code{tm} as returned by a call to \code{tm}}
}
\value{
    A vector representing the distribution of the time to coalescent event for
    \code{ttc}, or the time to mutation for \code{ttm}.
}
\examples{
# Load example dataset
data(pdm)

transitions <- full.transitions(pdm$unitary.transitions, pdm$loci)
pi <- stationary.dist(transitions)
n <- 10
mu <- 1
samples <- 10
est.res <- tm(transitions, pi, pdm$population, n, mu, samples)

# Plot distributions of time to coalescent event and time to mutation
par(mfrow=c(1, 2))
plot(ttc(est.res), type="b", pch=20, lty=3, xlab="SEN", ylab="",
     main="Time to coalescent event")
plot(ttm(est.res), type="b", pch=20, lty=3, xlab="SEN", ylab="",
     main="Time to mutation")
}

