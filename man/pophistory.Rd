\name{pophistory}
\alias{pophistory}
\alias{pophistory.tm}
\title{History of Population Size}
\description{Estimates the sequential event number (SEN) quantiles for each
             simulated population size.}
\usage{
    \method{pophistory}{tm}(x, probs=seq(0, 1, 0.25))
}
\arguments{
    \item{x}{An object of class \code{tm} as returned by a call to \code{tm}}
    \item{probs}{Probabilities passed to \code{quantile}}
}
\value{
    A matrix of column vectors containing the estimated SEN quantiles, with row
    names specifying the population size.
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

# Plot minimum, median and maximum SEN as a function of population size
h <- pophistory(est.res, probs=c(0, 0.5, 1))
matplot(h, rownames(h), type="s", lty=c(3, 1, 3), col="black",
        xlab="Population size", ylab="SEN")
}

