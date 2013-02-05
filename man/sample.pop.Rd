\name{sample.pop}
\alias{sample.pop}
\title{Sample Populations from the Wright-Fisher Model}
\description{Samples a population of the desired size from the Wright-Fisher
             model.}
\usage{sample.pop(transitions, pi=NULL, pop.size, mu)}
\arguments{
    \item{transitions}{Transition matrix between types}
    \item{pi}{Stationary distribution associated with the transition matrix, or
              \code{NULL} to compute it}
    \item{pop.size}{Desired population size}
    \item{mu}{Mutation rate}
}
\value{A vector representing the sampled population.}

