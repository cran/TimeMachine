\name{full.transitions}
\alias{full.transitions}
\title{Full Transition Matrix Computation}
\description{Computes the full transition matrix from a given unitary transition
             matrix.}
\usage{full.transitions(unitary.transitions, loci)}
\arguments{
    \item{unitary.transitions}{Square matrix of order 2 specifying per-locus
                               transition probabilities}
    \item{loci}{Number of loci}
}
\value{A square matrix of order 2^loci specifying all transition probabilities.}

