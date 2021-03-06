% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{MicroBias.test}
\alias{MicroBias.test}
\title{Hypothesis testing for the log-linear model}
\usage{
MicroBias.test(mod, C.matrix = NULL, n.perm = 1000, tol = 1e-08, C.list = NULL)
}
\arguments{
\item{mod}{The fitted log-linear model from LLBias.fit}

\item{C.matrix}{The contrast matrix to be tested. The number of columns should be n.otu \eqn{\times} n.variable (number of variable in the design matrix \eqn{x} for estimation). The algorithm will assess the testablity of the contrasts.}

\item{n.perm}{Number of permutations to be conducted. Default = 1000.}

\item{tol}{Tolerance level for matrix inversion. Default is 10^-8. Should be a small number.}

\item{C.list}{Another way to specify the contrast matrix. It consists a list of matrices each with the same dimensionality of the \eqn{\beta} coefficient matrix. The function will decipher this list and make it into a C.matrix. Only one of C.matrix and C.list needs to be specified.}
}
\value{
p-value from the hypothesis testing
}
\description{
This function tests a contrast of the beta coefficients in the log-linear model
}
