\name{avalue-package}
\alias{avalue}
\alias{avalue-package}
\docType{package}
\title{Abstention Values for Optimal Follow-Up Experimentation}
\description{
The \pkg{avalue} package computes abstention-values (A-values), the break-even
welfare cost at which a decision-maker switches from collecting follow-up
evidence to making an immediate policy decision.
}
\details{
The main function, \code{\link{calculate_a_values}}, accepts vectors of point
estimates and either standard errors or exact two-sided normal p-values. It
returns the inputs and calculated A-values together with a logarithmic
histogram and rank-comparison plots.

The argument \code{w} is the follow-up standard error divided by the initial
standard error, \eqn{s_2 / s_1}. Thus, smaller values of \code{w} correspond
to more precise follow-up experiments.
}
\references{
Epanomeritakis, A. and Viviano, D. (2026).
\dQuote{When Is the Statistical Evidence Strong Enough? Using Hypothesis
Tests to Value Data Collection.} Unpublished manuscript.
}
\author{
Aristotelis Epanomeritakis and Davide Viviano
}
\seealso{
\code{\link{calculate_a_values}}
}
\keyword{package}
