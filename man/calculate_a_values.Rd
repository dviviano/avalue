\name{calculate_a_values}
\alias{calculate_a_values}
\title{Calculate A-Values and Create Rank-Comparison Plots}
\usage{
calculate_a_values(
  point_estimate,
  std_error = NULL,
  p_value = NULL,
  w = 1,
  top_n = NULL,
  estimate_id = NULL,
  histogram_bins = 30L,
  rel.tol = 1e-09
)
}
\arguments{
\item{point_estimate}{A nonempty numeric vector of initial point estimates.
A named vector may be used to provide estimate IDs.}

\item{std_error}{An optional numeric vector of positive standard errors.
Supply exactly one of \code{std_error} and \code{p_value}. A scalar is
recycled.}

\item{p_value}{An optional numeric vector of exact two-sided normal p-values.
When supplied instead of \code{std_error}, the standard error is recovered
from the point estimate and p-value. A scalar is recycled.}

\item{w}{A nonnegative numeric scalar or vector equal to \eqn{s_2 / s_1}, the
follow-up standard error divided by the initial standard error. The default is
1, corresponding to an equal-precision follow-up.}

\item{top_n}{An optional positive integer giving the exact number of estimates
selected by each ranking rule. Cutoff ties are broken by estimate ID.}

\item{estimate_id}{An optional vector of unique estimate IDs. If omitted,
names on \code{point_estimate} are used; otherwise sequential IDs are
generated.}

\item{histogram_bins}{A positive integer specifying the number of logarithmic
histogram bins.}

\item{rel.tol}{A positive numerical tolerance used by integration and root
finding.}
}
\value{
An object of class \code{a_value_result}, which is a list with:
\itemize{
  \item \code{data}: a data frame containing the inputs, inferred quantities,
  A-values, ranks, and selection indicators.
  \item \code{plots}: a list of \pkg{ggplot2} objects named
  \code{histogram}, \code{p_value_ranks}, and
  \code{standard_error_ranks}.
  \item \code{settings}: the calculation and plotting settings.
}
}
\description{
Calculates abstention-values (A-values) for one or more initial estimates and
creates the distribution and ranking plots used to prioritize follow-up
experiments.
}
\details{
Under the conditions in the paper, the A-value is
\deqn{
A(X_1; s_1, w) =
s_1 \inf\{\gamma \geq 0: |X_1|/s_1 \geq \tau^\star(w,\gamma)\}.
}
The implementation evaluates the equivalent saddle-point first-order
conditions at \eqn{\tau = |X_1|/s_1}. The unconstrained-parameter result
requires the standardized parameter bound to satisfy
\eqn{\overline{\Delta}/s_1 \geq 0.76}.

Larger values receive smaller numerical ranks: rank 1 is the largest A-value,
p-value, or standard error. This is the follow-up experimentation convention
used in the empirical application, not the usual ranking of the smallest
p-values as most statistically significant.

When \code{p_value} is supplied, values must be strictly between zero and one
and their corresponding point estimates must be nonzero. The standard error
is recovered as
\deqn{
s_1 = |X_1| / \Phi^{-1}(1-p/2).
}
Because this calculation treats the p-values as exact, standard errors should
be supplied when reported p-values have been rounded.
}
\examples{
estimates <- c(study_A = 0.20, study_B = -0.10, study_C = 0.55)

result <- calculate_a_values(
  point_estimate = estimates,
  std_error = c(0.10, 0.30, 0.20),
  w = 1,
  top_n = 1
)

result$data
result$plots$histogram

# Exact two-sided normal p-values may be supplied instead.
result_from_p <- calculate_a_values(
  point_estimate = c(a = 0.2, b = -0.3),
  p_value = c(0.10, 0.25)
)
}
\references{
Epanomeritakis, A. and Viviano, D. (2026).
\dQuote{When Is the Statistical Evidence Strong Enough? Using Hypothesis
Tests to Value Data Collection.} Unpublished manuscript.
}
\seealso{
\code{\link{plot.a_value_result}}, \code{\link{print.a_value_result}}
}
