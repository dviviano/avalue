\name{a_value_result}
\alias{a_value_result}
\alias{print.a_value_result}
\alias{plot.a_value_result}
\title{Print and Plot A-Value Results}
\usage{
\method{print}{a_value_result}(x, ...)

\method{plot}{a_value_result}(
  x,
  y = NULL,
  which = c("all", "histogram", "p_value_ranks", "standard_error_ranks"),
  ...
)
}
\arguments{
\item{x}{An object returned by \code{\link{calculate_a_values}}.}

\item{...}{Additional arguments. For \code{print()}, these are passed to the
data-frame print method.}

\item{y}{Unused; included for compatibility with the \code{plot()} generic.}

\item{which}{The plot to display, or \code{"all"} to display all three plots
in sequence.}
}
\value{
Both methods return \code{x} invisibly.
}
\description{
The print method displays the augmented data frame and lists the available
plot objects. The plot method displays one or all of the paper-style plots.
}
\examples{
result <- calculate_a_values(
  point_estimate = c(a = 0.1, b = 0.3),
  std_error = c(0.2, 0.2),
  top_n = 1
)

print(result)
plot(result, which = "p_value_ranks")
}
