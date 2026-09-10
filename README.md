# avalue

`avalue` calculates abstention-values (A-values) for optimal follow-up
experimentation. An A-value is the break-even welfare cost below which a
decision-maker prefers to collect follow-up evidence instead of making an
immediate policy decision.

The package implements the framework in:

> Epanomeritakis, Aristotelis, and Davide Viviano. 2026. "When Is the
> Statistical Evidence Strong Enough? Using Hypothesis Tests to Value Data
> Collection." Unpublished manuscript.

## Installation

Install the source package from a local archive with:

```r
install.packages("avalue_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

```r
library(avalue)

estimates <- c(study_A = 0.20, study_B = -0.10, study_C = 0.55)

result <- calculate_a_values(
  point_estimate = estimates,
  std_error = c(0.10, 0.30, 0.20),
  w = 1,
  top_n = 1
)

result$data
result$plots$histogram
result$plots$p_value_ranks
result$plots$standard_error_ranks
```

The returned data frame contains the original inputs, z-statistics, p-values,
normalized and raw A-values, ranks, and selection indicators. The three plot
objects are ordinary `ggplot2` objects and can be modified or saved with
`ggplot2::ggsave()`.

## Input options

Supply exactly one of `std_error` and `p_value`.

```r
from_standard_errors <- calculate_a_values(
  point_estimate = c(a = 0.2, b = -0.3),
  std_error = c(0.1, 0.25)
)

from_p_values <- calculate_a_values(
  point_estimate = c(a = 0.2, b = -0.3),
  p_value = c(0.0455, 0.2301)
)
```

Direct p-value input assumes exact two-sided normal p-values. The package
recovers the standard error from the point estimate and p-value, so rounded
p-values can produce inaccurate A-values. Supply standard errors whenever they
are available.

## Follow-up precision

The paper defines

```text
w = follow-up standard error / initial standard error = s2 / s1.
```

Therefore, `w = 1` is an equal-precision follow-up and smaller values represent
more precise follow-up experiments. `w` may be one number or a vector with one
value per estimate.

## Ranking convention

Rank 1 is the largest value for A-values, p-values, and standard errors. If
`top_n` is supplied, each rule selects exactly that many estimates. Cutoff ties
are broken by `estimate_id`, following the empirical application.

## Interpreting units

A-values have the same welfare units as the point estimates and standard
errors. Estimates must be expressed in comparable welfare units before their
A-values are ranked across studies.


## Documentation and citation

After installation, open the function reference with
`?calculate_a_values`, read the worked guide with
`vignette("introduction", package = "avalue")`, and obtain the package and
methodology citations with `citation("avalue")`.
