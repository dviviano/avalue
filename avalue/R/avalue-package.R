#' avalue: Abstention Values for Optimal Follow-Up Experimentation
#'
#' The avalue package computes abstention-values (A-values), the break-even
#' welfare cost at which a decision-maker switches from collecting follow-up
#' evidence to making an immediate policy decision. It also produces the
#' logarithmic distribution and rank-comparison plots used in the empirical
#' application of Epanomeritakis and Viviano (2026).
#'
#' @section Main function:
#' Use [calculate_a_values()] to calculate A-values from point estimates and
#' either standard errors or exact two-sided normal p-values.
#'
#' @section Important convention:
#' The argument `w` is `s2 / s1`, the follow-up standard error divided by the
#' initial standard error. Smaller values therefore describe more precise
#' follow-up experiments.
#'
#' @references
#' Epanomeritakis, A. and Viviano, D. (2026). "When Is the Statistical
#' Evidence Strong Enough? Using Hypothesis Tests to Value Data Collection."
#' Unpublished manuscript.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "plot_a_value",
  "comparison_rank",
  "a_value_rank",
  "selection_group"
))
