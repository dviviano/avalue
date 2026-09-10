# A-values for optimal follow-up experimentation
#
# This file implements the abstention-value (A-value) in
# "Optimal Follow-up Experimentation" and returns the three plots used in the
# paper's empirical application.  The paper previously called this quantity a
# C-value.
#
# Paper notation:
#   X1 = point estimate from the initial study
#   s1 = standard error of X1
#   w  = s2 / s1, the follow-up standard error divided by the initial standard
#        error.  Thus, a smaller w means a more precise follow-up study.
#
# Under the conditions in Corollary 1, the A-value is
#   s1 * inf{gamma >= 0 : abs(X1) / s1 >= tau_star(w, gamma)}.
#
# The numerical implementation evaluates the saddle-point first-order
# conditions at tau = abs(X1) / s1.  It uses base R for the numerical work and
# ggplot2 for the returned plot objects.

.avalue_integral_terms <- function(tau, d, w, rel.tol) {
  # The second-period error-probability integral in Theorem 1 and its
  # derivative with respect to the standardized effect d.
  if (tau == 0 || w <= sqrt(.Machine$double.eps)) {
    return(c(I = 0, I_d = 0))
  }

  i_fun <- function(x) {
    a <- -w * x - d / w
    pnorm(a) * dnorm(x - d)
  }
  id_fun <- function(x) {
    a <- -w * x - d / w
    (-dnorm(a) / w + pnorm(a) * (x - d)) * dnorm(x - d)
  }

  i_val <- integrate(
    i_fun,
    lower = -tau,
    upper = tau,
    rel.tol = rel.tol,
    abs.tol = 0,
    subdivisions = 200L,
    stop.on.error = TRUE
  )$value
  id_val <- integrate(
    id_fun,
    lower = -tau,
    upper = tau,
    rel.tol = rel.tol,
    abs.tol = 0,
    subdivisions = 200L,
    stop.on.error = TRUE
  )$value

  c(I = i_val, I_d = id_val)
}

.avalue_gamma_from_tau_d <- function(tau, d, w) {
  # At an interior saddle point, the derivative of worst-case regret with
  # respect to tau is zero.  Solving that equation gives gamma = c / s1.
  denominator <- dnorm(tau - d) + dnorm(tau + d)
  if (!is.finite(denominator) || denominator == 0) return(0)

  if (w <= sqrt(.Machine$double.eps)) {
    f_plus <- 0
    f_minus <- 0
  } else {
    f_plus <- pnorm(-w * tau - d / w) * dnorm(tau - d)
    f_minus <- pnorm(w * tau - d / w) * dnorm(tau + d)
  }

  gamma <- d * (dnorm(tau + d) - f_plus - f_minus) / denominator
  max(0, gamma)
}

.avalue_d_first_order_condition <- function(d, tau, w, rel.tol) {
  terms <- .avalue_integral_terms(tau, d, w, rel.tol)
  gamma <- .avalue_gamma_from_tau_d(tau, d, w)

  p <- pnorm(-tau - d) + unname(terms["I"])
  p_d <- -dnorm(tau + d) + unname(terms["I_d"])
  q_d <- -dnorm(tau - d) + dnorm(tau + d)

  p + d * p_d + gamma * q_d
}

.normalized_a_value <- function(z, w, rel.tol = 1e-9) {
  tau <- abs(z)
  if (is.infinite(tau)) return(0)

  # Corollary 1 shows that the least-favorable standardized effect is below
  # 0.76.  It is interior in the cases covered by the paper.
  d_lo <- 1e-10
  d_hi <- 0.76
  objective <- function(d) {
    .avalue_d_first_order_condition(d, tau, w, rel.tol)
  }

  f_lo <- objective(d_lo)
  f_hi <- objective(d_hi)
  if (!is.finite(f_lo) || !is.finite(f_hi)) {
    stop(
      "Numerical integration failed while locating the least-favorable ",
      "standardized effect.",
      call. = FALSE
    )
  }

  if (f_lo * f_hi > 0) {
    # Defensive fallback for extreme inputs: locate a sign-changing interval.
    d_grid <- seq(d_lo, d_hi, length.out = 101L)
    f_grid <- vapply(d_grid, objective, numeric(1))
    left <- f_grid[-length(f_grid)]
    right <- f_grid[-1L]
    changes <- which(is.finite(left) & is.finite(right) & left * right <= 0)
    if (length(changes) == 0L) {
      stop(
        "Could not bracket the least-favorable standardized effect. ",
        "Try a less extreme estimate/standard-error ratio or a looser rel.tol.",
        call. = FALSE
      )
    }
    d_lo <- d_grid[changes[1L]]
    d_hi <- d_grid[changes[1L] + 1L]
  }

  d_star <- uniroot(
    objective,
    lower = d_lo,
    upper = d_hi,
    tol = max(rel.tol, sqrt(.Machine$double.eps))
  )$root

  .avalue_gamma_from_tau_d(tau, d_star, w)
}

.av_recycle <- function(x, n, label) {
  if (!(length(x) %in% c(1L, n))) {
    stop(label, " must have length 1 or length ", n, ".", call. = FALSE)
  }
  rep(x, length.out = n)
}

.av_bad_ids <- function(ids, rows) {
  paste(utils::head(ids[rows], 8L), collapse = ", ")
}

.av_rank <- function(values, ids, eligible) {
  out <- rep(NA_integer_, length(values))
  if (length(eligible) == 0L) return(out)

  # The paper breaks cutoff ties by program ID.  Character ordering gives the
  # same deterministic behavior for general estimate IDs.
  ordered <- eligible[order(-values[eligible], ids[eligible], method = "radix")]
  out[ordered] <- seq_along(ordered)
  out
}

.av_selection_group <- function(selected_a, selected_other) {
  label <- ifelse(
    selected_a & selected_other,
    "Both",
    ifelse(
      selected_a,
      "A-value only",
      ifelse(selected_other, "Other measure only", "Neither")
    )
  )
  factor(
    label,
    levels = c("A-value only", "Other measure only", "Both", "Neither")
  )
}

.av_rank_breaks <- function(n) {
  if (n <= 1L) return(1)
  candidates <- unique(c(1L, as.integer(pretty(c(1, n), n = 4)), n))
  candidates[candidates >= 1L & candidates <= n]
}

.av_plot_theme <- function() {
  ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "plain"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "plain"),
      legend.box = "horizontal"
    )
}

.av_histogram <- function(data, top_n, histogram_bins) {
  plot_data <- data[is.finite(data$a_value), , drop = FALSE]
  positive <- plot_data$a_value[plot_data$a_value > 0]

  zero_note <- NULL
  if (length(positive) == 0L) {
    floor_value <- 1e-300
  } else {
    floor_value <- max(min(positive) / 10, .Machine$double.xmin)
  }
  nonpositive <- plot_data$a_value <= 0
  plot_data$plot_a_value <- plot_data$a_value
  plot_data$plot_a_value[nonpositive] <- floor_value
  if (any(nonpositive)) {
    zero_note <- paste0(
      sum(nonpositive),
      " nonpositive A-value(s) are shown at ",
      format(floor_value, scientific = TRUE, digits = 3),
      " so they can appear on the logarithmic axis."
    )
  }

  caption <- zero_note
  if (!is.null(top_n)) {
    cutoff <- min(plot_data$plot_a_value[plot_data$selected_by_a_value])
    cutoff_note <- paste0(
      "Red dotted line: smallest A-value in the top-", top_n, " selection."
    )
    caption <- paste(c(cutoff_note, zero_note), collapse = " ")
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = plot_a_value)) +
    ggplot2::geom_histogram(
      bins = histogram_bins,
      fill = "#23858A",
      colour = "white",
      linewidth = 0.2
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(
      title = "Distribution of A-values",
      x = "A-value (logarithmic scale)",
      y = "Number of estimates",
      caption = caption
    ) +
    .av_plot_theme()

  if (!is.null(top_n)) {
    p <- p + ggplot2::geom_vline(
      xintercept = cutoff,
      colour = "#D9433F",
      linetype = "dotted",
      linewidth = 0.8
    )
  }
  p
}

.av_rank_scatter <- function(data, comparison, top_n) {
  if (comparison == "p_value") {
    plot_data <- data[
      is.finite(data$a_value_rank) & is.finite(data$p_value_rank),
      ,
      drop = FALSE
    ]
    plot_data$comparison_rank <- plot_data$p_value_rank
    plot_data$selection_group <- plot_data$p_value_selection_group
    comparison_title <- "p-value"
  } else {
    plot_data <- data[
      is.finite(data$a_value_rank) & is.finite(data$standard_error_rank),
      ,
      drop = FALSE
    ]
    plot_data$comparison_rank <- plot_data$standard_error_rank
    plot_data$selection_group <- plot_data$standard_error_selection_group
    comparison_title <- "standard error"
  }

  n <- nrow(plot_data)
  limits <- if (n == 1L) c(0.5, 1.5) else c(1, n)
  rank_breaks <- .av_rank_breaks(n)
  palette_fill <- c(
    "A-value only" = "#D9433F",
    "Other measure only" = "#2478B8",
    "Both" = "#D9433F",
    "Neither" = "#AEB6BC"
  )
  palette_colour <- c(
    "A-value only" = "#D9433F",
    "Other measure only" = "#2478B8",
    "Both" = "#2478B8",
    "Neither" = "#AEB6BC"
  )
  legend_title <- if (is.null(top_n)) NULL else paste0("Top ", top_n, " selection")

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = comparison_rank,
      y = a_value_rank,
      fill = selection_group,
      colour = selection_group
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      colour = "#9DA5AA",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(shape = 21, size = 2.8, stroke = 0.85) +
    ggplot2::scale_fill_manual(
      values = palette_fill,
      breaks = names(palette_fill),
      drop = FALSE,
      name = legend_title
    ) +
    ggplot2::scale_colour_manual(
      values = palette_colour,
      breaks = names(palette_colour),
      drop = FALSE,
      name = legend_title
    ) +
    ggplot2::scale_x_continuous(limits = limits, breaks = rank_breaks) +
    ggplot2::scale_y_continuous(limits = limits, breaks = rank_breaks) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0("A-value ranks versus ", comparison_title, " ranks"),
      x = paste0(
        if (comparison == "p_value") "p-value" else "Standard error",
        " rank (1 = largest)"
      ),
      y = "A-value rank (1 = largest)"
    ) +
    .av_plot_theme()

  if (is.null(top_n)) {
    p <- p + ggplot2::guides(fill = "none", colour = "none")
  }
  p
}

#' Calculate A-values and reproduce the paper-style empirical plots
#'
#' @param point_estimate Numeric vector of initial estimates.  A named vector
#'   may be used to provide estimate IDs.
#' @param std_error Optional numeric vector of positive standard errors.  Supply
#'   exactly one of std_error and p_value.
#' @param p_value Optional numeric vector of two-sided normal p-values.  When
#'   supplied instead of std_error, the function recovers s1 as
#'   abs(point_estimate) / qnorm(1 - p_value / 2).  Consequently, exact p-values
#'   strictly between zero and one and nonzero point estimates are required.
#'   Supplying standard errors is preferable when p-values have been rounded.
#' @param w Nonnegative scalar or vector equal to s2 / s1, the follow-up
#'   standard error divided by the initial standard error.  Defaults to 1.
#' @param top_n Optional positive integer giving the exact number of estimates
#'   selected by each ranking rule.  Cutoff ties are broken by estimate ID.
#' @param estimate_id Optional vector of unique estimate IDs.  If omitted, names
#'   on point_estimate are used; otherwise sequential IDs are generated.
#' @param histogram_bins Positive integer number of logarithmic histogram bins.
#' @param rel.tol Positive numerical tolerance for integration and root finding.
#'
#' @return An object of class a_value_result.  Its data element contains the
#'   inputs, A-values, ranks, and selection indicators.  Its plots element
#'   contains ggplot2 objects named histogram, p_value_ranks, and
#'   standard_error_ranks.
#'
#' @details The implementation uses the paper's unconstrained-parameter result,
#'   which requires the standardized parameter bound Delta_bar / s1 >= 0.76.
#'   Larger p-values, larger standard errors, and larger A-values all receive
#'   smaller numerical ranks; rank 1 is the largest value, as in Figure 4.
#'
#' @examples
#' estimates <- c(study_A = 0.20, study_B = -0.10, study_C = 0.55)
#' result <- calculate_a_values(
#'   point_estimate = estimates,
#'   std_error = c(0.10, 0.30, 0.20),
#'   w = 1,
#'   top_n = 1
#' )
#' result$data
#' result$plots$histogram
#' result$plots$p_value_ranks
#' result$plots$standard_error_ranks
#'
#' # If only exact two-sided normal p-values are available:
#' result_from_p <- calculate_a_values(
#'   point_estimate = c(a = 0.2, b = -0.3),
#'   p_value = c(0.10, 0.25)
#' )
#' @export
calculate_a_values <- function(
    point_estimate,
    std_error = NULL,
    p_value = NULL,
    w = 1,
    top_n = NULL,
    estimate_id = NULL,
    histogram_bins = 30L,
    rel.tol = 1e-9) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required to create the returned plots. ",
      "Install it with install.packages('ggplot2').",
      call. = FALSE
    )
  }
  if (!is.numeric(point_estimate) || length(point_estimate) == 0L) {
    stop("point_estimate must be a nonempty numeric vector.", call. = FALSE)
  }
  n <- length(point_estimate)

  if (is.null(estimate_id)) {
    estimate_names <- names(point_estimate)
    if (!is.null(estimate_names) && all(!is.na(estimate_names)) &&
        all(nzchar(estimate_names))) {
      estimate_id <- estimate_names
    } else {
      estimate_id <- seq_len(n)
    }
  }
  if (length(estimate_id) != n) {
    stop("estimate_id must have the same length as point_estimate.", call. = FALSE)
  }
  if (!is.atomic(estimate_id)) {
    stop("estimate_id must be an atomic vector.", call. = FALSE)
  }
  id_tie_key <- if (is.numeric(estimate_id)) {
    estimate_id
  } else {
    as.character(estimate_id)
  }
  estimate_id <- as.character(estimate_id)
  if (anyNA(estimate_id) || any(!nzchar(estimate_id))) {
    stop("estimate_id cannot contain missing or empty values.", call. = FALSE)
  }
  if (anyDuplicated(estimate_id)) {
    stop("estimate_id values must be unique.", call. = FALSE)
  }

  point_estimate <- as.numeric(point_estimate)
  bad_estimate <- !is.na(point_estimate) & !is.finite(point_estimate)
  if (any(bad_estimate)) {
    stop(
      "Every non-missing point_estimate must be finite; invalid ID(s): ",
      .av_bad_ids(estimate_id, which(bad_estimate)),
      call. = FALSE
    )
  }
  if (is.null(std_error) == is.null(p_value)) {
    stop("Supply exactly one of std_error and p_value.", call. = FALSE)
  }
  if (!is.numeric(w)) stop("w must be numeric.", call. = FALSE)
  w <- .av_recycle(as.numeric(w), n, "w")
  bad_w <- !is.na(w) & (!is.finite(w) | w < 0)
  if (any(bad_w)) {
    stop(
      "Every non-missing w must be finite and nonnegative; invalid ID(s): ",
      .av_bad_ids(estimate_id, which(bad_w)),
      call. = FALSE
    )
  }

  if (!is.numeric(rel.tol) || length(rel.tol) != 1L ||
      !is.finite(rel.tol) || rel.tol <= 0) {
    stop("rel.tol must be one finite positive number.", call. = FALSE)
  }
  if (!is.numeric(histogram_bins) || length(histogram_bins) != 1L ||
      !is.finite(histogram_bins) || histogram_bins < 1 ||
      histogram_bins != floor(histogram_bins)) {
    stop("histogram_bins must be one positive integer.", call. = FALSE)
  }
  histogram_bins <- as.integer(histogram_bins)

  if (!is.null(std_error)) {
    if (!is.numeric(std_error)) stop("std_error must be numeric.", call. = FALSE)
    std_error <- .av_recycle(as.numeric(std_error), n, "std_error")
    bad_se <- !is.na(std_error) & (!is.finite(std_error) | std_error <= 0)
    if (any(bad_se)) {
      stop(
        "Every non-missing std_error must be finite and positive; invalid ID(s): ",
        .av_bad_ids(estimate_id, which(bad_se)),
        call. = FALSE
      )
    }
    z_statistic <- point_estimate / std_error
    p_value <- 2 * pnorm(abs(z_statistic), lower.tail = FALSE)
    standard_error_source <- rep("supplied", n)
    p_value_source <- rep("computed from standard error", n)
  } else {
    if (!is.numeric(p_value)) stop("p_value must be numeric.", call. = FALSE)
    p_value <- .av_recycle(as.numeric(p_value), n, "p_value")
    bad_p <- !is.na(p_value) & (!is.finite(p_value) | p_value < 0 | p_value > 1)
    if (any(bad_p)) {
      stop(
        "Every non-missing p_value must be finite and between zero and one; ",
        "invalid ID(s): ", .av_bad_ids(estimate_id, which(bad_p)),
        call. = FALSE
      )
    }
    infer <- !is.na(point_estimate) & !is.na(p_value)
    cannot_infer <- infer &
      (point_estimate == 0 | p_value <= 0 | p_value >= 1)
    if (any(cannot_infer)) {
      stop(
        "A positive finite standard error cannot be recovered from these ",
        "point-estimate/p-value pairs; provide std_error for ID(s): ",
        .av_bad_ids(estimate_id, which(cannot_infer)),
        call. = FALSE
      )
    }
    z_abs <- rep(NA_real_, n)
    z_abs[infer] <- qnorm(p_value[infer] / 2, lower.tail = FALSE)
    std_error <- rep(NA_real_, n)
    std_error[infer] <- abs(point_estimate[infer]) / z_abs[infer]
    z_statistic <- sign(point_estimate) * z_abs
    standard_error_source <- ifelse(
      infer,
      "inferred from p-value",
      NA_character_
    )
    p_value_source <- ifelse(infer, "supplied", NA_character_)
  }

  a_value <- rep(NA_real_, n)
  normalized_a_value <- rep(NA_real_, n)
  complete <- is.finite(z_statistic) & is.finite(std_error) & is.finite(w)

  # Reuse solutions for repeated (absolute z, w) pairs.
  if (any(complete)) {
    complete_rows <- which(complete)
    keys <- paste(
      format(abs(z_statistic[complete]), digits = 17, scientific = TRUE),
      format(w[complete], digits = 17, scientific = TRUE),
      sep = "|"
    )
    unique_keys <- unique(keys)
    solved <- setNames(numeric(length(unique_keys)), unique_keys)
    for (key in unique_keys) {
      local_index <- which(keys == key)[1L]
      row_index <- complete_rows[local_index]
      solved[key] <- .normalized_a_value(
        z_statistic[row_index],
        w[row_index],
        rel.tol = rel.tol
      )
    }
    normalized_a_value[complete] <- unname(solved[keys])
    a_value[complete] <- std_error[complete] * normalized_a_value[complete]
  }

  eligible <- which(
    is.finite(a_value) & is.finite(p_value) & is.finite(std_error)
  )
  if (length(eligible) == 0L) {
    stop("At least one complete estimate is required.", call. = FALSE)
  }

  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1L || !is.finite(top_n) ||
        top_n < 1 || top_n != floor(top_n)) {
      stop("top_n must be one positive integer.", call. = FALSE)
    }
    if (top_n > length(eligible)) {
      stop(
        "top_n cannot exceed the number of complete estimates (",
        length(eligible), ").",
        call. = FALSE
      )
    }
    top_n <- as.integer(top_n)
  }

  a_value_rank <- .av_rank(a_value, id_tie_key, eligible)
  p_value_rank <- .av_rank(p_value, id_tie_key, eligible)
  standard_error_rank <- .av_rank(std_error, id_tie_key, eligible)

  selected_by_a_value <- rep(FALSE, n)
  selected_by_p_value <- rep(FALSE, n)
  selected_by_standard_error <- rep(FALSE, n)
  if (!is.null(top_n)) {
    selected_by_a_value <- !is.na(a_value_rank) & a_value_rank <= top_n
    selected_by_p_value <- !is.na(p_value_rank) & p_value_rank <= top_n
    selected_by_standard_error <-
      !is.na(standard_error_rank) & standard_error_rank <= top_n
  }

  result_data <- data.frame(
    estimate_id = estimate_id,
    point_estimate = point_estimate,
    standard_error = std_error,
    standard_error_source = standard_error_source,
    z_statistic = z_statistic,
    p_value = p_value,
    p_value_source = p_value_source,
    w = w,
    normalized_a_value = normalized_a_value,
    a_value = a_value,
    a_value_rank = a_value_rank,
    p_value_rank = p_value_rank,
    standard_error_rank = standard_error_rank,
    selected_by_a_value = selected_by_a_value,
    selected_by_p_value = selected_by_p_value,
    selected_by_standard_error = selected_by_standard_error,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result_data$p_value_selection_group <- .av_selection_group(
    selected_by_a_value,
    selected_by_p_value
  )
  result_data$standard_error_selection_group <- .av_selection_group(
    selected_by_a_value,
    selected_by_standard_error
  )

  plots <- list(
    histogram = .av_histogram(result_data, top_n, histogram_bins),
    p_value_ranks = .av_rank_scatter(result_data, "p_value", top_n),
    standard_error_ranks = .av_rank_scatter(
      result_data,
      "standard_error",
      top_n
    )
  )

  structure(
    list(
      data = result_data,
      plots = plots,
      settings = list(
        top_n = top_n,
        histogram_bins = histogram_bins,
        rel.tol = rel.tol,
        w_definition = "s2 / s1 (follow-up SE / initial SE)"
      )
    ),
    class = "a_value_result"
  )
}

#' @export
print.a_value_result <- function(x, ...) {
  print(x$data, row.names = FALSE, ...)
  cat(
    "\nPlot objects: $plots$histogram, $plots$p_value_ranks, ",
    "$plots$standard_error_ranks\n",
    sep = ""
  )
  invisible(x)
}

#' @export
plot.a_value_result <- function(
    x,
    y = NULL,
    which = c(
      "all",
      "histogram",
      "p_value_ranks",
      "standard_error_ranks"
    ),
    ...) {
  which <- match.arg(which)
  plot_names <- if (which == "all") names(x$plots) else which
  for (plot_name in plot_names) print(x$plots[[plot_name]])
  invisible(x)
}
