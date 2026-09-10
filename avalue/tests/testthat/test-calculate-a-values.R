test_that("paper benchmark values are reproduced", {
  equal_precision <- calculate_a_values(
    point_estimate = 1.96,
    std_error = 1,
    w = 1
  )
  perfect_follow_up <- calculate_a_values(
    point_estimate = 1.96,
    std_error = 1,
    w = 0
  )

  expect_equal(
    equal_precision$data$normalized_a_value,
    0.001499735,
    tolerance = 1e-7
  )
  expect_equal(
    perfect_follow_up$data$normalized_a_value,
    0.0691755,
    tolerance = 1e-7
  )
})

test_that("paper p-value comparison ratios are reproduced", {
  p_values <- c(0.05, 0.10, 0.15, 0.20, 0.25)
  z_values <- qnorm(p_values / 2, lower.tail = FALSE)
  result <- calculate_a_values(
    point_estimate = z_values,
    std_error = 1,
    w = 1
  )
  ratios <- result$data$a_value / result$data$a_value[1]

  expect_equal(
    ratios,
    c(1.000, 2.520, 4.380, 6.540, 8.980),
    tolerance = 0.004
  )
})

test_that("A-values scale with the standard error", {
  result <- calculate_a_values(
    point_estimate = c(a = 1.96, b = 3.92),
    std_error = c(1, 2),
    w = 1
  )

  expect_equal(result$data$normalized_a_value[1],
               result$data$normalized_a_value[2])
  expect_equal(result$data$a_value[2], 2 * result$data$a_value[1])
})

test_that("exact p-value input agrees with standard-error input", {
  estimates <- c(a = 0.2, b = -0.3)
  standard_errors <- c(0.1, 0.25)
  exact_p <- 2 * pnorm(abs(estimates / standard_errors), lower.tail = FALSE)

  from_se <- calculate_a_values(estimates, std_error = standard_errors)
  from_p <- calculate_a_values(estimates, p_value = exact_p)

  expect_equal(from_p$data$standard_error, standard_errors, tolerance = 1e-12)
  expect_equal(from_p$data$a_value, from_se$data$a_value, tolerance = 1e-10)
})

test_that("default and scalar w are recycled", {
  estimates <- c(a = 0.1, b = 0.2)
  standard_errors <- c(0.3, 0.4)

  default <- calculate_a_values(estimates, std_error = standard_errors)
  explicit <- calculate_a_values(
    estimates,
    std_error = standard_errors,
    w = c(1, 1)
  )

  expect_equal(default$data$w, c(1, 1))
  expect_equal(default$data$a_value, explicit$data$a_value)
})

test_that("selection is exact and cutoff ties use estimate IDs", {
  result <- calculate_a_values(
    point_estimate = c(b = 0.2, a = 0.2, c = 0.3),
    std_error = c(0.2, 0.2, 0.2),
    top_n = 1
  )

  expect_equal(sum(result$data$selected_by_a_value), 1)
  tied <- result$data$estimate_id %in% c("a", "b")
  expect_lt(
    result$data$a_value_rank[result$data$estimate_id == "a"],
    result$data$a_value_rank[result$data$estimate_id == "b"]
  )
  expect_equal(length(result$data$a_value_rank[tied]), 2)
})

test_that("missing rows are retained but excluded from rankings", {
  result <- calculate_a_values(
    point_estimate = c(a = 0.1, b = NA_real_, c = 0.3),
    std_error = c(0.2, 0.2, NA_real_)
  )

  expect_equal(nrow(result$data), 3)
  expect_true(is.finite(result$data$a_value[1]))
  expect_true(all(is.na(result$data$a_value[-1])))
  expect_true(all(is.na(result$data$a_value_rank[-1])))
})

test_that("input validation catches ambiguous or invalid inputs", {
  expect_error(
    calculate_a_values(0.2, std_error = 0.1, p_value = 0.05),
    "exactly one"
  )
  expect_error(calculate_a_values(0.2), "exactly one")
  expect_error(calculate_a_values(0.2, std_error = 0), "positive")
  expect_error(calculate_a_values(0.2, std_error = 0.1, w = -1),
               "nonnegative")
  expect_error(
    calculate_a_values(c(a = 0.2, a = 0.3), std_error = c(0.1, 0.2)),
    "unique"
  )
  expect_error(
    calculate_a_values(0, p_value = 0.5),
    "provide std_error"
  )
})

test_that("result and all plots have their documented classes", {
  result <- calculate_a_values(
    point_estimate = c(a = 0.1, b = 0.3, c = -0.2),
    std_error = c(0.2, 0.2, 0.4),
    top_n = 1
  )

  expect_s3_class(result, "a_value_result")
  expect_named(result$plots,
               c("histogram", "p_value_ranks", "standard_error_ranks"))
  lapply(result$plots, expect_s3_class, class = "ggplot")
  lapply(result$plots, function(plot) {
    expect_no_warning(ggplot2::ggplot_build(plot))
  })
  expect_invisible(print(result))
})
