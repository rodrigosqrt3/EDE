#' @keywords internal
#' @noRd
new_ede_estimate <- function(estimate, lower = NA_real_, upper = NA_real_,
                              method, alpha) {
  structure(
    list(
      estimate = estimate,
      lower = lower,
      upper = upper,
      method = method,
      alpha = alpha
    ),
    class = "ede_estimate"
  )
}

#' @export
print.ede_estimate <- function(x, ...) {
  cat(sprintf("<%s>\n", x$method))
  cat(sprintf("  estimate: %s\n", format(round(x$estimate, 2))))
  if (!is.na(x$lower) || !is.na(x$upper)) {
    cat(sprintf(
      "  %.0f%% CI: [%s, %s]\n",
      100 * (1 - x$alpha), format(round(x$lower, 2)), format(round(x$upper, 2))
    ))
  }
  invisible(x)
}
