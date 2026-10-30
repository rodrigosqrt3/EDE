#' Apply every applicable method to one sighting record
#'
#' Runs the estimators and persistence tests that need only a sighting record
#' and returns their results side by side. Methods that need additional input
#' (such as sighting reliabilities) are not included. A method that cannot be
#' applied to the record, or that warns, is kept in the table with the message
#' in `note`.
#'
#' The methods rest on different models, so their results answer different
#' questions and should not be averaged; see the vignette for guidance on
#' choosing a method.
#'
#' @param sd A [sighting_data] object.
#' @param test_year Latest time to test for the persistence tests. Must be
#'   later than the last sighting.
#' @param alpha Significance level, in (0, 1).
#'
#' @return A data frame with one row per method and columns `method`,
#'   `function_name`, `estimate`, `lower`, `upper`, `interval_type` and
#'   `note`. Persistence tests report the one-sided upper bound in `upper`.
#'
#' @examples
#' ede_compare(sighting_data(monk_seal), test_year = 2050)
#'
#' @export
ede_compare <- function(sd, test_year, alpha = 0.05) {
  stopifnot(inherits(sd, "sighting_data"))
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).", call. = FALSE)
  if (any(sd$count > 0)) {
    check_test_year(test_year, max(sd$time[sd$count > 0]))
  }

  calls <- list(
    ole = function() ole(sd, alpha = alpha),
    robson1964 = function() robson1964(sd, alpha = alpha),
    strauss1989 = function() strauss1989(sd, alpha = alpha),
    solow1993 = function() solow1993(sd, alpha = alpha, test_year = test_year),
    solow1993b = function() solow1993b(sd, alpha = alpha, test_year = test_year),
    solow2005 = function() solow2005(sd, alpha = alpha, test_year = test_year),
    mcinerny2006 = function() mcinerny2006(sd, alpha = alpha, test_year = test_year),
    burgman1995 = function() burgman1995(sd, alpha = alpha, test_year = test_year),
    solow_roberts2003 = function() solow_roberts2003(sd, alpha = alpha, test_year = test_year),
    jaric2010 = function() jaric2010(sd, alpha = alpha, test_year = test_year)
  )

  rows <- lapply(names(calls), function(name) {
    notes <- character(0)
    result <- tryCatch(
      withCallingHandlers(
        calls[[name]](),
        warning = function(w) {
          notes <<- c(notes, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        notes <<- c(notes, conditionMessage(e))
        NULL
      }
    )
    if (is.null(result)) {
      return(data.frame(method = name, function_name = name, estimate = NA_real_,
                        lower = NA_real_, upper = NA_real_,
                        interval_type = NA_character_,
                        note = paste(notes, collapse = " ")))
    }
    data.frame(method = result$method, function_name = name,
               estimate = result$estimate, lower = result$lower,
               upper = result$upper, interval_type = result$interval_type,
               note = paste(notes, collapse = " "))
  })
  do.call(rbind, rows)
}
