## Input validation helpers for ClockSHAP

.validate_named_numeric <- function(x, name) {
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }
  if (is.null(names(x))) {
    stop(sprintf("`%s` must be a named vector.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_same_names <- function(...) {
  objs <- list(...)
  name_lists <- lapply(objs, names)

  ref <- name_lists[[1]]
  for (i in seq_along(name_lists)) {
    if (!identical(ref, name_lists[[i]])) {
      stop("All inputs must have identical names.", call. = FALSE)
    }
  }
  invisible(TRUE)
}
