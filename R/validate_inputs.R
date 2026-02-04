## Input validation helpers for ClockSHAP

.validate_named_numeric <- function(x, name) {
  if (!is.numeric(x)) {
    stop(sprintf("`%s` must be numeric.", name), call. = FALSE)
  }
  if (length(x) < 1) {
    stop(sprintf("`%s` must have length >= 1.", name), call. = FALSE)
  }
  nm <- names(x)
  if (is.null(nm)) {
    stop(sprintf("`%s` must be a named vector.", name), call. = FALSE)
  }
  if (anyNA(nm) || any(nm == "")) {
    stop(sprintf("`%s` must have non-empty, non-missing names.", name), call. = FALSE)
  }
  if (anyDuplicated(nm)) {
    stop(sprintf("`%s` must not contain duplicated names.", name), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }
  invisible(TRUE)
}

.validate_same_names <- function(...) {
  objs <- list(...)
  if (length(objs) <= 1) return(invisible(TRUE))

  name_lists <- lapply(objs, names)
  ref <- name_lists[[1]]

  if (is.null(ref)) {
    stop("All inputs must be named.", call. = FALSE)
  }

  ref_set <- sort(ref)

  for (i in seq_along(name_lists)) {
    nm <- name_lists[[i]]
    if (is.null(nm)) {
      stop("All inputs must be named.", call. = FALSE)
    }

    if (!identical(sort(nm), ref_set)) {
      missing <- setdiff(ref, nm)
      extra   <- setdiff(nm, ref)

      msg <- "All inputs must have identical names."
      if (length(missing) > 0) {
        msg <- paste0(msg, " Missing: ", paste(missing, collapse = ", "), ".")
      }
      if (length(extra) > 0) {
        msg <- paste0(msg, " Extra: ", paste(extra, collapse = ", "), ".")
      }
      stop(msg, call. = FALSE)
    }
  }

  invisible(TRUE)
}

.reorder_named <- function(x, target_names, name = "x") {
  if (is.null(names(x))) {
    stop(sprintf("`%s` must be named.", name), call. = FALSE)
  }
  missing <- setdiff(target_names, names(x))
  extra   <- setdiff(names(x), target_names)

  if (length(missing) > 0 || length(extra) > 0) {
    msg <- sprintf("`%s` must have exactly the target names.", name)
    if (length(missing) > 0) {
      msg <- paste0(msg, " Missing: ", paste(missing, collapse = ", "), ".")
    }
    if (length(extra) > 0) {
      msg <- paste0(msg, " Extra: ", paste(extra, collapse = ", "), ".")
    }
    stop(msg, call. = FALSE)
  }

  x[target_names]
}
