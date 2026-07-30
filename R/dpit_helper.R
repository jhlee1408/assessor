#' @rawNamespace S3method(dpit,default)
dpit.default <- function(model) {
  cls <- paste(class(model), collapse = ", ")
  stop(sprintf("Unsupported model class for dpit(): %s", cls), call. = FALSE)
}

.new_dpit <- function(res, method = NULL) {
  eps <- .Machine$double.eps
  res <- pmin(pmax(res, eps), 1 - eps)

  structure(
    list(
      uniform = as.numeric(res),
      method = method
    ),
    class = "dpit"
  )
}

.attach_model_call <- function(x, model) {
  model_call <- tryCatch(
    stats::getCall(model),
    error = function(e) NULL
  )

  if (is.null(model_call) && isS4(model)) {
    model_call <- tryCatch(
      model@call,
      error = function(e) NULL
    )
  }

  x$call <- model_call
  x
}

#' Methods for DPIT residual objects
#'
#' Methods for printing, summarizing, plotting, and extracting values from
#' objects returned by `dpit()`, `dpit_2pm()`, and the distribution-specific
#' DPIT calculators.
#' The print method displays the fitted model call or calls, when available,
#' and the sample size.
#' The summary method reports quantiles, mean, and standard deviation for the
#' selected scale.
#'
#' @param x A `dpit` object.
#' @param object A `dpit` object.
#' @param scale You can choose the scale of the residuals among `normal` and `uniform`.
#' The sample quantiles of the residuals are plotted against
#' the theoretical quantiles of a standard normal distribution under the normal scale,
#' and against the theoretical quantiles of a uniform (0,1) distribution under the uniform scale.
#' The default scale is `normal`.
#' @param line_args A named list of graphical parameters passed to
#'   `graphics::abline()` to modify the reference (red) 45° line
#'   in the QQ plot. If left empty, a default red dashed line is drawn.
#' @param ... Additional arguments passed to or from methods. For `plot.dpit()`,
#'   graphical arguments are passed to `stats::qqplot()` for customizing the QQ
#'   plot (e.g., `pch`, `col`, `cex`, `xlab`, `ylab`).
#'
#' @returns `residuals.dpit()` returns a numeric vector of DPIT residuals.
#' `summary.dpit()` returns an object of class `summary.dpit`.
#' The print and plot methods return their input invisibly.
#'
#' @name dpit-methods
NULL

.print_dpit_calls <- function(calls) {
  if (is.null(calls)) {
    return(invisible(NULL))
  }

  if (is.list(calls) && !is.call(calls)) {
    cat("Model calls:\n")
    call_names <- names(calls)
    for (i in seq_along(calls)) {
      if (!is.null(call_names) && nzchar(call_names[[i]])) {
        cat(sprintf("%s:\n", call_names[[i]]))
      }
      print(calls[[i]])
    }
  } else {
    cat("Model call:\n")
    print(calls)
  }
  cat("\n")
  invisible(NULL)
}

#' @rdname dpit-methods
#' @export
print.dpit <- function(x, ...) {
  cat("DPIT residual object\n\n")
  .print_dpit_calls(x$call)
  cat(sprintf("Sample size: %d\n", length(x$uniform)))
  cat("Use residuals() to extract the residual values.\n")
  invisible(x)
}

#' @rdname dpit-methods
#' @export
residuals.dpit <- function(object,
                           scale = c("normal", "uniform"),
                           ...) {
  scale <- match.arg(scale)
  uniform_values <- object$uniform

  if (is.null(uniform_values) || !is.numeric(uniform_values)) {
    stop(
      "The dpit object does not contain numeric uniform-scale values.",
      call. = FALSE
    )
  }

  if (any(uniform_values < 0 | uniform_values > 1, na.rm = TRUE)) {
    stop(
      "Uniform-scale DPIT values must lie between 0 and 1.",
      call. = FALSE
    )
  }

  switch(
    scale,
    normal = stats::qnorm(uniform_values),
    uniform = uniform_values
  )
}

#' @rdname dpit-methods
#' @export
summary.dpit <- function(object,
                         scale = c("normal", "uniform"),
                         ...) {
  scale <- match.arg(scale)
  values <- residuals(object, scale = scale)
  keep <- is.finite(values)
  finite_values <- values[keep]

  statistics <- if (length(finite_values) == 0L) {
    stats <- rep(NA_real_, 7L)
    names(stats) <- c(
      "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "Std. Dev."
    )
    stats
  } else {
    c(
      "Min." = min(finite_values),
      "1st Qu." = stats::quantile(finite_values, 0.25, names = FALSE),
      "Median" = stats::median(finite_values),
      "Mean" = mean(finite_values),
      "3rd Qu." = stats::quantile(finite_values, 0.75, names = FALSE),
      "Max." = max(finite_values),
      "Std. Dev." = stats::sd(finite_values)
    )
  }

  structure(
    list(
      call = object$call,
      scale = scale,
      n = length(values),
      statistics = statistics
    ),
    class = "summary.dpit"
  )
}

#' @rdname dpit-methods
#' @export
print.summary.dpit <- function(x, ...) {
  cat("Summary of DPIT residuals\n\n")
  .print_dpit_calls(x$call)
  cat(sprintf("Residual scale: %s\n", x$scale))
  cat(sprintf("Sample size: %d\n\n", x$n))
  print(x$statistics, ...)
  invisible(x)
}

#' @rdname dpit-methods
#' @export
plot.dpit <- function(x,
                      scale = c("normal", "uniform"),
                      line_args = list(),
                      ...) {
  scale <- match.arg(scale)
  qqplot.resid(x$uniform, scale = scale, line_args = line_args, ...)
  invisible(x)
}

#' @keywords internal
qqplot.resid <- function(empcdf, scale, line_args, ...) {
  keep <- is.finite(empcdf)
  u <- empcdf[keep]
  n <- length(u)

  if (n == 0L) {
    stop("No finite residuals available for QQ plot.", call. = FALSE)
  }

  if (scale == "normal") {
    y <- stats::qnorm(u)
    x <- stats::qnorm(stats::ppoints(n))
  } else {
    y <- u
    x <- stats::ppoints(n)
  }

  qq_defaults <- list(
    main = "QQ plot",
    xlab = "Theoretical Quantiles",
    ylab = "Sample Quantiles",
    cex.lab = 1, cex.axis = 1, cex.main = 1.5, lwd = 1.5
  )
  qq_args <- utils::modifyList(qq_defaults, list(...))
  do.call(stats::qqplot, c(list(x = x, y = y), qq_args))

  abline_defaults <- list(a = 0, b = 1, col = "red", lty = 5, lwd = 1.5)
  ab_args <- utils::modifyList(abline_defaults, line_args)
  do.call(graphics::abline, ab_args)
}
