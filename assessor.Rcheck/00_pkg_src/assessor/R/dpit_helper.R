#' @rawNamespace S3method(dpit,default)
dpit.default <- function(model, ...) {
  cls <- paste(class(model), collapse = ", ")
  stop(sprintf("Unsupported model class for dpit(): %s", cls), call. = FALSE)
}

.dpit_finalize <- function(res,
                           plot,
                           scale,
                           line_args,
                           ...) {
  scale <- match.arg(scale, c("normal", "uniform"))
  eps <- .Machine$double.eps
  res <- pmin(pmax(res, eps), 1 - eps)
  if (isTRUE(plot)) {
    qqplot.resid(res, scale = scale, line_args = line_args, ...)
  }
  out <- if (scale == "normal") stats::qnorm(res) else res
  return(out)
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
