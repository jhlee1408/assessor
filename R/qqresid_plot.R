#' @keywords internal
qqplot.resid <- function(empcdf, scale, line_args, ...){
    n <- length(empcdf)
    if (scale == "normal") {
      empcdf <- qnorm(empcdf)
      x <- qnorm(ppoints(n))
      y <- empcdf[is.finite(empcdf)]
    } else {
      x <- ppoints(n)
      y <- empcdf[is.finite(empcdf)]
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
    do.call(abline, ab_args)
}
