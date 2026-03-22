#' @rawNamespace S3method(dpit,zeroinfl)
dpit.zeroinfl <- function(model,
                          plot = TRUE,
                          scale = "normal",
                          line_args = list(),
                          ...) {
  res_u <- .dpit_zeroinfl(model, ...)
  .dpit_finalize(res_u, plot = plot, scale = scale, line_args = line_args, ...)
}

.dpit_zeroinfl <- function(model, ...) {
  key <- list(model = model)
  dist <- NULL
  if (!is.null(model$dist)) {
    dist <- model$dist
  } else if (!is.null(model$call$dist)) {
    dist <- eval(model$call$dist, parent.frame())
  } else if (!is.null(model$call$distribution)) {
    dist <- eval(model$call$distribution, parent.frame())
  } else if (!is.null(model$call$family)) {
    dist <- eval(model$call$family, parent.frame())
  }

  dist_chr <- if (length(dist) == 1L && is.character(dist)) tolower(dist) else NA_character_

  class(key) <- if (!is.na(dist_chr) && identical(dist_chr, "poisson")) {
    "zeroinfl_poisson"
  } else if (!is.na(dist_chr) && dist_chr %in% c("negbin", "negbinom", "negative binomial", "nbinom", "nb")) {
    "zeroinfl_negbin"
  } else {
    "zeroinfl_default"
  }

  .dpit_zeroinfl_key(key, ...)
}

.dpit_zeroinfl_key <- function(key, ...) {
  UseMethod(".dpit_zeroinfl_key")
}

#' @rawNamespace S3method(.dpit_zeroinfl_key,zeroinfl_poisson)
.dpit_zeroinfl_key.zeroinfl_poisson <- function(key, ...) {
  model <- key$model

  # response y
  mf <- try(stats::model.frame(model), silent = TRUE)
  if (inherits(mf, "try-error")) {
    # fallback: older objects may store model matrix in model$model
    if (!is.null(model$model)) {
      y <- model$model[, 1]
    } else {
      stop("Cannot extract response y from zeroinfl model.", call. = FALSE)
    }
  } else {
    y <- stats::model.response(mf)
  }

  mu    <- stats::predict(model, type = "count")
  pzero <- stats::predict(model, type = "zero")

  dpit_zpois(y = y, pzero = pzero, mu = mu)
}

#' @rawNamespace S3method(.dpit_zeroinfl_key,zeroinfl_negbin)
.dpit_zeroinfl_key.zeroinfl_negbin <- function(key, ...) {
  model <- key$model
  mf <- try(stats::model.frame(model), silent = TRUE)
  if (inherits(mf, "try-error")) {
    if (!is.null(model$model)) {
      y <- model$model[, 1]
    } else {
      stop("Cannot extract response y from zeroinfl model.", call. = FALSE)
    }
  } else {
    y <- stats::model.response(mf)
  }

  mu    <- stats::predict(model, type = "count")
  pzero <- stats::predict(model, type = "zero")
  size <- model$theta
  if (is.null(size)) {
    stop("Cannot extract theta/size from zeroinfl(negbin) model (model$theta is NULL).", call. = FALSE)
  }
  dpit_znb(y = y, pzero = pzero, mu = mu, size = size)
}

#' @rawNamespace S3method(.dpit_zeroinfl_key,zeroinfl_default)
.dpit_zeroinfl_key.zeroinfl_default <- function(key, ...) {
  model <- key$model
  dist <- if (!is.null(model$dist)) model$dist else NULL
  msg <- if (is.null(dist)) {
    "Unsupported zeroinfl distribution for dpit() (expected poisson or negbin)."
  } else {
    sprintf("Unsupported zeroinfl distribution for dpit(): %s (expected poisson or negbin).", as.character(dist))
  }
  stop(msg, call. = FALSE)
}
