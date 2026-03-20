#' @rawNamespace S3method(dpit,glm)
dpit.glm <- function(model,
                     plot = TRUE,
                     scale = "normal",
                     line_args = list(),
                     ...) {
  .dpit_glm(model, ...)
}

.dpit_glm <- function(model, ...) {
  fam_obj <- stats::family(model)
  fam_raw  <- fam_obj$family
  link_raw <- fam_obj$link

  fam_key <- fam_raw
  if (startsWith(fam_raw, "Negative Binomial")) {
    fam_key <- "negbin"
  }

  key <- list(
    model  = model,
    family = fam_raw,
    link   = link_raw
  )

  class(key) <- fam_key
  .dpit_glm_key(key, ...)
}



.dpit_glm_key <- function(key, ...) {
  UseMethod(".dpit_glm_key")
}

#' @rawNamespace S3method(.dpit_glm_key,poisson)
.dpit_glm_key.poisson <- function(key, ...) {
  model <- key$model
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted.values(model)

  if (!exists("dpit_pois", mode = "function", inherits = TRUE)) {
    stop("dpit_pois() not found. Please define it in the package namespace.", call. = FALSE)
  }
  dpit_pois(y, mu)
}

#' @rawNamespace S3method(.dpit_glm_key,binomial)
.dpit_glm_key.binomial <- function(key, ...) {
  model <- key$model
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted.values(model)
  if (!exists("dpit_bin", mode = "function", inherits = TRUE)) {
    stop("dpit_bin() not found. Please define it in the package namespace.", call. = FALSE)
  }
  dpit_bin(y = y, prob =mu)
}


#' @rawNamespace S3method(.dpit_glm_key,negbin)
.dpit_glm_key.negbin <- function(key, ...) {
  model <- key$model
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted.values(model)

  size <- if (!is.null(model$theta)) model$theta else summary(model)$theta
  if (!exists("dpit_nb", mode = "function", inherits = TRUE)) {
    stop("dpit_nb() not found. Please define it in the package namespace.", call. = FALSE)
  }
  dpit_nb(y, mu, size)
}


#' @rawNamespace S3method(.dpit_glm_key,Tweedie)
.dpit_glm_key.Tweedie <- function(key, ...) {
  model <- key$model
  y1 <- model$y
  p.max <- get("p", envir = environment(model$family$variance))
  lambda1f <- model$fitted.values
  phi1f <- summary(model)$dis
  dpit_tweedie(y=y1, mu=lambda1f, xi= p.max, phi=phi1f)
}

#' @rawNamespace S3method(.dpit_glm_key,default)
.dpit_glm_key.default <- function(key, ...) {
  fam  <- key$family
  link <- key$link
  stop(sprintf("Unsupported GLM family/link for dpit(): %s / %s", fam, link),
       call. = FALSE)
}

#' @rawNamespace S3method(dpit,polr)
dpit.polr <- function(model,
                      plot = TRUE,
                      scale = "normal",
                      line_args = list(),
                      ...) {
  if (!exists("dpit_ordi", mode = "function", inherits = TRUE)) {
    stop("dpit_ordi() not found. Please define it in the package namespace.", call. = FALSE)
  }
  y <- model$model[,1]
  lev <- model$lev
  p1f <- fitted(model)
  dpit_ordi(y=y, level=lev, fitprob=p1f)
}





