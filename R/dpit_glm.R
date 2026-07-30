#' @rawNamespace S3method(dpit,glm)
dpit.glm <- function(model) {
  .attach_model_call(.dpit_glm(model), model)
}



.dpit_glm <- function(model) {
  fam_obj <- stats::family(model)
  fam_raw <- fam_obj$family
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

  .dpit_glm_key(key)
}


.dpit_glm_key <- function(key) {
  UseMethod(".dpit_glm_key")
}

#' @rawNamespace S3method(.dpit_glm_key,poisson)
.dpit_glm_key.poisson <- function(key) {
  model <- key$model
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted.values(model)

  dpit_pois(
    y = y,
    mu = mu
  )
}

#' @rawNamespace S3method(.dpit_glm_key,binomial)
.dpit_glm_key.binomial <- function(key) {
  model <- key$model
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted.values(model)


  dpit_bin(
    y = y,
    prob = mu
  )
}

#' @rawNamespace S3method(.dpit_glm_key,negbin)
.dpit_glm_key.negbin <- function(key) {
  model <- key$model
  y <- stats::model.response(stats::model.frame(model))
  mu <- stats::fitted.values(model)

  size <- if (!is.null(model$theta)) model$theta else summary(model)$theta

  dpit_nb(
    y = y,
    mu = mu,
    size = size
  )
}

#' @rawNamespace S3method(.dpit_glm_key,Tweedie)
.dpit_glm_key.Tweedie <- function(key) {
  model <- key$model
  y1 <- model$y
  p.max <- get("p", envir = environment(model$family$variance))
  lambda1f <- model$fitted.values
  phi1f <- summary(model)$dis

  dpit_tweedie(
    y = y1,
    mu = lambda1f,
    xi = p.max,
    phi = phi1f
  )
}

#' @rawNamespace S3method(.dpit_glm_key,default)
.dpit_glm_key.default <- function(key) {
  fam  <- key$family
  link <- key$link
  stop(sprintf("Unsupported GLM family/link for dpit(): %s / %s", fam, link),
       call. = FALSE)
}

#' @rawNamespace S3method(dpit,polr)
dpit.polr <- function(model) {
  y <- model$model[,1]
  lev <- model$lev
  p1f <- fitted(model)
  out <- dpit_ordi(y = y, level = lev, fitprob = p1f)
  .attach_model_call(out, model)
}



