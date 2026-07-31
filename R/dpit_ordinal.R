#' Residuals for regression models with ordinal outcomes
#'
#' Computes DPIT residuals for regression models with ordinal outcomes
#' using observed outcomes (`y`), ordinal outcome levels (`level`) and their fitted category
#' probabilities (`fitprob`).
#'
#'
#' @usage dpit_ordi(y, level, fitprob)
#' @param y An observed ordinal outcome vector.
#' @param level The response levels in their ordinal order. For instance,
#'   `c(0, 1, 2)` or `c("low", "medium", "high")`.
#' @param fitprob A matrix of fitted category probabilities. Each row
#'   corresponds to an observation, and the columns must follow the order in
#'   `level`. Each row must sum to one.
#' @returns A `dpit` object containing DPIT residuals.
#'
#' @details
#' For formulation details on discrete outcomes, see \code{\link{dpit}}.
#'
#' @examples
#' ## Ordinal example
#' library(MASS)
#' n <- 500
#' x1 <- rnorm(n, mean = 2)
#' beta1 <- 3
#' # True model
#' p0 <- plogis(1, location = beta1 * x1)
#' p1 <- plogis(4, location = beta1 * x1) - p0
#' p2 <- 1 - p0 - p1
#' genemult <- function(p) {
#'  rmultinom(1, size = 1, prob = c(p[1], p[2], p[3]))
#' }
#' test <- apply(cbind(p0, p1, p2), 1, genemult)
#' y1 <- rep(0, n)
#' y1[which(test[1, ] == 1)] <- 0
#' y1[which(test[2, ] == 1)] <- 1
#' y1[which(test[3, ] == 1)] <- 2
#' multimodel <- polr(as.factor(y1) ~ x1, method = "logistic")
#'
#' y1 <- multimodel$model[,1]
#' lev1 <- multimodel$lev
#' fitprob1 <- fitted(multimodel)
#'
#' dpit.ord <- dpit_ordi(y=y1, level=lev1, fitprob=fitprob1)
#' resid.ord <- residuals(dpit.ord)
#' plot(dpit.ord)
#' @export
dpit_ordi <- function(y, level, fitprob) {
  fitprob <- as.matrix(fitprob)
  k <- length(level)

  if (k < 2L) {
    stop("level must contain at least two ordinal levels.", call. = FALSE)
  }
  if (anyDuplicated(level)) {
    stop("level must not contain duplicates.", call. = FALSE)
  }
  if (length(y) != nrow(fitprob)) {
    stop("nrow(fitprob) must equal length(y).", call. = FALSE)
  }
  if (k != ncol(fitprob)) {
    stop("length(level) must equal ncol(fitprob).", call. = FALSE)
  }
  if (any(!is.finite(fitprob)) || any(fitprob < 0) || any(fitprob > 1)) {
    stop("fitprob must contain finite probabilities between 0 and 1.", call. = FALSE)
  }
  if (any(abs(rowSums(fitprob) - 1) > 1e-6)) {
    stop("Each row of fitprob must sum to 1.", call. = FALSE)
  }

  out <- match(as.character(y), as.character(level))
  if (anyNA(out)) {
    stop("All values of y must occur in level.", call. = FALSE)
  }

  n <- length(out)
  q <- t(apply(fitprob, 1, cumsum))
  inde <- cbind(seq_len(n), out)
  res <- q[inde]

  empcdf <- rep(NA, n)
  for(i in seq_len(n)){
    if(i %in% which(out==k)) next
    note <- matrix(NA, ncol=k, nrow=n)
    for(p in seq_len(k)){
        note[,p] <- fitprob[, p] * (res[i] > q[, p])
    }
    note.sum <- rowSums(note)
    note.sum[i] <- 0
    empcdf[i] <- sum(note.sum)/(n-1)
  }

  # for loop with max values
  ses <- ifelse(out == k, q[, 1], 0)
  for(i in seq_len(n)){
    if(i %in% which(out != k)) next
    pses <- (ses[i] < q[,1])*q[,k-1]
    pses[pses==0] <- 1
    pses[i] <- 0
    empcdf[i] <- sum(pses)/(n-1)
  }
  .new_dpit(empcdf, method = "Ordinal")
}
