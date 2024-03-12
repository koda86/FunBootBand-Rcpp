#' @title FunBootBand
#'
#' @author Daniel Koska
#'
#' @description Creates Functional Bootstraped (statistical) Bands.
#' The 'band' function rests on two basic ASSUMPTIONS:
#'
#' - Assumption 1 (A1): All curves are of equal length. If necessary, any interpolation must be
#' performed externally.
#'
#' - Assumption 2 (A2): Curves originate from stationary processes.
#'
#' If these assumptions are not met, the function will return an error (A1) or
#' result in potentially erroneous outputs (A2).
#'
#' @usage band(data, type, alpha, iid = TRUE, k.coef = 50,
#' B = 400)
#'
#' @param data A data frame of dimensions [t, n], where 'n' is the number of
#' curves and 't' is the length of the curves. All curves need to be of equal
#' length.
#' @param type Type of the statistical band (c("confidence", "prediction")).
#' @param alpha Desired type I error probability.
#' @param iid Dependency of the curves (iid = c(TRUE, FALSE)). Setting iid=TRUE
#' runs an ordinary (naive) bootstrap, where all curves are assumed to be
#' independent. When setting iid=FALSE, a two-stage bootstrap is run, where
#' curve clusters (comprising all of their curves) are resampled with
#' replacement in the initial stage, and one curve per cluster is sampled
#' without replacement in the second stage. If iid is set to FALSE, curves are
#' assumed to be nested in curve clusters. The cluster structure needs to be
#' specified in the colnames of the data frame using letters to indicate
#' clusters (see 'Format').
#' @param k.coef Number of Fourier coefficients (e.g., k = 50). Determines the
#' smoothness of the curve approximation.
#' @param B Number of bootstrap iterations (e.g., B = 1000). Default is 400.
#'
#' @export
#'
#' @return A data frame object that contains upper and lower band boundaries,
#' alongside a mean curve.
#'
#' @examples
#' library(FunBootBand)
#' band.limits <- band(data, type = "prediction", alpha = 0.05, iid = TRUE, B = 50)
#' plot(band.limits[1,], # upper band limit
#'      xlim = c(0, 101),
#'      ylim = range(band.limits),
#'      type = "l")
#' lines(band.limits[2, ]) # mean curve
#' lines(band.limits[3, ]) # lower band limit

band <- function(data, type, alpha, iid = TRUE, k.coef = 50, B = 400) {

  # Argument checking
  if (!inherits(type, "character")) {
    stop("'type' must be a variable of type 'character'.")
  } else if (!(type %in% c("confidence", "prediction"))) {
    stop("'type' must be either 'confidence' or 'prediction'.")
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a numeric value between 0 and 1.")
  }
  if (!is.logical(iid)) {
    stop("'iid' must be a logical value (TRUE or FALSE).")
  }
  if (!is.numeric(k.coef) || k.coef <= 0) {
    stop("'k.coef' must be a positive integer.")
  }
  if (!is.numeric(B) || B <= 0) {
    stop("'B' must be a positive integer.")
  }
  if (any(is.na(data))) {
    stop("Function stopped due to NA's in the input data.")
  }

  if (is.data.frame(data) == FALSE) {stop("Input data is not a data frame.")}

  # Check if all data elements are numeric
  if (!all(sapply(data, is.numeric))) {
    stop("Non-numeric data found in input.")
  }

  n.curves  <- dim(data)[2]
  n.time <- dim(data)[1]
  time <- seq(0, (n.time - 1))

  if (iid != TRUE) {
    # Check colnames to make sure the nested structure is correctly identified
    if (colnames(data)[1] != colnames(data)[2]) {
      # This is necessary to make sure curves are not detected as iid
      new_colnames <- substr(colnames(data), 1, 1)
      colnames(data) <- new_colnames
    }

    if (colnames(data)[1] != colnames(data)[2]) {
      stop("Header does not indicate a nested structure even though 'iid' is set to 'FALSE'.")
      # \n Make sure curves within the same cluster all have the exact same column label.")
    }

    # --------------------------------------------------------------------------
    # This implements a more robust approach that allows the detection of the
    # number and size (aka number of curves per cluster) of cluster in the case
    # where different clusters contain a different amount of curves per cluster.
    # Function to determine the base name (or cluster identifier) of a column
    # --------------------------------------------------------------------------
    get_base_name <- function(name) {
      # Split the name at any non-alphanumeric character (e.g., '.', '-')
      parts <- strsplit(name, "[^[:alnum:]]")[[1]]
      return(parts[1])
    }

    clusters <- list()
    for (name in colnames(data)) {
      base_name <- get_base_name(name)
      # Check if the base name already exists in the clusters
      if (!base_name %in% names(clusters)) {
        # If not, create a new entry in clusters with this base name
        clusters[[base_name]] <- 0
      }
      clusters[[base_name]] <- clusters[[base_name]] + 1
    }

    cluster_boundaries <- list()
    start_idx <- 1
    for (cluster_id in names(clusters)) {
      end_idx <- start_idx + clusters[[cluster_id]] - 1
      cluster_boundaries[[cluster_id]] <- c(start = start_idx, end = end_idx)
      start_idx <- end_idx + 1
    }

    # Function to get the indices for a single cluster
    get_cluster_indices <- function(cluster_id, cluster_boundaries) {
      if (cluster_id %in% names(cluster_boundaries)) {
        boundaries <- cluster_boundaries[[cluster_id]]
        return(seq(from = boundaries["start"], to = boundaries["end"]))
      } else {
        stop("Cluster ID not found.")
      }
    }

    n.cluster <- length(clusters)
    if (n.cluster < 2 | n.cluster == ncol(data)) {
      stop("Header does not indicate a nested structure even though 'iid' is set to 'FALSE'.")
      }
  }

  # Approximate curves using Fourier functions ---------------------------------
  fourier.koeffi  <- array(data = 0, dim = c(k.coef*2 + 1, n.curves))
  fourier.real    <- array(data = 0, dim = c(n.time, n.curves))
  fourier.mean    <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1))
  fourier.real_mw <- array(data = 0, c(n.time, 1))
  fourier.std1    <- array(data = 0, c(k.coef*2 + 1, k.coef*2 + 1, n.curves))
  fourier.cov     <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1))
  fourier.std_all <- array(data = 0, c(n.time, n.time))
  fourier.std     <- array(data = 0, dim = c(n.time, 1))

  # Construct Fourier series
  # General: f(t) = mu + sum(alpha cos(2pi*k*t/T) + beta sin(2pi*k*t/T))
  fourier.s = rep(1, times = n.time)
  for (k in seq(1, k.coef*2, 2)) {
    fourier.s <- cbind(fourier.s, cos(2*pi*(k/2)*time / (n.time-1)))
    fourier.s <- cbind(fourier.s, sin(2*pi*(k/2)*time / (n.time-1)))
    # '-1' to match the equations in Lenhoff Appendix A ('T')
  }

  # Helper function to calculate the pseudoinverse (Moore-Penrose)
  pseudo_inverse <- function(A, tol = .Machine$double.eps^(2/3)) {
    stopifnot(is.numeric(A) || is.complex(A), is.matrix(A))

    svd_result <- svd(A)
    U <- svd_result$u
    S <- svd_result$d
    V <- svd_result$v
    if (is.complex(A)) {U <- Conj(U)} # Adjust for complex matrices
    # Calculate the pseudoinverse
    threshold <- max(tol * S[1], 0)
    non_zero_indices <- S > threshold

    if (all(non_zero_indices)) {
      inverse <- V %*% diag(1/S) %*% t(U)
    } else if (any(non_zero_indices)) {
      V_filtered <- V[, non_zero_indices, drop = FALSE]
      S_filtered <- S[non_zero_indices]
      U_filtered <- U[, non_zero_indices, drop = FALSE]
      inverse <- V_filtered %*% diag(1/S_filtered) %*% t(U_filtered)
    } else {
      inverse <- matrix(0, nrow = ncol(A), ncol = nrow(A))
    }
    return(inverse)
  }

  for (i in 1:n.curves) {
    # Least squares Regression
    fourier.koeffi[, i] = pseudo_inverse(t(fourier.s) %*% fourier.s) %*%
      t(fourier.s) %*% data[, i]
    # Fourier curve
    fourier.real[, i] = fourier.s %*% fourier.koeffi[, i]
  }

  # Mean Fourier curve
  fourier.mean[, 1] = rowMeans(fourier.koeffi)
  fourier.real_mw[, 1] = fourier.s %*% fourier.mean[, 1]

  # Standard deviation of the Fourier curve
  for (i in 1:n.curves) {
    # Variance-covariance matrix
    fourier.std1[, , i] <- (fourier.koeffi[, i] - fourier.mean[, 1]) %*%
                           t(fourier.koeffi[, i] - fourier.mean[, 1])
  }

  fourier.cov <- apply(fourier.std1, c(1, 2), mean)
  # Lenhoff, Appendix A, Eq. (0.5)
  fourier.std_all <- suppressWarnings(sqrt(fourier.s %*% fourier.cov %*%
                     t(fourier.s))
                     )

  for (i in 1:n.time) {
    # Values are on the diagonal of the square matrix fourier.std_all
    fourier.std[i, 1] = fourier.std_all[i, i]
  }

  # Bootstrap ------------------------------------------------------------------
  bootstrap_sample        <- array(data = 0, dim = c(n.time, 4))
  bootstrap.mean          <- array(data = 0, dim = c(k.coef*2 + 1, B))
  bootstrap.real_mw       <- array(data = 0, dim = c(n.time, B))
  bootstrap.zz            <- array(data = 0, dim = c(n.curves, B))
  bootstrap.pseudo_koeffi <- array(data = 0, dim = c(k.coef*2 + 1, n.curves, B))
  bootstrap.real          <- array(data = 0, dim = c(n.time, n.curves, B))
  bootstrap.std1          <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1, n.curves))
  bootstrap.cov           <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1, B))
  bootstrap.std_all       <- array(data = 0, dim = c(n.time, n.time, B))
  bootstrap.std           <- array(data = 0, dim = c(n.time, B))

  # Run two-stage (cluster) bootstrap ------------------------------------------
  #
  # Quote from the Journal of Biomechanics paper (Koska et al., 2023)
  # [...] we implemented a second, modified version of the BOOT method, in
  # which multiple curves per subject are accounted for (BOOTrep).
  # Therefore, BOOTrep includes the two-stage bootstrap process described
  # in Davison and Hinkley (1997), in which subjects (including all of
  # their curves) are sampled with replacement in the first stage, and one
  # curve per subject is drawn without replacement in the second
  # stage.

  for (i in 1:B) {
    if (iid == FALSE) {
      for (k in 1:n.curves) {
        # STAGE 1: Sample curve clusters (including all curves) with replacement
        stage.1.idx <- sample(1:n.cluster, size = n.cluster, replace = TRUE)
        curves.stage.2 <- c()
        for (curve.idx in stage.1.idx) {
          # Here, the indices of a single cluster (drawn with replacement) are selected
          curve.idx.clustername <- names(cluster_boundaries)[curve.idx]
          curve.numbers.stage.1 <- get_cluster_indices(curve.idx.clustername, cluster_boundaries)

          # STAGE 2: Sample within stage clusters without replacement
          sample.curve.index <- sample(curve.numbers.stage.1, size = 1, replace = FALSE)
          # while (sample.curve.index %in% curves.stage.2) { # Assure drawing without replacement
          #   sample.curve.index <- sample(curve.numbers.stage.1, size = 1)
          # }
          curves.stage.2 <- c(curves.stage.2, sample.curve.index)
        }

        for (clust.idx in 1:n.cluster) {
          bootstrap.zz[k, i] = curves.stage.2[clust.idx] # Old version: curves[k] # Hier liegt der Hase im Pfeffer! Ab k=12 wirft es NA's
          bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
          bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
        }
      }
    # If iid == TRUE: Run ordinary (naive) bootstrap
    } else {
      for (k in 1:n.curves) {
        bootstrap.zz[k, i] = sample(n.curves, size=1)
        bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
        bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
      }
    }

    # Mean bootstrap curve and standard deviation
    bootstrap.mean[, i] <- rowMeans(bootstrap.pseudo_koeffi[, , 1])
    bootstrap.real_mw[, i] <- fourier.s %*% bootstrap.mean[, i]

    for (k in 1:n.curves) {
      bootstrap.std1[, , k] <- (bootstrap.pseudo_koeffi[, k, i] -
                                  bootstrap.mean[, i]) %*%
                        t(bootstrap.pseudo_koeffi[, k, i] - bootstrap.mean[, i])
    }

    bootstrap.cov[, , i] <- apply(bootstrap.std1, c(1, 2), mean)
    bootstrap.std_all[, , i] <- suppressWarnings(sqrt(fourier.s %*%
                                bootstrap.cov[, , i] %*% t(fourier.s))
                                )

    for (k in 1:n.time) {
      bootstrap.std[k, i] <- bootstrap.std_all[k, k, i]
    }
  }

  # Construct bands ------------------------------------------------------------
  band.mean <- rowMeans(bootstrap.real_mw)
  band.sd   <- rowMeans(bootstrap.std)

  if (type == "prediction") {
    cp.data   <- array(data = 0, dim = c(n.curves, B))
    cp.data_i <- array(data = 0, dim = c(n.curves, B))

    cp.mean <- 0
    cp.bound <- 0
    while (cp.mean < (1-alpha)) {
      for (i in 1:B) {
        for (k in 1:n.curves) {
          # Lenhoff et al., Appendix A, Eq. (0.6)
          cp.data[k, i] <- max(abs(fourier.real[, k] - bootstrap.real_mw[, i]) / bootstrap.std[, i])
          cp.data_i[k, i] <- cp.data[k, i] < cp.bound
        }
      }
      cp.mean <- mean(cp.data_i)
      cp.bound <- cp.bound + 0.05
    }

    band.boot <- rbind(band.mean + cp.bound * band.sd,
                       band.mean,
                       band.mean - cp.bound * band.sd
    )
  } else if (type == "confidence") {
    cc.data <- array(data = 0, dim = c(n.curves, B))

    for (i in 1:B) {
      for (k in 1:n.curves) {
        # Lenhoff, Appendix A, Eq. (0.8)
        cc.data[k, i] <- max(abs(bootstrap.real_mw[, i] - fourier.real_mw) /
                               bootstrap.std[, i])
      }
    }
    cc <- quantile(cc.data, probs = 1-alpha)

    band.boot <- rbind(band.mean + cc * band.sd,
                       band.mean,
                       band.mean - cc * band.sd
    )
  }

  row.names(band.boot) <- c("upper", "mean", "lower")

  return(band.boot)
}



################################################################################
# The functions quantile() and format_perc() (the latter is used in quantile())
# are included at the end of this script. Otherwise, they would have to be im-
# ported from the 'stats' package (via: @importFrom stats quantile)
################################################################################

quantile <- function (x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE,
                      type = 7, digits = 7, ...)
{
  if (is.factor(x)) {
    if (is.ordered(x)) {
      if (!any(type == c(1L, 3L)))
        stop("'type' must be 1 or 3 for ordered factors")
    }
    else stop("(unordered) factors are not allowed")
    lx <- levels(x)
    x <- as.integer(x)
  }
  else {
    if (is.null(x))
      x <- numeric()
    lx <- NULL
  }
  if (na.rm)
    x <- x[!is.na(x)]
  else if (anyNA(x))
    stop("missing values and NaN's not allowed if 'na.rm' is FALSE")
  eps <- 100 * .Machine$double.eps
  if (any((p.ok <- !is.na(probs)) & (probs < -eps | probs >
                                     1 + eps)))
    stop("'probs' outside [0,1]")
  n <- length(x)
  probs <- pmax(0, pmin(1, probs))
  np <- length(probs)
  {
    if (type == 7) {
      index <- 1 + max(n - 1, 0) * probs
      lo <- floor(index)
      hi <- ceiling(index)
      x <- sort(x, partial = if (n == 0)
        numeric()
        else unique(c(lo, hi)[p.ok]))
      qs <- x[lo]
      i <- which(!p.ok | (index > lo & x[hi] != qs))
      h <- (index - lo)[i]
      qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
    }
    else {
      if (type <= 3) {
        nppm <- if (type == 3)
          n * probs - 0.5
        else n * probs
        j <- floor(nppm)
        h <- switch(type, !p.ok | (nppm > j), ((nppm >
                                                  j) + 1)/2, !p.ok | (nppm != j) | ((j%%2L) ==
                                                                                      1L))
      }
      else {
        switch(type - 3, {
          a <- 0
          b <- 1
        }, a <- b <- 0.5, a <- b <- 0, a <- b <- 1, a <- b <- 1/3,
        a <- b <- 3/8)
        fuzz <- 4 * .Machine$double.eps
        nppm <- a + probs * (n + 1 - a - b)
        j <- floor(nppm + fuzz)
        h <- nppm - j
        if (any(sml <- abs(h) < fuzz, na.rm = TRUE))
          h[sml] <- 0
      }
      x <- sort(x, partial = if (n == 0)
        numeric()
        else unique(c(1, j[p.ok & j > 0L & j <= n], (j +
                                                       1)[p.ok & j > 0L & j < n], n)))
      x <- c(x[1L], x[1L], x, x[n], x[n])
      qs <- x[j + 2L]
      qs[!is.na(h) & h == 1] <- x[j + 3L][!is.na(h) & h ==
                                            1]
      other <- (0 < h) & (h < 1) & (x[j + 2L] != x[j +
                                                     3L])
      other[is.na(other)] <- TRUE
      if (any(other))
        qs[other] <- ((1 - h) * x[j + 2L] + h * x[j +
                                                    3L])[other]
    }
  }
  qs[!p.ok] <- probs[!p.ok]
  if (is.character(lx))
    qs <- factor(qs, levels = seq_along(lx), labels = lx,
                 ordered = TRUE)
  if (names && np > 0L) {
    stopifnot(is.numeric(digits), digits >= 1)
    names(qs) <- format_perc(probs, digits = digits)
  }
  qs
}



format_perc <- function (x, digits = max(2L, getOption("digits")), probability = TRUE,
                         use.fC = length(x) < 100, ...)
{
  if (length(x)) {
    if (probability)
      x <- 100 * x
    ans <- paste0(if (use.fC)
      formatC(x, format = "fg", width = 1, digits = digits)
      else format(x, trim = TRUE, digits = digits, ...), "%")
    ans[is.na(x)] <- ""
    ans
  }
  else character(0)
}
