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
#' @importFrom stats quantile
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

# TODO:
# - Preprint des JoB papers auf der Github Seite adden
# - test coverage (covr) ... chapter 13 https://r-pkgs.org/testing-basics.html
# - Checken ob die Funktion Fehler auswirft wenn Kurven unterschiedlich lang sind
# - Peer review, e.g. https://ropensci.org/software-review/
# - Vignette schreiben

library(Rcpp)

# Source the C++ functions
sourceCpp('~/FunBootBand-Rcpp/FunBootBand_cpp.cpp')

# Dummy Rcpp function to check if implementation works
hello_world()

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

  n.curves  <- dim(data)[2]
  n.time <- dim(data)[1]
  time <- seq(0, (n.time - 1))

  if (iid == FALSE) {
    if (is.data.frame(data) == FALSE) {stop("Input data is not a data frame.")}
    n.cluster <- length(unique(colnames(data)))
    curves.per.cluster <- n.curves / n.cluster
    if (n.cluster < 2 | n.cluster == ncol(data)) {stop("The header does not\n
        indicate a nested structure even though 'iid' is set to 'FALSE'.")}
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
  # fouriers.s <- constructFourierSeries(100, 50) # Use Cpp function
  
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
  bootstrap.std1          <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1,
                                                     n.curves))
  bootstrap.cov           <- array(data = 0, dim = c(k.coef*2 + 1, k.coef*2 + 1,
                                                     B))
  bootstrap.std_all       <- array(data = 0, dim = c(n.time, n.time, B))
  bootstrap.std           <- array(data = 0, dim = c(n.time, B))

  for (i in 1:B) {
    if (iid == FALSE) { # Run two-stage (cluster) bootstrap
      for (k in 1:curves.per.cluster) {
        # STAGE 1: Sample curve clusters with replacement
        stage.1.idx <- sample(1:n.cluster, size = n.cluster, replace = TRUE)
        # STAGE 2: Sample within stage clusters without replacement
        curves <- c()
        for (curve.idx in stage.1.idx) {
          curve.numbers.stage.1 <- seq(from = curve.idx*curves.per.cluster -
                                          curves.per.cluster + 1,
                                        to = curve.idx*curves.per.cluster)
          tmp <- sample(curve.numbers.stage.1, size = 1, replace = FALSE)
          while (tmp %in% curves) { # Assure drawing without replacement
            tmp <- sample(curve.numbers.stage.1, size = 1)
          }
          curves <- c(curves, tmp)
        }
        bootstrap.zz[k, i] = curves[k]
        bootstrap.pseudo_koeffi[, k, i] = fourier.koeffi[, bootstrap.zz[k, i]]
        bootstrap.real[, k, i] = fourier.s %*% bootstrap.pseudo_koeffi[, k, i]
    }
    } else { # Run 'ordinary' (naive) bootstrap
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
    cp.bound <- 0 # cp.begin
    while (cp.mean < (1-alpha)) {
      for (i in 1:B) {
        for (k in 1:n.curves) {
          # Lenhoff et al., Appendix A, Eq. (0.6)
          cp.data[k, i] <- max(abs(fourier.real[, k] - bootstrap.real_mw[, i]) /
                                 bootstrap.std[, i])
          cp.data_i[k, i] <- cp.data[k, i] < cp.bound
        }
      }
      cp.mean <- mean(cp.data_i)
      cp.bound <- cp.bound + 0.05
    }
    cp_out <- cp.bound

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

