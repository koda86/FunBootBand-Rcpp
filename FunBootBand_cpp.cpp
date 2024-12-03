// #include <Rcpp.h>
#include <RcppArmadillo.h> // Armadillo: C++ library for linear algebra and scientific computing

using namespace Rcpp;

// [[Rcpp::export]]
void hello_world () {
  std::cout << "Hello, world!";
}

// [[Rcpp::export]]
List getDimensions(DataFrame data) {
  int n_curves = data.ncol();
  int n_time = data.nrow();
  
  return List::create(Named("n_time") = n_time,
                      Named("n_curves") = n_curves);
}

// [[Rcpp::export]]
NumericMatrix constructFourierSeries(int n_time, int k_coef) {
  NumericMatrix fourier_s(n_time, k_coef * 2 + 1);
  
  // Fill the first column with ones
  for (int i = 0; i < n_time; ++i) {
    fourier_s(i, 0) = 1;
  }
  
  for (int k = 1; k <= k_coef * 2; k += 2) {
    for (int t = 0; t < n_time; ++t) {
      double normalized_time = static_cast<double>(t) / (n_time - 1);
      fourier_s(t, k) = cos(2 * M_PI * (k / 2) * normalized_time);
      fourier_s(t, k + 1) = sin(2 * M_PI * (k / 2) * normalized_time);
    }
  }
  return fourier_s;
}

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function to calculate the pseudoinverse (Moore-Penrose)
// [[Rcpp::export]]
arma::mat pseudo_inverse(const arma::mat& A, Rcpp::Nullable<double> tol = R_NilValue) {
  // Set default value for tol if not provided
  double tolerance = tol.isNotNull() ? Rcpp::as<double>(tol) : std::pow(std::numeric_limits<double>::epsilon(), 2.0 / 3.0);
  
  // Perform singular value decomposition
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, A);
  
  // Set the threshold for small singular values
  double threshold = tolerance * S(0); // S(0) is the largest singular value
  arma::uvec non_zero_indices = arma::find(S > threshold);
  
  if (non_zero_indices.n_elem == 0) {
    // All singular values are below the threshold
    return arma::zeros(A.n_cols, A.n_rows);
  } else {
    // Filter out small singular values and compute the pseudoinverse
    arma::vec S_inv = arma::zeros(S.n_elem);
    S_inv(non_zero_indices) = 1 / S(non_zero_indices);
    return V * arma::diagmat(S_inv) * U.t();
  }
}

/*
// [[Rcpp::export]]
NumericMatrix pseudoInverse(const NumericMatrix& A) {
  SVD svd(A);
  NumericMatrix U = svd.u;
  NumericVector S = svd.d;
  NumericMatrix V = svd.v;
  
  // Compute the pseudo-inverse
  double tolerance = std::max(A.nrow(), A.ncol()) * S[0] * std::numeric_limits<double>::epsilon();
  NumericMatrix S_inv = diag(1 / S);
  for (int i = 0; i < S.length(); ++i) {
    if (S[i] < tolerance) S_inv(i, i) = 0;
  }
  return V * S_inv * transpose(U);
}
 */