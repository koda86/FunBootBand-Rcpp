#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void hello_world () {
  std::cout << "Hello, world!";
}

/*
// [[Rcpp::export]]
NumericMatrix constructFourierSeries(int n_time, int k_coef) {
  NumericMatrix fourier_s(n_time, k_coef * 2 + 1);
  fourier_s(_, 0) = 1; // First column as ones
  
  for (int k = 1; k <= k_coef * 2; k += 2) {
    fourier_s(_, k) = cos(2 * M_PI * (k / 2) * seq_len(n_time) / (n_time - 1));
    fourier_s(_, k + 1) = sin(2 * M_PI * (k / 2) * seq_len(n_time) / (n_time - 1));
  }
  return fourier_s;
}

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