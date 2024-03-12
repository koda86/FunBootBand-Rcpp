#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void hello_world () {
  std::cout << "Hello, world!";
}

// [[Rcpp::export]]
List getDimensions(DataFrame data) {
  int nColumns = data.ncol();
  int nRows = data.nrow();
  
  return List::create(Named("nRows") = nRows,
                      Named("nColumns") = nColumns);
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