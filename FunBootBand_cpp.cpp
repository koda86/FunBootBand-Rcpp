// #include <Rcpp.h>
#include <RcppArmadillo.h> // Armadillo: C++ library for linear algebra and scientific computing
#include <iostream> // For std::cout and std::endl
#include <string>   // For std::string
#include <regex>    // For std::regex
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
void validate_band_arguments(
    const DataFrame& data,    // Input data
    const std::string& type,  // Type: "confidence" or "prediction"
    double alpha,             // Alpha: between 0 and 1
    bool iid,                 // IID: true or false
    int k_coef,               // k.coef: positive integer
    int B                     // B: positive integer
) {
  // Check 'type'
  if (type.empty() || (type != "confidence" && type != "prediction")) {
    stop("'type' must be either 'confidence' or 'prediction'.");
  }
  
  // Check 'alpha'
  if (alpha <= 0.0 || alpha >= 1.0) {
    stop("'alpha' must be a numeric value between 0 and 1.");
  }
  
  // Check 'iid' (logical in R maps to bool in C++)
  if (!LogicalVector::is_na(iid) && (iid != true && iid != false)) {
    stop("'iid' must be a logical value (TRUE or FALSE).");
  }
  
  // Check 'k.coef'
  if (k_coef <= 0) {
    stop("'k.coef' must be a positive integer.");
  }
  
  // Check 'B'
  if (B <= 0) {
    stop("'B' must be a positive integer.");
  }
  
  // Check for NA values in the DataFrame
  LogicalVector na_check = is_na(data);
  if (is_true(any(na_check))) {
    stop("Function stopped due to NA's in the input data.");
  }
  
  // Check if the input data is a data frame (already enforced by Rcpp::DataFrame)
  // Additional validation could be added if necessary.
  
  // Check if all elements of the data frame are numeric
  for (int i = 0; i < data.size(); i++) {
    if (!Rf_isNumeric(data[i])) {
      stop("Non-numeric data found in input.");
    }
  }
}
 
// [[Rcpp::export]]
List getDimensions(DataFrame data) {
 int n_curves = data.ncol();
 int n_time = data.nrow();
 
 return List::create(Named("n_time") = n_time,
                     Named("n_curves") = n_curves);
}

// [[Rcpp::export]]

 
// --------------------------------------------------------------------------
// This implements a more robust approach that allows the detection of the
// number and size (i.e., the number of curves per cluster) in cases where 
// different clusters contain a varying number of curves per cluster.
//
// Function to determine the base name (or cluster identifier) of a column
// --------------------------------------------------------------------------

// Helper function to split a string using a regular expression
std::vector<std::string> split(const std::string& str, const std::string& regex_pattern) {
  std::regex regex_delim(regex_pattern);
  std::sregex_token_iterator iter(str.begin(), str.end(), regex_delim, -1);
  std::sregex_token_iterator end;
  return {iter, end};
}

// Function to get the base name
// [[Rcpp::export]]
std::string get_base_name(const std::string& name) {
  // Split the string at any non-alphanumeric character
  std::vector<std::string> parts = split(name, R"([^a-zA-Z0-9])");
  // Return the first part if it exists
  return !parts.empty() ? parts[0] : "";
}

// [[Rcpp::export]]
List process_clusters(CharacterVector colnames) {
 List clusters; // Equivalent to R's list()
 
 for (auto& name : colnames) {
   // Extract base name
   std::string base_name = get_base_name(as<std::string>(name));
   
   // Check if the base name exists in the clusters
   if (!clusters.containsElementNamed(base_name.c_str())) {
     // If not, initialize it with 0
     clusters[base_name] = 0;
   }
   
   // Increment the count for this base name
   clusters[base_name] = as<int>(clusters[base_name]) + 1;
 }
 
 return clusters;
}

 // [[Rcpp::export]]
 List calculate_cluster_boundaries(List clusters) {
   List cluster_boundaries; // Equivalent to R's list()
   int start_idx = 1;
   
   // Iterate over the names of the clusters
   CharacterVector cluster_ids = clusters.names();
   for (const auto& cluster_id : cluster_ids) {
     std::string cluster_name = as<std::string>(cluster_id);
     int cluster_size = as<int>(clusters[cluster_name]);
     int end_idx = start_idx + cluster_size - 1;
     
     // Create a named vector for start and end indices
     NumericVector boundary = NumericVector::create(
       _["start"] = start_idx,
       _["end"] = end_idx
     );
     
     // Assign the boundary to the cluster ID
     cluster_boundaries[cluster_name] = boundary;
     
     // Update start index for the next cluster
     start_idx = end_idx + 1;
   }
   
   return cluster_boundaries;
 }
 
// [[Rcpp::export]]
IntegerVector get_cluster_indices(const std::string& cluster_id, List cluster_boundaries) {
 // Check if the cluster ID exists in the cluster boundaries
 if (cluster_boundaries.containsElementNamed(cluster_id.c_str())) {
   // Retrieve the boundaries for the given cluster ID
   NumericVector boundaries = cluster_boundaries[cluster_id];
   
   // Extract start and end values
   int start = boundaries["start"];
   int end = boundaries["end"];
   
   // Generate and return the sequence of indices
   return seq(start, end);
 } else {
   // Throw an error if the cluster ID is not found
   stop("Cluster ID not found.");
 }
}

// [[Rcpp::export]]
void validate_nested_structure(List clusters, int ncol_data, bool iid) {
 int n_cluster = clusters.size(); // Equivalent to length(clusters) in R
 
 // Check the condition for nested structure
 if (n_cluster < 2 || n_cluster == ncol_data) {
   if (!iid) {
     // Throw an error if 'iid' is set to FALSE
     stop("Header does not indicate a nested structure even though 'iid' is set to 'FALSE'.");
   }
 }
}



// Approximate curves using Fourier functions ---------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List initialize_fourier_arrays(int k_coef, int n_curves, int n_time) {
  // Create and initialize arrays
  NumericMatrix fourier_koeffi((k_coef * 2) + 1, n_curves); // Create the matrix
  std::fill(fourier_koeffi.begin(), fourier_koeffi.end(), 0.0); // Initialize to 0
  
  NumericMatrix fourier_real(n_time, n_curves);
  std::fill(fourier_real.begin(), fourier_real.end(), 0.0);
  
  NumericMatrix fourier_mean((k_coef * 2) + 1, (k_coef * 2) + 1);
  std::fill(fourier_mean.begin(), fourier_mean.end(), 0.0);
  
  NumericMatrix fourier_real_mw(n_time, 1);
  std::fill(fourier_real_mw.begin(), fourier_real_mw.end(), 0.0);
  
  NumericVector fourier_std1((k_coef * 2 + 1) * (k_coef * 2 + 1) * n_curves);
  std::fill(fourier_std1.begin(), fourier_std1.end(), 0.0);
  
  NumericMatrix fourier_cov((k_coef * 2) + 1, (k_coef * 2) + 1);
  std::fill(fourier_cov.begin(), fourier_cov.end(), 0.0);
  
  NumericMatrix fourier_std_all(n_time, n_time);
  std::fill(fourier_std_all.begin(), fourier_std_all.end(), 0.0);
  
  NumericMatrix fourier_std(n_time, 1);
  std::fill(fourier_std.begin(), fourier_std.end(), 0.0);
  
  // Return the arrays as a named list
  return List::create(
    _["fourier_koeffi"] = fourier_koeffi,
    _["fourier_real"] = fourier_real,
    _["fourier_mean"] = fourier_mean,
    _["fourier_real_mw"] = fourier_real_mw,
    _["fourier_std1"] = fourier_std1,
    _["fourier_cov"] = fourier_cov,
    _["fourier_std_all"] = fourier_std_all,
    _["fourier_std"] = fourier_std
  );
}


// Construct Fourier series
// General: f(t) = mu + sum(alpha cos(2pi*k*t/T) + beta sin(2pi*k*t/T))
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

// Helper function to calculate the pseudoinverse (Moore-Penrose)
// [[Rcpp::depends(RcppArmadillo)]]
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

// [[Rcpp::export]]
void calculate_fourier_coefficients_and_curves(
    arma::mat& fourier_koeffi,   // Matrix to store Fourier coefficients
    arma::mat& fourier_real,     // Matrix to store Fourier real curves
    const arma::mat& fourier_s,  // Fourier basis matrix
    const arma::mat& data        // Data matrix
) {
  int n_curves = data.n_cols;
  
  for (int i = 0; i < n_curves; i++) {
    // Least squares regression to calculate Fourier coefficients
    arma::mat pseudo_inv = arma::pinv(fourier_s.t() * fourier_s); // Pseudo-inverse
    fourier_koeffi.col(i) = pseudo_inv * fourier_s.t() * data.col(i);
    
    // Calculate Fourier curve
    fourier_real.col(i) = fourier_s * fourier_koeffi.col(i);
  }
}

// [[Rcpp::export]]
void calculate_mean_fourier_curve(
    arma::mat& fourier_mean,       // Matrix to store Fourier means
    arma::mat& fourier_real_mw,    // Matrix to store mean Fourier curve
    const arma::mat& fourier_koeffi, // Matrix of Fourier coefficients
    const arma::mat& fourier_s     // Fourier basis matrix
) {
  // Compute the mean of Fourier coefficients across columns (equivalent to rowMeans)
  fourier_mean.col(0) = arma::mean(fourier_koeffi, 1); // Mean along rows
  
  // Compute the mean Fourier curve
  fourier_real_mw.col(0) = fourier_s * fourier_mean.col(0); // Matrix multiplication
}

// [[Rcpp::export]]
void calculate_fourier_std1(
    arma::cube& fourier_std1,          // 3D array to store variance-covariance matrices
    const arma::mat& fourier_koeffi,   // Matrix of Fourier coefficients
    const arma::vec& fourier_mean      // Vector of mean Fourier coefficients
) {
  int n_curves = fourier_koeffi.n_cols;
  
  for (int i = 0; i < n_curves; i++) {
    // Compute the difference (fourier.koeffi[, i] - fourier.mean[, 1])
    arma::vec diff = fourier_koeffi.col(i) - fourier_mean;
    
    // Compute the variance-covariance matrix and assign to fourier.std1[, , i]
    fourier_std1.slice(i) = diff * diff.t();
  }
}


// [[Rcpp::export]]
void calculate_fourier_cov_and_std(
    arma::mat& fourier_cov,          // Matrix to store the covariance matrix
    arma::mat& fourier_std_all,      // Matrix to store the standard deviation results
    const arma::cube& fourier_std1, // 3D array of variance-covariance matrices
    const arma::mat& fourier_s       // Fourier basis matrix
) {
  int n_rows = fourier_std1.n_rows;
  int n_cols = fourier_std1.n_cols;
  int n_slices = fourier_std1.n_slices;
  
  // Calculate the mean covariance matrix (apply mean across slices)
  fourier_cov.zeros(n_rows, n_cols); // Initialize to zero
  for (int i = 0; i < n_slices; i++) {
    fourier_cov += fourier_std1.slice(i);
  }
  fourier_cov /= n_slices; // Divide by the number of slices to compute the mean
  
  // Compute the standard deviation matrix
  // fourier_std_all = sqrt(fourier_s %*% fourier_cov %*% t(fourier_s))
  arma::mat temp = fourier_s * fourier_cov * fourier_s.t();
  fourier_std_all = arma::sqrt(temp); // Element-wise square root
}

// [[Rcpp::export]]
void extract_diagonal_to_fourier_std(
    arma::mat& fourier_std,       // Matrix to store extracted diagonal values
    const arma::mat& fourier_std_all // Square matrix (fourier.std_all)
) {
  int n_time = fourier_std_all.n_rows;
  
  for (int i = 0; i < n_time; i++) {
    // Extract diagonal element and store in fourier_std[i, 1]
    fourier_std(i, 0) = fourier_std_all(i, i);
  }
}



