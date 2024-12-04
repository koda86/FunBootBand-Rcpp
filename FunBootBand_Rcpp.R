# Rcpp implementation of FunBootBand
library(Rcpp)

data <- read.csv("~/FunBootBand-Rcpp/example_data.csv")

# Load the functions created in the R-Cpp source file (FunBootBand_cpp.cpp)
sourceCpp("~/FunBootBand-Rcpp/FunBootBand_cpp.cpp")

# Testweise Laden der Function getDimensions
getDimensions(data)

# Example usage of the Moore-Penrose pseudoinverse function
A <- matrix(c(1, 2, 3, 4), nrow = 2)
pseudo_inverse(A)

# Split string function
get_base_name("cluster.1.curve") 


