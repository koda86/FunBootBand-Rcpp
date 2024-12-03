# Rcpp implementation of FunBootBand
library(Rcpp)

data <- read.csv("~/FunBootBand-Rcpp/example_data.csv")

# Load the functions created in the R-Cpp source file (FunBootBand_cpp.cpp)
sourceCpp("~/FunBootBand-Rcpp/FunBootBand_cpp.cpp")

# Testweise Laden der Funtion getDimensions
getDimensions(data)
