# Load the experiment: compile the C++ expectile solver (cached by Rcpp) and
# source the R side. Every script starts with source("R/setup.R").
Rcpp::sourceCpp("src/gev_expectile.cpp", cacheDir = "src/.rcpp-cache")
for (f in c("gev", "dgp", "weight", "quadrature", "estimators"))
  source(file.path("R", paste0(f, ".R")))
