# Load the experiment: compile the C++ expectile solver (cached by Rcpp) and
# source the R side. Every script starts with source("R/setup.R").
Rcpp::sourceCpp("src/gev_expectile.cpp", cacheDir = "src/.rcpp-cache")
Rcpp::sourceCpp("src/gev_elastile.cpp", cacheDir = "src/.rcpp-cache")
for (f in c("gev", "dgp", "graft", "weight", "quadrature", "estimators"))
  source(file.path("R", paste0(f, ".R")))
