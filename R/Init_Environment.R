# Initialisation file to load all functions into the environment
if (!require("Rcpp")) install.packages("Rcpp")
library("Rcpp")

path <- getwd()

# Load R functions
suppressWarnings(eval(parse(file = paste0(path,"/R/get_loss.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/process_group.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/process_lambda.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/process_penalty.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/process_X.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/process_y.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/process_Z.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/sgp.R"))))
suppressWarnings(eval(parse(file = paste0(path,"/R/sgp_cv.R"))))

# Load C++ functions
suppressWarnings(Rcpp::sourceCpp(paste0(path,"/CPP/penalty_functions.cpp")))
suppressWarnings(Rcpp::sourceCpp(paste0(path,"/CPP/get_functions.cpp")))
suppressWarnings(Rcpp::sourceCpp(paste0(path,"/CPP/max_cor.cpp")))
suppressWarnings(Rcpp::sourceCpp(paste0(path,"/CPP/linear_fit.cpp")))
suppressWarnings(Rcpp::sourceCpp(paste0(path,"/CPP/log_fit.cpp")))

