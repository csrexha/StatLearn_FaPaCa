#' data simulator
#'
#' @param N no of samples
#' @param P no of predictors
#' @param p_ref no of relevant predictors
#' @param tau bias
#' @param sigma variance
#' @param rho correlation
#' @param link link function
#'
#' @returns list of y and x
#' @export

binary_data_generator <- function(N, P, p_ref, tau, sigma, rho, link = "logit"){
    
  linkinv <- make.link(link)$linkinv
  
  f <- rnorm(N, 0, 1)
  y <- rbinom(N, 1, linkinv(tau + sigma*f))
  
  # regenerate if total number of minor group smaller than 4
  while(sum(y) < 4){
    y <- rbinom(N, 1, linkinv(tau + sigma*f))
  }
  
  x <- cbind(sqrt(1 - rho) * matrix(rnorm(N*p_ref, 0, 1), nrow = N) +
               sqrt(rho) * matrix(rep(f, p_ref), nrow = N),
             matrix(rnorm(N * (P - p_ref), 0, 1), nrow = N))
  
  return(list(y = y, x = x))
}
