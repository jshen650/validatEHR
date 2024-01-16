#' Induce misclassification in a binary outcome
#'
#' @param yVec Vector of outcomes
#' @param p11 Sensitivity, P(Y*=1 | Y=1)
#' @param p10 Value for P(Y*=1 | Y=0), otherwise known as (1-Specificity)
#'
#' @return
#' @export
#'
#' @examples
induceME <- function(yVec, p11, p10){
  
  y.ABcopy = yVec
  Y1 = which(y.ABcopy==1)
  Y0 = which(y.ABcopy==0)
  
  Ystar = y.ABcopy
  Ystar[Y1] = rbinom(length(Y1), 1,p11)
  Ystar[Y0] = rbinom(length(Y0), 1,p10)
  
  return(Ystar)
  
}