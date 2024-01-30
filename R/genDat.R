#' Generate a simulated probability sample with a nested validation sample
#'
#' @param seedNum Seed used for random generation
#' @param p11 Sensitivity, P(Y*=1 | Y=1)
#' @param p10 Value for P(Y*=1 | Y=0), otherwise known as (1-Specificity)
#' @param valType Type of validation sample (SRS or non-probability)
#' @param n Size of entire probability sample
#' @param alphaInt Intercept value relevant for size of a non-probability validation sample
#' 
#' @import stats
#'
#' @return Data frame containing error-free outcomes, misclassified outcomes, 
#' and indices for individuals belonging to the validation sample
#' @export
genDat <- function(seedNum, p11, p10, valType, n=NULL, alphaInt= NULL){
  set.seed(seedNum)
  
  ## for n total individuals, consider 5 covariates
  ## and binry exposure variable
  if(is.null(n)){
    n <- 5000
  }
  
  p <- 5
  
  x.A <- matrix( rnorm(n*(p),0,1),n,(p))
  pt.A <- expit(0.8+0.3*rowSums(x.A)) 
  
  t.A <- rbinom(n, 1, pt.A)
  py.A <- expit(-3.9 + t.A+rowSums(x.A))
  
  y <- rbinom(n, 1, py.A) 
  
  
  partDat = data.frame(t.A, x.A)
  estAlpha <- glm(t.A ~ . , data=partDat, family="binomial")
  estAlpha_coef <- coef(estAlpha)
  estE <- expit(cbind(1,x.A)%*%estAlpha_coef)
  
  ## generate outcomes subject to misclassification
  yStar <- induceME(y, p11=p11, p10=p10)
  y.Both = data.frame(y, yStar)
  
  if(is.null(alphaInt)){
    alphaInt <- -2.8
  }
  
  ## generate validation sample
  if(valType=="SRS"){
    inVal <- rbinom(n, size=1, prob=0.17)
  } else{
    alpha0 <- c(alphaInt, 0.5, 1, 1, 1, 1, rep(0, p-4))
    probB <- (1+exp(-cbind(1,t.A, x.A)%*%alpha0))^(-1)
    inVal <- rbinom(n,size = 1,prob = (probB) )
  }
  
  B.loc <- which(inVal == 1)
  
  
  allDat = data.frame(y.Both, inVal, t.A, x.A)
  
  return(allDat)
  
}