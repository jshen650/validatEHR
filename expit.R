#' Expit function
#'
#' @param x Input value(s)
#'
#' @return Resulting values from applying function
expit<-function(x){
  exp(x)/(1+exp(x))
}
