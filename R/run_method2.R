#' Run approach from Equation 2 of paper that uses gold outcomes from validation data 
#' and silver standard outcomes from non-validation individuals
#'
#' @param df Data frame containing exposure status, gold-standard outcomes, silver-standard outcomes, 
#' and all relevant variables for analyses
#' @param inVal Vector that indicates whether an individual was included in the validation data (=1 for included)
#' @param varTrt Variable name for treatment/exposure
#' @param varGold Variable name for gold-standard outcome
#' @param varSilver Variable name for silver-standard outcome
#' @param varTrtMod Variable names for modeling treatment propensity model
#' @param varSelectMod Variable names for modeling validation sample selection model
#'
#' @return Data frame containing weight used, final estimate of ATE, final variance estimate of ATE,
#' estimate of ATE using validation data, variance estimate of ATE using validation data,
#' estimate of ATE using silver standard outcomes, and variance estimate of ATE using silver standard outcomes
#' @export
run_method2 <- function(df, inVal, varTrt, varGold, varSilver, varTrtMod, varSelectMod){
  
  ## gather validation data
  valDat <- df[which(inVal==1),]
  
  
  t.A = df[,varTrt]
  y.A = df[,varGold]
  nB <- nrow(valDat)
  inVal.loc = which(inVal==1)
  
  ## xCols must be columns of data frame pertaining to main effects
  x.A = df[, varTrtMod]
  
  
  t.B = valDat[,varTrt]
  y.B = valDat[,varGold]
  x.B = x.A[inVal.loc,]
  yStar = df[,varSilver]
  y.AB=yStar
  
  nA = nrow(df)
  nB = nrow(valDat)
  
  ## estimate propensity scores for treatment
  partDat = data.frame(t.A, x.A)
  estAlpha <- glm(t.A ~ ., data=partDat, family="binomial")
  estE = predict.glm(estAlpha, partDat, type="response")
  
  partDat2 = data.frame(inVal, t.A, df[,c(varSelectMod)])
  estAlphaB <- glm(inVal ~ . , data=partDat2, family="binomial")
  estAlphaB_coef <- coef(estAlphaB)
  estB <- predict.glm(estAlphaB, partDat2, type="response")
  
  gold1 = which(valDat[,varGold]==1)
  gold0 = which(valDat[,varGold]==0)
  y.Both_Y1 = valDat[gold1,]
  y.Both_Y0 = valDat[gold0,]
  
  silver1 = which(y.Both_Y1[,varSilver]==1)
  silver1b = which(y.Both_Y0[,varSilver]==1)
  p11 = length(silver1)/nrow(y.Both_Y1)
  p10 = length(silver1b)/nrow(y.Both_Y0)
  
  ## estimate ATE using gold standard outcomes from validation individuals
  w1 = estE*estB
  w2 = estB*(1-estE)
  tauV_sum1 = mean (inVal*t.A*y.A/w1)
  tauV_sum2 = mean (inVal*(1-t.A)*y.A/w2)
  estATE = tauV_sum1 - tauV_sum2
  
  ## estimate ATE using silver standard outcomes from non-validation individuals
  wN1 = estE*(1-estB)
  wN2 = (1-estE)*(1-estB)
  tauN_sum1 = (1/(sum((1-inVal)*t.A/wN1)))*sum((1-inVal)*t.A*yStar/wN1)
  tauN_sum2 = (1/(sum((1-inVal)*(1-t.A)/wN2)))*sum((1-inVal)*(1-t.A)*yStar/wN2)
  tauN_B1_other = (1/(p11-p10))*(tauN_sum1 - tauN_sum2)
  
  ## design matrix relevant for treatment propensity model
  xDesign = model.matrix(estAlpha)
  xM = xDesign
  nXM = ncol(xM) 
  
  ## design matrix relevant for validation sample selection propensity model
  xDesign_Samp = model.matrix(estAlphaB)
  xtilde = xDesign_Samp
  nXtilde = ncol(xtilde) 
  
  ## weighted regression equivalent for estimating Hajek form of the ATE
  R = ((1-inVal)*t.A)/wN1 + ((1-inVal)*(1-t.A))/wN2
  regCoef = coef(glm(yStar ~ t.A, weights=R))
  alpha = regCoef[1]
  beta = regCoef[2]
  
  yScore = yStar - alpha - beta*(p11-p10)*t.A ## estimating equation involving Y*
  tScore = (t.A - estE)*estB ## estimating equation involving T
  V_est = inVal-estB
  dConst = ((1-inVal)*(1-t.A)/wN2)*estE - ((1-inVal)*(t.A)/wN1)*(1-estE)
  zConst = ((1-inVal)*t.A/wN1)*estB + ((1-inVal)*(1-t.A)/wN2)*estB
  g1 =(y.A*yStar - p11*y.A )*inVal*(nA/nB)
  g2 = ((1-y.A)*yStar - p10*(1-y.A))*inVal*(nA/nB)
  
  
  nVS = nXM+nXtilde+1 ## number of cols/rows relevant for estimating variance for validation data
  nAS = 4 ## number of cols/rows meant for non-val; to account for alpha, beta, p11, and p10 params
  
  matDim = nVS +nAS
  
  ## estimate sandwich variance
  
  ### the bread A
  AA = matrix(0,matDim, matDim)
  
  AA[1,1] <- 1
  
  for(i in 2:(nXM+1)){
    for(j in 2:(nXM+1)){
      AA[i, j] = mean(xM[,i-1]*xM[,j-1]*estE*(1-estE)*estB)
    }
  }
  
  for(i in (nXM+2):(nVS) ){
    for(j in (nXM+2):(nVS)){
      AA[i, j] = mean(xtilde[,(i-nXM-1)]*xtilde[,(j-nXM-1)]*estB*(1-estB))
    }
  }
  
  for(j in 2:(nXM+1)){
    AA[1, j] = mean( ( (inVal*t.A*y.A*(1-estE))/(estB*(estE)) + (inVal*(1-t.A)*y.A*estE)/(estB*(1-estE)))*xM[,j-1]  )
  }
  
  for(j in (nXM+2):(nVS)){
    AA[1, j] = mean( ( (inVal*t.A*y.A*(1-estB))/(estB*(estE)) - (inVal*(1-t.A)*y.A*(1-estB))/(estB*(1-estE)))*xtilde[,j-nXM-1]  )
  }
  
  for(i in 2:(nXM+1)){
    for(j in (nXM+2):nVS){
      AA[i, j] = mean(-1*(t.A-estE)*estB*(1-estB)*xM[,(i-1)]*xtilde[,(j-nXM-1)])
    }
  }
  
  AA[(nVS+1),(nVS+1)] <- mean(R)
  AA[(nVS+2),(nVS+1)] <- mean(R*t.A)
  AA[(nVS+1),(nVS+2)] <- mean(R*t.A*(p11-p10))
  AA[(nVS+2),(nVS+2)] <- mean(R*(t.A^2)*(p11-p10) )
  
  for(j in 3:(nXM+2)){
    AA[(nVS+1), (j-1)] = mean(-dConst*xM[,j-2]*yScore)
  }
  
  for(j in 3:(nXM+2)){
    AA[(nVS+2), (j-1)] = mean(-dConst*xM[,j-2]*t.A*yScore)
  }
  
  AA[(nVS+1),(nVS+3)] <- mean(R*beta*t.A)
  AA[(nVS+1), (nVS+4)] <- mean(-R*beta*t.A)
  AA[(nVS+2), (nVS+3)] <- mean(R*beta*t.A^2)
  AA[(nVS+2), (nVS+4)] <- mean(-R*beta*t.A^2)
  AA[(nVS+3), (nVS+3)] <- mean(y.A*inVal*(nA/nB) )
  AA[(nVS+4), (nVS+4)] <- mean((1-y.A)*inVal*(nA/nB))
  
  for(j in (nXM+2):(nXM+1+nXtilde)){
    AA[(nVS+1), (j-1)] = mean(-zConst*xtilde[,j-nXM-1]*yScore)
  }
  
  for(j in (nXM+2):(nXM+1+nXtilde)){
    AA[(nVS+2), (j-1)] = mean(-zConst*t.A*xtilde[,j-nXM-1]*yScore)
  }
  
  
  
  
  AA_inv <- solve(AA)
  
  ##### the meat B
  BB = matrix(0,matDim, matDim)
  
  getEst = ( inVal*t.A*y.A/w1 ) - (inVal*(1-t.A)*y.A/w2) - estATE
  BB[1,1] = mean(getEst^2) 
  
  for(j in 2:(nXM+1)){
    BB[1, j] = mean(getEst*(t.A-estE)*xM[,j-1]*estB)
  }
  
  for(i in 2:(nXM+1)){
    BB[i,1] = mean(getEst*(t.A-estE)*xM[,i-1]*estB)
  }
  
  for(j in (nXM+2):nVS){
    BB[1, j] = mean(getEst*((inVal-estB)*xtilde[,j-nXM-1]))
  }
  
  for(i in (nXM+2):nVS){
    BB[i, 1] = mean(getEst*((inVal-estB)*xtilde[,i-nXM-1]))
  }
  
  for(i in 2:(nXM+1)){
    for(j in 2:(nXM+1)){
      BB[i, j] = mean( ( (t.A-estE)*xM[,i-1])*estB*(t.A-estE)*xM[,j-1]*estB) 
    }
  }
  
  for(i in (nXM+2):nVS){
    for(j in (nXM+2):nVS){
      BB[i, j] = mean( ( (inVal-estB)*xtilde[,(i-nXM-1)])* ((inVal-estB)*xtilde[,(j-nXM-1)]))
    }
  }
  
  for(i in (nXM+2):nVS){
    for(j in 2:(nXM+1)){
      BB[i, j] = mean( ((inVal-estB)*xtilde[,(i-nXM-1)]) *  ( (t.A-estE)*xM[,j-1]*estB) )  
    }
  }
  

  for(i in 2:(nXM+1)){
    for(j in (nXM+2):nVS){
      BB[i, j] = mean( ( (t.A-estE)*xM[,i-1]*estB)*((inVal-estB)*xtilde[,(j-nXM-1)] ))
    }
  }
  
  
  BB[(nVS+1), 1] <- BB[1, (nVS+1)] <- mean(R*yScore*getEst)
  BB[(nVS+2), 1] <- BB[1, (nVS+2)] <- mean(R*t.A*yScore*getEst)
  BB[(nVS+3),1] <- BB[1, (nVS+3)] <- mean(g1*getEst)
  BB[(nVS+4), 1] <- BB[1, (nVS+4)] <- mean(g2*getEst)
  
  for(j in 3:(nXM+2)){
    BB[(nVS+1), (j-1)] = mean(xM[,j-2]*R*yScore*tScore)
  }
  
  for(i in 3:(nXM+2)){
    BB[(i-1), (nVS+1)] = mean(xM[,i-2]*R*yScore*tScore)
  }
  
  
  for(j in 3:(nXM+2)){
    BB[(nVS+2), (j-1)] =mean(xM[,j-2]*R*t.A*yScore*tScore)
  }
  
  for(i in 3:(nXM+2)){
    BB[(i-1), (nVS+2)] = mean(xM[,i-2]*R*t.A*yScore*tScore)
  }
  
  for(j in 3:(nXM+2)){
    BB[(nVS+3), (j-1)] =mean(g1*tScore*xM[,j-2])
  }
  
  for(i in 3:(nXM+2)){
    BB[(i-1), (nVS+3)] = mean(g1*tScore*xM[,i-2])
  }
  
  for(j in 3:(nXM+2)){
    BB[(nVS+4), (j-1)] =mean(g2*tScore*xM[,j-2])
  }
  
  for(i in 3:(nXM+2)){
    BB[(i-1), (nVS+4)] = mean(g2*tScore*xM[,i-2])
  }
  
  for(j in (nXM+2):(nVS)){
    BB[(nVS+1), j] = mean(xtilde[,j-nXM-1]*V_est*R*yScore)
    
  }
  
  for(i in (nXM+2):(nVS)){
    BB[i, (nVS+1)] = mean(xtilde[,i-nXM-1]*V_est*R*yScore)
  }
  
  for(j in (nXM+2):(nVS)){
    BB[(nVS+2), j] = mean(xtilde[,j-nXM-1]*V_est*R*t.A*yScore)
    
  }
  
  for(i in (nXM+2):(nVS)){
    BB[i, (nVS+2)] = mean(xtilde[,i-nXM-1]*V_est*R*t.A*yScore)
  }
  
  for(j in (nXM+2):(nVS)){
    BB[(nVS+3), j] = mean(xtilde[,j-nXM-1]*V_est*g1)
    
  }
  
  for(i in (nXM+2):(nVS)){
    BB[i, (nVS+3)] = mean(xtilde[,i-nXM-1]*V_est*g1)
  }
  
  for(j in (nXM+2):(nVS)){
    BB[(nVS+4), j] = mean(xtilde[,j-nXM-1]*V_est*g2)
    
  }
  
  for(i in (nXM+2):(nVS)){
    BB[i, (nVS+4)] = mean(xtilde[,i-nXM-1]*V_est*g2)
  }
  
  
  BB[(nVS+1), (nVS+2)] = BB[(nVS+2), (nVS+1)] = mean((R^2)*(yScore^2)*t.A)
  BB[(nVS+1), (nVS+3)] = BB[(nVS+3), (nVS+1)] = mean((R)*yScore*g1)
  BB[(nVS+1), (nVS+4)] = BB[(nVS+4), (nVS+1)] = mean((R)*yScore*g2)
  BB[(nVS+2), (nVS+3)] = BB[(nVS+3), (nVS+2)] = mean(R*t.A*yScore*g1)
  BB[(nVS+2), (nVS+4)] = BB[(nVS+4), (nVS+2)] = mean(R*t.A*yScore*g2)
  BB[(nVS+3), (nVS+4)] = BB[(nVS+4), (nVS+3)] = mean(g1*g2)
  BB[(nVS+1),(nVS+1)] <- mean(R^2 * yScore^2)
  BB[(nVS+2),(nVS+2)] <- mean((R*t.A*yScore)^2)
  BB[(nVS+3), (nVS+3)] <- mean(g1*g1)
  BB[(nVS+4), (nVS+4)] <- mean(g2*g2)
  
  sv_mat = (AA_inv%*%BB)%*%t(AA_inv)/nA
  
  VtauN = sv_mat[(nVS+2),(nVS+2)]
  VtauV2 = sv_mat[1,1]
  Cov2 = sv_mat[(nVS+2),1]
  w = nB/nA
  estNew = w*estATE + (1-w)*tauN_B1_other
  VNew = w^2*VtauV2+(1-w)^2*VtauN+2*w*(1-w)*Cov2
  
  res_df = data.frame(w, estNew, VNew, estATE,VtauV2, tauN_B1_other, VtauN)
  colnames(res_df) = c("w", "estATE", "varATE", "estATE_gold", "varATE_gold", "estATE_silver",
                       "varATE_silver")
  return(res_df)
  
}