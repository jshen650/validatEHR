#' Run approach from Equations 3 and 4 that uses gold standard outcomes from validation data
#' and all silver standard outcomes
#'
#' @param df Data frame containing exposure status, gold-standard outcomes, silver-standard outcomes, 
#' and all relevant variables for analyses
#' @param B.index Vector that indicates whether an individual was included in the validation data (=1 for included)
#' @param varExp Variable name for exposure
#' @param varGold Variable name for gold-standard outcome
#' @param varSilver Variable name for silver-standard outcome
#' @param varX Variable names for modeling treatment propensity model
#' @param varSelect Variable names for modeling validation sample selection model
#' @param opt Whether to use the optimal weight; TRUE will use the optimal weight 
#' while FALSE will weight by proportion of validation sample size
#'
#' @return Data frame containing weight used, final estimate of ATE, final variance estimate of ATE,
#' estimate of ATE using validation data, variance estimate of ATE using validation data,
#' estimate of ATE using silver standard outcomes, and variance estimate of ATE using silver standard outcomes
#' @export
run_method4 <- function(df,B.index, varExp, varGold, varSilver, varX, varSelect, opt){
  ## gather validation data
  valDat <- df[which(B.index==1),]
  
  
  t.A = df[,varExp]
  y.A = df[,varGold]
  nB <- nrow(valDat)
  B.loc = which(B.index==1)
  
  ## xCols must be columns of data frame pertaining to main effects
  x.A = df[, varX]
  
  
  t.B = valDat[,varExp]
  y.B = valDat[,varGold]
  x.B = x.A[B.loc,]
  yStar = df[,varSilver]
  y.AB=yStar
  
  nA = nrow(df)
  nB = nrow(valDat)
  
  ## estimate propensity scores for treatment
  partDat = data.frame(t.A, x.A)
  estAlpha <- glm(t.A ~ ., data=partDat, family="binomial")
  estE = predict.glm(estAlpha, partDat, type="response")
  
  partDat2 = data.frame(B.index, t.A, df[,c(varSelect)])
  estAlphaB <- glm(B.index ~ . , data=partDat2, family="binomial")
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
  
  
  w1 = estE*estB
  w2 = estB*(1-estE)
  tauV_sum1 = mean (B.index*t.A*y.A/w1)
  tauV_sum2 = mean (B.index*(1-t.A)*y.A/w2)
  estATE = tauV_sum1 - tauV_sum2
  
  tauN_mod = (1/(p11-p10))*( (1/(sum(t.A/(estE) )))*((sum(t.A*y.AB/(estE ) ))) - 
                               (1/(sum((1-t.A)/(1-estE) ))*((sum( (1-t.A)*y.AB/(1-estE ) )))))
  
  ## design matrix relevant for treatment propensity model
  xDesign = model.matrix(estAlpha)
  xM = xDesign
  nXM = ncol(xM) 
  
  ## design matrix relevant for validation sample selection propensity model
  xDesign_Samp = model.matrix(estAlphaB)
  xtilde = xDesign_Samp
  nXtilde = ncol(xtilde) 
  
  ## weighted regression equivalent for estimating Hajek form of the ATE
  R = t.A/estE + (1-t.A)/(1-estE)
  regCoef = coef(glm(yStar ~ t.A, weights=R)) 
  alpha = regCoef[1]
  beta = regCoef[2]
  
  ## estimating equation involving Y*
  yScore = yStar - alpha - beta*(p11-p10)*t.A 
  ## estimating equations involving T
  tScore = (t.A - estE)*estB
  tScore2 = t.A - estE ## estimating equation involving T
  V_est = B.index-estB
  dConst = ( (1-t.A)/(1-estE))*estE - (t.A/estE)*(1-estE)  ## derivative of R but missing X component (will include later)
  g1 =(y.A*yStar - p11*y.A )*B.index*(nA/nB)
  g2 = ((1-y.A)*yStar - p10*(1-y.A))*B.index*(nA/nB)
  
  
  nVS = nXM+nXtilde+1 ## number of cols/rows relevant for estimating variance for validation data
  nAS = 4+nXM ## number of cols/rows relevant for estimating variance for all silver standard outcomes
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
    AA[1, j] = mean( ( (B.index*t.A*y.A*(1-estE))/(estB*(estE)) + (B.index*(1-t.A)*y.A*estE)/(estB*(1-estE)))*xM[,j-1]  )
  }
  
  for(j in (nXM+2):(nVS)){
    AA[1, j] = mean( ( (B.index*t.A*y.A*(1-estB))/(estB*(estE)) + (-1*B.index*(1-t.A)*y.A*(1-estB))/(estB*(1-estE)))*xtilde[,j-nXM-1]  )
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
    AA[(nVS+1), (j+nVS+2)] = mean(-dConst*xM[,j-2]*yScore)
  }
  
  for(j in 3:(nXM+2)){
    AA[(nVS+2), (j+nVS+2)] = mean(-dConst*xM[,j-2]*t.A*yScore)
  }
  
  AA[(nVS+1),(nVS+3)] <- mean(R*beta*t.A)
  AA[(nVS+1), (nVS+4)] <- mean(-R*beta*t.A)
  AA[(nVS+2), (nVS+3)] <- mean(R*beta*t.A^2)
  AA[(nVS+2), (nVS+4)] <- mean(-R*beta*t.A^2)
  AA[(nVS+3), (nVS+3)] <- mean(y.A*B.index*(nA/nB) )
  AA[(nVS+4), (nVS+4)] <- mean((1-y.A)*B.index*(nA/nB))
  
  for(i in 3:(nXM+2)){
    for(j in 3:(nXM+2)){
      AA[(i+nVS+2), (j+nVS+2)] = mean(xM[,i-2]*xM[,j-2]*estE*(1-estE))
    }
  }
  
  AA_inv <- solve(AA)
  
  ##### the meat B
  BB = matrix(0,matDim, matDim)
  
  getEst = ( B.index*t.A*y.A/w1 ) - (B.index*(1-t.A)*y.A/w2) - estATE
  BB[1,1] = mean(getEst^2) 

  for(j in 2:(nXM+1)){
    BB[1, j] = mean(getEst*(t.A-estE)*xM[,j-1]*estB)
  }
  
  for(i in 2:(nXM+1)){
    BB[i,1] = mean(getEst*(t.A-estE)*xM[,i-1]*estB)
  }
  
  for(j in (nXM+2):nVS){
    BB[1, j] = mean(getEst*((B.index-estB)*xtilde[,j-nXM-1]))
  }
  
  for(i in (nXM+2):nVS){
    BB[i, 1] = mean(getEst*((B.index-estB)*xtilde[,i-nXM-1]))
  }
  
  for(i in 2:(nXM+1)){
    for(j in 2:(nXM+1)){
      BB[i, j] = mean( ( (t.A-estE)*xM[,i-1])*estB*(t.A-estE)*xM[,j-1]*estB) 
    }
  }
  
  for(i in (nXM+2):nVS){
    for(j in (nXM+2):nVS){
      BB[i, j] = mean( ( (B.index-estB)*xtilde[,(i-nXM-1)])* ((B.index-estB)*xtilde[,(j-nXM-1)]))
    }
  }
  
  for(i in (nXM+2):nVS){
    for(j in 2:(nXM+1)){
      BB[i, j] = mean( ((B.index-estB)*xtilde[,(i-nXM-1)]) *  ( (t.A-estE)*xM[,j-1]*estB) )  
    }
  }
  
  for(i in 2:(nXM+1)){
    for(j in (nXM+2):nVS){
      BB[i, j] = mean( ( (t.A-estE)*xM[,i-1]*estB)*((B.index-estB)*xtilde[,(j-nXM-1)] ))
    }
  }
  
  BB[(nVS+1), 1] <- BB[1, (nVS+1)] <- mean(R*yScore*getEst)
  BB[(nVS+2), 1] <- BB[1, (nVS+2)] <- mean(R*t.A*yScore*getEst)
  BB[(nVS+3),1] <- BB[1, (nVS+3)] <- mean(g1*getEst)
  BB[(nVS+4), 1] <- BB[1, (nVS+4)] <- mean(g2*getEst)
  
  for(j in 2:(nAS-3)){
    BB[(nVS+1), j] = mean(xM[,j-1]*R*yScore*tScore)
  }
  
  for(i in 2:(nAS-3)){
    BB[i, (nVS+1)] = mean(xM[,i-1]*R*yScore*tScore)
  }
  
  for(j in 2:(nAS-3)){
    BB[(nVS+2), j] =mean(xM[,j-1]*R*t.A*yScore*tScore)
  }
  
  for(i in 2:(nAS-3)){
    BB[i, (nVS+2)] = mean(xM[,i-1]*R*t.A*yScore*tScore)
  }
  
  for(j in 2:(nAS-3)){
    BB[(nVS+3), j] =mean(g1*tScore*xM[,j-1])
  }
  
  for(i in 2:(nAS-3)){
    BB[i, (nVS+3)] = mean(g1*tScore*xM[,i-1])
  }
  
  for(j in 2:(nAS-3)){
    BB[(nVS+4), j] =mean(g2*tScore*xM[,j-1])
  }
  
  for(i in 2:(nAS-3)){
    BB[i, (nVS+4)] = mean(g2*tScore*xM[,i-1])
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
  
  for(j in 2:(nAS-3)){
    BB[1,(nVS+4+j-1)] = mean(xM[,j-1]*tScore2*getEst)
  }
  
  for(i in 2:(nAS-3)){
    BB[(nVS+4+i-1),1] = mean(xM[,i-1]*tScore2*getEst)
  }
  
  for(i in 3:(nXM+2)){
    for(j in 3:(nXM+2)){
      BB[(i+nVS+2), (j+nVS+2)] = mean(xM[,i-2]*tScore2*xM[,j-2]*tScore2)
    }
  }
  
  for(i in 3:(nXM+2)){
    for(j in 3:(nXM+2)){
      BB[(i+nVS+2), (j-1)] = mean(xM[,i-2]*tScore2*xM[,j-2]*tScore)
    }
  }
  
  for(i in 3:(nXM+2)){
    for(j in 3:(nXM+2)){
      BB[(i-1), (j+nVS+2)] = mean(xM[,i-2]*tScore2*xM[,j-2]*tScore)
    }
  }
  
  for(i in 3:(nXM+2)){
    for(j in 3:(nXM+3)){
      BB[(i+nVS+2), (j-1+nXM)] = mean(xM[,i-2]*tScore2*V_est*xtilde[,j-2])
    }
  }
  
  for(i in 3:(nXM+3)){
    for(j in 3:(nXM+2)){
      BB[(i-1+nXM), (j+nVS+2)] = mean(xM[,j-2]*tScore2*V_est*xtilde[,i-2])
    }
  }
  
  for(i in 3:(nXM+2)){
    BB[(i+nVS+2),(nVS+1)] = mean(xM[,i-2] * tScore2*R*yScore)
  }
  
  for(j in 3:(nXM+2)){
    BB[(nVS+1),(j+nVS+2)] = mean(xM[,j-2] * tScore2*R*yScore)
  }
  
  for(i in 3:(nXM+2)){
    BB[(i+nVS+2),(nVS+2)] = mean(xM[,i-2] * tScore2*R*t.A*yScore)
  }
  
  for(j in 3:(nXM+2)){
    BB[(nVS+2),(j+nVS+2)] = mean(xM[,j-2] * tScore2*R*t.A*yScore)
  }
  
  for(i in 3:(nXM+2)){
    BB[(i+nVS+2),(nVS+3)] = mean(xM[,i-2] * tScore2*g1)
  }
  
  for(j in 3:(nXM+2)){
    BB[(nVS+3),(j+nVS+2)] = mean(xM[,j-2] * tScore2*g1)
  }
  
  for(i in 3:(nXM+2)){
    BB[(i+nVS+2),(nVS+4)] = mean(xM[,i-2] * tScore2*g2)
  }
  
  for(j in 3:(nXM+2)){
    BB[(nVS+4),(j+nVS+2)] = mean(xM[,j-2] * tScore2*g2)
  }
  
  sv_mat = (AA_inv%*%BB)%*%t(AA_inv)/nA
  VtauN = sv_mat[nVS+2,nVS+2]
  VtauV2 = sv_mat[1,1]
  Cov2 = sv_mat[nVS+2,1]
  
  if(opt==FALSE){
    w = nB/nA
    estNew = w*estATE + (1-w)*tauN_mod
    VNew = w^2*VtauV2+(1-w)^2*VtauN+2*w*(1-w)*Cov2
    
  }
  
  
  if(opt==TRUE){
    w = ifelse(((VtauN-Cov2)/(VtauN+VtauV2-2*Cov2))*(1-(VtauN-Cov2)/(VtauN+VtauV2-2*Cov2))>=0
                   &VtauN+VtauV2-2*Cov2>0,
                   (VtauN-Cov2)/(VtauN+VtauV2-2*Cov2),as.integer(VtauV2<VtauN))
    estNew = w*estATE + (1-w)*tauN_mod
    VNew = w^2*VtauV2+(1-w)^2*VtauN+2*w*(1-w)*Cov2
  }
  
  estVal = estATE
  estAS = tauN_mod
  
  res_df = data.frame(w, estNew, VNew, estVal,VtauV2, estAS, VtauN)
  colnames(res_df) = c("w", "estATE", "varATE", "estATE_gold", "varATE_gold", "estATE_silver",
                       "varATE_silver")
  return(res_df)
  
}
