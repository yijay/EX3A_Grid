
 ##########################################################################################
 # Beyeisan single structural model by using Grid search
 # EX3A in Punt and Hilborn (2001) BAYES-SA, Bayesian Stock Assessment Methods in Fisheries
 # Yi-Jay Chang
 # 4/21/2016 10:35:19 PM
 ##########################################################################################
 
 setwd("C:\\Users\\Yi-Jay.Chang\\Desktop\\New folder\\")
 DATA = read.table("EX3A_Grid.csv",header=T,sep=",")
 
 n_year = nrow(DATA)
 
 # initial parameters 
 S = 0.7
 R_bar = 315
 q = 0.00065
 sigma = 0.4

 #-------------------------------------------------------------------------------
 # likelihood section
 #------------------------------------------------------------------------------- 
 # population dynamics model
 pop_model <- function(q,R_bar){
  
 Predict_B = NULL
 Predict_B[1] = R_bar/(1-S)
 for(t in 2:n_year){
 Predict_B[t] = max(1,S*Predict_B[t-1] + R_bar - DATA$Catch[t-1])  
 }
 
 Predict_CPUE = q*Predict_B
 SSQ = sum((log(DATA$CPUE)-log(Predict_CPUE))^2)
 neg_likelihood =  exp(-SSQ/(2*sigma^2))
 return(neg_likelihood)
 }
 #------------------------------------------------------------------------------- 
 pop_model(q=0.00065, R_bar=315)
 
 #-------------------------------------------------------------------------------
 # Prior section
 #-------------------------------------------------------------------------------
 # define the grid of each parameter
 q_grid = seq(0.00005,0.002,0.00001)
 R_bar_grid = seq(155,500,2.5)
 
 # prior porbability
 pR_bar_grid = dnorm(R_bar_grid, mean=300, sd=70)
 pq_grid = rep(1,length(q_grid))
 pTheta_matrix = expand.grid(q = pq_grid, R_bar = pR_bar_grid)
  
 # expand the combination of grid 
 Theta_matrix = expand.grid(q = q_grid, R_bar = R_bar_grid)
 #-------------------------------------------------------------------------------
 
 #-------------------------------------------------------------------------------
 # Posterior section
 #------------------------------------------------------------------------------- 
 # negtive likelihood values of each grid
 pDataGivenTheta = mapply(pop_model, Theta_matrix[,1], Theta_matrix[,2])
 
 # Compute the posterior density
 pData = sum(pDataGivenTheta*pTheta_matrix[,1]*pTheta_matrix[,2])
 pThetaGivenData = (pDataGivenTheta*pTheta_matrix[,1]*pTheta_matrix[,2])/pData
 
 # arrange the results
 RESULTs = cbind(Theta_matrix,pThetaGivenData)
 #-------------------------------------------------------------------------------
 
 par(mfrow=c(2,1),mar=c(4,4,2,2))
 
 # calculate the marginal probability
 marginal_q = with(RESULTs, tapply(pThetaGivenData, q, sum))
 plot(q_grid,marginal_q,type="h",lwd=3,xlab="Catchability",ylab=bquote(paste("p(",theta,"|D)")))
 
 # calculate the marginal probability
 marginal_R_bar = with(RESULTs, tapply(pThetaGivenData, R_bar, sum))
 plot(R_bar_grid,marginal_R_bar,type="h",lwd=3,xlab="Average recruitment",ylab=bquote(paste("p(",theta,"|D)"))) 
 
 