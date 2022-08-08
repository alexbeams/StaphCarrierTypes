# This is code for an ensemble Metropolis Hastings sampler.

#It uses the "stretch move".
#At each timestep, the ensemble moves as follows:
#Loop through ensemble. Given Xk, select another ensemble member Xj at random
#Move Xk along that direction.
#Xk_proposed=Xj + Z*(Xk-Xj), Z from some distribution
#Accept or reject based on the likelihood value.

dimTheta = length(thetaCurrent)

ensemble <- array( dim=c(dimTheta+1,dimTime,dimEnsemble) )

#initialize the ensemble - Gaussian around the initial point
for(i in 1:dimEnsemble){
 ensemble[1:dimTheta,1,i] <- rnorm(dimTheta,thetaCurrent,1)
}
#calculate the initial likelihoods
for(i in 1:dimEnsemble){
 ensemble[dimTheta+1,1,i] <- likelihood(ensemble[1:dimTheta,1,i])
}



#generate all the independent stretches initially
y <- runif(2*dimTime*dimEnsemble)
a <- 2 #tuning parameter; stretches at most by a factor of a or compresses down by a factor of 1/a
Z <- (y*(a-1)+1)^2/a #normalized the g(z) from Goodman and Weare (2010), generated unif(0,1)'s and inverted to get stretches
#Z <- rnorm(2*dimTime*dimEnsemble,1,1)
Zind <- 1

for(i in 1:(dimTime-1)){
 for(k in 1:dimEnsemble){
  Xk = ensemble[1:dimTheta,i,k]
  likelihood_Current = ensemble[1+dimTheta,i,k]
  j = sample(setdiff(1:dimEnsemble,k),1) 
  Xj = ensemble[1:dimTheta,i,j]
  Y = Xj + (Xk-Xj)*Z[Zind]
  likelihood_Proposed = likelihood(Y)
  if(is.na(likelihood_Proposed)){likelihood_Proposed=-Inf}
  proposal_Probability = Z[Zind]*exp( likelihood_Proposed - likelihood_Current ) #metropolis-hastings acceptance probability from Goodman and Weare (2010)
  if(runif(1)<proposal_Probability){
   ensemble[,i+1,k] <- c(Y,likelihood_Proposed)
  }else{
   ensemble[,i+1,k] <- c(Xk,likelihood_Current)
  }
  Zind <- Zind + 1
 }
}



