
getLogLikelihood.2cr.importation <- function(theta, X){

 alpha= exp(theta[1])
 gamma1 = exp(theta[2])
 gamma2 = gamma1/(1+ exp(-theta[3]) )
 sigma = 1/(1+exp(-theta[4]))
 phi1 = 1/(1+exp(-theta[5]))
 psi1 = 1/(1+exp(-theta[6]))
 phi2 = 1/(1+exp(-theta[7]))
 psi2 = 1/(1+exp(-theta[8]))

 Amat <- function(tm,alpha,gamma1,gamma2){
 tau.1 <- alpha+gamma1
 pi.X1 <- 1/(1+exp(-theta[9])) #gamma1/(alpha+gamma1)
 pi.Y1 <- 1- pi.X1 #alpha/(alpha+gamma1)
 
 tau.2 <- alpha+gamma2
 pi.X2 <- 1/(1+exp(-theta[10])) #gamma2/(alpha+gamma2)
 pi.Y2 <- 1-pi.X2 #alpha/(alpha+gamma2)
 
 probMat = matrix(byrow=T, ncol=4,
   c(
    pi.X1 + pi.Y1*exp(-tau.1*tm),pi.X1*(1-exp(-tau.1*tm)),0,0,
    pi.Y1*(1-exp(-tau.1*tm)),pi.Y1 + pi.X1*exp(-tau.1*tm),0,0,
    0,0,pi.X2 + pi.Y2*exp(-tau.2*tm),pi.X2*(1-exp(-tau.2*tm)),
    0,0,pi.Y2*(1-exp(-tau.2*tm)),pi.Y2 + pi.X2*exp(-tau.2*tm)
    )
   )
   return(probMat)
 }

 Adim=4

 b <- function(o_t,phi1,psi1,phi2,psi2){
  diag(
   c(ifelse(o_t==1,psi1,1-psi1),
     ifelse(o_t==1,1-phi1,phi1),
     ifelse(o_t==1,psi2,1-psi2),
     ifelse(o_t==1,1-phi2,phi2)
   )
  )
 }


 Am = function(tm) Amat(tm, alpha, gamma1, gamma2)

 #pre-calculate the transition probs for all of the unique test times
 times.unique <- unique(unlist(lapply(X, function(x) x$times)))
 Aji <- lapply(times.unique,Am)


 pi_1 = sigma*gamma1/(alpha+gamma1)
 pi_2 = sigma*alpha/(alpha+gamma1)
 pi_3 = (1-sigma)*gamma2/(alpha+gamma2)
 pi_4 = (1-sigma)*alpha/(alpha+gamma2)

 Pi = c(pi_1, pi_2, pi_3, pi_4)
 
 bb = function(o_t) b(o_t,phi1,psi1,phi2,psi2)

 getLj <- function(testData){#uses the forward method
  observation <- testData$test
  timevec <- testData$times
  #to look up the pre-calculated transition probs:
  timevec <- sapply(timevec, function(x) which(times.unique==x))

  bmatrix <- lapply(observation,bb)

  alpha_0 <- bmatrix[[1]] %*% Pi  #initialize the forward method

  Aji_t <- Aji[timevec] #At =  {a_ij}

  prod = Map('%*%',bmatrix[-1],Aji_t) #forward iterations
  l = Reduce('%*%',rev(prod) ) %*% alpha_0 #multiplied by Init. Cond.
  return( sum( l ) )      #Probability(observing a sequence)
 }

 L = unlist(lapply(X, getLj)) #apply(X, 2, getLk )
 return( sum(log(L)) )
}


#alpha,gamma1,gamma2,sigma,phi1,psi1,phi2,psi2
thetaCurrent.2cr.importation <-  c(-4.9690386, -2.5080328, -2.5712596, -0.1944147,  2.8535757,  3.8000464,  2.8535757,  3.8000464, 
                                   log(gamma1.mle.2cr/alpha.mle.2cr), log(gamma2.mle.2cr/alpha.mle.2cr) ) 

