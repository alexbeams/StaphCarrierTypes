require(expm)

getLogLikelihood.2s.importation <- function(theta, X){

 alpha1 = exp(theta[4])
 alpha2 = exp(theta[5])
 gamma1 = exp(theta[6]) 
 gamma2 = exp(theta[7])
 psiU = 1/(1+exp(-theta[1]))
 phiC1 = 1/(1+exp(-theta[2]))
 phiC2 = 1/(1+exp(-theta[3]))
 
 pc1 = 1/(1+exp(-theta[8]))
 pc2 = 1/(1+exp(-theta[9]))
 
 
 Adim = 3

 model <- function(x_0, t, alpha1, alpha2, gamma1, gamma2){
  A = matrix(c(-alpha1-alpha2,gamma1,gamma2,
               alpha1, -gamma1,0,
               alpha2,0,-gamma2),3,byrow=TRUE)
  expAt = expm(A * t)
  return(t(expAt %*% x_0))
 }

 A <- function(t,alpha1, alpha2, gamma1, gamma2){
  c(model(c(1,0,0),t,alpha1, alpha2, gamma1, gamma2),
    model(c(0,1,0),t,alpha1, alpha2, gamma1, gamma2),
    model(c(0,0,1),t,alpha1, alpha2, gamma1, gamma2) )
 }

 Amat <- function(t,alpha1, alpha2, gamma1, gamma2){
   matrix(A(t,alpha1, alpha2, gamma1, gamma2),Adim,byrow=F)
 }

 b <- function(o_t,psiU,phiC1,phiC2){
  diag(
   c(ifelse(o_t==1,psiU,1-psiU),
     ifelse(o_t==1,1-phiC1,phiC1),
     ifelse(o_t==1,1-phiC2, phiC2)
   )
  )
 }

 
 Am = function(t) Amat(t,alpha1, alpha2, gamma1, gamma2)

 #pre-calculate the transition probs for all of the unique test times
 times.unique <- unique(unlist(lapply(X, function(x) x$times)))
 Aji <- lapply(times.unique,Am)

 pi_1 = 1-pc1-pc2 #gamma1*gamma2/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)
 pi_2 = pc1 #alpha1*gamma2/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)
 pi_3 = pc2 #alpha2*gamma1/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)

 Pi = c(pi_1, pi_2, pi_3)
 
 bb = function(o_t) b(o_t,psiU,phiC1,phiC2)

 getLj <- function(testData){#uses the forward method
  observation <- testData$test
  timevec <- testData$times
  #to look up the pre-calculated transition probs:
  timevec <- sapply(timevec, function(x) which(times.unique==x))

  bmatrix <- lapply(observation,bb)

  a_0 <- bmatrix[[1]] %*% Pi  #initialize the forward method

  Aji_t <- Aji[timevec] #At =  {a_ij}

  prod = Map('%*%',bmatrix[-1],Aji_t) #forward iterations
  l = Reduce('%*%',rev(prod) ) %*% a_0 #multiplied by Init. Cond.
  return( sum( l ) )      #Probability(observing a sequence)
 }


 L = unlist(lapply(X, getLj ) ) #apply to all subjects

 logLikelihood = sum(log(L)) 
 return(logLikelihood)
}

thetaCurrent.2s.importation <- log(c(.975/(1-.975), .891/(1-.891), .999/.001, .00482, .00499,.0157, .00249))
pivals <- c(.34-.294,.294)
pivals <- log(pivals/(1-pivals))

thetaCurrent.2s.importation <- c(thetaCurrent.2s.importation,pivals)
