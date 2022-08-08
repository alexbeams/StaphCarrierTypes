
getLogLikelihood.ih.importation = function(theta, X){

 alpha= exp(theta[1])
 gam = exp(theta[2])
 phi = 1/(1+exp(-theta[3]) )
 psi = 1/(1+exp(-theta[4]) )
 pi.X <- 1/(1+exp(-theta[5]))
 pi.Y <- 1-pi.X
 
 Amat <- function(tm,alpha,gama){
 Tau <- alpha+gama
 
 probMat = matrix(byrow=T, ncol=2,
  c(
   pi.X + pi.Y*exp(-Tau*tm),pi.X*(1-exp(-Tau*tm)),
   pi.Y*(1-exp(-Tau*tm)),pi.Y + pi.X*exp(-Tau*tm)
   )
  )
  return(probMat)
 }

 b <- function(o_t,phi,psi){
   diag(c(ifelse(o_t==1,psi,1-psi), ifelse(o_t==1,1-phi,phi) ) )
 }

 Adim = 2

 Am = function(tm) Amat(tm,alpha,gam)
   
 #pre-calculate the transition probs for all of the unique test times
 times.unique <- unique(unlist(lapply(X, function(x) x$times)))
 Aji <- lapply(times.unique,Am)

 pi_1 = gam/(alpha+gam)
 pi_2 = alpha/(alpha+gam)

 Pi = c(pi_1, pi_2)
 
 bb = function(o_t) b(o_t,phi,psi)

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

 L = unlist(lapply(X, getLj ) ) #apply to all subjects
 return( sum(log(L)) )
}

thetaCurrent.ih.importation <- log(c(.00489, .00983, .972/.028, .944/.056, .6677989/.3322011))
