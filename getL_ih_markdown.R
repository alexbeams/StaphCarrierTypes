
getLogLikelihood.ih = function(theta, X){

 alpha= exp(theta[1])
 gam = exp(theta[2])
 phi = 1/(1+exp(-theta[3]) )
 psi = 1/(1+exp(-theta[4]) )

 Amat <- function(tm,alpha,gama){
 Tau <- alpha+gama
 pi.X <- gama/(alpha+gama)
 pi.Y <- alpha/(alpha+gama)
 
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

thetaCurrent.ih <- log(c(.00489, .00983, .944/.056, .972/.028))

argmax <- function(x) which(x==max(x))

getTrueState.j.ih <- function(theta, data){
 
 #use the Viterbi algorithm to produce the most likely underlying state(s)
 observation = data$test 
 timevec = data$times

 alpha= exp(theta[1])
 gam = exp(theta[2])
 phi = 1/(1+exp(-theta[3]) )
 psi = 1/(1+exp(-theta[4]) )

 testInds <- length(observation)
 T1 <- matrix(nrow=Adim, ncol=testInds)
 T2 <- T1

 
 Amat <- function(tm,alpha,gama){
   Tau <- alpha+gama
   pi.X <- gama/(alpha+gama)
   pi.Y <- alpha/(alpha+gama)
   
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
 times.unique <- unique(timevec)
 Aji <- lapply(times.unique,Am)

 pi_1 = gam/(alpha+gam)
 pi_2 = alpha/(alpha+gam)
 Pi = c(pi_1, pi_2)
 bb = function(o_t) b(o_t,phi,psi)

 bmatrix <- lapply(observation,bb)
 timevec <- sapply(timevec, function(x) which(times.unique==x))

 Aji_t <- Aji[timevec] #At =  {a_ij}

 T1[,1] <-  bmatrix[[1]] %*% Pi
 T2[,1] <- 0

 for(j in 2:length(observation) ){
  Aki = t(Aji_t[[ j-1 ]])
  crud <- T1[,j-1] * Aki %*% bmatrix[[j-1]]
  T1[,j] <- apply( crud, 2, max)
  T2[,j] <- apply(crud, 2, argmax )
 }

 Z <- vector(length=length(observation))

 Z[j] = argmax(T1[,j])
 
 for(j in length(observation):2){
  Z[j-1] = T2[ Z[j], j ]
 }

 observation = Z
 return( list(test=observation,times=data$times) )
}

getTrueState.ih <- function(data) getTrueState.j.ih(thetaCurrent.ih,data)
