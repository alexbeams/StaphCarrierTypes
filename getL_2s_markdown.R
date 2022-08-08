require(expm)

getLogLikelihood.2s <- function(theta, X){

 alpha1 = exp(theta[4])
 alpha2 = exp(theta[5])
 gamma1 = exp(theta[6]) 
 gamma2 = exp(theta[7])
 psiU = 1/(1+exp(-theta[1]))
 phiC1 = 1/(1+exp(-theta[2]))
 phiC2 = 1/(1+exp(-theta[3]))


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

 pi_1 = gamma1*gamma2/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)
 pi_2 = alpha1*gamma2/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)
 pi_3 = alpha2*gamma1/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)

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

theta1Current.2s <- log(c(.979/(1-.979), .805/(1-.805), .999/.001))
theta2Current.2s <- log(c( .00612, .468*.0531,(1-.468)*.0531, .0173))

dimTheta1.2s <- length(theta1Current.2s)
dimTheta2.2s <- length(theta2Current.2s)

thetaCurrent.2s <- c(theta1Current.2s,theta2Current.2s)
thetaCurrent.2s <- thetaCurrent.2s

argmax <- function(x) which(x==max(x))

getTrueState.j.2s <- function(theta, data){
 
 #use the Viterbi algorithm to produce the most likely underlying state(s)
 observation = data$test 
 timevec = data$times

 alpha1 = exp(theta[4])
 alpha2 = exp(theta[5])
 gamma1 = exp(theta[6]) 
 gamma2 = exp(theta[7])
 psiU = 1/(1+exp(-theta[1]))
 phiC1 = 1/(1+exp(-theta[2]))
 phiC2 = 1/(1+exp(-theta[3]))
 
 
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
 times.unique <- unique(timevec)
 Aji <- lapply(times.unique,Am)

 pi_1 = gamma1*gamma2/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)
 pi_2 = alpha1*gamma2/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)
 pi_3 = alpha2*gamma1/(gamma1*gamma2+alpha1*gamma2+alpha2*gamma1)

 Pi = c(pi_1, pi_2, pi_3)
 
 bb = function(o_t) b(o_t,psiU,phiC1,phiC2)


 testInds <- length(observation)
 T1 <- matrix(nrow=Adim, ncol=testInds)
 T2 <- T1

 
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

getTrueState.2s <- function(data) getTrueState.j.2s(thetaCurrent.2s,data)