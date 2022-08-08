require(expm)

getLogLikelihood.ucc <- function(theta, X){

 rhoUE = exp(theta[4]) #0.006134599
 rhoEU = exp(theta[5]) #0.028223768
 rhoEC = exp(theta[6]) #0.024990619
 rhoCE = exp(theta[7]) #0.01711582
 psiU = 1/(1+exp(-theta[1]))
 phiE = 1/(1+exp(-theta[2]))
 phiC = 1/(1+exp(-theta[3]))


 Adim = 3

 model <- function(x_0, t, rhoUE, rhoEU, rhoEC, rhoCE){
  A = matrix(c(-rhoUE, rhoEU, 0,
               rhoUE, -(rhoEU+rhoEC), rhoCE,
               0, rhoEC, -rhoCE), 3, byrow=TRUE)
  expAt = expm(A * t)
  return(t(expAt %*% x_0))
 }

 A <- function(t,rhoUE,rhoEU,rhoEC,rhoCE){
  c(model(c(1,0,0),t,rhoUE,rhoEU,rhoEC,rhoCE),
    model(c(0,1,0),t,rhoUE,rhoEU,rhoEC,rhoCE),
    model(c(0,0,1),t,rhoUE,rhoEU,rhoEC,rhoCE) )
 }

 Amat <- function(t,rhoUE,rhoEU,rhoEC,rhoCE){
   matrix(A(t,rhoUE,rhoEU,rhoEC,rhoCE),Adim,byrow=F)
 }

 b <- function(o_t,psiU,phiE,phiC){
  diag(
   c(ifelse(o_t==1,psiU,1-psiU),
     ifelse(o_t==1,1-phiE,phiE),
     ifelse(o_t==1,1-phiC, phiC)
   )
  )
 }

 
 Am = function(t) Amat(t,rhoUE,rhoEU,rhoEC,rhoCE)

 #pre-calculate the transition probs for all of the unique test times
 times.unique <- unique(unlist(lapply(X, function(x) x$times)))
 Aji <- lapply(times.unique,Am)

 pi_1 = rhoEU*rhoCE/(rhoEU*rhoCE+rhoUE*rhoCE+rhoUE*rhoEC)
 pi_2 = rhoUE*rhoCE/(rhoEU*rhoCE+rhoUE*rhoCE+rhoUE*rhoEC)
 pi_3 = rhoUE*rhoEC/(rhoEU*rhoCE+rhoUE*rhoCE+rhoUE*rhoEC)

 Pi = c(pi_1, pi_2, pi_3)
 
 bb = function(o_t) b(o_t,psiU,phiE,phiC)

 getLj <- function(testData){#uses the forward method
  observation <- testData$test
  timevec <- testData$times
  #to look up the pre-calculated transition probs:
  timevec <- sapply(timevec, function(x) which(times.unique==x))

  bmatrix <- lapply(observation,bb)

  rhoUE_0 <- bmatrix[[1]] %*% Pi  #initialize the forward method

  Aji_t <- Aji[timevec] #At =  {a_ij}

  prod = Map('%*%',bmatrix[-1],Aji_t) #forward iterations
  l = Reduce('%*%',rev(prod) ) %*% rhoUE_0 #multiplied by Init. Cond.
  return( sum( l ) )      #Probability(observing a sequence)
 }


 L = unlist(lapply(X, getLj ) ) #apply to all subjects

 logLikelihood = sum(log(L)) 
 return(logLikelihood)
}

theta1Current.ucc <- log(c(.979/(1-.979), .805/(1-.805), .999/.001))
theta2Current.ucc <- log(c( .00612, .468*.0531,(1-.468)*.0531, .0173))

dimTheta1.ucc <- length(theta1Current.ucc)
dimTheta2.ucc <- length(theta2Current.ucc)

thetaCurrent.ucc <- c(theta1Current.ucc,theta2Current.ucc)
thetaCurrent.ucc <- thetaCurrent.ucc

argmax <- function(x) which(x==max(x))

getTrueState.j.ucc <- function(theta, data){
 
 #use the Viterbi algorithm to produce the most likely underlying state(s)
 observation = data$test 
 timevec = data$times

 rhoUE = exp(theta[4]) #0.006134599 #exp(theta[1])
 rhoEU = exp(theta[5]) #0.028223768 #exp(theta[2])
 rhoEC = exp(theta[6]) #0.024990619 #exp(theta[3])
 rhoCE = exp(theta[7]) #0.01711582 
 psiU = 1/(1+exp(-theta[1]))
 phiE = 1/(1+exp(-theta[2]))
 phiC = 1/(1+exp(-theta[3]))
 
 
 Adim = 3
 
 model <- function(x_0, t, rhoUE, rhoEU, rhoEC, rhoCE){
   A = matrix(c(-rhoUE, rhoEU, 0,
                rhoUE, -(rhoEU+rhoEC), rhoCE,
                0, rhoEC, -rhoCE), 3, byrow=TRUE)
   expAt = expm(A * t)
   return(t(expAt %*% x_0))
 }
 
 A <- function(t,rhoUE,rhoEU,rhoEC,rhoCE){
   c(model(c(1,0,0),t,rhoUE,rhoEU,rhoEC,rhoCE),
     model(c(0,1,0),t,rhoUE,rhoEU,rhoEC,rhoCE),
     model(c(0,0,1),t,rhoUE,rhoEU,rhoEC,rhoCE) )
 }
 
 Amat <- function(t,rhoUE,rhoEU,rhoEC,rhoCE){
   matrix(A(t,rhoUE,rhoEU,rhoEC,rhoCE),Adim,byrow=F)
 }
 
 b <- function(o_t,psiU,phiE,phiC){
   diag(
     c(ifelse(o_t==1,psiU,1-psiU),
       ifelse(o_t==1,1-phiE,phiE),
       ifelse(o_t==1,1-phiC, phiC)
     )
   )
 }
 
 
 Am = function(t) Amat(t,rhoUE,rhoEU,rhoEC,rhoCE)

 #pre-calculate the transition probs for all of the unique test times
 times.unique <- unique(timevec)
 Aji <- lapply(times.unique,Am)

 pi_1 = rhoEU*rhoCE/(rhoEU*rhoCE+rhoUE*rhoCE+rhoUE*rhoEC)
 pi_2 = rhoUE*rhoCE/(rhoEU*rhoCE+rhoUE*rhoCE+rhoUE*rhoEC)
 pi_3 = rhoUE*rhoEC/(rhoEU*rhoCE+rhoUE*rhoCE+rhoUE*rhoEC)

 Pi = c(pi_1, pi_2, pi_3)
 
 bb = function(o_t) b(o_t,psiU,phiE,phiC)


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

getTrueState.ucc <- function(data) getTrueState.j.ucc(thetaCurrent.ucc,data)