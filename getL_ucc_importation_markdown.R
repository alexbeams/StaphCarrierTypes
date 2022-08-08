require(expm)

getLogLikelihood.ucc.importation <- function(theta, X){

 rhoUE = exp(theta[4]) #0.006134599
 rhoEU = exp(theta[5]) #0.028223768
 rhoEC = exp(theta[6]) #0.024990619
 rhoCE = exp(theta[7]) #0.01711582
 psiU = 1/(1+exp(-theta[1]))
 phiE = 1/(1+exp(-theta[2]))
 phiC = 1/(1+exp(-theta[3]))
 
 pc <- 1/(1+exp(-theta[8]))
 pcc <- 1/(1+exp(-theta[9]))
 
 pi_1 = 1-pc
 pi_2 = pc*(1-pcc)
 pi_3 = pc*pcc
 
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

thetaCurrent.ucc.importation <- c( log(c(.979/(1-.979), .805/(1-.805), .999/.001)),
                                   log(c( .00612, .468*.0531,(1-.468)*.0531, .0173)),
                                   log(.34/.66), log(.6/.4))



