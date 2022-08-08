
getLogLikelihood.2cr <- function(theta, X){

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
 pi.X1 <- gamma1/(alpha+gamma1)
 pi.Y1 <- alpha/(alpha+gamma1)
 
 tau.2 <- alpha+gamma2
 pi.X2 <- gamma2/(alpha+gamma2)
 pi.Y2 <- alpha/(alpha+gamma2)
 
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
thetaCurrent.2cr <-  c(-4.9690386, -2.5080328, -2.5712596, -0.1944147,  2.8535757,  3.8000464,  2.8535757,  3.8000464)
#alpha,gamma1,sigma,phi1,psi1,phi2,psi2
thetaCurrent.7 <-  thetaCurrent.2cr[-3]
#alpha,gamma1,gamma2,sigma,phi1,psi1,phi2,
thetaCurrent.6 <-  thetaCurrent.7[-7]


argmax <- function(x) which(x==max(x))

getTrueState.j.2cr <- function(theta, data){
 
 #use the Viterbi algorithm to produce the most likely underlying state(s)
 observation = data$test
 timevec = data$times 

 alpha= exp(theta[1])
 gamma1 = exp(theta[2])
 gamma2 = gamma1*1/(1+ exp(-theta[3]) )
 sigma = 1/(1+exp(-theta[4]))
 phi1 = 1/(1+exp(-theta[5]))
 psi1 = 1/(1+exp(-theta[6]))
 phi2 = phi1 #1/(1+exp(-theta[7]))
 psi2 = psi1 #1/(1+exp(-theta[8]))

 testInds <- length(observation)
 T1 <- matrix(nrow=Adim, ncol=testInds)
 T2 <- T1

 Amat <- function(tm,alpha,gamma1,gamma2){
   tau.1 <- alpha+gamma1
   pi.X1 <- gamma1/(alpha+gamma1)
   pi.Y1 <- alpha/(alpha+gamma1)
   
   tau.2 <- alpha+gamma2
   pi.X2 <- gamma2/(alpha+gamma2)
   pi.Y2 <- alpha/(alpha+gamma2)
   
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
 times.unique <- unique(timevec)
 Aji <- lapply(times.unique,Am)

 pi_1 = sigma*gamma1/(alpha+gamma1)
 pi_2 = sigma*alpha/(alpha+gamma1)
 pi_3 = (1-sigma)*gamma2/(alpha+gamma2)
 pi_4 = (1-sigma)*alpha/(alpha+gamma2)

 Pi = c(pi_1, pi_2, pi_3, pi_4)
 
 bb = function(o_t) b(o_t,phi1,psi1,phi2,psi2)
 
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

 # to just determine colonized or uncolonized
 #x[x==4] =2
 #x[x==3] =1
 return( list(test=observation,times=data$times) )
}

getTrueState.2cr <- function(data) getTrueState.j.2cr(thetaCurrent.2cr,data)
