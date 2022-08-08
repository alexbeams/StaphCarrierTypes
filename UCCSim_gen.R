# load transition probabilities #
require(expm)
model <- function(x_0, tm, alpha, tau, p, gamma){
 A = matrix(c(-alpha, (1-p)*tau, 0,
              alpha, -tau, gamma,
              0, p*tau, -gamma), 3, byrow=TRUE)
 expAt = expm(A * tm)
 return(t(expAt %*% x_0))
}



getNewState <- function(state, time, alpha, tau, p, gam){                           #function updates status in 'time' from now based on 'state'
	probs = model(diag(1,3)[state,], time, alpha, tau, p, gam)
	newState = which(rmultinom(1,1,probs)>0)
	return(newState)
}

getSubsequentStates = function(state,time, alpha, tau, p, gam){                     #given a vector of times to test and initial state, generate chain of states
		if(length(time) == 1){
			subsequentState = getNewState( tail(state, n=1), time, alpha, tau, p, gam)
			return(c(state, subsequentState) )
		}else{
			subsequentState = getNewState( tail(state, n=1), time[1], alpha, tau, p, gam)
			state = c(state, subsequentState)
			return( getSubsequentStates(state,time[-1],alpha, tau, p, gam) )
		}
}

getStateSeq = function(time,alpha,tau,p,gam){                                    #given a vector of times to test, generates initial state and chain of subsequent states
	state  = which(rmultinom(1,1, c(tau*gam*(1-p), alpha*gam, alpha*p*tau ) ) >0 )
	return( getSubsequentStates(state, time, alpha, tau, p, gam)  )
}

getStateSeq.time = function(time) getStateSeq(time, alpha.mle.ucc, tau.mle.ucc, p.mle.ucc, gamma.mle.ucc)
StateSequenceData. = lapply(times,getStateSeq.time)

getSimTestDat = function(psi, phi1, phi2, StateSequenceData){
 simTestDat <- StateSequenceData
 simTestDat.true <- list()
 for(i in 1:length(StateSequenceData)){
  simTestDat.true[[i]]<- list(test=simTestDat[[i]], times=times[[i]])
 }


 getIncorrectTests <- function(y){
  x=y$test
  results <- x

  results[sample(which(x==2), 
   rbinom(1,length(which(x==2)), 1-phi1) )] <- 1 #false negs

  results[sample(which(x==3), 
   rbinom(1,length(which(x==3)), 1-phi2) )] <- 1 #false negs

  results[sample(which(x==1), 
   rbinom(1,length(which(x==1)), 1-psi) )] <- 2 #false pos

  #for feeding into simCohorts:
  results[results==3] <- 2

  return(results)
 }

 simTestDat.obs.crud <- lapply(simTestDat.true, getIncorrectTests )
 simTestDat.obs <- list()
 for(i in 1:length(simTestDat.obs.crud)){
  simTestDat.obs[[i]]<- list(test=simTestDat.obs.crud[[i]], times=times[[i]])
 }

 return(list(simTestDat.obs,simTestDat.true))
}

simTestDat <- getSimTestDat(psiU.mle.ucc, phiE.mle.ucc, phiC.mle.ucc, StateSequenceData.)

simTestDat.obs <- simTestDat[[1]]
simTestDat.true <- simTestDat[[2]]