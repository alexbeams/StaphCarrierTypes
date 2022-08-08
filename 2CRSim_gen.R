# load transition probabilities #
pCtoC <- function(tm,alpha,gam) alpha/(alpha+gam) + gam/(alpha+gam)*exp(-(alpha+gam)*tm)
pStoC <- function(tm,alpha,gam) alpha/(alpha+gam)*(1-exp(-(alpha+gam)*tm))


getNewState <- function(state, time, alpha, gam){                           #function updates status in 'time' from now based on 'state'
	if(state == 1){
		newState  = rbinom(1,1,pCtoC(time, alpha, gam) )    #stay colonized at time with prob pCtoC
	}else{newState = rbinom(1,1,pStoC(time, alpha, gam))}     #become colonized at time t with prob pStoC
	return(newState)
}

getSubsequentStates = function(state,time,alpha,gam){                     #given a vector of times to test and initial state, generate chain of states
		if(length(time) == 1){
			subsequentState = getNewState( tail(state, n=1), time, alpha, gam)
			return(c(state, subsequentState) )
		}else{
			subsequentState = getNewState( tail(state, n=1), time[1], alpha, gam)
			state = c(state, subsequentState)
			return( getSubsequentStates(state,time[-1],alpha,gam) )
		}
}

getStateSeq = function(time,alpha,gamma1,gamma2,sigma){                                    #given a vector of times to test, generates initial state and chain of subsequent states
	state  = which(rmultinom(1,1, c(sigma*gamma1/(alpha+gamma1), sigma*alpha/(alpha+gamma1), (1-sigma)*gamma2/(alpha+gamma2), (1-sigma)*alpha/(alpha+gamma2) ) ) >0 )
	if(state == 1 | state==2){
		return( getSubsequentStates(state-1, time, alpha, gamma1) + 1 )
	}else{
		return( getSubsequentStates(state-3, time, alpha, gamma2) + 3)
	}
}

getStateSeq.time = function(time) getStateSeq(time, alpha.mle.2cr, gamma1.mle.2cr, gamma2.mle.2cr, sigma.mle.2cr)
StateSequenceData. = lapply(times,getStateSeq.time)

getSimTestDat = function(phi1, psi1, phi2, psi2, StateSequenceData){
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

  results[sample(which(x==4), 
   rbinom(1,length(which(x==4)), 1-phi2) )] <- 3 #false negs

  results[sample(which(x==1), 
   rbinom(1,length(which(x==1)), 1-psi1) )] <- 2 #false pos

  results[sample(which(x==3), 
   rbinom(1,length(which(x==3)), 1-psi2) )] <- 4 #false pos

  #for feeding into simCohorts:
  results[results==3] <- 1
  results[results==4] <- 2

  return(results)
 }

 simTestDat.obs.crud <- lapply(simTestDat.true, getIncorrectTests )
 simTestDat.obs <- list()
 for(i in 1:length(simTestDat.obs.crud)){
  simTestDat.obs[[i]]<- list(test=simTestDat.obs.crud[[i]], times=times[[i]])
 }

 return(list(simTestDat.obs,simTestDat.true))
}

simTestDat<- getSimTestDat(phi1.mle.2cr, psi1.mle.2cr, phi2.mle.2cr, psi2.mle.2cr, StateSequenceData.)

simTestDat.obs <- simTestDat[[1]]
simTestDat.true <- simTestDat[[2]]
