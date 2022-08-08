#perfect specificity:

#alpha = 0.0118 
#gam = 0.0220
#phi = 0.944
#psi = 1

#full model:

alpha = alpha.mle.ih #.00489
gam = gam.mle.ih #.00983
phi = phi.mle.ih #.944
psi = psi.mle.ih #.972

PA = alpha/(alpha+gam)

# load transition probabilities #
pCtoC <- function(tm,alpha,gam) alpha/(alpha+gam) + gam/(alpha+gam)*exp(-(alpha+gam)*tm)
pStoC <- function(tm,alpha,gam) alpha/(alpha+gam)*(1-exp(-(alpha+gam)*tm))


getNewState <- function(state, time){                           #function updates status in 'time' from now based on 'state'
	if(state == 1){
		newState  = rbinom(1,1,pCtoC(time, alpha, gam) )    #stay colonized at time with prob pCtoC
	}else{newState = rbinom(1,1,pStoC(time, alpha, gam))}     #become colonized at time t with prob pStoC
	return(newState)
}

getSubsequentStates = function(state,time){                     #given a vector of times to test and initial state, generate chain of states
		if(length(time) == 1){
			subsequentState = getNewState( tail(state, n=1), time)
			return(c(state, subsequentState) )
		}else{
			subsequentState = getNewState( tail(state, n=1), time[1] )
			state = c(state, subsequentState)
			return( getSubsequentStates(state,time[-1]) )
		}
}

getStateSeq = function(time){                                    #given a vector of times to test, generates initial state and chain of subsequent states
	state  = rbinom(1,1,PA)
	return( getSubsequentStates(state, time) )
}

StateSequenceData = lapply(times,getStateSeq)

simTestDat =  lapply(StateSequenceData, function(x) x+1) #1=MSSA negative, 2=MSSA positive



simTestDat.true <- list()
for(i in 1:length(StateSequenceData)){
 simTestDat.true[[i]]<- list(test=simTestDat[[i]], times=times[[i]])
}

getIncorrectTests <- function(y){
 x=y$test
 results <- x
 results[sample(which(x==2), 
  rbinom(1,length(which(x==2)), 1-phi) )] <- 1 #false negs

 results[sample(which(x==1), 
  rbinom(1,length(which(x==1)), 1-psi) )] <- 2 #false pos
 return(results)
}


simTestDat.obs.crud <- lapply(simTestDat.true, getIncorrectTests )
simTestDat.obs <- list()
for(i in 1:length(simTestDat.obs.crud)){
 simTestDat.obs[[i]]<- list(test=simTestDat.obs.crud[[i]], times=times[[i]])
}
