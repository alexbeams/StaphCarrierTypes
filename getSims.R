for(iter in 1:maxIters){
 source(simNm)
 source('simCohorts.R')

 #use the viterbi algorithm to deduce underlying states, based on test test results
 simTestDat.estimate <- lapply(simTestDat.obs,getTrueState)
 estimates <- unlist(lapply(simTestDat.estimate, function(x) x$test))
 truths <- unlist(lapply(simTestDat.true, function(x) x$test))
 
 #how reliable is the viterbi algorithm?
 psi_model <- sum(estimates[truths==1]==1)/sum(estimates==1)
 phi_model <- sum(estimates[truths==2]==2)/sum(estimates==2)

 #compute the carrier indices (# positive tests)
 ci <- data.frame(table(unlist(lapply(simTestDat.obs, function(x) sum(x$test==2) ))))
 names(ci) <- c('numPos','Freq')
 ci$numPos <- as.numeric(levels(ci$numPos))
 crud <- data.frame('numPos'=0:numTests, 'Freq'=rep(0,numTests+1))
 crud$Freq[crud$numPos %in% intersect(crud$numPos,ci$numPos)] <- ci$Freq
 ci <- crud
 rm(crud)
 
 
 #uncomment these lines to calculate likelihoods/MLE for the simulated data sets
 #loglikelihood <- getLogLikelihood(thetaCurrent,simTestDat.obs)
 #obj <- function(theta) -getLogLikelihood(theta,simTestDat.obs)

 #uncoment this line when you aren't calculating likelihoods for the simulated data sets
 loglikelihood <- 0
 
 #calculate the runlengths of consecutive positive tests and negative tests
 getConsecutivePositives <- function(x){
   x <- x$test
   vls <- rle(x)$values
   lths <- rle(x)$lengths
   ConsecutivePositives <- sort(lths[vls==2])
   return(ConsecutivePositives)
 } 
 
 getConsecutiveNegatives <- function(x){
   x <- x$test
   vls <- rle(x)$values
   lths <- rle(x)$lengths
   ConsecutiveNegatives <- sort(lths[vls==1])
   return(ConsecutiveNegatives)
 } 

 consecPos <- data.frame(table(unlist(lapply(simTestDat.obs, getConsecutivePositives ))))
 names(consecPos) <- c('lengths','Freq')
 consecPos$lengths <- as.numeric(levels(consecPos$lengths))
 crud <- data.frame('lengths'=1:numTests, 'Freq'=rep(0,numTests))
 crud$Freq[crud$lengths %in% intersect(crud$length,consecPos$lengths)] <- consecPos$Freq
 consecPos <- crud
 rm(crud)
 
 consecNeg <- data.frame(table(unlist(lapply(simTestDat.obs, getConsecutiveNegatives ))))
 names(consecNeg) <- c('lengths','Freq')
 consecNeg$lengths <- as.numeric(levels(consecNeg$lengths))
 crud <- data.frame('lengths'=1:numTests, 'Freq'=rep(0,numTests))
 crud$Freq[crud$lengths %in% intersect(crud$length,consecNeg$lengths)] <- consecNeg$Freq
 consecNeg <- crud
 rm(crud)
 
 #uncomment this line to get MLEs from the simulated data sets
 #opt <- optim(thetaCurrent,obj,control=list(maxit=2000))
 
 #uncomment these lines when not estimating MLEs from the simulated data sets
 opt <- list()
 opt$par <- thetaCurrent
 opt$convergence <- 0
 


 dts[iter,] <- c(cohortDatSim, psi_model, phi_model, loglikelihood, opt$convergence,
                 ci$Freq,consecPos$Freq,consecNeg$Freq,opt$par)
}

