simCohorts = lapply(simTestDat.obs, function(x) if(2 %in% x$test){if(1 %in% x$test){'Intermittently'}else{'Persistently'}}else{'never'} )
neverColonizedSim = length(which(simCohorts == 'never'))
intermittentlyColonizedSim = length(which(simCohorts == 'Intermittently'))
persistentlyColonizedSim = length(which(simCohorts == 'Persistently'))
cohortDatSim = c(never = neverColonizedSim, intermittently = intermittentlyColonizedSim,
 persitently = persistentlyColonizedSim)
