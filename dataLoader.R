# load and clean data #

priceData = read.csv("hcwDat05022019.csv") #negative test result -1, MSSA 1, MRSA 2, both 3, no test 0
dat.mrsa = apply(priceData , c(1,2) , function(x) if(x > 0){x <- 2}else{x^2} )  #change -1 -> 1, positive result -> 2

#subset out the MRSA
inds1 = apply(priceData, 2, function(x) ifelse(2 %in% x, 0,1) )
inds2 = apply(priceData, 2, function(x) ifelse(3 %in% x, 0,1) )
inds = intersect( as.vector(which(inds1 == 1) ), as.vector(which(inds2==1)) )

dat.mssa = apply(priceData[,inds] , c(1,2) , function(x) if(x > 0){x <- 2}else{x^2} )  #change -1 -> 1, positive 

# generate times to simulate data points #
testDays.mrsa = apply(dat.mrsa , 2, function(x) which(x != 0))  #list of the days in the study HCW were tested
intervalsBetweenTests.mrsa = lapply(testDays.mrsa, diff)   #list of 4-wk periods between tests while enrolled; 1 = good attendance, 2=missed the next test
times.mrsa = lapply(intervalsBetweenTests.mrsa, function(x) x*4 ) #time in weeks

testDays.mssa = apply(dat.mssa , 2, function(x) which(x != 0))  #list of the days in the study HCW were tested
intervalsBetweenTests.mssa = lapply(testDays.mssa, diff)   #list of 4-wk periods between tests while enrolled; 1 = good attendance, 2=missed the next test
times.mssa = lapply(intervalsBetweenTests.mssa, function(x) x*4 ) #time in weeks


list.mrsa <- list()
for(i in 1:dim(dat.mrsa)[2]){
 list.mrsa[[i]] <- list(test = dat.mrsa[dat.mrsa[,i]>0,i],times = times.mrsa[[i]])
}


list.mssa <- list()
for(i in 1:dim(dat.mssa)[2]){
 list.mssa[[i]] <- list(test = dat.mssa[dat.mssa[,i]>0,i],times = times.mssa[[i]])
}