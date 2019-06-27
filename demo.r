source("lmcensor.r")
source("funs_censor.r")
library(fda);
library(dplyr)
library(reshape2)
library(data.table)



meanfit <- readRDS("mean.rds")
pc1s <- readRDS("pc1.rds")
pc2s <- readRDS("pc2.rds")
timegrid = seq(1,10,by=0.5)

ssize=100


observed = lapply(1:ssize,function(i){
mu = eval.fd(timegrid, meanfit$pc_fit)
pc1add=eval.fd(timegrid, pc1s$pc_fit)*rnorm(1,0,50)
pc2add= eval.fd(timegrid, pc2s$pc_fit)*rnorm(1,0,25)
err=rnorm(1,mean=0,sd=10)
as.numeric(mu+pc1add+pc2add+err)
})

timepoints = lapply(1:ssize,function(i){
timegrid
})


indexcensored = lapply(1:ssize,function(i){

temp = observed[[i]]
idxcensor = which(temp<45)
#idxcensor = (length(timegrid)-7):(length(timegrid)-1)

idxcensor
})

observed_c = lapply(1:ssize,function(i){
temp = observed[[i]]
#temp[indexcensored[[i]]] = max(temp[indexcensored[[i]]])
temp[indexcensored[[i]]] = 45
temp
})

delta = lapply(1:ssize,function(i){
temp = rep(0,length(timegrid))
temp[indexcensored[[i]]] = 1
temp
})

sigma=1
spline_basis=create.bspline.basis(rangeval=c(1,10),nbasis=6,norder=4)

meanfit = findmean(observed=observed_c, timepoints=timepoints, delta=delta,minit=6,threshold=1e-5)
mu = eval.fd(seq(1,10,by=0.5), meanfit$pc_fit)

observedcenter = lapply(1:length(observed_c),function(i){
	 (observed_c[[i]]-eval.fd(timepoints[[i]], meanfit$pc_fit))[,1]
})


### first fpc
pc1s = first_FPC(rnorm(6),observed=observedcenter, timepoints=timepoints, delta=delta,threshold=1e-3)


##second fpc
previous_beta = list()
previous_beta[[1]] = pc1s$beta
beta3=pc1s$beta
 
pc2s=third_FPC_conditional(rnorm(6), 2, observedcenter, timepoints, delta,betalist =previous_beta , threshold=1e-3)


