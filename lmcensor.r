#install.packages("nloptr")
library(nloptr)


if(FALSE){
a=3
a2=1
b=3


numofcoef = 5
numofsample=10000

x1=rnorm(numofsample)
x2=rnorm(numofsample)
x3=rnorm(numofsample)
x4=rnorm(numofsample)
x5=rnorm(numofsample)

y=rep(0,numofsample)
delta=rep(0,numofsample)

for(i in 1:numofsample){
y[i]=rnorm(1,a*x1[i]+a2*x2[i]+a2*x3[i]+a2*x4[i]+a2*x5[i],3)
if (y[i]< -1){
	y[i]= -1
	delta[i]=1
}
}
X=cbind(x1,x2,x3,x4,x5)

sigma=3




cmat = rbind( c(1,2,0,0,0),
					  c(3,-0.5,1,0,0)
					)
lmcensor(y,X,delta)
lmcensor(y,X,delta,cmat)	





}



findsigma<-function(y,X,delta){

sigma = sqrt( mean(lm.fit(as.matrix(X[which(delta==0),]),y[which(delta==0)])$residual^2))
return(sigma)
}


lmcensor<-function(y,X,delta, cmat=0){
numofcoef  = ncol(X)
numofsample = length(y)

eval_f <- function(beta) {
	temp= log(sigma^2) * numofsample / 2
	for (i in 1:numofsample){

	mu= X[i,] %*% beta

	if(delta[i]==1){
	temp = temp + log(pnorm(y[i]-mu,mean=0,sd=sigma) )
	}else if(delta[i]==0){
	temp=temp + (-1)*(y[i]- mu)^2 / 2 / sigma^2
	}

	}

	return(-temp)
}



eval_grad_f<-function(beta){

	temp=matrix(0,nrow=numofcoef)
	for (i in 1:numofsample){

	if(delta[i]==1){
	temp = temp + c(dnorm(y[i] - X[i,] %*% beta,mean=0,sd=sigma)/pnorm(y[i] - X[i,] %*% beta,mean=0,sd=sigma)) * X[i,]
	}else if(delta[i]==0){
	temp = temp - c(2*(y[i] - X[i,] %*% beta))/ 2 / sigma^2 * X[i,]  
	}

	}
	return(temp)

}
					

			
eval_g_eq <- function( betax ) {
  constr <- cmat %*% betax
  grad <- cmat
  
  return( list( "constraints"=constr, "jacobian"=grad ) )
}




x0 <- rep(0.001,ncol(X))

#sigma = findsigma(y,X,delta)
#if (sigma==0)
#sigma = 37.33
#if (length(which(delta == 0))<=4)sigma = 37.37

local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                    "xtol_rel" = 1.0e-7 )
opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-7,
              "maxeval" = 1000,
              "local_opts" = local_opts)#, "print_level"= 3 ,"check_derivatives"=TRUE )
			  


if(is.vector(cmat) == FALSE && cmat != 0){
res <- nloptr( x0=x0,
eval_f=eval_f,eval_grad_f=eval_grad_f,eval_g_eq =eval_g_eq,
opts=opts )

return(res$solution)
}else{
x0= as.matrix(x0,ncol=1)


res <- nloptr( x0=x0,
eval_f=eval_f,eval_grad_f=eval_grad_f,
opts=opts )


return(res$solution)
}

}

