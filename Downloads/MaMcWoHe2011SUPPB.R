#########################################################################
#########################################################################
#### SupplementaryMaterialsModelFitting.R	     
#### David S. Matteson						
#### 2010 - 11 - 18							
#### Annals of Applied Statistics 			
#### Forecasting Emergency Medical Service Call Arrival Rates 
#### By D.S. Matteson, M.W. McLean, D.B. Woodard, S.G. Henderson 
#########################################################################
#########################################################################
#### Implemented using R version 2.10.1
#### The period is assumed to be 1 hour and 
#### no missing observations within a day.
#### User must supply the following data:
#### y : observed counts per hour
#### hour : 24 hour of the day 1-24 (same length as y)
#### day : day of the year 1-365 (same length as y)
#### dofw : day of the week 1-7 (same length as y)
#### week : week of the year 1-53 (same length as y)
#########################################################################
#########################################################################

# rm(list = ls()) # clear workspace

library(mgcv) # Main library for Multiple Smoothing Parameter Estimation
# library(chron) # Helpful for manipulating chronological objects
# library(lattice) # Trellis graphics for R

#########################################################################
#### Simulate 12 weeks of hourly observations T = 24*7*12 = 2016
#### with a few of the stylized facts observed in the call 
#### arrival data which is not available for distribution
#########################################################################
T = 24*7*12
hour = rep(1:24, 7*12)
day = rep(1:(7*12), each = 24)
dofw = rep(1:7, each = 24, 12)
week = rep(1:12, each = 24*7)

err = arima.sim(n = T, list(ar = 0.95, ma = -0.8), rand.gen = function(n, ...) rnorm(n, 0, 1))
rate = 20 + 10*sin(2*pi*hour/24) + 2*sin(2*pi*week/26) + err
y = rpois(T,rate)

ts.plot(y)
lines(rate, col = 2)
#########################################################################
#########################################################################


N = 24 # periods per cycle, i.e. hours per day
D = length(y)/(N)  ; D # number of "days"
ND = length(y) # total number of observations

dofwindex = as.factor(dofw)
weekindex = as.factor(week)
Y = t(matrix(y,N,D))
DoW = t(matrix(dofw,N,D))
WEEK = t(matrix(week,N,D))

########################################
########################################
#### The main estimation algorithm for 
#### fitting the K-factor model 
#### using constrainsts and 
#### smothing splines
########################################
########################################
K.max = 4
muhat = matrix(0,N*D,K.max)
Max.iter = 40
# Set exit level for relative reduction in deviance
dev.exit = 0.0001

for(k in 1:K.max){
########################################
# Initialization:
dim(Y);	min(Y);	min(ifelse(Y==0,0.01,Y))
gY = log(ifelse(Y==0,0.01,Y))
gYsvd = svd(gY)

# coefs
B.new = matrix(0,D,k)
for(i in 1:k){  B.new[,i] = gYsvd$d[i]*gYsvd$u[,i]  }

#factors
F.new = matrix(0,N,k)
for(i in 1:k){  F.new[,i] = gYsvd$v[,i]  }

########################################
# Begin iterative algorithm
iter = 1
dev.new = Inf

while(iter < Max.iter){
	 tic = proc.time()[3] 

dev.old = dev.new
F.old = F.new
B.old = B.new

########################################
X.temp = matrix(0,ND,k)
for(kk in 1:k){ X.temp[,kk] = rep(F.old[,kk],D) }

xnam <- paste(paste("s(as.numeric(weekindex),by = X.temp[,", 1:k, sep=""), "],bs='cc')", sep = "")
fmla <- as.formula(paste("y ~ -1 + X.temp:dofwindex +", paste(xnam, collapse= "+")))

fit6 = gam(fmla, family = poisson)
B.tempD = matrix(as.vector(fit6$coefficients[1:(7*k)]), 7, k, byrow=TRUE)

# Extracting fitted values
n = 52 # number of weeks in the year
S = NULL  

for(s in 1:k){  
  raw <- fit6$model[fit6$smooth[[s]]$term]
  xx <- seq(min(raw), max(raw), length = n)
  by <- rep(1, n)
  dat <- data.frame(x = xx, by = by)
	names(dat) <- c(fit6$smooth[[s]]$term, fit6$smooth[[s]]$by)
	Xmat <- PredictMat(fit6$smooth[[s]], dat)
	first <- fit6$smooth[[s]]$first.para
	last <- fit6$smooth[[s]]$last.para
	p <- fit6$coefficients[first:last]
	S.temp <- Xmat %*% p
	S = c(S,S.temp)
}

B.tempW = matrix(as.vector(S), 52, k, byrow=FALSE)

########################################
B.temp = matrix(0, D, k, byrow=TRUE)
rm(fit6) 

# 7 days in the week
for(j in 1:7){  
	for(ell in 1:k){
 		B.temp[which(DoW[1:D,1] == levels(dofwindex)[j]),ell] = B.temp[which(DoW[1:D,1] == levels(dofwindex)[j]),ell] + as.numeric(B.tempD[j,ell])
	}  
}

for(j in 2:53){  
	for(ell in 1:k){
 		B.temp[which(WEEK[1:D,1] == levels(weekindex)[j]),ell] = B.temp[which(WEEK[1:D,1] == levels(weekindex)[j]),ell] + as.numeric(B.tempW[(j-1),ell])
	}  
}

Z.temp = matrix(0,ND,k)
for(kk in 1:k){ Z.temp[,kk] = rep(B.temp[,kk],each=N) }

########################################

znam <- paste(paste("s(hour,by = Z.temp[,", 1:k, sep=""), "])", sep = "")
fmla <- as.formula(paste("y ~ -1 +", paste(znam, collapse= "+")))

fit4 = gam(fmla, family = poisson)

# Extracting fitted values
n = 24 # 24 hours per day
S = NULL  

for(s in 1:k){  
  raw <- fit4$model[fit4$smooth[[s]]$term]
  xx <- seq(min(raw), max(raw), length = n)
  by <- rep(1, n)
  dat <- data.frame(x = xx, by = by)
	names(dat) <- c(fit4$smooth[[s]]$term, fit4$smooth[[s]]$by)
	Xmat <- PredictMat(fit4$smooth[[s]], dat)
	first <- fit4$smooth[[s]]$first.para
	last <- fit4$smooth[[s]]$last.para
	p <- fit4$coefficients[first:last]
	S.temp <- Xmat %*% p
	S = c(S,S.temp)
}

F.temp = matrix(as.vector(S),N,k, byrow=FALSE)

# Save most recent fit before orthogonalization
fit.final = fit4
rm(fit4)

########################################
# Orthogonalize Factors F
G.temp = B.temp %*% t(F.temp)
Gsvd = svd(G.temp)

B.new = matrix(0, D, k)
	for(i in 1:k){  B.new[,i] = Gsvd$d[i]*Gsvd$u[,i]  }
F.new = matrix(0,N,k)
	for(i in 1:k){  F.new[,i] = Gsvd$v[,i]  }

dev.new = fit.final$deviance
if(0 < dev.old - dev.new & dev.old - dev.new < dev.exit) iter = Inf

 toc = proc.time()[3] - tic ; toc 
# optional print statements 
print(c(iter, toc/60))
iter = iter + 1
print(fit.final$deviance)
flush.console()

}

muhat[,k] = fit.final$fitted

# optional print statements
print(k)
#print(summary(F.old - F.new)) ; 
#print(max(abs(F.old - F.new)))
#print(summary(B.old - B.new)) ; 
#print(max(abs(B.old - B.new)))
#print(round(crossprod(F.old,F.new),4))
#print(diag(round(crossprod(F.old,F.new),4)))

}

# fitted values in vector form (same length as y) for k = K.max
index = seq(1,24,by=1)
mu.hat = numeric(ND)

for(i in 1:D){
	mu.hat[((i-1)*N+1):(i*N)] = as.vector(exp(F.new[index,]%*%B.new[i,]))
}

# multiplicative residual
Et = y/mu.hat

# a couple residual plots
par(mfrow=c(3,1))
ts.plot(Et[1:1000]) ; abline(h = 1)
acf(Et, ylim=c(-0.01,0.1), lag.max = 96*2+16)
abline(v = c(96.6, 192.6), lty = 2, col = 2)
acf(Et, ylim=c(-0.02,0.1), lag.max = 96*2+16, type = "partial")
abline(v = c(96.6, 192.6), lty = 2, col = 2)

###############################################################
# if some missing days were removed use 'misshour' below
# to reinitilize the conditional likelihoods below
misshour = c(1, ifelse(diff(day) > 1 , 1, 0))
sum(misshour)

###############################################################
######## For conditional ML estiation of Equation 6     #######
######## Int-GARCH(1,1)                                 #######
###############################################################
"condPoissonInt11" = function(parms, y, mu, misshour, llik){
	alpha = parms[1]
	beta  = parms[2]
	omega = 1 - alpha - beta
	N = length(y)
	lambda = numeric(N) 
	eta = numeric(N)
	epsilon = y/mu
	eta[1] = 1
	lambda[1] = 1	 
	loglik = 0 # -sum(lfactorial(y))
	for(i in 2:N){
		eta[i] = omega + alpha*epsilon[(i-1)] + beta*ifelse(misshour[(i-1)] == 1, 1, eta[(i-1)])
		lambda[i] = mu[i]*eta[i]
#		if(lambda[i] <= 0){print(c(i,lambda[i],alpha,beta))}
		temp = -lambda[i] + y[i]*log(lambda[i]) - (lfactorial(y[i]))
		loglik = loglik + ifelse(misshour[i] == 1, 0, temp)
	}
	if(llik==TRUE){-loglik}
	else{eta}
}
###############################################################
theta.0 = c(0.05, 0.5) 
           
condPoissonInt11(parms = theta.0, y = y, mu = mu.hat, misshour = misshour, llik = TRUE)

outInt11 = optim(par=theta.0, fn = condPoissonInt11, y = y, mu=mu.hat, llik=TRUE, 
	misshour= misshour, method = "L-BFGS-B", lower = c(0.000,0.000), 
	upper = c(0.2,0.9), hessian=T, control = 
	list(trace = TRUE, ndeps = rep.int(0.000001, 2), 
		  maxit = 200L, factr = 1e+31, pgtol = 0))

# parameter estimates
igparInt11 = outInt11$par ; igparInt11 ; 1 - sum(igparInt11)

# approximate SEs
igseInt11 = sqrt(diag(solve(outInt11$hessian))) ; igseInt11

# CIIR
etaInt11 = condPoissonInt11(parms = outInt11$par,y= y, mu=mu.hat, misshour= misshour, llik=FALSE)

# Mltiplicative residuals
e = y/mu.hat

# Fitted values
lambdaInt11 = mu.hat*etaInt11

# Standardized multiplicative residuals
epsilonInt11 = y/lambdaInt11

# a couple standardized residual plots
par(mfrow=c(2,1))
acf(epsilonInt11, ylim=c(-0.05,0.1), lag.max = 96*14)
abline(v = c(seq(96,96*14,96)+.6), lty = 2, col = 2)
acf(epsilonInt11, ylim=c(-0.05,0.1), lag.max = 96*2+16)
abline(v = c(96.6, 192.6), lty = 2, col = 2)

###############################################################
######## For conditional ML estiation of Equation 7     #######
######## Int-ExpAR(1,1)  Note: this is quite slow       #######
###############################################################
"IntExpAR11" = function(parms, y, mu, misshour, llik){
	alpha = parms[1]
	beta  = parms[2]
	omega = parms[3]
	gamma = parms[4]
	N = length(y)
	lambda = numeric(N) 
	eta = numeric(N)
	epsilon = y/mu
	eta[1] = 1
	lambda[1] = 1	 
	loglik = 0 # -sum(lfactorial(y))
	for(i in 2:N){
#		eta[i] = omega + alpha*epsilon[(i-1)] + beta*ifelse(misshour[(i-1)] == 1, 1, eta[(i-1)])
       etatm1 = ifelse(misshour[(i-1)] == 1, 1, eta[(i-1)])
		eta[i] = (omega + beta*exp(-gamma*etatm1^2))*etatm1 + alpha*epsilon[(i-1)]
		lambda[i] = mu[i]*eta[i]
	if(lambda[i] <= 0){print(c(i,lambda[i],alpha,beta,omega,gamma))}
		temp = -lambda[i] + y[i]*log(lambda[i]) - (lfactorial(y[i]))
		loglik = loglik + ifelse(misshour[i] == 1, 0, temp)
	}
	if(llik==TRUE){-loglik}
	else{eta}
}
###############################################################
theta.0 = c(0.113, 0.0539, 0.8427, 1)

IntExpAR11(parms = theta.0, y = y, mu = mu.hat, misshour = misshour, llik = TRUE)

outExpAR11 = nlminb(start = theta.0, objective = IntExpAR11,  y = y, mu=mu.hat, llik=TRUE, misshour= misshour, lower = c(0.0,0.0,0.0,0.0), upper = c(0.1,Inf,0.95,Inf))

outExpAR11 = optim(par=theta.0, fn = IntExpAR11, y = y, mu=mu.hat, llik=TRUE, 
	misshour= misshour, method = "L-BFGS-B", lower = c(0.03,0.15,0.60,0.100), 
	upper = c(0.1,Inf,0.95,Inf), hessian=T, control = 
	list(trace = 5, ndeps = rep.int(0.000001, 4), 
		  maxit = 10L, factr = 1e+31, pgtol = 0))

# parameter estimates		  
igparExpAR11 = outExpAR11$par ; igparExpAR11 

# approximate SEs
igseExpAR11 = sqrt(diag(solve(outExpAR11$hessian))) ; igseExpAR11

# CIIR
etaExpAR11 = IntExpAR11(parms = outExpAR11$par, y=y, mu=mu.hat, misshour= misshour, llik=FALSE)

# Fitted values
lambdaExpAR11 = mu.hat*etaExpAR11

# Standardized multiplicative residuals
epsilonExpAR11 = y/lambdaExpAR11

###############################################################
######## For conditional ML estiation of Equation 8     #######
######## Int-ThreshGARCH(1,1)  Note: this is quite slow #######                           
###############################################################
"IntThresh11" = function(parms, y, mu, misshour, llik){
	alpha = parms[1]
	beta  = parms[2]
	omega = parms[3]
	gamma = parms[4]
	delta = parms[5]
	C = 1.15 #parms[6]
	N = length(y)
	lambda = numeric(N) 
	eta = numeric(N)
	epsilon = y/mu
	eta[1] = 1
	lambda[1] = 1	 
	loglik = 0 # -sum(lfactorial(y))
	for(i in 2:N){
#		eta[i] = omega + alpha*epsilon[(i-1)] + beta*ifelse(misshour[(i-1)] == 1, 1, eta[(i-1)])
       etatm1 = ifelse(misshour[(i-1)] == 1, 1, eta[(i-1)])
       epsBig = ifelse( 1/C < epsilon[(i-1)] && epsilon[(i-1)] < C, 0, 1)
		eta[i] = omega + (alpha + gamma*epsBig)*epsilon[(i-1)] + (beta + delta*epsBig)* etatm1
		lambda[i] = mu[i]*eta[i]
	if(lambda[i] <= 0){print(c(i,lambda[i],alpha,beta,omega,gamma,delta,C))}
		temp = -lambda[i] + y[i]*log(lambda[i]) - (lfactorial(y[i]))
		loglik = loglik + ifelse(misshour[i] == 1, 0, temp)
	}
	if(llik==TRUE){-loglik}
	else{eta}
}
###############################################################

theta.0 =  c(0.027, 0.858, 0.115, 0.025, -0.028)

IntThresh11(parms = theta.0, y = y, mu = mu.hat, misshour = misshour, llik = TRUE)

outThresh11 = nlminb(start = theta.0, objective = IntThresh11, y = y, mu=mu.hat, llik=TRUE, misshour= misshour, 			   lower = c(0.0, 0.80, 0.0001, 0.0,-0.10), upper = c(0.1,0.9,1,0.08,0.0))

outThresh11 = optim(par=theta.0, fn = IntThresh11, y = y, mu=mu.hat, llik=TRUE, 
	misshour= misshour, method = "L-BFGS-B", lower = c(0.0, 0.80, 0.0001, 0.0,-0.10), 
	upper = c(0.1,0.9,Inf,0.05,0.05), hessian=T, control = 
	list(trace = 5, ndeps = rep.int(0.000001, 5), 
		  maxit = 10L, factr = 1e+31, pgtol = 0))

# parameter estimates
igparThresh11 = outThresh11$par ; igparThresh11

# approximate SEs
igseThresh11 = sqrt(diag(solve(outThresh11$hessian))) ; igseThresh11

# CIIR
etaThresh11 = IntThresh11(parms = outThresh11$par, y=y, mu=mu.hat, misshour= misshour, llik=FALSE)

# Fitted values
lambdaThresh11 = mu.hat*etaThresh11

# Standardized multiplicative residuals
epsilonThresh11 = y/lambdaThresh11

###############################################################
######## For conditional ML estiation of Equation 9     #######
########  Int-rsmGARCH(1,1)  Note: this is quite slow   #######
###############################################################
"IntRsm11" = function(parms, y, mu, t1,t2, misshour, llik, hour){
	#if wish to estimate omega in both regimes
	#parms vector contains: c(omega1,alpha1,beta1,omega2,alpha2,beta2)
	#omega=c(parms[1],parms[4])
	#alpha= c(parms[2],parms[5])
	#beta= c(parms[3],parms[6])

	#OR
	#if wish to have omega_i = 1-alpha_i-beta_i
	#parms vector contains: c(alpha1,beta1,alpha2,beta2)
	alpha=c(parms[1],parms[3])
	beta=c(parms[2],parms[4])
	omega=c(1-alpha[1]-beta[1],1-alpha[2]-beta[2])

	N = length(y)
	lambda = numeric(N) 
	eta = numeric(N)
	regind=ifelse(hour>=t1 & hour <t2,1,2)
	epsilon = y/mu
	eta[1] = 1
	lambda[1] = 1	 
	loglik = 0 # -sum(lfactorial(y))
	for(i in 2:N){
		r=regind[i]
		eta[i] = omega[r] + alpha[r]*epsilon[(i-1)] + beta[r]*ifelse(misshour[(i-1)] == 1, 1, eta[(i-1)])
		lambda[i] = mu[i]*eta[i]
#		if(lambda[i] <= 0){print(c(i,lambda[i],alpha,beta))}
		temp = -lambda[i] + y[i]*log(lambda[i]) - (lfactorial(y[i]))
		loglik = loglik + ifelse(misshour[i] == 1, 0, temp)
	}
	if(llik==TRUE){-loglik}
	else{eta}
}
###############################################################

theta.0 = c(0.074, 0.7901, 0.0366, 0.906) 

IntRsm11(parms = theta.0, y = y, mu = mu.hat, t1=10, t2=16, misshour = misshour, llik = TRUE, hour=hour)

outRsm11 = nlminb(start = theta.0, objective = IntRsm11, y = y, mu=mu.hat, t1=10, t2=16, llik=TRUE, misshour= misshour, hour=hour, lower = c(0.0,0.0,0.0,0.0), upper = c(0.1,0.9,0.1,0.94))

outRsm11 = optim(par=theta.0, fn = IntRsm11, y = y, mu=mu.hat, t1=10, t2=16, llik=TRUE, misshour= misshour, hour=hour, method = "L-BFGS-B", lower = c(0.0,0.0,0.0,0.0), 
	upper = c(0.1,0.9,0.1,0.94), hessian=T, control = 
	list(trace = 5, ndeps = rep.int(0.000001, 4), 
		  maxit = 10L, factr = 1e+31, pgtol = 0))

# parameter estimates
igparRsm11 = outRsm11$par ; igparRsm11

# approximate SEs
igseRsm11 = sqrt(diag(solve(outRsm11$hessian))) ; igseRsm11

# CIIR
etaRsm11 = IntRsm11(parms = outRsm11$par,y = y, mu = mu.hat, t1=10, t2=16, misshour = misshour, llik = FALSE, hour=hour)

# Fitted values
lambdaRsm11 = mu.hat*etaRsm11

# Standardized multiplicative residuals
epsilonRsm11 = y/lambdaRsm11

###############################################################
###############################################################
#### 1-step forecast evaluation
#### In-sample example below, out-of-sample was used
#### in the final analysis
###############################################################
###############################################################

# RMSME (multiplicitive)
rM0 = y/mu.hat-1
rM1 = y/lambdaInt11-1
rM2 = y/lambdaExpAR11-1
rM3 = y/lambdaThresh11-1
rM4 = y/lambdaRsm11-1

sum(rM0^2) ; sqrt(mean(rM0^2))
sum(rM1^2) ; sqrt(mean(rM1^2))
sum(rM2^2) ; sqrt(mean(rM2^2))
sum(rM3^2) ; sqrt(mean(rM3^2))
sum(rM4^2) ; sqrt(mean(rM4^2))

# RMSPE (Pearson)
rP0 = (y-mu.hat)/sqrt(mu.hat)
rP1 = (y-lambdaInt11)/sqrt(lambdaInt11)
rP2 = (y-lambdaExpAR11)/sqrt(lambdaExpAR11)
rP3 = (y-lambdaThresh11)/sqrt(lambdaThresh11)
rP4 = (y-lambdaRsm11)/sqrt(lambdaRsm11)

sum(rP0^2) ; sqrt(mean(rP0^2))
sum(rP1^2) ; sqrt(mean(rP1^2))
sum(rP2^2) ; sqrt(mean(rP2^2))
sum(rP3^2) ; sqrt(mean(rP3^2))
sum(rP4^2) ; sqrt(mean(rP4^2))

# RMSAE (Anscombe)
rA0 = (3/2)*(y^(2/3) - mu.hat^(2/3))/mu.hat^(1/6)
rA1 = (3/2)*(y^(2/3) - lambdaInt11^(2/3))/lambdaInt11^(1/6)
rA2 = (3/2)*(y^(2/3) - lambdaExpAR11^(2/3))/lambdaExpAR11^(1/6)
rA3 = (3/2)*(y^(2/3) - lambdaThresh11^(2/3))/lambdaThresh11^(1/6)
rA4 = (3/2)*(y^(2/3) - lambdaRsm11^(2/3))/lambdaRsm11^(1/6)

sum(rA0^2) ; sqrt(mean(rA0^2))
sum(rA1^2) ; sqrt(mean(rA1^2))
sum(rA2^2) ; sqrt(mean(rA2^2))
sum(rA3^2) ; sqrt(mean(rA3^2))
sum(rA4^2) ; sqrt(mean(rA4^2))

