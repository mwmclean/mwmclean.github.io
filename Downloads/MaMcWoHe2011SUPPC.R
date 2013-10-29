#########################################################################
#########################################################################
#### QueueSimulations.R	     
#### Mathew W. McLean						
#### 2010 - 11 - 18							
#### Annals of Applied Statistics 			
#### Forecasting Emergency Medical Service Call Arrival Rates 
#### By D.S. Matteson, M.W. McLean, D.B. Woodard, S.G. Henderson 
#########################################################################
#########################################################################
#### Implemented using R version 2.09.1
#### The period is assumed to be 1 hour and 
#### no missing observations within a day.
#### User must supply the following data:
#### y : observed counts per hour
#### nsim : number of simulations (of size sum(y) calls) to run
#### preds : predicted arrival rate for each period
#### nu : service rate
#### theta : targeted proportion of calls served immediately
#### q : cost of not serving a customer immediately	
#### simdatadir : directory for storing simulated data
#########################################################################
#########################################################################


######################################################################
########        Functions for choosing number of servers       #######
########            for staffing and evaluating cost           #######
######################################################################

# compute p_0, the formula for which can be found in paragraph 3 on p. 23
# r: ratio of service rate to arrival rate
# c: number of servers

"SpecialSumP0" = function(r,c){
	rho = r/c
	if(rho >= 1) stop("r/c = rho >= 1 not allowed")
	n = 0:c
	numer = r^n
	n2 = c(1,n[-1]) 
	denom = cumprod(n2)
	denom[(c+1)] = denom[(c+1)] * (1 - rho)
	1/(sum(numer/denom))
}

# compute g, the formula for which can be found in paragraph 3 on p. 23
# lambda : arrival rate
# s : number of servers
# nu : service rate

"g" = function(lambda, s, nu = 1){
	if(is.na(lambda)) return(NA)
	r = lambda/nu
	rho = r/s
	p0 = SpecialSumP0(r,s) # p.70 Gross&Harris 3rd Ed.
 	1 - (r^s*p0)/(factorial(s)*(1-rho))
}

# compute Pen function given in paragraph 3 on p. 22
# N.t : number of callers served immediately in period t
# Y.t : number of callers in period t
# theta : targeted proportion of call served immediately
# q : cost of not serving a customer immediately		

"Pen" = function(N.t, Y.t, theta = 0.9, q = 10){
	thetaTimesY = theta*Y.t
	ifelse(N.t >= thetaTimesY, 0, q*(thetaTimesY - N.t))
}

# compute C whose formula is given in the last equation on p. 22
# N.t : number of callers served immediately in each period
# Y.t : number of callers in each period (same length as N.t)
# s.t : number of servers in use for each period (same length as N.t)
# theta : targeted proportion of call served immediately
# q : cost of not serving a customer immediately		

"TotalCostOverPeriods" = function(N.t, Y.t, s.t, theta = 0.9, q = 10){
	sum(Pen(N.t,Y.t, theta, q) + s.t)
}

# compute cost for the period (first formula in para. 3, p. 22)
# N.sim : simulated (binomial) number of callers served immediately
#	 in period t (length equal to J)
# Y.t : simulated number of callers in period t (same length as N.sim)
# s.t : number of servers in use in period t 
# theta : targeted proportion of call served immediately
# q : cost of not serving a customer immediately

"periodCost" = function(s.t, Y.t, N.sim, theta = 0.8, q = 10){
       Pen(N.sim, Y.t, theta, q) + s.t
}

# approximate the expected period cost using the method discussed in
#	paragraphs 1 and 2 on p. 23
# s.t : number of servers in use in period t
# J : number of realizations to generate of number of calls served 
#	immediately for the period (see para. 1 and 2, p.23
# lambda : arrival rate
# nu : service rate
# theta : targeted proportion of call served immediately
# q : cost of not serving a customer immediately		 

"ExpectedPeriodCost" = function(s.t, J = 25000, lambda, nu = 1, 
	theta= 0.8, q = 10){

       if(is.na(lambda)|lambda<=0) stop("issues with lambda!")
       g.temp = g(lambda, s.t, nu)
       Y.temp = rpois(J, lambda)
       N.temp = rbinom(J, Y.temp, g.temp)

       C.temp = periodCost(s.t, Y.t = Y.temp, N.sim = N.temp, theta, q)
       mean(C.temp)
}

# compute \hat{s}_t, found in paragraph 1, Eq. 11, p.23
# lambda : arrival rate
# J : number of realizations to generate of number of calls served 
#	immediately for the period (see para. 1 and 2, p.23
# theta : targeted proportion of call served immediately
# nu : service rate
# q : cost of not serving a customer immediately

"min.s.hat" = function(lambda, J = 25000, theta = 0.8, nu = 1, q = 10){

       if(is.na(lambda)) return(NA)
       r = floor(lambda/nu)
       s.min = max(r + 1, 1)
       C.temp = rep(Inf,100)
       i = s.min
       diff.temp = -1
       while(diff.temp < 0){

               C.temp[i] = ExpectedPeriodCost(s.t = i, J,  lambda, nu, theta, q)

               diff.temp = C.temp[i] - C.temp[(i-1)]
               i = i+1
       }
#       print(C.temp)
       (i-2)
}

######################################################################
########  Functions for simulating queueing model from paper  ########
######################################################################


# simulate arrival time and service times for every caller that will arrive
# in the system based on the inputted number of calls for each period
# y : vector of observed number of calls for each period
# nu : constant service rate
# the output is a matrix whose first column is random uniform(0,1) arrival
# 	times for each call received.  The column is sorted from first to last
#     arrival within each hour.  For example, if the input y=c(3,2), then
#     the first outputted column to 3 decimal places could be
#     c(.112,.454,.822,.398,.745)
#	the second column of the output matrix is the random Exponential(serv_time)
#	service time attached to each caller

sim_data=function(y,nu)
{

	total=sum(y)
	unis=runif(total,0,1)
	counter=1
	for(i in 1:length(y)){
		unis[counter:(y[i]+counter-1)]=sort(unis[counter:(y[i]+
			counter-1)])
		counter=counter+y[i]
	}
	cbind(unis,rexp(total,nu))
}

# simulate the queueing system outlined in section 4.4 of the paper
# data : output from running sim_data(y,nu) (service+arrival times for each caller)
# s : vector of number of servers to use for each period
# y : vector of observed number of callers for each period
# misshour : vector of 0's and 1's that is equal to 1 if the previous
#       period is missing
# iservtimes : vector of service times for each customer starting in the
#	system. length of the vector is the initial # of customers in the system
# the output is a vector of 0's and 1's with ith component equal to 1 if
#       the ith caller was served immediately

queuesim=function(data,s,y,misshour=numeric(length(y)),iservtimes)
{


	MAX_SERVERS=100 # most servers that could conceivably be in queue 
			   # needed to define dimension of the server list
	MAX_QUEUE_LENGTH=1000 # upper bound for number of callers in queue
	T=nrow(data) # total number of calls for the year
	hourgone = 0 # proportion of current hour that has past
	h=1      	 # current hour in the year
	c=1      	 # current call
	servedimmind = numeric(T) #indicator if caller was served immediately
	servers=rep(NA,MAX_SERVERS) # -Inf indicates server doesn't exist
						    # -1 indicates server is idle
						    #positive number is remaining service
						    # time
	servers[1:s[1]]=-1	#initialize num. of active servers for first hour
	queue=rep(NA,MAX_QUEUE_LENGTH)

	#initialize number of customers to start in queue and their service times

	if(!missing(iservtimes)){
		if(length(iservtimes)!=0){
			initinsystem=length(iservtimes)
			for(i in 1:initinsystem){
				if(!allbusy(servers)){
					servers=sadd(servers,iservtimes[i])
				}
				else{
					queue=qadd(queue,iservtimes[i])
				}
			}
		}
	}
	nextcaller=data[1,]    # next customer arriving for service
	totalcallsforhour=y[1] # total number of callers to be seen for the
				     #    current hour
	callssofarforhour=0    # current number of callers seen for the hour
	
	# start simulation

	while(c < T){

		#arrival before end of hour and all servers idle 

		if(sempty(servers) & nextcaller[1] <= 1-hourgone){
			servers=sadd(servers,nextcaller[2])
			hourgone=hourgone+nextcaller[1]
			servedimmind[c]=1
			callssofarforhour=callssofarforhour+1
			c=c+1
			nextcaller=data[c,]-c(hourgone,0)
		}
		# arrival before service finishing and before end of hour

		else if( nextcaller[1] <= servers[1] & 
			nextcaller[1] <= 1-hourgone){
			hourgone=hourgone+nextcaller[1]
			servers=update(servers,nextcaller[1])

			#there is an idle server
			
			if(!allbusy(servers)){
				servedimmind[c]=1
				servers=sadd(servers,nextcaller[2])
			}

			#no idle server, add to queue
	
			else{ 
				queue=qadd(queue,nextcaller[2])
			}
			callssofarforhour=callssofarforhour+1
			c=c+1
			nextcaller=data[c,]-c(hourgone,0)
		}

		#service finishes before next arrival and before hour finishes
		
		else if( nextcaller[1] > servers[1] &
			servers[1] <= 1-hourgone){
			hourgone=hourgone+servers[1]
			nextcaller[1]=nextcaller[1]-servers[1]
			temp=servers[1]
			servers=sremove(servers)
			servers=update(servers,temp)
			if(!qempty(queue)){
				templist=qremove(queue)
				queue=templist[[1]]
				servers=sadd(servers,templist[[2]])
			}
		}

		#seen all calls for that hour

		if(callssofarforhour==totalcallsforhour){
			callssofarforhour=0
			h=h+1
			totalcallsforhour=y[h]
			nextcaller=data[c,]

			#if still have services finishing before end of hour

			while(!sempty(servers) & servers[1] <= 1- hourgone){
				hourgone=hourgone+servers[1]
				temp=servers[1]
				servers=sremove(servers)
				servers=update(servers,temp)
				if(!qempty(queue)){
					templist=qremove(queue)
					queue=templist[[1]]
					servers=sadd(servers,templist[[2]])
				}
			}

			#day of missing data next so restart queue

			if(misshour[h]==1){
				queue=rep(NA,MAX_QUEUE_LENGTH)
				servers=rep(NA,MAX_SERVERS) 
				servers[1:s[h]]=-1
				hourgone=0
			}

			#shift change

			else{
				servers=update(servers,1-hourgone)
				hourgone=0
				servers=hourlyupdate(servers,s[h],s[h-1])
				while(!allbusy(servers) & !qempty(queue)){
					templist=qremove(queue)
					queue=templist[[1]]
					servers=sadd(servers,templist[[2]])
				}
				
			}

			#for case where no callers for the hour still need to
			#	have servers serve

			while(totalcallsforhour==0){
				hourgone=0		 
				while(!sempty(servers) & servers[1] <= 1- hourgone){
					hourgone=hourgone+servers[1]
					temp=servers[1]
					servers=sremove(servers)
					servers=update(servers,temp)
					if(!qempty(queue)){
						templist=qremove(queue)
						queue=templist[[1]]
						servers=sadd(servers,templist[[2]])
					}
				}
				h=h+1
				totalcallsforhour=y[h]
				hourgone=0
				
				#perform same code as just above
				#need to because have gone through another period
				#day of missing data next so "start over"

				if(misshour[h]==1){
					queue=rep(NA,MAX_QUEUE_LENGTH)
					servers=rep(NA,MAX_SERVERS) 
					servers[1:s[h]]=-1
					hourgone=0
				}

				#shift change

				else{
					servers=update(servers,1-hourgone)
					hourgone=0
					servers=hourlyupdate(servers,s[h],s[h-1])
					while(!allbusy(servers) & !qempty(queue)){
						templist=qremove(queue)
						queue=templist[[1]]
						servers=sadd(servers,templist[[2]])
					}
				
				}
			}
#		print(paste('period= ',h,', %done= ',100*h/length(y)))
#		flush.console()
		}
	}

	# determine if very last caller served immediately

	# there is an empty server he is served immediately whenever he arrives

	if(!allbusy(servers)){
		servedimmind[T]=1
	}	

	# else see if enough servers finish and the queue empties by the time
	# he arrives
	
	else if( nextcaller[1] > servers[1] ){
		while(nextcaller[1] > servers[1] & !sempty(servers)){
			nextcaller[1]=nextcaller[1]-servers[1]
			temp=servers[1]
			servers=sremove(servers)
			servers=update(servers,temp)
			if(!qempty(queue)){
				templist=qremove(queue)
				queue=templist[[1]]
				servers=sadd(servers,templist[[2]])
			}
		}
		if(!allbusy(servers)){
			servedimmind[T]=1
		}
	}
	return(servedimmind)
}

# check if the queue is empty
# datavec : vector of service times of everyone in queue

qempty=function(datavec){
	sum(is.na(datavec))==length(datavec)
}

#check if all the servers are idle
# datavec : vector of service times of everyone in service

sempty=function(datavec){
	sum(is.na(datavec) | datavec==-1)==length(datavec)
}

# add caller to empty server
# datavec : vector of service times of everyone in service
# newservice : servicet time of caller being added to service

sadd=function(datavec,newservice){
	index=which(datavec==-1)[1]
	if(is.na(index)){
		print('error')
		cat(datavec)	
		stop('trying to add to full server')
	}
	datavec[index]=newservice
	index=!is.na(datavec) & datavec!=-1
	datavec[index]=sort(datavec[index]) #servers always stay in order of
	return(datavec)				#least service time remaining
}

# add caller to queue
# datavec : vector of service times of everyone in queue
# newservice : servicet time of caller being added to queue

qadd=function(datavec,newservice){
	index=which(is.na(datavec))[1]

	if(is.na(index))stop('trying to add to full queue')
	datavec[index]=newservice
	return(datavec)
}

# update remaining service time to account for elapsed time for customers
#	 still in service
# datavec : vector of service times of everyone in service
# elapsedtime : amount of time that has elapsed since last update, to be
#		    subtracted from service time of all remaining in service

update=function(datavec,elapsedtime){
	indices=!is.na(datavec) & datavec!=-1
	datavec[indices]=datavec[indices]-elapsedtime
	return(datavec)
}

#check if all servers busy
# datavec : vector of service times of everyone in service

allbusy=function(datavec){
	sum(datavec==-1,na.rm=TRUE)==0
}

# remove caller that has finished service
# datavec : vector of service times of everyone in service
# pos specified means to remove server at position pos in vector
# pos not specified means to remove server that has least service time left
# pos specified is for removing server from queue if decreasing servers
#    at start of hour

sremove=function(datavec,pos){
	if(missing(pos)){
		if(allbusy(datavec)){#index is location of last busy server
			index=which(is.na(datavec))[1]-1
		}
		else{
			index=which(datavec==-1)[1]-1
		}
		if(is.na(index))stop('trying to remove from empty server')
		if(index==1){#only one server was in use
			datavec[1]=-1
		}
		else{
			datavec[1:(index-1)]=datavec[2:index]
			datavec[index]=-1
		}
	}
	else{
		datavec[pos:(length(datavec)-1)]=datavec[(pos+1):length(datavec)]
	}
	return(datavec)
}

# remove caller from queue to join empty server
# datavec : vector of service times of everyone in queue

qremove=function(datavec){
	temp=datavec[1]
	index=which(is.na(datavec))[1]-1
	if(is.na(index))stop('trying to remove from empty queue')
	if(index==1){#queue only had one caller
		
		datavec[1]=NA
	}
	else{
		datavec[1:(index-1)]=datavec[2:index]
		datavec[index]=NA
	}
	return(list(datavec,temp))
}

# add or remove servers for shift change at start of each period
# datavec : vector of service times of everyone in service
# numservers : number of servers that will be in use for the period
# prevnumservers : number of servers that were in use for the previce period

hourlyupdate=function(datavec,numservers,prevnumservers){
	diff=numservers-prevnumservers

	# need to increase number of servers in use
	
	if(diff > 0){
		index=which(is.na(datavec))[1]
		datavec[index:(index+diff-1)]=-1
	}

	#need to decrease number of servers in use. first drop idle servers
	#if necessary, drop busy servers in order of least
	#service time remaining

	else if(diff < 0){
		diff=abs(diff)
		if(allbusy(datavec)){
			while(diff > 0){
				datavec=sremove(datavec,1)
				diff=diff-1
			}
		}
		else{
			temp=which(datavec==-1)
			ltemp=length(temp)
			if(diff <= ltemp){
				datavec[temp[(ltemp-diff+1):ltemp]]=NA
			}
			else{
				datavec[temp]=NA
				diff=diff-ltemp
				while(diff > 0){
					datavec=sremove(datavec,1)
					diff=diff-1
				}
			}
		}
		#alternative method of randomly choosing servers for removal
		#for(j in 1:abs(diff)){
		#	server_to_remove=sample(prevnumservers,1)
		#	datavec=sremove(datavec,server_to_remove)
		#	prevnumservers=prevnumservers-1
		#}
	}
	return(datavec)
}

############################################################################
# once the above functions are loaded, the below code will run the queue
# and output the performance measures used in the paper
# the code is written assuming no missing days
# to account for missing days simply repeatedly run this code
# for each block that does not contain a missing day,
# reinitializing the queue each time, and then aggregate the results
############################################################################


# determine the random poisson number of callers starting in the queue
# for each simulation

initcallers = rpois(1,y[1])
	
# now get their service times
initservtimes = vector('list',nsim)
for(i in 1:nsim){
	initservtimes[[i]]=rexp(initcallers,nu)
}

# write the initial service data to a file for later use

setwd(simdatadir)
for(i in 1:nsim){
	filename=paste('initservtimes_sim',i,'.txt',sep='')
	write.table(initservtimes[[i]],file=filename,
		row.names=FALSE, col.names=FALSE,sep=' ')
}

# create arrival and service times for each call in y for each simulation
#	and write them to a file

for( i in 1:nsim){
	tdata=sim_data(y,nu)
	filename=paste('queuesimdata',i,'.txt',sep='')
	write.table(tdata,file=filename,row.names=FALSE,col.names=FALSE,
		sep=' ')
}

# if necessary, read in initservtimes	
	
#initservtimes=vector('list',length(temp))
#for(i in 1:nsim){
#		filename=paste('initservtimes_sim',i,'.txt',sep='')
#		initservtimes[[i]]=scan(filename)
#}

# determine number of servers to staff using approach outlined on p. 23

numserv=rep(0,length(y))
for(i in 1:length(y)){
	numserv[i]=min.s.hat(lambda = preds[i], J = 25000, theta = theta, 
		nu = nu, q = q)
}

siimat=matrix(0,sum(y),nsim) # matrix for storing whether a caller was
					    #    served immediately for every period
					    #    of every simulation
for(j in 1:nsim){
	filename=paste('queuesimdata',j,'.txt',sep='')
	data=as.matrix(read.csv(filename,sep=' ',header=FALSE))
	siimat[,j]=queuesim(data,numserv,y,numeric(length(y)),
		initservtimes[[j]])
}

# compute total number of callers served immediately for each simulation

totalservimm=colSums(siimat)

# mean and standard deviatioin % served immediately for this simulation

100*mean(totalservimm/nrow(siimat))
100*sd(totalservimm/nrow(siimat))

# compute sum total cost for each period for each simulation
# (i.e. the last equation on p. 22 for each simulation)

tcop=numeric(nsim)
for(l in 1:nsim){
	
	# first compute the number of callers served immediately in each period
	#	(i.e. observed value of N_t for each simulation)

	totalservimm.t=numeric(length(y))
	ind=1
	for(k in 1:length(y)){
		totalservimm.t[k]=sum(siimat[ind:(ind+y[k]-1),l])
		ind=ind+y[k]
	}

	tcop[l]=TotalCostOverPeriods(totalservimm.t, y, s.t=numserv, 
		theta = theta, q = q)
	
	#print(paste('%done=',100*l/nsim))
	#flush.console()
}

# mean and standard deviation of total per period cost over nsim simulations

mean(tcop)
sd(tcop)

#########################################################################
#     a simple simulation of hourly call totals for 12 weeks     ########
#                  to test the above code                        ########
#########################################################################

T = 24*7*12
hour = rep(1:24, 7*12)
day = rep(1:(7*12), each = 24)
dofw = rep(1:7, each = 24, 12)
week = rep(1:12, each = 24*7)

err = arima.sim(n = T, list(ar = 0.95, ma = -0.8), rand.gen =
function(n, ...) rnorm(n, 0, 1))
rate = 20 + 10*sin(2*pi*hour/24) + 2*sin(2*pi*week/26) + err
y = rpois(T,rate)

ts.plot(y)
lines(rate, col = 2)

preds=rnorm(T,mean(y),sd(y))
preds=ifelse(preds <10,10,preds) #ensure that the predicted lambda is not <= 0
lines(preds,col='green')

nsim=10
q=10
theta=.8 #.9
nu=1     #2/3
simdatadir=getwd()

