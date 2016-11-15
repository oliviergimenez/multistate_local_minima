# Local minima in multistate capture-recapture models
# R analyses using Mark and RMark
# Cooch & White, Gentle introduction to Mark, chapter 9

############################
########### DATA ###########
############################

# fake data from J. Dupuis used in Gimenez et al. (2005) 
# true parameters are:
# survival = 1 ; detection = 0.6 ; psi12 = 0.6 ; psi21 = 0.85

his <- matrix(c(2, 0, 2, 1, 2, 0, 2, 2, 0, 2, 1, 2, 0, 2,
2, 0, 2, 1, 2, 0, 2,
2, 0, 2, 1, 2, 0, 2,
1, 1, 1, 0, 1, 0, 1,
1, 1, 1, 0, 1, 0, 1,
1, 1, 1, 0, 1, 0, 1,
1, 1, 1, 0, 1, 0, 1,
2, 0, 2, 0, 2, 0, 1,
2, 0, 2, 0, 2, 0, 1,
2, 0, 2, 0, 2, 0, 1,
2, 0, 2, 0, 2, 0, 1,
1, 0, 1, 0, 1, 0, 1,1, 0, 1, 0, 1, 0, 1,
1, 0, 1, 0, 1, 0, 1,
1, 0, 1, 0, 1, 0, 1,
2, 0, 2, 0, 2, 0, 2,2, 0, 2, 0, 2, 0, 2,
2, 0, 2, 0, 2, 0, 2,
2, 0, 2, 0, 2, 0, 2,
1, 0, 1, 0, 1, 0, 2,
1, 0, 1, 0, 1, 0, 2,
1, 0, 1, 0, 1, 0, 2,
1, 0, 1, 0, 1, 0, 2,
2, 2, 0, 1, 0, 2, 1,
2, 2, 0, 1, 0, 2, 1,
2, 2, 0, 1, 0, 2, 1,
2, 2, 0, 1, 0, 2, 1,
2, 1, 0, 2, 0, 1, 1,
2, 1, 0, 2, 0, 1, 1,
2, 1, 0, 2, 0, 1, 1,
2, 1, 0, 2, 0, 1, 1),byrow=T,ncol=7)

# format data for analysis in RMark
k = ncol(his) # nb of capture occasions	
n = nrow(his) # nb of individuals	
out = array(dim=n)	
for (i in 1:n){	
	y = (his[i,] > 0) * his[i,]	
	out[i] = paste(y,collapse="")	
}	
capt.hist = data.frame(ch = out)	

#############################
########### MODEL ###########
#############################

# load RMark package	
library(RMark)	
	
# Process data	
mstrata.processed=process.data(capt.hist,model="Multistrata")	
	
# Create default design data	
mstrata.ddl=make.design.data(mstrata.processed)	
	
#---- fit with sin link function
	
# Define survival probability	
S.stratum=list(formula=~1) # survival depends on states	
	
#  Define detection probability	
p.dot=list(formula=~1) # constant over time, does not depend on states	
	
# Define transition probs	
Psi.s=list(formula=~-1+stratum:tostratum,link='sin')	
	
# Run model with state effect on survival	
mstrata.mod = mark(mstrata.processed,mstrata.ddl,model.parameters=list(S=S.stratum,p=p.dot,Psi=Psi.s),output = FALSE,delete=T)	
mstrata.mod$results$real

#                            estimate        se       lcl       ucl fixed    note
# S s1 g1 c1 a0 o1 t1       1.0000000 0.0000000 1.0000000 1.0000000              
# p s1 g1 c1 a1 o1 t2       0.5833333 0.0355797 0.5123869 0.6509884              
# Psi s1 to2 g1 c1 a0 o1 t1 0.2480010 0.0666486 0.1406682 0.3991872              
# Psi s2 to1 g1 c1 a0 o1 t1 0.3419952 0.0753736 0.2123359 0.5005178

# idem with simulated annealing
mstrata.mod = mark(mstrata.processed,mstrata.ddl,model.parameters=list(S=S.stratum,p=p.dot,Psi=Psi.s),output = FALSE,delete=T,options="SIMANNEAL")	
mstrata.mod$results$real

#                            estimate           se       lcl       ucl fixed    note
# S s1 g1 c1 a0 o1 t1       1.0000000 3.299808e-06 0.9999935 1.0000065              
# p s1 g1 c1 a1 o1 t2       0.5833309 3.557970e-02 0.5123844 0.6509861              
# Psi s1 to2 g1 c1 a0 o1 t1 0.5993189 7.381930e-02 0.4501934 0.7320690              
# Psi s2 to1 g1 c1 a0 o1 t1 0.8442048 7.056110e-02 0.6543753 0.9394246

#---- fit with logit link function
	
# Define transition probs	
Psi.s=list(formula=~-1+stratum:tostratum,link='mlogit')	
	
# Run model with state effect on survival	
mstrata.mod = mark(mstrata.processed,mstrata.ddl,model.parameters=list(S=S.stratum,p=p.dot,Psi=Psi.s),output = FALSE,delete=T)	
mstrata.mod$results$real[c(1:3,24),1:4]	


#                            estimate        se       lcl       ucl
# S s1 g1 c1 a0 o1 t1       1.0000000 0.0000000 1.0000000 1.0000000
# p s1 g1 c1 a1 o1 t2       0.5833333 0.0355797 0.5123869 0.6509885
# Psi s1 to2 g1 c1 a0 o1 t1 0.5993209 0.0738200 0.4501940 0.7320719
# Psi s2 to1 g1 c1 a0 o1 t1 0.8441961 0.0705638 0.6543621 0.9394203


###############################
########### PROFILE ###########
###############################

grid = seq(.2,0.99,0.01)
dev = rep(NA,length(grid))
ind = 1
for (i in grid){
	Psi.s=list(formula=~-1+stratum:tostratum,fixed=list(index=22:42,value=rep(i,21)))
	res = mark(mstrata.processed,mstrata.ddl,model.parameters=list(S=S.stratum,p=p.dot,Psi=Psi.s),output = FALSE,delete=T)	
	dev[ind] = res$results$deviance
	ind = ind + 1
}
plot(grid,dev,type='l',xlab=expression(psi^{21}),ylab='deviance',col='black',lwd=3)

# add deviance value when initial value for psi21 is changed (other parameters to estimated value)
initial = mstrata.mod$results$beta$estimate
initial_back_transformed = 1/(1+exp(-initial))
#initial_back_transformed = rep(.5,4)
dev2 = rep(NA,length(grid))
ind = 1
Psi.s=list(formula=~-1+stratum:tostratum) # to unfix the psi's
for (i in grid){
	initial_back_transformed[3] = i # fix psi21 to diff values
	init = log(initial_back_transformed/(1-initial_back_transformed))
	res = mark(mstrata.processed,mstrata.ddl,model.parameters=list(S=S.stratum,p=p.dot,Psi=Psi.s),output = FALSE,delete=T,initial= init)	
	dev2[ind] = res$results$deviance
	ind = ind + 1
}
# initial values in red area lead to local minimum
polygon(c(min(grid[1:24]), max(grid[1:24]), max(grid[1:24]), min(grid[1:24]), min(grid[1:24])), c(min(dev2),min(dev2),max(dev), max(dev),min(dev2)),col=rgb(1, 0, 0,0.5), border=NA)
# initial values in green area lead to local minimum
polygon(c(min(grid[25:80]), max(grid[25:80]), max(grid[25:80]), min(grid[25:80]), min(grid[25:80])), c(min(dev2),min(dev2),max(dev), max(dev),min(dev2)),col=rgb(0, 1, 0,0.5), border=NA)
lines(grid,dev,col='black',lwd=3)

# save results
save(dev,dev2,file='loc_min.RData')
#load('loc_min.RData')



###############################
########### MCMC ###########
###############################

# Let us get the occasion of first capture for each individual:
get.first <- function(x) min(which(x!=0))
f <- apply(his, 1, get.first)
f

# Recode the data such that 1 = seen alive in A, 2 = seen alive in B, 3 = not seen:
his_recoded <- his
his_recoded[his_recoded==0] <- 3

# Then we fit a multistate model:
sink("state_on_survival.jags")
cat("
model {

# -------------------------------------------------
# Parameters:
# phi: survival probability
# psiAB: movement probability from site A to site B
# psiBA: movement probability from site B to site A
# p: recapture probability
# -------------------------------------------------
# States (S):
# 1 alive at A
# 2 alive at B
# 3 dead
# Observations (O):  
# 1 seen at A 
# 2 seen at B
# 3 not seen
# -------------------------------------------------

# Priors
   phi ~ dunif(0, 1)
   psiAB ~ dunif(0, 1)
   psiBA ~ dunif(0, 1)
   p ~ dunif(0, 1)

# Define state-transition and observation matrices
for (i in 1:nind){  
   # Define probabilities of state S(t+1) given S(t)
   for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phi * (1-psiAB)
      ps[1,i,t,2] <- phi * psiAB
      ps[1,i,t,3] <- 1-phi
      ps[2,i,t,1] <- phi * psiBA
      ps[2,i,t,2] <- phi * (1-psiBA)
      ps[2,i,t,3] <- 1-phi
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- p
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 1-p
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- p
      po[2,i,t,3] <- 1-p
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- 1
      } #t
   } #i

# Likelihood 
for (i in 1:nind){
   # Define latent state at first capture
   z[i,f[i]] <- y[i,f[i]]
   for (t in (f[i]+1):n.occasions){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
   } #i
}
",fill = TRUE)
sink()

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
   # notseen: label for ënot seení
   state <- ms
   state[state==notseen] <- NA
   for (i in 1:dim(ms)[1]){
      m <- min(which(!is.na(state[i,])))
      state[i,m] <- NA
      }
   return(state)
   }

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
   for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
   states <- max(ch, na.rm = TRUE)
   known.states <- 1:(states-1)
   v <- which(ch==states)
   ch[-v] <- NA
   ch[v] <- sample(known.states, length(v), replace = TRUE)
   return(ch)
   }

# Bundle data
jags.data <- list(y = his_recoded, f = f, n.occasions = dim(his_recoded)[2], nind = dim(his_recoded)[1], z = known.state.ms(his_recoded, 3))

# Initial values
inits <- function(){list(phi = runif(1, 0, 1), psiAB = runif(1, 0, 1), p = runif(1, 0, 1), z = ms.init.z(his_recoded, f))}  

# Parameters monitored
parameters <- c("phi", "psiAB", "psiBA", "p")

# MCMC settings
ni <- 5000
nt <- 6
nb <- 2000
nc <- 2

# Call JAGS from R
library(R2jags)
ms <- jags(jags.data, inits, parameters, "state_on_survival.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
print(ms, digits = 3)

par(mfrow=c(1,2))
hist(ms$BUGSoutput$sims.matrix[,2],main='detection',xlab='')
hist(ms$BUGSoutput$sims.matrix[,5],main='transition 2->1',xlab='')