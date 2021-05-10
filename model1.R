## Multimodel bayesian species richness analysis
# Script for analyzing camera-trap data from Jamari National Forest, Brazil

# Based on template from: 
# https://github.com/zipkinlab/Community_model_examples-covariate_model/blob/master/covariate%20model%20code.r
# and Yamaura et al. 2011. Journal of Applied Ecology, doi: 10.1111/j.1365-2664.2010.01922.x
# with modifications by Elildo @ICMBio/CENAP

##----- 1 - Load libraries -----
library(R2jags)
library(here)


##----- 2 - Write the model -----
# Specify the model in JAGS language
sink(here("src", "model_foreco_LvsU.txt"))
cat("
   model{

   ## Prior distributions for community-level model parameters
   
   omega ~ dunif(0,1)  # inclusion probability

   for (i in 1:(n+nzeroes)) {
     for (b in 1:nblock){
       a.psi.block[i,b] ~ dnorm(mu.block, tau.block)  # random site effects
     } #j
   } #i
   mu.block ~ dnorm(0, 0.5)  # hyperparameter
   tau.block <- 1/(sd.block*sd.block)
   sd.block ~ dunif(0,5)
   
   for (i in 1:(n+nzeroes)) {
     for (y in 1:nyear){
       a.psi.year[i,y] ~ dnorm(mu.year, tau.year)  # random year effects
    } #j
   } #i
   mu.year ~ dnorm(0, 0.5)  # hyperparameter
   tau.year <- 1/(sd.year*sd.year)
   sd.year ~ dunif(0,5)


   ## Coefficients
   
   # mean value
   # parameter related to occupancy
   mu.psi0 ~ dnorm(0, 0.1)	# overall intercept
   mu.psi1 ~ dnorm(0, 0.1)	# treatment

   # parameter related to detection
   mu.p0 ~ dnorm(0, 0.1) # overall intercept
   mu.p1 ~ dnorm(0, 0.1)  # treatment
   mu.p2 ~ dnorm(0, 0.1)  # dates

   # precision
   tau.psi0 ~ dgamma(0.1, 0.1)
   tau.psi1 ~ dgamma(0.1, 0.1)

   tau.p0 ~ dgamma(0.1, 0.1)
   tau.p1 ~ dgamma(0.1, 0.1)
   tau.p2 ~ dgamma(0.1, 0.1)

   for (i in 1:(n+nzeroes)) {
   
   # Create priors for species i from the community level prior distributions
   
   w[i] ~ dbern(omega)  # inclusion indicators
   
   a.psi0[i] ~ dnorm(mu.psi0, tau.psi0)  # intercept
   a.psi1[i] ~ dnorm(mu.psi1, tau.psi1)  # treatment on psi

   a.p0[i] ~ dnorm(mu.p0, tau.p0)
   a.p1[i] ~ dnorm(mu.p1, tau.p1)  # treatment on p
   a.p2[i] ~ dnorm(mu.p1, tau.p1)  # dates on p

   ## Likelihood

   # Estimate the Z matrix (true occurrence for species i at point j)
   for (j in 1:J) {
     logit(psi[j,i]) <-  a.psi0[i] + a.psi.block[i,block[j]] + a.psi.year[i,year[j]]
     + a.psi1[i]*treatment[j]
     
     mu.psi[j,i] <- psi[j,i]*w[i]
     Z[j,i] ~ dbern(mu.psi[j,i])
       
   # Estimate detection for species i at point j during sampling period k
   for (k in 1:K[j]) {
     logit(p[j,k,i]) <-  a.p0[i]
     + a.p1[i]*treatment[j] + a.p2[i]*dates[j,k]
     
     mu.p[j,k,i] <- p[j,k,i]*Z[j,i]  # can only be detected if Z=1
     X[j,k,i] ~ dbern(mu.p[j,k,i])
     
     }#k
    }#j
   }#i
  
  # Derived quantities:
  
  # Sum al observed (n) and unobserved species (n0) to find estimated richness
  n0 <- sum(w[(n+1):(n+nzeroes)])
  N <- n + n0
  
  # Determine point level richness
  for(j in 1:J){
    Nsite[j] <- inprod(Z[j,1:(n+nzeroes)],w[1:(n+nzeroes)])
   }
}", fill=TRUE)
sink()


##----- 3 - Create the necessary arguments to run the jags.model() command -----

# Load all the data
data <- readRDS(here("data", "data_FORECO.rds"))
attach(data)
names(data)


treatment <- recovery.original
treatment[treatment < 4] <- 1
treatment[treatment == 4] <- 0
treatment


# define which data will be used 
jags_data = list(n=n, nzeroes=nzeroes, J=J, K=effort, X=X, block=block, nblock=nblock, year=year, nyear=nyear,
                 treatment=treatment, dates=dates)

# Specify the parameters to be monitored
params = c("N", "Nsite",
           "a.psi0", "a.psi1",  
           "a.p0", "a.p1", "a.p2", 
           "a.psi.block", "a.psi.year",
           "Z", "mu.psi")


# Specify the initial values
inits = function() {
   omegaGuess = runif(1, 0, (n/(n+nzeroes)))
   psi.meanGuess = runif(1, 0.1, .5)
   list(omega = omegaGuess,
        #w = c(rep(1, n), rbinom(nzeroes, size = 1, prob=omegaGuess)),
        w = rep(1,n+nzeroes),
        #Z = cbind(Zobs, matrix(rbinom((nzeroes)*J, size = 1, prob = psi.meanGuess), nrow = J, ncol = (nzeroes))),
        #Z = matrix(rbinom((n+nzeroes)*J, size = 1, prob = psi.meanGuess), nrow = J, ncol = (n+nzeroes)),
        Z = matrix(1, nrow=J, ncol=(n+nzeroes)),
        a.psi0=rnorm(n+nzeroes), a.psi1=rnorm(n+nzeroes),
        a.psi.block = matrix(rnorm((n+nzeroes)*nblock, 0, 0.5), nrow = n+nzeroes, ncol = nblock),
        a.psi.year = matrix(rnorm((n+nzeroes)*nyear, 0, 0.5), nrow = n+nzeroes, ncol = nyear),
        a.p0 = rnorm(n+nzeroes), a.p1=rnorm(n+nzeroes), a.p2=rnorm(n+nzeroes) )
}

##----- 4 - Run the model and save the results -----

# Run the model in jags
#fit <- jags(data=jags_data, inits=inits, parameters.to.save=params, n.chains=3, n.iter=1000, n.burnin=50, n.thin=10, model.file=here("src","model_foreco_LvsU.txt"))
fit <- jags(data=jags_data, inits=inits, parameters.to.save=params, n.chains=3, n.iter=100000, n.burnin=50000, n.thin=100, model.file=here("src","model_foreco_LvsU.txt"))

# Save results  
saveRDS(fit, here("results", "model_FORECO_LvsU.rds"))

