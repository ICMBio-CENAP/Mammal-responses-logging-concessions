## Multimodel bayesian species richness analysis
# Code used in Carvalho Jr et al. paper "Mammal responses to reduced-impact logging in Amazonian forest concessions"

# The Bayesian model is based on a template written by Elise Zipkin, available at: 
# https://github.com/zipkinlab/Community_model_examples-covariate_model/blob/master/covariate%20model%20code.r
# with modifications by Elildo Carvalho Jr @ICMBio/CENAP

##----- 1 - Load libraries -----
library(R2jags)
library(here)


##----- 2 - Write the model -----
# Specify the model in JAGS language
sink(here("model.txt"))
cat("
   model{

## Define prior distributions for community-level model parameters

omega ~ dunif(0,1)                                  # inclusion probability


for (i in 1:(n+nzeroes)) {
  for (b in 1:nblock){                              
  alpha.psi.block[i,b] ~ dnorm(mu.block, tau.block) # random block effects on psi
    } #b
  } #i
mu.block ~ dnorm(0, 0.5)                            # Hyperparameter for random block effects on psi
tau.block <- 1/(sd.block*sd.block)
sd.block ~ dunif(0,5)


for (i in 1:(n+nzeroes)) {
  for (y in 1:nyear){                              
  alpha.psi.year[i,y] ~ dnorm(mu.year, tau.year)    # random year effects on psi
    } #y
  } #i
mu.year ~ dnorm(0, 0.5)                             # Hyperparameter for random year effects on psi
tau.year <- 1/(sd.year*sd.year)
sd.year ~ dunif(0,5)


# coefficients
mu.alpha.psi ~ dnorm(0, 0.1)                      # grand psi mean
mu.beta1.psi ~ dnorm(0, 0.1)                      # logging intensity
mu.beta2.psi ~ dnorm(0, 0.1)                      # recovery time
mu.beta3.psi ~ dnorm(0, 0.1)                      # road length

mu.alpha.p ~ dnorm(0, 0.1)                        # grand p mean 
mu.beta1.p ~ dnorm(0, 0.1)                        # Intensity
mu.beta2.p ~ dnorm(0, 0.1)                        # Recovery
mu.beta3.p ~ dnorm(0, 0.1)                        # road length
mu.beta4.p ~ dnorm(0, 0.1)                        # Dates


tau.alpha.psi ~ dgamma(0.1, 0.1)
tau.beta1.psi ~ dgamma(0.1, 0.1)
tau.beta2.psi ~ dgamma(0.1, 0.1)
tau.beta3.psi ~ dgamma(0.1, 0.1)

tau.alpha.p ~ dgamma(0.1, 0.1)
tau.beta1.p ~ dgamma(0.1, 0.1)
tau.beta2.p ~ dgamma(0.1, 0.1)
tau.beta3.p ~ dgamma(0.1, 0.1)
tau.beta4.p ~ dgamma(0.1, 0.1)


for (i in 1:(n+nzeroes)) {

# Create priors for species i from the community level prior distributions

    w[i] ~ dbern(omega)                                 # inclusion indicators
    
    alpha.psi[i] ~ dnorm(mu.alpha.psi, tau.alpha.psi)
    beta1.psi[i] ~ dnorm(mu.beta1.psi, tau.beta1.psi)   # logging intensity on psi
    beta2.psi[i] ~ dnorm(mu.beta2.psi, tau.beta2.psi)   # recovery time on psi
    beta3.psi[i] ~ dnorm(mu.beta3.psi, tau.beta3.psi)   # road dens on psi
    
    alpha.p[i] ~ dnorm(mu.alpha.p, tau.alpha.p)
    beta1.p[i] ~ dnorm(mu.beta1.p, tau.beta1.p)        # intensity on p
    beta2.p[i] ~ dnorm(mu.beta2.p, tau.beta2.p)        # recovery on p
    beta3.p[i] ~ dnorm(mu.beta3.p, tau.beta3.p)        # road dens on p
    beta4.p[i] ~ dnorm(mu.beta4.p, tau.beta4.p)        # dates on p


# Create a loop to estimate the Z matrix (true occurrence for species i at point j)
   for (j in 1:J) {
       logit(psi[j,i]) <-  alpha.psi[i] + alpha.psi.block[i,block[j]] + alpha.psi.year[i,year[j]] 
       + beta1.psi[i]*intensity[j] + beta2.psi[i]*recovery[j] + beta3.psi[i]*road.dens[j]

      mu.psi[j,i] <- psi[j,i]*w[i]
      Z[j,i] ~ dbern(mu.psi[j,i])

# Create a loop to estimate detection for species i at point j during sampling period k.      
   for (k in 1:K[j]) {  
      logit(p[j,k,i]) <-  alpha.p[i]
      + beta1.p[i]*intensity[j] + beta2.p[i]*recovery[j] + beta3.p[i]*road.dens[j] + beta4.p[i]*dates[j,k]
      
       mu.p[j,k,i] <- p[j,k,i]*Z[j,i]              # can only be detected if Z=1
       X[j,k,i] ~ dbern(mu.p[j,k,i])
       
    }#k
  }#j
}#i

# Derived quantities:

# Sum all species observed (n) and unobserved species (n0) to find the total estimated richness
n0 <- sum(w[(n+1):(n+nzeroes)])
N <- n + n0
#N <- sum(w[])

# Create a loop to determine point level richness for the whole community and for subsets of interest
for(j in 1:J){
Nsite[j] <- inprod(Z[j,1:(n+nzeroes)],w[1:(n+nzeroes)])
  }
}", fill=TRUE); sink()


##----- 3 - Create the necessary arguments to run the jags.model() command -----

# Load the data
sp.data <- readRDS(here("data_jamari.rds"))
attach(sp.data)
names(sp.data)

# define which data will be used 
jags.data = list(n=n, nzeroes=nzeroes, J=J, K=K, X=X, block=block, nblock=nblock, year=year, nyear=nyear,
                 intensity=intensity, recovery=recovery, road.dens=road.dens, dates=dates)

# Specify the parameters to be monitored
sp.params = c('alpha.psi', 'alpha.psi.block', 'alpha.psi.year',
              'beta1.psi', 'beta2.psi', 'beta3.psi',
              'alpha.p',
              'beta1.p', 'beta2.p', 'beta3.p', 'beta4.p',
              'omega', 'N', 'Nsite')


# Specify the initial values
sp.inits = function() {
  omegaGuess = runif(1, 0, (n/(n+nzeroes)))
  psi.meanGuess = runif(1, 0.1, .5)
  list(omega = omegaGuess,
       w = rep(1,n+nzeroes),
       Z = matrix(1, nrow=J, ncol=(n+nzeroes)),
       alpha.psi = rnorm(n+nzeroes), beta1.psi=rnorm(n+nzeroes), beta2.psi=rnorm(n+nzeroes), beta3.psi=rnorm(n+nzeroes),
       alpha.psi.block = matrix(rnorm((n+nzeroes)*nblock, 0, 0.5), nrow = n+nzeroes, ncol = nblock),
       alpha.psi.year = matrix(rnorm((n+nzeroes)*nyear, 0, 0.5), nrow = n+nzeroes, ncol = nyear),
       alpha.p = rnorm(n+nzeroes),
       beta1.p = rnorm(n+nzeroes), beta2.p = rnorm(n+nzeroes), beta3.p = rnorm(n+nzeroes), beta4.p = rnorm(n+nzeroes) )
  }

##----- 4 - Run the model and save the results -----

# Run the model in jags
fit <- jags(data=jags.data, inits=sp.inits, parameters.to.save=sp.params, n.chains=1, n.iter=100, n.burnin=10, n.thin=10, model.file="model.txt")
fit <- jags(data=jags.data, inits=sp.inits, parameters.to.save=sp.params, n.chains=3, n.iter=100000, n.burnin=50000, n.thin=100, model.file="model.txt")


# Save results  
saveRDS(fit, here("results", "model_results.rds"))

