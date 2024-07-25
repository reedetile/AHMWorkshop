# R script for abundance-mediated interaction models between coyotes and bobcats from  the MesoCarnivores dataset #

setwd('D:/gitrepos/AHMWorkshop/Day4')

# load packages
library(unmarked)
library(AHMbook)
library(nimble)
library(coda)
library(ggplot2)
library(runjags)
library(ggfortify)

# load the data
data("MesoCarnivores")

# look at  data
lapply(MesoCarnivores, head)

# get detection/non-detection data
y_b <- MesoCarnivores$bobcat
y_c <- MesoCarnivores$coyote

# Get sample sizes
nsites <- dim(y_b)[1]       # 1437 sites
nsurveys <- dim(y_b)[2]     # 3 surveys

# Prepare site covariates for analysis
siteCovs <- MesoCarnivores$sitecovs

dist_cov <- standardize(siteCovs$Dist_5km)
hdens_cov <- standardize(siteCovs$HDens_5km)
lat_cov <- standardize(siteCovs$Latitude)
long_cov <- standardize(siteCovs$Longitude)
people_cov <- standardize(siteCovs$People_site)

# Prepare observation covariate for detection
table(Trail <- MesoCarnivores$sitecovs[,'Trail'])    # summarize
Trail <- matrix(cbind(Trail, Trail, Trail), ncol = 3)
Trail_cov <- Trail

# Set number of K subsamples in each occasion j
K = 1 

# Single chain run

# specify model
NimModel <- nimbleCode({
  #Model for dominant species (D)
  # Likelihood
  # State model for abundance of dominant species (D)
  for (i in 1:nsites) {
    nD[i] ~ dpois(lambdaD[i])
    log(lambdaD[i]) <-  beta0D + betaD[1] * dist[i] + betaD[2] * longitude[i] + betaD[3] * latitude[i] + betaD[4] * Hdens[i]
  }
  
  #Observation model for detection/non-detection data of dominant species (D)
  for (i in 1:nsites) {
    for (j in 1:nsurveys){
      logit(rD[i,j]) <- alpha0D + alphaD[1] * trail[i,j] 
      pD[i,j] <- 1 - (1 - rD[i,j])^nD[i] #prob detecting at least 1 ind per day
      yD.new[i,j] ~ dbinom(pD[i,j],size=K)
    }
  }
  
  #Model for subordinate species (S)
  # Likelihood
  # State model for abundance of subordinate species (S)
  for (i in 1:nsites) {
    nS[i] ~ dpois(lambdaS[i])
    log(lambdaS[i]) <-  beta0S + betaS[1] * dist[i] + betaS[2] * longitude[i] + betaS[3] * latitude[i] + betaS [4] * Hdens[i] + gamma0DS * nD[i] 
  }
  
  #Observation model for detection/non-detection data of subordinate species (S)
  for (i in 1:nsites) {
    for (j in 1:nsurveys){
      logit(rS[i,j]) <- alpha0S + alphaS[1] * trail[i,j] 
      pS[i,j] <- 1 - (1 - rS[i,j])^nS[i] #prob detecting at least 1 ind per day
      yS.new[i,j] ~ dbinom(pS[i,j],size=K)
    }
  }
  # --- Priors ---
  #Priors for dominant species (D)
  beta0D ~ dnorm(0, 0.1) 
  alpha0D ~ dnorm(0, 0.1) 
  for (k in 1:1) {
    alphaD[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:4) {
    betaD[k] ~ dnorm(0, 0.1)
  }
  #Priors for subordinate species (s)
  beta0S ~ dnorm(0, 0.1) 
  alpha0S ~ dnorm(0, 0.1) 
  for (k in 1:1) {
    alphaS[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:4) {
    betaS[k] ~ dnorm(0, 0.1)
  }
  gamma0DS ~ dnorm(0, 0.1)
  
} )

# example of plotting prior (if you wondering the difference between the values in the normal distribution in the model above and here are that nimble takes
# mean and precision for its dnorm function, while , ggdistribution takes mean and sd. Where precision is the inverse of the variance, and variance = sd^2
ggdistribution(dnorm, seq(-5, 5, 0.1), mean = 0, sd = 3.16) + theme_classic()

# Parameters monitored
parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0S", "alphaS", "beta0S", "betaS", "gamma0DS")

#data and constants
Nimdata <- list(yD.new = y_c, yS.new = y_b)
constants=list(nsites = nsites, nsurveys = nsurveys, K = K, trail= Trail_cov, dist = dist_cov, Hdens = hdens_cov, latitude = lat_cov, longitude = long_cov)

#inits list
Niminits=list(beta0D=0,beta0S=0, nD=apply(y_c,1,max), nS =apply(y_b,1,max), alpha0S=0, alpha0D=0,
              alphaD=0,alphaS=0, betaD=c(0,0,0,0), betaS=c(0,0,0,0), gamma0DS=0)

#thinning rate
nt = 4

# Build the model, configure the mcmc, and compile
start.time<-Sys.time()
Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=nt,useConjugacy = TRUE)

# Build and compile
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#mini run iter
n.iter <- 10000

# Run the model.
start.time2<-Sys.time()
Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
end.time<-Sys.time()
time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
time2=end.time-start.time2 # post-com

#get the samples
mvSamples = as.matrix(Cmcmc$mvSamples)

#plot the chain
n.iter=2500
n.burn = 500
plot(mcmc(mvSamples[n.burn:nrow(mvSamples),]))

# package it up and run 3 chains together
# This is commented out as it will take several hours to run - load the prerun chains below this chunk of code
# n.chains = 3
# chains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
# for(chain in 1:n.chains){
# 
#   # Parameters monitored
#   parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0S", "alphaS", "beta0S", "betaS", "gamma0DS")
# 
#   #data and constants
#   Nimdata <- list(yD.new = y_c, yS.new = y_b)
#   constants=list(nsites = nsites, nsurveys = nsurveys, trail= Trail_cov, dist = dist_cov, Hdens = hdens_cov, latitude = lat_cov, longitude = long_cov, K = K)
# 
#   #inits list
#   Niminits=list(beta0D=0,beta0S=0, nD=apply(y_c,1,max), nS =apply(y_b,1,max), alpha0S=0, alpha0D=0, alphaD=0,alphaS=0, betaD=c(0,0,0,0), betaS=c(0,0,0,0), gamma0DS=0)
# 
#   #thinning rate
#   nt = 2
# 
#   # Build the model, configure the mcmc, and compile
#   start.time<-Sys.time()
#   Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
#                         inits=Niminits)
#   conf <- configureMCMC(Rmodel,monitors=parameters,thin=nt, useConjugacy = TRUE, enableWAIC=TRUE)
# 
#   # Build and compile
#   Rmcmc <- buildMCMC(conf)
#   Cmodel <- compileNimble(Rmodel)
#   Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
# 
#   #full run iter
#   n.iter <- 50000
# 
#   # Run the model.
#   start.time2<-Sys.time()
#   Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
#   end.time<-Sys.time()
#   time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
#   time2=end.time-start.time2 # post-com
# 
#   #get the chains
#   mvSamples = as.matrix(Cmcmc$mvSamples)
#   mvSamples2 = as.matrix(Cmcmc$mvSamples2)
# 
#   chains[[chain]]=mvSamples
#   chains2[[chain]]=mvSamples2
# }

# load the pre run chains
load('Two-species_abundance-abundance_model_Mesocarnivores__coyote_bobcat_50k_chains_default_samplers.RData')

# combine the chains and burn
n.iter=25000
n.burn = 20000
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 
# Have these converged??

# Gelman-rubin diagnostics
gelman_a <- gelman.diag(a)
gelman_a

# plot the chains
plot(a)

# combine the 3 chains into a single chain
a=runjags::combine.mcmc(a)

# get parameter names
colnames(mvSamples)

# check for posterior correlation between the parameters
cor <- cor(a[, c("alpha0D",   "alpha0S",   "alphaD[1]", "alphaS[1]", "beta0D",    "beta0S",    "betaD[1]",  "betaD[2]",  "betaD[3]",  "betaD[4]",  "betaS[1]", "betaS[2]",  "betaS[3]",  "betaS[4]",  "gamma0DS")])

# look at it, hard to spot the correlated parameters
cor

# turn to absolute values and identify which parameters have >0.5 posterior correlation
abs <- abs(cor)>0.5

# look at it, much nicer.
abs

## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 
# CHALLENGE 02- Identify the correlated parameters  #

# package up 3 chains with RW block updates on parameters showing posterior correlation, this is commented out as takes a few hours to run, load chains below
# n.chains = 3
# chains = vector("list", n.chains)
# chains2 = vector("list", n.chains)
# for(chain in 1:n.chains){
# 
#     # Parameters monitored
#   parameters <-  c("alpha0D", "alphaD", "beta0D", "betaD","alpha0S", "alphaS", "beta0S", "betaS", "gamma0DS") 
#   
#   #data and constants
#   Nimdata <- list(yD.new = y_c, yS.new = y_b)
#   
#   # Nimdata <- list(yD.new = yF.new)
#   
#   constants=list(nsites = nsites, nsurveys = nsurveys, trail= Trail_cov, dist = dist_cov, Hdens = hdens_cov, latitude = lat_cov, longitude = long_cov, K = K)
# 
#   #inits list
#   
#   Niminits=list(beta0D=0,beta0S=0, nD=apply(y_c,1,max), nS =apply(y_b,1,max), alpha0S=0, alpha0D=0, alphaD=0,alphaS=0, betaD=c(0,0,0,0), betaS=c(0,0,0,0), gamma0DS=0)
#   
#   #thinning rate
#   nt = 8
#   nt2 = 400
#   # 
#   # nt = 2
#   # Build the model, configure the mcmc, and compile
#   start.time<-Sys.time()
#   Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
#                         inits=Niminits)
#   conf <- configureMCMC(Rmodel,monitors=parameters, monitors2 = parameters, thin=nt, thin2 = nt2,useConjugacy = TRUE, enableWAIC=TRUE)
#   
#   conf$addSampler(target = c("beta0S", "betaS[4]"),
#                   type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
#   conf$addSampler(target = c( "beta0D", "alpha0D", "alphaD[1]"),
#                   type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
#   conf$addSampler(target = c("betaS[2]","betaS[3]"),
#                   type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
#   conf$addSampler(target = c("betaD[2]","betaD[3]"),
#                   type = 'RW_block',control = list(adaptive=TRUE),silent = TRUE)
#   
#   # Build and compile
#   Rmcmc <- buildMCMC(conf)
#   Cmodel <- compileNimble(Rmodel)
#   Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#   
#   #full run iter
#   n.iter <- 200000
#   
#   # Run the model.
#   start.time2<-Sys.time()
#   Cmcmc$run(n.iter,reset=FALSE) #can keep running this line to get more sample
#   end.time<-Sys.time()
#   time1=end.time-start.time  # total time for compilation, replacing samplers, and fitting
#   time2=end.time-start.time2 # post-com
#   
#   #get the chains
#   mvSamples = as.matrix(Cmcmc$mvSamples)
#   mvSamples2 = as.matrix(Cmcmc$mvSamples2)
#   
#   chains[[chain]]=mvSamples
#   chains2[[chain]]=mvSamples2
#   
# }

# load pre-run chains!
load('Two-species_abundance-abundance_model_Mesocarnivores__coyote_bobcat_200k_chains_custom_RWblocksamplers.RData')
load('Two-species_abundance-abundance_model_Mesocarnivores__coyote_bobcat_200k_chains2_custom_RWblocksamplers.RData')

## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 
# CHALLENGE 03 - Assess if these model runs with custom RW block samplers have converged (5 minutes)

#combine the chains and burn
n.iter=25000
n.burn = 20000
a=mcmc.list(mcmc(chains[[1]][n.burn:n.iter,]),
            mcmc(chains[[2]][n.burn:n.iter,]),
            mcmc(chains[[3]][n.burn:n.iter,]))

# gelman diagnostics
gelman.diag(a)

# traceplots
plot(a)

# how are we looking??? can we proceed?

#They look decent, although the density of betaS[4] looks a little skewed

# Once we are happy that our chains have converged we can proceed. 

# outputs
summary(a)

# name summary
sum <- summary(a)

# create a matrix and take parameter names, means, credible intervals etc.
esti <- matrix(NA, nrow = ncol(mvSamples), ncol = 6)
esti[,1] <- colnames(mvSamples)
esti[,2] <- sum$statistics[,1]
esti[,3] <- sum$quantiles[,1]
esti[,4] <- sum$quantiles[,5]
esti[,5] <- c("coyote", "bobcat", "coyote", "bobcat", "coyote", "bobcat", "coyote", "coyote", "coyote","coyote", "bobcat", "bobcat", "bobcat","bobcat","bobcat")
esti[,6] <- c("observation","observation", "observation", "observation", "state", "state", "state", "state", "state", "state", "state", "state", "state", "state", "state")

# give the columns names
colnames(esti)<- c("parameter", "mean", "LCI", "UCI", "species", "submodel")

# remove the intercepts so we can plot the slopes
coefs_esti<- esti[-c(1,2,5,6),]

# change parameter names to covariate names
coefs_esti[,1] <- c("trail", "trail", "disturbance", "longitude", "latitude", "housing density","disturbance", "longitude", "latitude", "housing density", "coyote abundance")

#colour blind friendly pallette
cbbPalette <- c("#E69F00", "#000000", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# make df
coefs_esti <- as.data.frame(coefs_esti)

# take just the bobcat estimates and reorder the parameter levels
bobcat <- subset(coefs_esti, species == "bobcat")
bobcat$parameter <- factor(bobcat$parameter, levels = c("coyote abundance", "disturbance", "longitude", "latitude", "housing density", "trail"))

# take the coyote estimates and reorder the parameter levels
coyote <- subset(coefs_esti, species == "coyote")
coyote$parameter <- factor(coyote$parameter, levels = c("disturbance", "longitude", "latitude", "housing density", "trail"))

# make everything numeric
bobcat$mean <-  as.numeric(bobcat$mean)
bobcat$LCI <- as.numeric(bobcat$LCI)
bobcat$UCI <- as.numeric(bobcat$UCI)
coyote$mean <-  as.numeric(coyote$mean)
coyote$LCI <- as.numeric(coyote$LCI)
coyote$UCI <- as.numeric(coyote$UCI)

# plot them

ggplot(data = bobcat, aes(x = mean, y = parameter, group = submodel, colour = submodel, fill = submodel))+geom_point(aes(), size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+ 
  ylab("parameter") + xlab("coefficient estimates")  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid

ggplot(data = coyote, aes(x = mean, y = parameter, group = submodel, colour = submodel, fill = submodel))+geom_point(aes(), size = 2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), width = 0.2, size = 0.75, position=position_dodge(width=0.5))+ theme_classic()+
  scale_colour_manual(values = cbbPalette) + scale_fill_manual(values = cbbPalette)+ guides(colour = "none", fill= "none")+
  ylab("parameter") + xlab("coefficient estimates")  + geom_vline(xintercept=0, linetype="dashed")# Make plot black and white with no background grid


######### Producing predictions occupancy estimates ##########
# Now create mcmc object for heavily thinned chains (from which we will do our posterior calcs and visualizations)
n.iter=500
n.burn = 250
b=mcmc.list(mcmc(chains2[[1]][n.burn:n.iter,]),
            mcmc(chains2[[2]][n.burn:n.iter,]),
            mcmc(chains2[[3]][n.burn:n.iter,]))

b=runjags::combine.mcmc(b)

# This example produces expected abundance ests for coyote and bobcat as function of housing density
# number of posterior samples
n.samples = nrow(b)

# range of housing density covariate
r<- range(hdens_cov)

# create sequence along range of housing density
hdens_data <- seq(r[1], r[2], length.out=500)

# set all other covs at mean
dist_mean = 0
latitude_mean = 0
longitude_mean = 0


# create matrices to stick estimates in
coyote_log_lambda = matrix(NA, n.samples, length(hdens_data))

# Lambda predictions for dominant species
# Sample from posterior for the sequence of values for dominant species 
for (i in 1:n.samples){
  for (j in 1:length(dist_data)){
    # create the linear predictors for dominant species
    coyote_log_lambda[i,j] = b[,'beta0D'][[i]] + b[,'betaD[1]'][[i]] * dist_mean + b[,'betaD[2]'][[i]] * latitude_mean + b[,'betaD[3]'][[i]] * longitude_mean + b[,'betaD[4]'][[i]] * hdens_data[j]
  }}

# create array 
coyote_lambda = matrix(NA, n.samples, length(hdens_data))

# transform lambda off log-scale
coyote_lambda <- exp(coyote_log_lambda)

# calculate means and credible intervals
coyote_hdens_lambda_means = colMeans(coyote_lambda)
coyote_hdens_lambda_CIs <- apply(coyote_lambda,2,quantile, c(0.025,0.975), na.rm=TRUE)

# stuff into df
coyote_hdens_preds<- data.frame(housing_density = hdens_data, 
                               predicted = coyote_hdens_lambda_means, 
                               species = "coyote",
                               LCI = coyote_hdens_lambda_CIs[1,],
                               UCI = coyote_hdens_lambda_CIs[2,])

# but now if we are interested in subordinate species expected abundance as function of housing density, then we need realized abundance for dominant species at
# the covariate values of interest (N_D)
# create an array to stick realized abundance of dominant species (N_D) in
n_D_samples = matrix(NA, n.samples, length(hdens_data))

# draw from poisson distribution using posterior lambda estimates
for (i in 1:n.samples){
  for (j in 1:length(dist_data)){
    n_D_samples[i,j] <- rpois(n = 1, lambda=coyote_lambda[i,j])
  }
}

# create matrix for estimates 
log_lambda_bobcat_w_coyote = matrix(NA, n.samples, length(hdens_data))

# sample from posterior for linear predictors for subordinate species including impacts of dominant species 
for (i in 1:n.samples){
  for (j in 1:length(hdens_data)){
    # create the linear predictors for first order natural parameters
    log_lambda_bobcat_w_coyote[i,j] = b[,'beta0S'][[i]]  + b[,'betaS[1]'][[i]] * dist_mean + b[,'betaS[2]'][[i]] * latitude_mean +
      b[,'betaS[3]'][[i]] * longitude_mean + b[,'betaS[4]'][[i]] * hdens_data[j] + b[,'gamma0DS'][[i]] * n_D_samples[i,j]
  }}

# transform lambda off log-scale
bobcat_hdens_lambda_w_coyote = exp(log_lambda_bobcat_w_coyote)

# calculate means and 95% credible intervals
bobcat_hdens_lambda_w_coyote_means = colMeans(bobcat_hdens_lambda_w_coyote)
bobcat_hdens_lambda_w_coyote_CIs <- apply(bobcat_hdens_lambda_w_coyote,2,quantile, c(0.025,0.975), na.rm=TRUE)


# package up the bobcat lambda preds with the sequence of housing density values
bobcat_hdens_preds_w_coyote <- data.frame(housing_density = hdens_data, 
                                         predicted = bobcat_hdens_lambda_w_coyote_means, 
                                         species = "bobcat",
                                         LCI = bobcat_hdens_lambda_w_coyote_CIs[1,],
                                         UCI = bobcat_hdens_lambda_w_coyote_CIs[2,])

# combine estimates for both species
bothspeices_dist_preds <- rbind(coyote_hdens_preds, bobcat_hdens_preds_w_coyote)

# plot their expected abundance
ggplot(bothspeices_dist_preds, aes(x = housing_density, y = predicted, colour = species, fill = species)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ ylab(bquote("Predicted expected abundance ("*lambda~")"))+ xlab("Housing density (scaled and centered)")+
  theme_classic()


## # # # # # # # # # # # # # # # # # # # # # # # #  WORKSHOP CHALLENGE  # # # # # # # # # # # # # # # # 

###### CHALLENGE 04 Produce predictions for coyote expected abundance as a function of disturbance

n.samples <- nrow(b)
r <- range(dist_cov)
dist_data <- seq(r[1],r[2], length.out = 500)

hdens_mean <- 0
latitude_mean <- 0
longitude_mean <- 0
for (i in 1:n.samples){
  for (j in 1:length(dist_data)){
    # create the linear predictors for dominant species
    coyote_log_lambda[i,j] <- b[,'beta0D'][[i]] + 
      b[,'betaD[1]'][[i]] * dist_data[j] + 
      b[,'betaD[2]'][[i]] * latitude_mean + 
      b[,'betaD[3]'][[i]] * longitude_mean + 
      b[,'betaD[4]'][[i]] * hdens_mean
  }}

# create array 
coyote_lambda = matrix(NA, n.samples, length(hdens_data))

# transform lambda off log-scale
coyote_lambda <- exp(coyote_log_lambda)

# calculate means and credible intervals
coyote_dist_lambda_means = colMeans(coyote_lambda)
coyote_dist_lambda_CIs <- apply(coyote_lambda,2,quantile, c(0.025,0.975), na.rm=TRUE)

# stuff into df
coyote_dist_preds<- data.frame(disturbance = dist_data, 
                                predicted = coyote_dist_lambda_means, 
                                species = "coyote",
                                LCI = coyote_dist_lambda_CIs[1,],
                                UCI = coyote_dist_lambda_CIs[2,])

ggplot(coyote_dist_preds, aes(x = disturbance, y = predicted, colour = species, fill = species)) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.1, linetype = "dashed")+
  geom_path(size=1)+ ylab(bquote("Predicted expected abundance ("*lambda~")"))+ xlab("Disturbance (scaled and centered)")+
  theme_classic()
###### CHALLENGE 05 Produce predictions for bobcat expected abundance as a function of disturbance





###### CHALLENGE 06 specify and run a three species abundance mediated interaction model with abundant mediated interaction terms between C->B, C->RF, and B->RF

# Single chain run

# specify model
NimModel_3spp <- nimbleCode({
  #Model for coyote
  # Likelihood
  # State model for abundance of Coyotes
  for (i in 1:nsites) {
    n_coy[i] ~ dpois(lambda_coy[i])
    log(lambda_coy[i]) <-  beta0_coy + 
      beta_coy[1] * dist[i] + 
      beta_coy[2] * longitude[i] + 
      beta_coy[3] * latitude[i] + 
      beta_coy[4] * Hdens[i]
  }
  
  #Observation model for detection/non-detection data of coyotes
  for (i in 1:nsites) {
    for (j in 1:nsurveys){
      logit(r_coy[i,j]) <- alpha0_coy + alpha_coy[1] * trail[i,j] 
      p_coy[i,j] <- 1 - (1 - r_coy[i,j])^n_coy[i] #prob detecting at least 1 ind per day
      y_coy.new[i,j] ~ dbinom(p_coy[i,j],size=K)
    }
  }
  
  #Model for subordinate species (S)
  # Likelihood
  # State model for abundance of bobcats
  for (i in 1:nsites) {
    nS[i] ~ dpois(lambda_bob[i])
    log(lambda_bob[i]) <-  beta0_bob + 
      beta_bob[1] * dist[i] + 
      beta_bob[2] * longitude[i] + 
      beta_bob[3] * latitude[i] + 
      beta_bob[4] * Hdens[i] + 
      gamma_CoyBob * n_coy[i]
  }
  
  #Observation model for detection/non-detection data of bobcats
  for (i in 1:nsites) {
    for (j in 1:nsurveys){
      logit(r_bob[i,j]) <- alpha0_bob + alpha_bob[1] * trail[i,j] 
      p_bob[i,j] <- 1 - (1 - r_bob[i,j])^n_bob[i] #prob detecting at least 1 ind per day
      y_bob.new[i,j] ~ dbinom(p_bob[i,j],size=K)
    }
  }
  
  #Model for subordinate species (S)
  # Likelihood
  # State model for abundance of bobcats
  for (i in 1:nsites) {
    nS[i] ~ dpois(lambda_bob[i])
    log(lambda_bob[i]) <-  beta0_bob + 
      beta_bob[1] * dist[i] + 
      beta_bob[2] * longitude[i] + 
      beta_bob[3] * latitude[i] + 
      beta_bob[4] * Hdens[i] + 
      gamma_CoyBob * n_coy[i]
  }
  
  #Observation model for detection/non-detection data of bobcats
  for (i in 1:nsites) {
    for (j in 1:nsurveys){
      logit(r_bob[i,j]) <- alpha0_bob + alpha_bob[1] * trail[i,j] 
      p_bob[i,j] <- 1 - (1 - r_bob[i,j])^n_bob[i] #prob detecting at least 1 ind per day
      y_bob.new[i,j] ~ dbinom(p_bob[i,j],size=K)
    }
  }
  # --- Priors ---
  #Priors for dominant species (D)
  beta0D ~ dnorm(0, 0.1) 
  alpha0D ~ dnorm(0, 0.1) 
  for (k in 1:1) {
    alphaD[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:4) {
    betaD[k] ~ dnorm(0, 0.1)
  }
  #Priors for subordinate species (s)
  beta0S ~ dnorm(0, 0.1) 
  alpha0S ~ dnorm(0, 0.1) 
  for (k in 1:1) {
    alphaS[k] ~ dnorm(0, 0.1)
  }
  for (k in 1:4) {
    betaS[k] ~ dnorm(0, 0.1)
  }
  gamma0DS ~ dnorm(0, 0.1)
  
} )



# Still need to keep working on this
#Also want to figure out how to add in lag effect


## Interpret the model to test the hypotheses - abundant bobcat populations limit red fox populations in eastern US ##