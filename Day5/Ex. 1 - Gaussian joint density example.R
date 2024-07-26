
# Montpellier short-course on integrated models, Kéry, Schaub & Strebel, April 2023

# Example 1: data assumed to be iid (independent and identically distributed)

# Joint density is the basic start of any parametric statistical model
# when doing maximum likelihood or Bayesian posterior inference

# Imagine measurements of observed avian species richness at five sites,
# Want to estimate mean and sd of species richness in the statistical population of sites
# from which these 5 sites are sampled. We ignore imperfect detection and discreteness of 
# these data and we decide to work with a Normal density.

# Simulate some data: 5 data points with imaginary species counts
set.seed(1)
(y <- round(rnorm(5, mean = 50, sd = 10)))
# [1] 44 52 42 66 53

# Save truth
truth <- c("mean" = 50, "sd" = 10)

library(ASMbook)

# Define (negative log-)likelihood function
# (= negative of log of joint density of data, viewed as a function of params)
NLL <- function(param, y) {
  mn <- param[1]      # First parameter is mean (mu)
  sd <- param[2]        # Second is SD (sigma)
  -1* sum(dnorm(y, mean = mn, sd = sd, log=TRUE))
}

# (Try to) Get MLEs for single datum (the first one)
out1 <- optim(par = c(mn = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1])
MLE <- out1$par                         # MLEs
ASE <- sqrt(diag(solve(out1$hessian)))  # Asymptotic standard errors
print(cbind(MLE, ASE), 4)
# Numerical error, doesn't work, can't estimate 2 parameters from a single datum

# (Try to) Get MLEs for two data points
out2 <- optim(par = c(mn = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:2])
MLE <- out2$par                         # MLEs
ASE <- sqrt(diag(solve(out2$hessian)))  # Asymptotic standard errors
print(cbind(MLE, ASE), 4)

## note 'get_MLE' function does same and also gives Wald 95% CIs
get_MLE(out2, 4)

# Get MLEs for three data points
out3 <- optim(par = c(mn = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:3])
get_MLE(out3, 3)

# Get MLEs for four data points
out4 <- optim(par = c(mn = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:4])
get_MLE(out4, 3)

# Get MLEs for all five data points
out5 <- optim(par = c(mn = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:5])
get_MLE(out5, 3)



# See already some main features of 'integrated models' (though note: this is NOT an IM)
# - Cannot even estimate parameters when we don't have enough data (i.e., here with 1 datum)
# - Get better estimates when we include more data in the analysis
# - point estimate closer to target, smaller SEs



# Now imagine that your colleague brings you another set of counts of species
# which you don't have reason a priori to believe they stem from a different 'process'
# (e.g., imagine they're from the same region as the first set of data)
# Let's simulate this 2nd data set and then include it in the analysis

set.seed(2)
(y.add <- round(rnorm(5, mean = 50, sd = 10)))
# [1] 41 52 66 39 49


# Get MLEs for both data sets together. Is this an integrated model ?
# Not really, perhaps, because data sets not disparate
out10 <- optim(par = c(mn = 0, sd = 1), fn = NLL, hessian = TRUE, y = c(y, y.add))
get_MLE(out10, 4)

# Remember the truth
truth

