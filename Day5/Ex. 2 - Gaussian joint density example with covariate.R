
# Montpellier short-course on integrated models, Kéry, Schaub & Strebel, April 2023

# Example 2: Assume now that species richness counts depend on a covariate,
#            taken at different elevations in a mountain range
#            (species richness usually declines with increasing elevation)


library(ASMbook)

# First simulate some data that look about right to illustrate joint density etc.
set.seed(1)
n <- 20              # sample size
(x <- round(runif(n, min = 200, max = 3000)))  # Imagine elevation in metres
# [1]  943 1242 1804 2743  765 2715 2845 2050 1962  373
(y <- round(rnorm(n, mean = 65 - 15 * (x/1000), sd = 10) ))
(truth <- c(65, -15, 10))

par(mar = c(6,6,5,4), cex.lab = 2, cex.axis = 2)
plot(x/1000, y, xlab = 'Elevation (in km)', ylab = 'Observed species richness', ylim = c(0, 100), 
  pch = 16, cex = 2)


# Define (negative log-)likelihood function of what is ordinary linear regression
# (= negative of log of joint density of data, viewed as a function of params)
NLL <- function(param, y, x) {
  alpha <- param[1]     # Intercept (expected SR at sea-level)
  beta <- param[2]      # Slope of SR on elevation
  sd <- param[3]        # 'Residual standard error'
  mu <- alpha + beta * (x / 1000)    # Linear predictor
  LL <- dnorm(y, mean = mu, sd = sd, log=TRUE) # numerically more stable
  # Next is same, but numerically unstable
  # L <- dnorm(y, mean = mu, sd = sd, log=FALSE) # Likelihood contribution of each datum
  # LL <- log(L)          # log-likelihood contribution of each datum
  negLL <- -sum(LL)       # neg. loglikelihood of dataset is minus the sum
  negLL
}



# (Try to) Get MLEs for single datum (the first one)
out1 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1], x = x[1])
get_MLE(out1, 4)
# Numerical error, doesn't work, can't estimate 3 parameters from a single datum

# (Try to) Get MLEs for two data points
out2 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:2], x = x[1:2])
get_MLE(out2, 4)
# Numerical error, doesn't work, can't estimate 3 parameters from two data points single datum

# Get MLEs for three data points
out3 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:3], x = x[1:3])
get_MLE(out3, 4)
# That works, numerically

# Get MLEs for 4, 8, 12, 16 and 20 points
out4 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:4], x = x[1:4])
get_MLE(out4, 4)

out8 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:8], x = x[1:8])
get_MLE(out8, 4)

out12 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:12], x = x[1:12])
get_MLE(out12, 4)

out16 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y[1:16], x = x[1:16])
get_MLE(out16, 4)

out20 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, hessian = TRUE, y = y, x = x)
get_MLE(out20, 4)


# now imagine again that your dear friend collected another data set with another 20 data points for you
set.seed(2)
(x.add <- round(runif(n, min = 200, max = 3000)))  # Imagine elevation in metres
(y.add <- round(rnorm(n, mean = 65 - 15 * (x.add/1000), sd = 10) ))
# Get one negative, but for sake of illustration this doesn't matter....

# As before, we decide to add this 2nd data set into our analysis
out20 <- optim(par = c(alpha = 50, beta = 0, sd = 1), fn = NLL, 
  hessian = TRUE, y = c(y, y.add), x = c(x, x.add))
get_MLE(out20, 4)

# Compare again with the truth
truth

# Although this is ALSO not an integrated model according to our definition (since the data sets are 
# NOT 'disparate'), it does nevertheless illustrate two important features that any joint density
# shares with the joint density of an integrated model
# - the density of each data point is different, but we formulate the joint density
#   as a function of parameters that are shared: these are the alpha and the beta and the sd
# - and then, as always: more data produce 'better' (more precise) estimates
# - in the extreme case where we have fewer data points than params, can't even estimate things

