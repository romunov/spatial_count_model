# a script with examples from the book (see README on my github page) using Metropolis-Hastings sampler
# from page 443

LRMH <- function(formula, mu0, sigma0, theta, delta, niter, data) {
  out <- rep(NA, niter)
  get.data <- all.vars(formula)
  y <- data[, get.data[1]]
  x <- data[, get.data[2]]
  
  for (i in 1:niter) {
    # this part should be generalized to accomodate more than one x, probably
    # something along the lines of %*% and counting the number of explanatory
    # variables
    theta.cand <- rnorm(1, mean = theta, sd = delta)
    z.cand <- theta.cand * x # specify linear predictor a la GLM
    z <- theta * x
    
    # calculate likelihoods and priors
    loglik.cand <- sum(dbinom(y, size = 1, prob = exp(z.cand)/(1 + exp(z.cand)), log = TRUE))
    logprior.cand <- dnorm(theta.cand, mean = mu0, sd = sigma0, log = TRUE)

    loglik <- sum(dbinom(y, size = 1, prob = exp(z)/(1 + exp(z)), log = TRUE))
    logprior <- dnorm(theta, mean = mu0, sd = sigma0, log = TRUE)
    
    ratio <- exp((loglik.cand + logprior.cand) - (loglik + logprior))
    
    if (runif(1) < ratio) {
      theta <- theta.cand
    }
    
    out[i] <- theta
  }
  
  out
}

# simulate some data for logistic regression where response variables
# takes 1 with P(theta)

set.seed(357)
x <- rnorm(100)
z <- 5 * x
theta <- 1/(1 + exp(-z))

y <- rbinom(100, size = 1, prob = theta)

sim.data <- data.frame(y = y, x = x)
plot(sim.data)

# analyze using a "classical" GLM
summary(glm(y ~ -1 + x, data = sim.data, family = binomial))

# analyze using a MH algo
mh.res <- LRMH(y ~ x, data = sim.data, mu0 = 0, sigma0 = 3, theta = 3, delta = 3, niter = 50000)
hist(mh.res, xlab = "theta")
