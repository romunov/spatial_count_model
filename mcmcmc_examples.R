# a script with examples from the book (see README on my github page) using Metropolis-Hastings sampler
# from page 443

LRMH <- function(y, mu0, sigma0, theta, delta, niter) {
  out <- rep(NA, niter)
  
  for (i in 1:niter) {
    theta.cand <- rnorm(1, mean = theta, sd = delta)
    
    loglik.cand <- sum(dbinom(y, size = 1, prob = exp(theta.cand)/(1 + exp(theta.cand)), log = TRUE))
    logprior.cand <- dnorm(theta.cand, mean = mu0, sd = sigma0, log = TRUE)

    loglik <- sum(dbinom(y, size = 1, prob = exp(theta)/(1 + exp(theta)), log = TRUE))
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
summary(glm(y ~ -1 + x, data = sim.data, family = binomial)) # works

# analyze using a MH algo
mh.res <- LRMH(y = sim.data$y, mu0 = 3, sigma0 = 2, theta = 3, delta = 3, niter = 5000)
hist(mh.res, xlab = "theta") # wrong, I missed something
