library(scrbook) # metropolis-hastings function
library(coda) # displaying diagnostics of the MCMC run

data(beardata)

M <- 700
trapmat <- beardata$trapmat
# Count the number of bears on a specific trap on all sampling
# occasions. 
bearmat <- apply(beardata$bearArray, 1:2, sum)

# Create an empty matrix and add actual data at the top. You can
# use image(t(Xaug)) to check if the data is incorporated in the 
# augmented dataset. If you used t(), the data should be displayed
# at the bottom of the image near the x axis.
Xaug <- matrix(0, nrow = M, ncol = nrow(trapmat))
Xaug[1:nrow(bearmat), ] <- bearmat

# Add buffer to the sampling grid.
xl <- min(trapmat[, 1]) - 20
yl <- min(trapmat[, 2]) - 20
xu <- max(trapmat[, 1]) + 20
yu <- max(trapmat[, 2]) + 20


set.seed(13)

mod0 <- SCR0binom.cl(y = Xaug, X = trapmat, M = M, xl = xl, xu = xu,
                      yl = yl, yu = yu, K = 8, delta = c(0.1, 0.05, 2), niter = 5000)


mod0.mcmc <- mcmc(mod0)
plot(window(mod0.mcmc, start = 2001))
effectiveSize(window(mod0.mcmc, start = 2001))
