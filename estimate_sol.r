## ----------------------------------------------------------------
# The falling body example from the book
## ----------------------------------------------------------------

## ----------------------------------------------------------------
# Parameters for simulation
A <- matrix(c(1,1,
              0,1), nrow=2, byrow=TRUE)
B <- matrix(c(-0.5,
                -1), nrow=2)
C <- matrix(c(1,0), nrow=1)
Sigma1 <- matrix(c(2,  0.8,
                   0.8,  1), nrow=2, byrow=TRUE)
Sigma2 <- 1000
# The state-vector init value
X0 <- c(10000, 0)
# Number of time points to simulate
n <- 100
# Vector for saving the output
y <- rep(NA, n)
# The input (i.e. gravity, hence constant input not changing over time)
g <- 9.82
## ----------------------------------------------------------------

## ----------------------------------------------------------------
# Run the simulation
X <- X0
for (i in 2:n){
    X <- A %*% X + B%*%g + chol(Sigma1) %*% matrix(rnorm(2),ncol=1)
    y[i] <- C %*% X + sqrt(Sigma2) %*% rnorm(1)
}

# Plot the simulated height
plot(y, type="b", xlab="time", ylab="Observed altitude [m]", ylim=range(y,na.rm=TRUE))
## ----------------------------------------------------------------


## ----------------------------------------------------------------
# For keeping the likelihoods
lik <- rep(NA,n)
# Initialize the state
X <- X0
SigmaX <- matrix(c(1000,   0,
                      0,1000), nrow=2, byrow=TRUE)
# Start on the iterations
i <- 1
#while(i < n){
# ---- PREDICT ----
X <- A %*% X + B %*% g
SigmaX <- A %*% SigmaX %*% t(A) + Sigma1
SigmaY <- C %*% SigmaX %*% t(C) + Sigma2
# ---- Step ----
i <- i + 1
# The likelihood
HERE CALCULATE the likelihood and save it to lik[i]
# ---- UPDATE ----
K <- SigmaX %*% t(C) %*% solve(SigmaY)
X <- X + K %*% (y[i] - C %*% X)
SigmaX <- SigmaX - K %*% SigmaY %*% t(K)
#}

# Finally, calculate the negative log-likelihood
CALCULATE HERE the negative log-likelihood and sum them
## ----------------------------------------------------------------


## ----------------------------------------------------------------
negloglik <- function(prm){
  # prm is a vector, so make it back to B
  B <- matrix(prm[c("B1","B2")], nrow=2)
  #
  PUT THE negative log-likelihood calculation here
}
# Test it
# Make B a vector
prm <- c(B=matrix(c(-0.5,-1), nrow=2))
prm
negloglik(prm)
## ----------------------------------------------------------------


## ----------------------------------------------------------------
# Estimate the B matrix coefficients
# We can now pass it to an optimizer to calculate the maximum likelihood estimates
val <- optim(prm, negloglik, control=list(trace=0), method="L-BFGS-B")
val
## ----------------------------------------------------------------