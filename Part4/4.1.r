# Load necessary libraries
library(Matrix)  # For sparse and diagonal matrix operations

# Load data from CSV files
X1 <- read.csv('data/rain1.csv')
X2 <- read.csv('data/rain2.csv')
X3 <- read.csv('data/rain3.csv')
X4 <- read.csv('data/rain4.csv')

# Combine data from all events into a single data frame
data <- rbind(X1)

K_states <- 4

N <- length(data$minutes) # dependent on what experiment we use

B <- matrix(c(1,0,0,0),nrow=K_states) # column
C <- matrix(c(0,0,0,1),nrow=1) # row 

# Initialize the SigmaX
SigmaX <- diag(K_states)*1000

negloglik <- function(prm) {
  a <- prm[1]
  sigma1 <- prm[2]
  sigma2 <- prm[3]
  
  # System matrices based on parameters
  A <- matrix(c(1-a, 0, 0, 0,
                a, 1-a, 0, 0,
                0, a, 1-a, 0,
                0, 0, a, 0.98), nrow=K_states, byrow=TRUE)
  Sigma1 <- diag(rep(sigma1, K_states))
  Sigma2 <- matrix(sigma2, nrow=1, ncol=1)
  
  # Initialize state and likelihood
  X <- matrix(0, nrow=K_states, ncol=N)
  X0 <- matrix(0, nrow=K_states, ncol=1)
  X[,1] <- X0
  lik <- numeric(N)
  
  i <- 1
  # Kalman filter loop
  while (i < N) {
    # Prediction step
    X_pred <- A %*% X[,i] + B %*% data$u[i] + chol(Sigma1) %*% rnorm(K_states, 0, 1)
    SigmaX <- A %*% SigmaX %*% t(A) + Sigma1
    SigmaY <- C %*% SigmaX %*% t(C) + Sigma2
    
    i <- i + 1
    # Likelihood
    lik[i] <- dnorm(data$y[i], mean = C %*% X_pred, sd = sqrt(SigmaY), log = TRUE)

    # Update step
    innov <- data$y[i] - C %*% X_pred
    SigmaY <- C %*% SigmaX %*% t(C) + Sigma2
    K <- SigmaX %*% t(C) %*% solve(SigmaY)
    X[,i] <- X_pred + K %*% innov
    SigmaX <- SigmaX - K %*% C %*% SigmaX
  }
  
  # Negative log-likelihood
  -sum((lik[!is.na(lik)]))
}

params <- list(a = 0.04, sigma1 = 0.1, sigma2 = 0.5)
lower.bound <- c(0.01, 0.001, 0.1)
upper.bound <- c(0.1, 1, 5)

# Test it
opt_res <- optim(par = c(params$a, params$sigma1, params$sigma2), 
                fn = negloglik, method = "L-BFGS-B",
                lower = lower.bound, upper = upper.bound)

print(opt_res)
