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
  X[,1] <- X0
  lik <- numeric(N)
  
  # Kalman filter loop
  for (i in 2:N) {
    # Prediction step
    X_pred <- A %*% X[,i-1]
    SigmaX_pred <- A %*% SigmaX %*% t(A) + Sigma1
    
    # Update step
    innov <- data$y[i] - C %*% X_pred
    S <- C %*% SigmaX_pred %*% t(C) + Sigma2
    K <- SigmaX_pred %*% t(C) %*% solve(S)
    X[,i] <- X_pred + K %*% innov
    SigmaX <- SigmaX_pred - K %*% C %*% SigmaX_pred
    
    # Calculate likelihood
    lik[i] <- dnorm(data$y[i], mean = C %*% X[,i], sd = sqrt(S))
  }
  
  # Negative log-likelihood
  -sum(log(lik[-1]))  # Skip the first likelihood as it's not initialized
}

params <- list(a = 0.04, sigma1 = 0.1, sigma2 = 0.5)
lower.bound <- c(0.01, 0.001, 0.1)
upper.bound <- c(0.1, 1, 5)

# Test it
opt_res <- optim(par = c(params$a, params$sigma1, params$sigma2), 
                fn = negloglik, method = "L-BFGS-B",
                lower = lower.bound, upper = upper.bound)

print(opt_res)