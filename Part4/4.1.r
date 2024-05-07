# Load and combine data
X1 <- read.csv('data/rain1.csv', header=FALSE)
X2 <- read.csv('data/rain2.csv', header=FALSE)
X3 <- read.csv('data/rain3.csv', header=FALSE)
X4 <- read.csv('data/rain4.csv', header=FALSE)
X <- rbind(X1, X2, X3, X4)

# Constants and parameters
a <- 0.1
sigma1 <- 1
sigma2 <- 1
n <- nrow(X)  # Number of observations
ut <- X[, 1]  # Assuming first column is ut
yt <- X[, 2]  # Assuming second column is yt

# Matrices for state-space model
A <- matrix(c(1-a, 0, 0, 0,
              a, 1-a, 0, 0,
              0, a, 1-a, 0,
              0, 0, a, 0.98), nrow=4, byrow=TRUE)
B <- matrix(c(1, 0, 0, 0), nrow=4, byrow=TRUE)
C <- matrix(c(0, 0, 0, 1), nrow=1, byrow=TRUE)
Sigma1 <- diag(c(sigma1, sigma1, sigma1, sigma1))
Sigma2 <- diag(sigma2)

# Initial state
X1 <- matrix(c(0, 0, 0, 0), ncol=1)
SigmaX <- diag(c(sqrt(abs((X1[1, 1] - 0)^2)), sqrt(abs((X1[2, 1] - 0)^2)),
               sqrt(abs((X1[3, 1] - 0)^2)), sqrt(abs((X1[4, 1] - 0)^2)))

get_G <- function(X) {
  sqrt_abs_X <- sqrt(abs(X))
  G_matrix <- diag(sqrt_abs_X)
  return(G_matrix)
}


# Arrays to store predictions and estimates
ypred <- numeric(n)
yhat <- numeric(n)

# Kalman filter loop
for (i in 1:n) {
  # Prediction
  X <- A %*% X + B %*% ut[i]
  SigmaX <- A %*% SigmaX %*% t(A) + Sigma1
  SigmaY <- C %*% SigmaX %*% t(C) + Sigma2
  ypred[i] <- as.numeric(C %*% X)  # Keep the prediction
  
  # Update
  if (i < n) {
    K <- SigmaX %*% t(C) %*% solve(SigmaY)
    X <- X + K %*% (yt[i] - C %*% X)
    SigmaX <- SigmaX - K %*% SigmaY %*% t(K)
    yhat[i] <- as.numeric(C %*% X)  # Keep the estimate
  }
}
