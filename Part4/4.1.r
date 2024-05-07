# Load necessary libraries
library(Matrix)  # For sparse and diagonal matrix operations

# Load data from CSV files
X1 <- read.csv('data/rain1.csv', header=FALSE)[-1,]
X2 <- read.csv('data/rain2.csv', header=FALSE)[-1,]
X3 <- read.csv('data/rain3.csv', header=FALSE)[-1,]
X4 <- read.csv('data/rain4.csv', header=FALSE)[-1,]

# Combine data from all events into a single data frame
data <- rbind(X1, X2, X3, X4)

# Extract rain inputs (ut) and water levels (yt)
ut <- as.numeric(data[, 1])  # Assuming first column is ut
yt <- as.numeric(data[, 2])  # Assuming second column is yt

# Parameters
a <- 0.05  # Transition rate
sigma1 <- 0.01  # Variance of system noise
sigma2 <- 1  # Variance of observation noise

# State-space matrices
A <- matrix(c(1-a, 0, 0, 0,
              a, 1-a, 0, 0,
              0, a, 1-a, 0,
              0, 0, a, 0.98), nrow=4, byrow=TRUE)
B <- matrix(c(1, 0, 0, 0), nrow=4, ncol=1)
C <- matrix(c(0, 0, 0, 1), nrow=1, ncol=4)

# Initial states and covariance matrices
X <- matrix(rep(0, 4), ncol=1)  # Initial state
SigmaX <- diag(4)  # Initial state covariance, large values to represent uncertainty

# Function to generate G(Xt) matrix
get_G <- function(X) {
  s <- sqrt(abs(X))
  diag(as.vector(s))
}

# Kalman filtering process
n <- length(ut)
predicted_levels <- numeric(n)
estimated_levels <- numeric(n)

for (i in 2:n) {
  # Predict step
  e1_t <- rnorm(4, mean = 0, sd = sqrt(sigma1))
  G_t_minus_1 <- get_G(X)
  system_noise <- G_t_minus_1 %*% e1_t
  X <- A %*% X + B %*% ut[i] + system_noise
  SigmaX <- A %*% SigmaX %*% t(A) + G_t_minus_1 %*% diag(rep(sigma1, 4)) %*% t(G_t_minus_1)
  predicted_levels[i] <- as.numeric(C %*% X)
  
  # Update step
  e2_t <- rnorm(1, mean = 0, sd = sqrt(sigma2))
  Y_t <- yt[i] + e2_t
  innov <- Y_t - C %*% X
  S <- C %*% SigmaX %*% t(C) + sigma2
  K <- SigmaX %*% t(C) %*% solve(S)
  X <- X + K %*% innov
  SigmaX <- SigmaX - K %*% C %*% SigmaX
  estimated_levels[i] <- as.numeric(C %*% X)
}

# Results
plot(predicted_levels, type = 'l', col = 'blue', ylim = range(c(predicted_levels, estimated_levels)), ylab = 'Water Level', main = 'Predicted vs Estimated Water Levels')
lines(estimated_levels, col = 'red')
legend("topright", legend=c("Predicted", "Estimated"), col=c("blue", "red"), lty=1)
