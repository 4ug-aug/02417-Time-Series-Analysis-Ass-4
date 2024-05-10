setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)));getwd() # automatically sets WD to repo-base
set.seed(69)

# Load data from CSV files
X1 <- read.csv('data/rain1.csv')
X2 <- read.csv('data/rain2.csv')
X3 <- read.csv('data/rain3.csv')
X4 <- read.csv('data/rain4.csv')

for (data in list(X1)) {
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
    Y_residuals <- numeric(N)
    
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

      Y_residuals[i] <- innov
    }
  
    # Negative log-likelihood
    return(Y_residuals)
  }

  # 2
  # 0.03942545 0.04444793 1.541325
  
    params <- list(a = 0.03942545, sigma1 = 0.04444793, sigma2 = 1.541325)
    Y_residuals <- negloglik(c(params$a, params$sigma1, params$sigma2))
}

# Plot of residuals
par(mfrow = c(2,2))
plot(Y_residuals, type = 'l', col = '#21378e', xlab = 'Time', ylab = 'Residuals', main = 'Residuals of Kalman Filter')
hist(Y_residuals, col = '#21378e', xlab = 'Residuals', main = 'Histogram of Residuals')
# QQ plot
qqnorm(Y_residuals, col = '#21378e', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles', main = 'QQ Plot of Residuals')
qqline(Y_residuals, col = 'red')

# ACF plot
acf(Y_residuals, col = '#21378e', xlab = 'Lag', ylab = 'ACF', main = 'ACF of Residuals')
