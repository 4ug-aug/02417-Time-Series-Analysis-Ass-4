setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)));getwd() # automatically sets WD to repo-base
set.seed(69)

# Load data from CSV files
X1 <- read.csv('data/rain1.csv')
X2 <- read.csv('data/rain2.csv')
X3 <- read.csv('data/rain3.csv')
X4 <- read.csv('data/rain4.csv')

# df for likelihood parameters
df <- data.frame(a = numeric(), sigma1 = numeric(), sigma2 = numeric())

for (data in list(X1, X2, X3, X4)) {
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
  
  params <- list(a = 0.04, sigma1 = 0.025, sigma2 = 0.5)
  lower.bound <- c(0.01, 0.001, 0.1)
  upper.bound <- c(0.1, 1, 5)

  # Test it
  opt_res <- optim(par = c(params$a, params$sigma1, params$sigma2), 
                  fn = negloglik, method = "L-BFGS-B",
                  lower = lower.bound, upper = upper.bound)

  df <- rbind(df, data.frame(a = opt_res$par[1], sigma1 = opt_res$par[2], sigma2 = opt_res$par[3], nll = opt_res$value))
}

df

kalman_plot <- function(data, simul, title, params, position, bottom, x_max, y_max){
  if (!(position %in% c("left", "right", "middle"))) {
    stop("Invalid side argument. Choose 'left', 'right' or 'middle'.")
  }
  
  plot(data$minutes, data$y, col="blue", type="l", 
       xlab="", ylab="", yaxt="n", xaxt="n",
       ylim=c(0, y_max),
       xlim=c(0, x_max))
  
  lw=3;lt=4
  lines(data$minutes,simul$X[1,], col="#946ca3", lwd=lw, lty=lt)
  lines(data$minutes,simul$X[2,], col="#da739e", lwd=lw, lty=lt)
  lines(data$minutes,simul$X[3,], col="#ffa760", lwd=lw, lty=lt)
  lines(data$minutes,simul$X[4,], col="#64aa3c", lwd=lw, lty=lt)
  lines(data$minutes,simul$Y, col="#573500", lwd=lw, lty=lt)
  
  if (position == "left"){
    axis(side=2, col.axis="black", las=1)
  }
  
  title <- bquote( ~ "a:" ~ .(params[1]) ~ sigma[1]: .(params[2]) ~ sigma[2]: .(params[3]))
  
  par(new=TRUE)
  plot(data$minutes, data$u, 
       col="red", lwd=1, type="l", 
       xlab="", ylab="", yaxt="n", xaxt="n", 
       main=title, 
       ylim=c(0, y_max/10))
  if (position == "right") {
    axis(side=4, col.axis="red", las=1)
  } else if (position == "middle") {
  }
  if (bottom == TRUE){
    axis(side=1, col.axis="black", las=1)
  }
}
