setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)));getwd() # automatically sets WD to repo-base

b1 <- read.csv("data/rain1.csv")
set.seed(69)

# parameters
a <- 0.04
sigma1 <- 0.1
sigma2 <- 0.1
K_states <- 4

N <- length(b1$minutes) # dependent on what experiment we use

# initial state
X <- matrix(nrow=K_states,ncol=N)
X0 <- c(0,0,0,0) # column vector
X[,1] <- X0
B <- matrix(c(1,0,0,0),nrow=K_states) # column
C <- matrix(c(0,0,0,1),nrow=1) # row 

kalman_filter <- function(input, a, sigma1, sigma2){
  A <- matrix(c(1-a,0,0,0,
                a,1-a,0,0,
                0,a,1-a,0,
                0,0,a,0.98),nrow=4, byrow=TRUE)
  Sigma1 <- diag(4)*sigma1
  Sigma2 <- matrix(sigma2)
  Y <- numeric(N)
  Y[1] <- C%*%X[,1] + sqrt(Sigma2) %*% rnorm(1) # scalar (the level of water in basin)
  # running N steps Kalman filtering with input u_t (rainfall)
  for (I in 2:N){
    X[,I] <- A %*% X[,I-1,drop=FALSE] + B%*%input[I-1] + chol(Sigma1) %*% matrix(rnorm(4),ncol=1)
    Y[I] <- C %*% X[,I] + sqrt(Sigma2) %*% rnorm(1)
  }
  return(list(X = X, Y = Y))
}

as <- c(0.04, 0.02, 0.06)
sigma1s <- c(0.1, 0.01, 0.5)
sigma2s <- c(0.1, 0.01, 1)

a_sims <- vector("list", 3)
sigma1_sims <- vector("list", 3)
sigma2_sims <- vector("list", 3)

for (i in 1:3) {
  a_sims[[i]] <- kalman_filter(b1$u, as[i], sigma1, sigma2)
} # seems to work. allocates all 4 states into $X and output into $Y for each a-value.
for (i in 1:3) {
  sigma1_sims[[i]] <- kalman_filter(b1$u, a, sigma1s[i], sigma2)
} 
for (i in 1:3) {
  sigma2_sims[[i]] <- kalman_filter(b1$u, a, sigma1, sigma2s[i])
}

### TODO: MAKE MANY PLOTS WITH DIFFERENT PARAMETER VALUES in kalman_filter. see 1_1 with more compact plotting
kalman_plot <- function(data, simul, params, position, bottom, y_max){
  if (!(position %in% c("left", "right", "middle"))) {
    stop("Invalid side argument. Choose 'left', 'right' or 'middle'.")
  }
  
  plot(data$minutes, data$y, col="blue", type="l", 
       xlab="", ylab="", yaxt="n", xaxt="n",
       ylim=c(0, y_max))
  
  lw=3;lt=4
  lines(data$minutes,simul$X[1,], col="#946ca3", lwd=lw, lty=lt)
  lines(data$minutes,simul$X[2,], col="#da739e", lwd=lw, lty=lt)
  lines(data$minutes,simul$X[3,], col="#ffa760", lwd=lw, lty=lt)
  lines(data$minutes,simul$X[4,], col="#64aa3c", lwd=lw, lty=lt)
  lines(data$minutes,simul$Y, col="#573500", lwd=lw, lty=lt)

  if (position == "left"){
    #mtext(expression(paste("(100 ", m^3, ")")), side=2, line=2.5, col="black")
    axis(side=2, col.axis="black", las=1)
  }
  
  title <- bquote(.(title) ~ "a:" ~ .(params[1]) ~ sigma[1]: .(params[2]) ~ sigma[2]: .(params[3]))
  
  par(new=TRUE)
  plot(data$minutes, data$u, 
       col="red", lwd=1, type="l",
       xlab="", ylab="", yaxt="n", xaxt="n",
       main=title, 
       ylim=c(0, y_max/10))
  if (position == "right") {
    #mtext(expression(paste("Rainfall 'u' (100 ", m^3, "/min)")), side=4, line=2.5, col="red")
    axis(side=4, col.axis="red", las=1)
  } else if (position == "middle") {
  }
  if (bottom == TRUE){
    axis(side=1, col.axis="black", las=1)
  }
}


### plotting
y_max <- max(b1$y)
old_mar <- par("mar")
par(mfrow=c(3,3), oma=c(0,5,0,3), mar=c(2.1, 0, 1.5, 1), mgp = c(5, 0.5, 0))

# right side will be the water level for, left will be input water rain
kalman_plot(b1, a_sims[[1]], params=c(as[1],sigma1, sigma2), 
            "left", bottom=FALSE, y_max = y_max)
kalman_plot(b1, a_sims[[2]], params=c(as[2],sigma1, sigma2),
            "middle", bottom=FALSE, y_max = y_max)
kalman_plot(b1, a_sims[[3]], params=c(as[3],sigma1, sigma2),
            "right", bottom=FALSE, y_max = y_max)

kalman_plot(b1, sigma1_sims[[1]], params=c(a,sigma1s[1], sigma2),
            "left", bottom=FALSE, y_max = y_max)
kalman_plot(b1, sigma1_sims[[2]], params=c(a,sigma1s[2], sigma2),
            "middle", bottom=FALSE, y_max = y_max)
kalman_plot(b1, sigma1_sims[[3]], params=c(a,sigma1s[3], sigma2),
            "right", bottom=FALSE, y_max = y_max)

kalman_plot(b1, sigma2_sims[[1]], params=c(a,sigma1, sigma2s[1]),
            "left", bottom=TRUE, y_max = y_max)
kalman_plot(b1, sigma2_sims[[2]], params=c(a,sigma1, sigma2s[2]),
            "middle", bottom=TRUE, y_max = y_max)
kalman_plot(b1, sigma2_sims[[3]], params=c(a,sigma1, sigma2s[3]),
            "right", bottom=TRUE, y_max = y_max)


mtext("Time (minutes)", side=1, outer=TRUE, line=2.5)
mtext(expression(paste("(100 ", m^3, ")")), side=2, outer=TRUE, line=2.5, col="black")
mtext(expression(paste("Rainfall 'u' (100 ", m^3, "/min)")), side=4, outer=TRUE, line=2.5, col="red")
lw=3;lt=4
legend("topright", 
       legend = c("rainfall", "water level", "State 1", "State 2", "State 3", "State 4", "Output"), 
       col = c("red","blue","#946ca3", "#da739e", "#ffa760", "#64aa3c", "#573500"), 
       lwd = c(1, 1, lw, lw, lw, lw, lw),
       lty = c(1, 1, lt, lt, lt, lt, lt),
       cex=0.7)

