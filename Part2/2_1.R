setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)));getwd() # automatically sets WD to repo-base

b1 <- read.csv("data/rain1.csv")
b2 <- read.csv("data/rain2.csv")
b3 <- read.csv("data/rain3.csv")
b4 <- read.csv("data/rain4.csv")

length(b1$minutes)


# parameters
a <- 0.05
sigma1 <- 0.01
sigma2 <- 0.01
K_states <- 4

N <- length(b1$minutes) # dependent on what experiment we use

# initial state
X <- matrix(nrow=K_states,ncol=N)
X0 <- c(0,0,0,0) # column vector
X[,1] <- X0

A <- matrix(c(1-a,0,0,0,
              a,1-a,0,0,
              0,a,1-a,0,
              0,0,a,0.98),nrow=K_states, byrow=TRUE)
B <- matrix(c(1,0,0,0),nrow=K_states) # column
C <- matrix(c(0,0,0,1),nrow=1) # row 
Sigma1 <- diag(K_states)*sigma1
Sigma2 <- matrix(sigma2)

Y <- numeric(N)
Y[1] <- C%*%X[,1] + sqrt(Sigma2) %*% rnorm(1) # scalar (the level of water in basin)

# running N steps Kalman filtering with input u_t (rainfall)
for (I in 2:N){
  X[,I] <- A %*% X[,I-1,drop=FALSE] + B%*%b1$u[I-1] + chol(Sigma1) %*% matrix(rnorm(K_states),ncol=1)
  Y[I] <- C %*% X[,I] + sqrt(Sigma2) %*% rnorm(1)
}

### plotting
y_max <- max(b1$y)

old_mar <- par("mar")
par(mfrow=c(2,2), oma=c(0,0,0,3), mar=c(4.1, 4, 2, 1))

plot(b1$minutes, b1$u, 
     col="red", lwd=2, type="l",
     xlab="Time (minutes)", ylab="", yaxt="n",
     main="", ylim=c(0, y_max/10))
axis(side=2, col.axis="red", las=1)
mtext(expression(paste("u (100 ", m^3, "/min)")), side=2, line=2.5, col="red")
par(new=TRUE)
plot(b1$minutes,b1$y, col="blue", type="l", axes=FALSE, xlab="", ylab="")
axis(side=4, col.axis="blue", las=1)
mtext(expression(paste("y (100 ", m^3, ")")), side=4, line=2.5, col="blue")

lw=3;lt=4
lines(b1$minutes,X[1,], col="#946ca3", lwd=lw, lty=lt)
lines(b1$minutes,X[2,], col="#da739e", lwd=lw, lty=lt)
lines(b1$minutes,X[3,], col="#ffa760", lwd=lw, lty=lt)
lines(b1$minutes,X[4,], col="#64aa3c", lwd=lw, lty=lt)
legend("topright", 
       legend = c("rainfall", "water level", "State 1", "State 2", "State 3", "State 4"), 
       col = c("red","blue","#946ca3", "#da739e", "#ffa760", "#64aa3c"), 
       lwd = c(2, 1, lw, lw, lw, lw),
       lty = c(1, 1, lt, lt, lt, lt),
       cex=0.8)


