# load R modules

b1 <- read.csv("data/rain1.csv")
b2 <- read.csv("data/rain2.csv")
b3 <- read.csv("data/rain3.csv")
b4 <- read.csv("data/rain4.csv")

# We plot u (measured incoming rain) and y (measured water level) for each basin
par(mfrow=c(2,2))
# We plot two lines in each plot
plot(b1$minutes, b1$u, type="l", col="red", 
        xlab="Time (minutes)", ylab="u (100m3/min)", 
        main="Basin 1", ylim=c(0, 10))
lines(b1$minutes, b1$y, col="blue")
legend("topright", legend=c("Rain (u)", "Water level (y)"), col=c("red", "blue"), lty=1, cex=0.8)

plot(b2$minutes, b2$u, type="l", col="red", 
        xlab="Time (minutes)", ylab="u (100m3/min)", 
        main="Basin 2", ylim=c(0, 10))
lines(b2$minutes, b2$y, col="blue")
legend("topright", legend=c("Rain (u)", "Water level (y)"), col=c("red", "blue"), lty=1, cex=0.8)

plot(b3$minutes, b3$u, type="l", col="red", xlab="Time (minutes)", ylab="u (100m3/min)", main="Basin 3")
lines(b3$minutes, b3$y, col="blue")
legend("topright", legend=c("Rain (u)", "Water level (y)"), col=c("red", "blue"), lty=1, cex=0.8)

plot(b4$minutes, b4$u, type="l", col="red", xlab="Time (minutes)", ylab="u (100m3/min)", main="Basin 4")
lines(b4$minutes, b4$y, col="blue")
legend("topright", legend=c("Rain (u)", "Water level (y)"), col=c("red", "blue"), lty=1, cex=0.8)


