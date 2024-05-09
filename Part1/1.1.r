setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)));getwd() # automatically sets WD to repo-base

# load R modules
b1 <- read.csv("data/rain1.csv")
b2 <- read.csv("data/rain2.csv")
b3 <- read.csv("data/rain3.csv")
b4 <- read.csv("data/rain4.csv")


# plot helpers
# limits
y_min <- 0
y_max <- max(c(b1$y, b2$y, b3$y, b4$y))
x_min <- 0
x_max <- max(c(b1$minutes, b2$minutes, b3$minutes, b4$minutes))

double_axis_plot <- function(data, title, side, x_max, y_max){
  if (side != "left" && side != "right") {
    stop("Invalid side argument. Choose 'left' or 'right'.")
  }
  
  plot(data$minutes, data$u, type="l", col="red",
       xlab="Time (minutes)", ylab="", yaxt="n",
       main=title, 
       ylim=c(0, y_max/10),
       xlim=c(0, x_max))
  axis(side=2, col.axis="red", las=1)
  par(new=TRUE)
  plot(data$minutes, data$y, type="l", col="blue", axes=FALSE, xlab="", ylab="", 
       ylim=c(0, y_max),
       xlim=c(0, x_max))
  axis(side=4, col.axis="blue", las=1)
  legend("topright", legend=c("Rain (u)", "Water level (y)"), col=c("red", "blue"), lty=1, cex=0.8)
  
  if (side == "left"){
    mtext(expression(paste("u (100 ", m^3, "/min)")), side=2, line=2.5, col="red")
  } else if (side == "right") {
    mtext(expression(paste("y (100 ", m^3, ")")), side=4, line=2.5, col="blue")
  } 
}

# We plot u (measured incoming rain) and y (measured water level) for each basin
old_mar <- par("mar")
par(mfrow=c(2,2), oma=c(0,2.5,0,2), mar=c(2.1, 2, 2, 2))
# We plot two lines in each plot, rainfall scaled by 10 visually
double_axis_plot(b1, "Event 1", "left", x_max, y_max)
double_axis_plot(b2, "Event 2", "right", x_max, y_max)
double_axis_plot(b3, "Event 3", "left", x_max, y_max)
double_axis_plot(b4, "Event 4", "right", x_max, y_max)

