# Team: Redrum
# Task 1

graphics.off()

moment <- function(x,plot=TRUE,...) {
  
  calc.center <- function(x,p=0,q=0,x.s=0,y.s=0) {
    sum(x*(col(x)-x.s)^p * (row(x)-y.s)^q) 
  }
  
  x.s <- calc.center(x,1,0)/calc.center(x)
  y.s <- calc.center(x,0,1)/calc.center(x)
  alpha <- 1/2 * atan2(2*calc.center(x,1,1,x.s,y.s),(calc.center(x,2,0,x.s,y.s)-calc.center(x,0,2,x.s,y.s)))
  
  if (plot) {
    x.s.p <- x.s-0.5
    y.s.p <- nrow(x) - (y.s-0.5)
    alpha <- -alpha
    if (alpha > 0) factor <- -90 else factor <- 90
    plot.array(x,...)
    mtext(round(factor +(alpha*180/pi),1),side = 1)
    abline(h = y.s.p, v = x.s.p, col = "navyblue")
    abline(a = y.s.p-x.s.p*tan(alpha),b=tan(alpha),col = "red")
  }
  
  c(x.center=x.s,y.center=y.s,alpha.grad=alpha*180/pi)
}

data <- load("../moment.rda")

moment.plot <- function(pic,...) {
  layout(matrix(1:4,2,2,byrow = T))
  moment(pic,...)
  moment(1-pic,...)
  moment(binarise(pic),...)
  moment(binarise(1-pic),...)
}

for (element in data) {
  if (is.array(get(element))) moment.plot(get(element),main= element)
}

#dev.off()