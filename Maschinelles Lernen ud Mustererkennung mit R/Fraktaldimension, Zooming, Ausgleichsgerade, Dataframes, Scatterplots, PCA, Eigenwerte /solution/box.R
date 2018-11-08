## Team: Redrum
## Task: 1

graphics.off()
## 1a
zoom <- function(x, delta = 1) {
	new_row <- (nrow(x) %/% delta) * delta
	new_col <- (ncol(x) %/% delta) * delta
	x <- x[1:new_row,1:new_col]
	v <- as.vector(mapply(1:(nrow(x)/delta),FUN=rep, MoreArgs=list(times=delta)))
	h <- as.vector(mapply(1:(ncol(x)/delta)-1, FUN=rep, MoreArgs=list(times=delta))) * (nrow(x)/delta)
	m <- outer(v, h, FUN='+')
	n <- tapply(x,m, FUN='all')
	matrix(n, nrow=(new_row/delta))
}

## 1b
data <- get(load("../box.rda"))
layout(matrix(1:9,3,3, byrow=TRUE))

for(i in 1:9) {
	plot.array(zoom(Couple, i), main=paste("Couple",i))
}

## 1c
slope <- function(x,y) {
	lm(y ~ x)$coeff[2]
}

x <- -5:5
y <- pi*x+sin(x)
invisible(slope(x, y))

## 1d/e
layout(matrix(1,1,1))
boxdim <- function(x, d, plot=FALSE) {
	delta <- 1:d
	vol <- sapply(delta, function(d, x) sum(!zoom(x,d)),x)
	dim <- -slope(log(delta), log(vol))
	if(plot) {
		plot(delta, vol, log="xy", col="blue", type="p", pch=16)
		mtext(dim,  3)
	}
	dim
}

invisible(boxdim(tonga, 8, T))

boxdim.plot <- function(pic,...) {
  plot.array(pic,...)
  boxdim(pic,8,T)
}
layout(matrix(1:4,2,2,byrow = T))

pics <- load("../box.rda")

for (element in pics) {
  if (is.array(get(element))) boxdim.plot(get(element),main = element)
}

dev.off()