# Team: Redrum
# Task 1
graphics.off()
# initializations
load("../filter2D.rda")

## 1a
translate <- function(x,dr,dc) {
  dim <- dim(x)
  row <- ((1:dim[1] - (dr+1)) %% dim[1]) +1
  col <- ((1:dim[2] - (dc+1)) %% dim[2]) +1
  x[row,col]
}

r.rot <- c(12,12,-12,-12)
c.rot <- c(25,-25,25,-25)

layout(matrix(1:4,2,2))

for (off in 1:length(r.rot)) {
    plot.array(translate(JFK,r.rot[off],c.rot[off]), main = paste("JFK Offset (",r.rot[off],", ",c.rot[off],")"))
}

## 1b
dc.seq <- rep(c(-1,0,1),3)
dr.seq <- c((c(1,1,1) %o% c(-1,0,1)))

filter.mean.3x3 <- function(x) {
   
   sum.translate.mean <-0
   for (i in 1:length(dr.seq)) {
     sum.translate.mean <- sum.translate.mean + translate(x,dr.seq[i],dc.seq[i])
   }
    1/9 * sum.translate.mean
}

filter.laplace4 <- function(x, norm = TRUE) {
  
  factor.seq <- c(0,1,0,1,-4,1,0,1,0)
  sum.translate.lp <- 0
  for (i in 1:length(dr.seq)) {
      sum.translate.lp <- sum.translate.lp + factor.seq[i] * translate(x,dr.seq[i],dc.seq[i])
  }
  result <- sum.translate.lp
  if (norm == T) result/8 + 1/2 else result
}

filter.prewitt.v <- function(x,norm = TRUE) {
  dr.seq.v <- rev(dc.seq)
  dc.seq.v <- rev(dr.seq)
  factor.seq <- c(as.vector(c(1,1,1) %o% c(-1,0,1)))
  sum.translate.p.v <- 0
  for (i in 1:length(dr.seq)) {
    sum.translate.p.v <- sum.translate.p.v + (factor.seq[i] * translate(x,dr.seq.v[i],dc.seq.v[i]))
  }
  result <- 1/6 * sum.translate.p.v
  if (norm == T) result + 1/2 else result
}

filter.prewitt.h <- function(x,norm = TRUE) {
  dr.seq.h <- rev(dr.seq)
  dc.seq.h <- rev(dc.seq)
  factor.seq <- dr.seq
  sum.translate.p.h <- 0
  for (i in 1:length(dr.seq)) {
    sum.translate.p.h <- sum.translate.p.h + factor.seq[i] * translate(x,dr.seq.h[i],dc.seq.h[i])
  }
  result <- 1/6 * sum.translate.p.h  
 if (norm == T) result + 1/2 else result
}

## 1c
filter.grad.mag <- function(x,norm = TRUE) {
	result <- sqrt((filter.prewitt.v(x, norm=F))^2 + (filter.prewitt.h(x, norm=F))^2)
	if (norm == TRUE) result  else result
}

filter.grad.angle <- function(x,norm = TRUE) {
  result <- atan2(filter.prewitt.v(x, norm=F),filter.prewitt.h(x, norm=F))
  if (norm == TRUE) result / (2 * pi) + 0.5 else result
}

## 1d
layout(matrix(1:6,3,2))
filter2D <- list(filter.mean.3x3(JFK),filter.laplace4(JFK),filter.prewitt.v(JFK),filter.prewitt.h(JFK),filter.grad.mag(JFK),filter.grad.angle(JFK))
filter2D.names <- c("filter.mean.3x3","filter.laplace4","filter.prewitt.v","filter.prewitt.h","filter.grad.mag","filter.grad.angle")
invisible(mapply(plot.array,filter2D,main = filter2D.names, sub = "pic: JFK"))

## 1e
layout(matrix(1:2,2,1))
grad.angle.JFK <- filter.grad.angle(JFK)*180/pi
breaks <- seq(-180,180,len = 45+1)
hist(grad.angle.JFK,breaks, main = "Kantenrichtungsbild JFK", xlab = "Phasenwinkel in Grad")
grad.angle.JFK.quantil80 <- quantile(grad.angle.JFK, 0.8)
hist(grad.angle.JFK[grad.angle.JFK > grad.angle.JFK.quantil80],breaks, main = "Richtungsbildpunkte > 80% Quantil des Gradientenbetrags",xlab = "Phasenwinkel in Grad")

## 1f

plot.gradientfilter <- function(picCollection,picsToDraw,mulaw.function,mulaw = F, ...) {
  
  for (i in 1:length(picCollection)) {
    if (picCollection[i] %in% picsToDraw) {
      pic <- eval(parse(text = picCollection[i]))
      
      stopifnot(is.array(pic),max(pic)<=1,min(pic)>=0)
      if (mulaw == F) plot.array(filter.grad.mag(pic),main = paste("filter.grad.mag(",picCollection[i],")"),sub = paste("Pic: ",picCollection[i]),...)
      else plot.array(mulaw.function(filter.grad.mag(pic)),main = paste("mulaw(filter.grad.mag(",picCollection[i],"))"),sub = paste("Pic: ",picCollection[i]),...)
    }
  }
}

data.objects <- c("algae", "cashmere", "muscle", "tonga", "ludwig", "xray")
data <- load("../filter2D.rda")
layout(matrix(1:6,2,3))
plot.gradientfilter(data,data.objects)

## 1g
mulaw <- function(x,mu=100) {
	x.min <- min(x)
	x.max <- max(x)
	x <- (x-x.min) / (x.max-x.min)
	log(1 + mu * x)/log(1+ mu)
}
layout(matrix(1:6,2,3))
plot.gradientfilter(data,data.objects,mulaw,T)


dev.off()