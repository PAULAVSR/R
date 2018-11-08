## Team: Redrum
## Task: 1

graphics.off()

## 1a/g
parzen <- function(x,sigma = NULL) {
  if (is.null(sigma)) {
    opt.L1o <- function(x,sd) sum(log(predict(parzen(x,sd))))
    sigma <- optimize(opt.L1o,x=x,lower = 0,7,tol = 0.001,maximum = T)[[1]]  
  }
  
  structure(.Data= list(support=x,sigma=sigma),class = "parzen")
}

## 1b
is.parzen <- function(o) inherits(o,"parzen")

predict.parzen <- function(o, newdata=NULL) {
  stopifnot(is.parzen(o))
  ## 1e
  if (is.null(newdata))
    sapply(1:length(o$support),function(i) mean(dnorm(o$support[i],o$support[-i],o$sigma)))
  else
    sapply(newdata,function(x) mean(dnorm(x,o$support,o$sigma)))
}

## 1c
plot.parzen <- function(o,xlim=range(o$support)+c(-5,5)*o$sigma,...) {
  stopifnot(is.parzen(o))
  curve(predict(o,x),xlim = xlim,n=1000,col = "blue",...)
  rug(o$support,-0.04,col = "red",lwd = 3)
  mtext(side = 3,paste("s=",signif(o$sigma,3)),cex=.8)
  ## 1f
  points(o$support,predict(o),pch=15,cex=.3)
}

## 1d/1f
parzen.data <- load("../parzen.rda")
layout(matrix(1:6,3,2,T))
invisible(sapply(c(7:-4),function(x) plot(parzen(samples,0.7^x))))

## 1h
layout(1)
plot(parzen(samples))

#dev.off()