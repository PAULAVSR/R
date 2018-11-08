# Team: Redrum
# Task 1

graphics.off()


## 1a
DFTmat <- function(N) {
  k <- 0:(N-1)
  X <- -2*pi*1i*k/N
  W <- 1/N * exp(X %o% k)
  return(W)
}

## 1b
W <- DFTmat(N = 5)
Conj(W) %*% W

## 1c
powspec <- function(x) {
  n <- length(x)
  W <- DFTmat(n)
  spec <-  W %*% x
  spec.norm <- 2 * abs(spec) / length(spec)
  return (20*log10(spec.norm[1:(n/2)]))
}

## 1d
plot.powspec <- function(x,ATF=10000,...) {
  y <- seq(0,(length(x)-1)*1000/ATF,by = 1000/ATF)
  plot(y,x,type = "l",col = "blue",xlab = "msec",ylab = "Amplitude",...)
  # db-spectre
  fft.x <- powspec(x)
  f.Nyquist <- 1/2*ATF
  f <- f.Nyquist * c(seq(length(y)/2))/(length(y)/2)
  #f <- f.Nyquist*c(seq(length(y)/2),rep(NA,length(y)/2))/(length(y)/2)
  plot(f,fft.x,type = "l",col = "red",xlab = "Hz",ylab = "dB")
  
}

#pdf("power.pdf")
## 1e
T <- 256
par(mfrow = c(3,2))
x1 <- rep(seq(0,1,len=16),T/16)
plot.powspec(x1 )
x2 <- rep(rep(0:1,each=8),T/16)
plot.powspec(x2 )
x3 <- rnorm(T,0,3)
plot.powspec(x3 )
x4 <- cos((1:T)*2)+cos((1:T)*1.2)+cos((1:T)*0.6)
plot.powspec(x4 )

## 1f
# hn = fn - fn-1
x.diff <- x4[1:(length(x4)-1)]
hn <- c(x4[1],x4[2:length(x4)] -x.diff)
plot.powspec(hn )

# hn = fn + hn-1
hn.discrete <- c(x4[1])
for(i in 2:length(x4)) {
  hn.discrete <- c(hn.discrete, x4[i] + hn.discrete[(i-1)])
}
plot.powspec(hn.discrete )

dev.off()