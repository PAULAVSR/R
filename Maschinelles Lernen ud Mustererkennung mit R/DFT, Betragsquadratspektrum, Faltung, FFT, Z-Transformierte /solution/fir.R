# Team: Redrum
# Task 2

graphics.off()


## 1a
sms.Gz <- function(g,n) {
  omega.freq <- seq(0,pi,length.out = n)
  z <- exp(1i*omega.freq)
  z.power <- 0:(-(length(g)-1))
  z.trans <- outer(z,z.power, "^")
  G <- g%*%t(z.trans)
  return(Mod(G)^2)
}

## 1b
sms.fft <- function(g, n){
  g<-c(g,rep(0,max((n-1)*2-length(g),0)))
  fourier <- Mod(fft(g))^2
  fourier[1:n]
}

## 1c
sms.conv <- function(g,n,f=seq(0,1,length.out = n)) {
  h.answer <- convolve(f,g,type = "o")
  H.ny <- sms.fft(h.answer,n)
  F.ny <- sms.fft(f,n)
  G.ny <- H.ny/F.ny
}

## 1d
smsplot <- function(g,n=24,...) {
  omega.freq <- seq(0,pi,length.out = n)
  plot(omega.freq,sms.conv(g,n),col = "red",pch = 23,...)
  points(omega.freq,sms.fft(g,n),pch = 20,col = "navyblue")
  lines(omega.freq,sms.Gz(g,n),col = "darkgreen")
}
#pdf("fir.pdf")
## 1e
### a
smsplot(c(1,-2,1),n = 24, main="Hochpass")
smsplot(c(1,2,1),n = 24, main="Tiefpass")
smsplot(c(1,0,-1),n = 24, main="Bandpass")
smsplot(c(1,0,1),n = 24, main="Kerbfilter")

### b
for (i in c(2,4,7,11)) {
  smsplot(rep(1/i,i), main = paste("Mittelwertfilter ",i))
}

### c
for (m in c(2,4,7,11)) {
  smsplot(dnorm(seq(-1,+1,length=2*m+1)),main = paste("gaußscher Koeff. für ",m))
}

### d
exp.filter <- function(lambda,n.max = 24) {
  (1-lambda)*lambda^(0:(n.max-1))
}

lambda <- c(1/2,3/4,7/8,15/16)
n.max <- c(24,30)
par(mfrow = c(4,2))

for (x in lambda) {
  for (nmax in n.max) {
    smsplot(exp.filter(x,nmax),nmax,main = paste("Exp.Filter nmax = ",nmax," für lambda ",x))
  }
}

## mit smsplot nmax+1 Aufruf 
#for (x in lambda) {
#  for (nmax in n.max) {
 #   smsplot(exp.filter(x,nmax),nmax+1,main = paste("Exp.Filter nmax = ",nmax," für lambda ",x))
  #}
#}


dev.off()
