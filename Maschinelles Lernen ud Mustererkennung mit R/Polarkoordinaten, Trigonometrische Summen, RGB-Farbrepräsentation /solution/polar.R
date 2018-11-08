######## WMM-Assignment 1 
### 
### A1
## a)
fcos <- function(x) { a * cos(x) }
fsin <- function(x) { b * sin(x) }
fsum <- function(x) { fcos(x) + fsin(x) }
fpol <- function(x) { A*cos(x+phy) }

## b)
a <- -3
b <- -4
#x <- c(seq(0,8,0.01))
phy.seq <- seq(0,2*pi,0.01)

A <- sqrt(a^2+b^2)
phy <- atan((-b/a))-pi

## c)
graphics.off()
pdf(file = "polar.pdf")
plot(fpol,0,8, type = "p", col = "orange", ylab = "fcos, fsin, fsum, fpol",add = F, main = "Task1")
plot(fcos,0,8, type = "l", col= "tomato", ylab = "",add = T)
plot(fsin,0,8, type = "l", col = "navyblue", ylab = "",add = T)
plot(fsum,0,8, type = "l", col = "darkgreen", ylab = "",add = T)
legend(0,5, legend=c("fcos", "fsin", "fsum", "fpol"),
       col=c("tomato", "navyblue", "darkgreen", "orange"), pch=c("-","-","-","°"), cex=1)
dev.off()
