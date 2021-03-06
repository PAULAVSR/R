# Team: Redrum
# Task 2

graphics.off()
pdf("zip.pdf")
## 1a
data <-get(load("../zip.rda"))

## 1b
bits <- function(x,compress = T) {
  if (compress) x <- memCompress(x,"gzip")
  sum(nchar(x,"bytes"))*8
}

## 1c
compress.factor <- sort(sapply(text,bits)/sapply(text,bits,F),F)
dotchart(compress.factor,main = "Kompressionsfaktor", color = "navyblue",cex = .7,pch = 16)

## 1d
entropy <- function(xp,xq) {
  (bits(c(xq,xp)) - bits(xq))/bits(xp,F)
}

## 1e
divergence <- function(xp,xq) {
  entropy(xp,xq) - entropy(xp,xp)
}

distance <- function(X) {
  S <- sapply(X,function (a) sapply(X,function(b) divergence(a,b)))
  as.dist(S+t(S))
}

## 1f
if (require("cluster") == F) install.packages("cluster")
require("cluster")
dist <- distance(text)
plot(agnes(dist),which.plots = 2, main = "Agnes")
plot(diana(dist),which.plots = 2, main = "Diana")


dev.off()