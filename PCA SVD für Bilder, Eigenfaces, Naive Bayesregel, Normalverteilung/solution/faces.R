## Team: Redum
## Task 1

graphics.off()

faces.data <- load("../faces.rda")

## 1a
data.arrays <- sapply(faces.data,function(x,log) if(is.array(get(x))) T else F)
length.data <- length(data.arrays[data.arrays == T])
pic.names <- faces.data[1:length.data]

faces.list <- lapply(1:length.data,function(x,faces.data) if(is.array(get(faces.data[x]))) get(faces.data[x]),faces.data)
names(faces.list) <- pic.names

layout(matrix(1:9,3,3))
invisible(mapply(plot.array,faces.list,main = names(faces.list)))

## 1b
impack <- function(x) {
  stopifnot(is.list(x))
  new <- matrix(1:length(x),nrow =length(x),byrow = T)
  t(apply(new,1,function(rows,x) as.vector(x[[rows]]),x))
}

## 1c
plot.flatpic <- function(x,dim,...) {
  stopifnot(is.vector(x))
  new.mat <- matrix(x,dim[1],dim[2])
  range.mat <- range(new.mat)
  norm.mat <- apply(new.mat,1:2,function(x,range) (x-range[1])/(range[2]-range[1]),range.mat)
  plot.array(as.array(norm.mat),...)
}

## 1d
layout(matrix(1:9,3,3))
faces.packed <- impack(faces.list)
invisible(lapply(1:27,function(i) plot.flatpic(faces.packed[i,],c(160,120),main=pic.names[i])))

## 1e
eigenface <- function(X,n=nrow(X)) {
  X <- sweep(X, 2, colMeans(X), "-")
  X.singulaer <- svd(X, n)
  names(X.singulaer)
  singval <- X.singulaer$d[1:n]
  PC <- X.singulaer$u %*% diag(X.singulaer$d[1:n]) / sqrt(n)
  PA <- X.singulaer$u %*% X / sqrt(n)
  list(PC=PC, PA=PA, singval=singval)
} 

## 1f
layout(matrix(1))
faces.matrix <- impack(faces.list)
faces.PCA <- eigenface(faces.matrix)
barplot(faces.PCA$singval,main = "Singulärwerte")


## 1g
plot.PC <- function(PC,toPlot,names="",cols,...) {
  stopifnot(is.vector(toPlot))
  stopifnot(length(toPlot) == 2)
  stopifnot(dim(PC)[2] >= toPlot[2])
  stopifnot(is.matrix(PC))
  
  plot(PC[,toPlot],col=cols[rep(c(1,10,19),each = 9)],type = "p", pch = 16,xlab = paste("PC #",toPlot[1],sep=""),ylab = paste("PC #",toPlot[2],sep=""), ...)
  text(PC[,toPlot],names,adj = c(0,1))
}

layout(matrix(1:9,3,3,byrow = T))
pc2 <- 2*1:(12)
pc1 <- pc2-1
cols <- palette(rainbow(dim(faces.PCA$PC)[1]))
invisible(apply(cbind(pc1,pc2),1,function(i) plot.PC(faces.PCA$PC,i,pic.names,cols)))

## 1h
invisible(apply(as.matrix(1:length.data),1,function(i) plot.flatpic(faces.PCA$PA[i,],c(160,120), main = pic.names[i])))

## 1i
invisible(apply(as.matrix(1:length.data),1,function(i) plot.flatpic(as.vector(faces.PCA$PC[,i]%*%faces.PCA$PA),c(160,120),main = faces.data[i])))
dev.off()