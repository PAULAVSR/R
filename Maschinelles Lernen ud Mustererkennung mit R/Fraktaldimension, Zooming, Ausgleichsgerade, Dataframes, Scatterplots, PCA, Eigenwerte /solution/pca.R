## Team: Redrum
## Task: 2

graphics.off()
## 2a
data(iris)
dataset <- load("../pca.rda")

## 2b
plot_lfd <- function(x, to_plot = 1:(ncol(x) - 1), ...) {
  palette(rainbow(nlevels(x[, ncol(x)])))
	plot(x[,to_plot],col ="black" ,bg=x[, ncol(x)],pch = 21, ...)
}

## 2c
for(j in 1:4) {
plot_lfd(iris[c(j:4,5)],main = paste("Attribute ",j," bis 4"))
}

## 2d
PCA <- function(x, n=(length(x))) {
  stopifnot(is.data.frame(x) == T)
	x <- x[-length(x)]
	mean <- colSums(x)/nrow(x)
	cov.mat <- cov(x)
	eigen.values <- eigen(cov.mat)
	eigen.values
	 value <- list(mean=mean, eigenval=eigen.values$val[1:n], eigenvec=matrix(eigen.values$vec,length(x)))
	class(value) <- "PCA"
	value
}

## 2e
predict.PCA <- function(o, newdata) {
  stopifnot(class(o) == "PCA")
	klassenfaktor <- newdata[length(newdata)]
	newdata <- newdata[-length(newdata)]
	newdata <- t(t(newdata) - o$mean)
	eigenvec <- o$eigenvec
	if(ncol(newdata) < ncol(eigenvec)) {
	 	eigenvec <- eigenvec[1:length(newdata)]
	 }
	new_data <- apply(eigenvec,2,function (eigenvec) t(eigenvec) %*% t(newdata))
	eigenval <- as.matrix(1/sqrt(o$eigenval))
	new_data <- apply(as.matrix(c(1:length(eigenval))),1,function(i,eigenval,new_data) eigenval[i,] %*% new_data[,i],eigenval,new_data)
	data.frame(cbind(new_data, klassenfaktor))
}

## 2f
for(j in 4:1) {
	plot_lfd (predict (PCA (iris,j),iris))
}

## 2g
layout(matrix(1:8,2,2))

ldata <- list(iris.original = iris,iris.part = iris[iris$Species == 'setosa',])
tdata <- ldata

for (i in 1:length(ldata))
  for (j in 1:length(tdata))
    plot_lfd(predict(PCA(ldata[[i]],2),tdata[[j]]),main = paste("tdata=",names(tdata[j])," : ldata=",names(ldata[i])),xlab = "PC.1", ylab = "PC.2")
  

## 2h
layout(matrix(1:4,2,2,byrow = T))
dataset <- load("../pca.rda")

for (data in dataset) {
  name <- data
  data <- get(data)
  stopifnot(is.data.frame(data) == T)
  last.two <- (ncol(data) - 2):(ncol(data))
  data.length <- length(data)
  plot_lfd(data[c(1,2,data.length)], main = "Ersten 2 Merkmale",sub = name)
  plot_lfd(data[last.two], main = "Letzten 2 Merkmale",sub = name)
  plot_lfd(predict(PCA(data,2),data), main = "Ersten 2 HK",sub = name)
  plot_lfd(predict(PCA(data,last.two[2]),data),last.two[1:2], main = "Letzten 2 HK",sub = name)
}

dev.off()