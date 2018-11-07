## Team: Redrum
## Task: 2

graphics.off()

## 2a
naivegauss <- function(x) {
  stopifnot(is.data.frame(x))
  stopifnot(is.factor(x[,ncol(x)]))
  
  x.splitted <- split(x[-(ncol(x))], x[,ncol(x)])
  
  x.prob <- tapply(x[,ncol(x)], x[,ncol(x)], length) / nrow(x)
  x.mean <- sapply(x.splitted,function(x) sapply(x,mean))
  x.variance <- sapply(x.splitted,function(x) sapply(x,var))
  
  classifier <- list(classes = levels(x[,ncol(x)]),class.prob = as.matrix(x.prob), class.means = as.matrix(x.mean), class.variance = as.matrix(x.variance))
  classifier <- structure(.Data = classifier,class = "naivegauss")
  classifier
}

## 2b
predict.naivegauss <- function(o,newdata) {
  stopifnot(class(o) == "naivegauss")
  stopifnot(is.factor(newdata[,ncol(newdata)]))
  
  newdata <- newdata[-ncol(newdata)]
  newdata <- as.matrix(newdata)

  u.kappa <- function(x,k,o) {
    num.feature <- length(o$class.variance)/length(o$classes)
    row.index <- 1:num.feature
    if (num.feature == 1) {row.index <- k; k <- num.feature}
    
    -2*log(o$class.prob[k]) + length(x)*log(2*pi) + sum(log(o$class.variance[row.index,k])) + sum(((x-o$class.means[row.index,k])/sqrt((o$class.variance[row.index,k])))^2)
  }
  u.mat <- matrix(nrow = nrow(newdata), ncol = length(o$classes))
  u.mat <- sapply(1:length(o$classes),function (i) apply(newdata,1,u.kappa,k=i,o=o))
  factor (o$classes[apply(u.mat,1, FUN = which.min)], levels=o$classes)
}

## 2c
heldout <- function(x,newdata=x,method,...) {
  stopifnot(is.factor(x[,ncol(x)]))
  stopifnot(is.factor(newdata[,ncol(newdata)]))
  
  classified <- predict.naivegauss(method(x,...),newdata)
  reclassification.rate <- sum(classified != newdata[,ncol(newdata)])/nrow(newdata)
  groups <- split(classified == newdata[,ncol(newdata)],newdata[,ncol(newdata)])
  attr(reclassification.rate,"confused") <- as.matrix(table(classified == newdata[,ncol(newdata)],newdata[,ncol(newdata)]))
  reclassification.rate
}

## 2d
heldout(iris,iris,method=naivegauss)
#write(paste("Fehlerraten\n2d:",error),"error-rate.txt")
#heldout(iris[4:5],method=naivegauss)

## 2e
data <- load("../vehicle.rda")
error.rates <- apply(combos <- expand.grid(data,data),1,function (c) heldout(get(c[1]),get(c[2]),naivegauss))
names(error.rates) <- paste(combos[,1],combos[,2],sep = ":")
error.rates
#dim(vehicle.lern)
#dim(vehicle.test)
#write("2e:","error-rate.txt",append = T)
#write.table(error.rates,"error-rate.txt",append = T,col.names = T)

#dev.off()