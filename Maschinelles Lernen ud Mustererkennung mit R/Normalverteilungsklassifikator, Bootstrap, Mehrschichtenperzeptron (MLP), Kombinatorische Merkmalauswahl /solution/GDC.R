## Team: Redrum
## Task: 1

graphics.off()
require(MASS)

## 1a/d
GDC <- function(x,linear=FALSE) {
  structure(list(q=ifelse(linear,lda, qda)(x[,-ncol(x)],x[,ncol(x)])),class="GDC")
}

## 1b
is.GDC <- function(o) inherits(o,"GDC")
predict.GDC <- function(o,newdata) {
  stopifnot(!is.factor(newdata[[ncol(newdata)]]))
  stopifnot(is.GDC(o))
  
  predict(o$q,newdata)$class
}

## 1c
heldout <- function(x,newdata=x,method,...) {
  GDC <- method(x,...)
  classified <- predict(GDC, newdata[-ncol(newdata)])
  reclassification.rate <- mean(classified != newdata[,ncol(newdata)])*100
  attr(reclassification.rate,"confused") <- table(classified == newdata[,ncol(newdata)],newdata[,ncol(newdata)])
  reclassification.rate
}

heldout.KxK <- function(data,...,plot=TRUE) {
  stopifnot(is.list(data))
  
  all <- do.call(rbind,data)
  comb.names <- c(paste(rev(names(data))),"all")
  comb.mat <- as.matrix(expand.grid(rev(comb.names),rev(comb.names)))
  rownames(comb.mat) <- rep(rev(comb.names),each=length(comb.names))
  heldout.caller <- function(x,...,name) {
    error <- heldout(x,...)
    names(error) <- name
    error
  }
  error.rates <- sapply(1:nrow(comb.mat),function(i) heldout.caller(get(comb.mat[i,1]),get(comb.mat[i,2]),GDC,name=rownames(comb.mat)[i],...),simplify = T)
  sol.mat <- matrix(error.rates,nrow=length(comb.names),ncol=length(comb.names),dimnames = list(rev(comb.names),rev(comb.names)))
  
  if (plot) {
    ifelse(nchar(comb.names[1])>5,dataset.name <- paste(substring(comb.names[1],first=1,last=nchar(comb.names[1])-5),"mit"),dataset.name <-"")
    barplot(sol.mat,beside = T,col = rainbow(length(comb.names)),xlab = "Abrufdatensatz",ylab="Fehlerrate[%]",main=paste(dataset.name,ifelse(length(substitute(...))==1 && substitute(...) == TRUE,"'lda'","'qda'")),ylim = c(0,round(max(sol.mat))))
  }
  
  sol.mat
  
}

satimage <- mget(load("../satimage.rda"))

satimage.qda <- heldout.KxK(satimage)
satimage.lda <- heldout.KxK(satimage,linear=TRUE)

## only for assignment
# write.table(satimage.qda,"satimage.qda_errors.txt")
# write.table(satimage.lda,"satimage.lda_errors.txt")

## 1e
bootstrap <- function(x,k,...) {
  heldout.caller <- function(x,...) {
    train <- sample(1:nrow(x),nrow(x),replace = T)
    test <- (1:nrow(x))[-train]
    heldout(x[train,],x[test,],GDC,...)
  }
  sapply(1:k,function(i) heldout.caller(x,...))
}

## 1f
plotter <- function(data,k=32) {
  layout(matrix(c(1,2,3,3),2,2,byrow=T))
  heldout.KxK(data,linear=FALSE)
  heldout.KxK(data,linear=TRUE)
  all.qda <- bootstrap(do.call(rbind,data),k,linear=FALSE)
  all.lda <- bootstrap(do.call(rbind,data),k,linear=TRUE)
  lern.lda <- bootstrap(data[[1]],k,linear=TRUE)
  test.lda <- bootstrap(data[[2]],k,linear=TRUE)
  boxplot(cbind(all.qda,all.lda,lern.lda,test.lda),col=c(2,2,3,4),main=paste("bootstrap errors:",k,"draws"),xlab="Abrufdatensatz",ylab="Fehlerrate[%]")
}

heart <- mget(load("../heart.rda"))
plotter(heart)
plotter(satimage)
vehicle <- mget(load("../vehicle.rda"))
plotter(vehicle)
diabetes <- mget(load("../diabetes.rda"))
plotter(diabetes)

dev.off()