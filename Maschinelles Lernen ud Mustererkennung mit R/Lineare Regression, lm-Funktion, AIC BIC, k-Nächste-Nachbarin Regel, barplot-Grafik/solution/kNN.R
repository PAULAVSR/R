## Team: Redrum
## Task

graphics.off()
set.seed (4711)	
if (require(class) == F) install.packages("class")
require(class)

## 2a
kNN <- function(x,neighbors=1) {
  stopifnot(is.data.frame(x))
  stopifnot(is.factor(x[,length(x)]))
  
  classifier <- list(train=x[,-length(x)],cl=x[,length(x)], neighbors=neighbors)
  classifier <- structure(.Data = classifier,class="kNN")
  classifier
}

## 2b/c
is.kNN <- function(o) inherits(o,"kNN")

predict.kNN <- function(o,newdata) {
  stopifnot(is.kNN(o))
  
  if (is.null(newdata)) 
      knn.cv(o$train, o$cl, o$neighbors)
  else 
  {
    stopifnot(is.matrix(newdata) || is.data.frame(newdata))
    stopifnot(is.factor(newdata[,length(newdata)])==F)
    
    knn(o$train,newdata,o$cl,o$neighbors)
  }
}

## 2d
heldout <- function(x, newdata=x, method, ...) {
  classified <- predict(method(x, ...),  newdata[-length(newdata)])
  #reclassification.rate <- 
    if(is.null(newdata)) {
      reclassification.rate <- mean(classified != x[,ncol(x)])
      attr(reclassification.rate,"confused") <- table(classified == x[,ncol(x)],x[,ncol(x)])
    } else {
      reclassification.rate <- mean(classified != newdata[,ncol(newdata)])
      attr(reclassification.rate,"confused") <- table(classified == newdata[,ncol(newdata)],newdata[,ncol(newdata)])
    }
  reclassification.rate
}

## 2e
run.1st <- function(x,y,choice=1+2*0:9) {
  to_calc <- function(j) {
    rbind(heldout(x,y,kNN,neighbors=j),heldout(y,x,kNN,neighbors=j),heldout(rbind(x,y),NULL,kNN,neighbors=j))
  }
  pe.1 <- sapply(choice,to_calc)  
  colnames(pe.1) <- choice
  rownames(pe.1) <- c("train/test","test/train","LO")
  pe.1
}

## 2f 
run.2nd <- function(x,y,choice=2^(0:13)) {
  pe.2 <- sapply(choice,function(j) heldout(x[sample(nrow(x),j,T),],y,kNN))
  names(pe.2) <- choice
  pe.2
}

## 2g
run.3rd <- function(x,choice=2:ncol(x)-1) {
  pe.3 <- sapply(choice,function(j) heldout(x[,-j],NULL,kNN))
  names(pe.3) <- choice
  pe.3
}

## 2h
diabetes <- load("../diabetes.rda")
letter <- load("../letter.rda")
germany <- load("../germany.rda")

pe.1 <- run.1st(diabetes.lern,diabetes.test)
pe.2 <- run.2nd(letter.lern,letter.test)
pe.3 <- run.3rd(rbind(germany.lern,germany.test))

#save(pe.1,pe.2,pe.3,file="kNN.rda")

## 2i
barplot(pe.1*100,beside = T,main = "diabetes",ylab = "error rate [%]",xlab = "number of nearest neighbour")
legend("bottomright",legend = rownames(pe.1),fill = gray.colors(3),cex = .65)
barplot(t(pe.1)*100,beside = T,main = "diabetes",ylab = "error rate [%]",xlab = "number of nearest neighbour")
barplot(pe.2*100,main = "letter",ylab = "error rate [%]",xlab = "training set size",col = "cyan")
barplot(pe.3*100,main = "germany",ylab = "error rate [%]",xlab = "knock-out feature",col="yellow")

dev.off()
