## Team: Redrum
## Task: 2

graphics.off()

## 2a
if (require(nnet) == F) install.packages("nnet")
require(nnet)

## 2e
p <- function(D, H, K) {
  (D + 1) * H + (H + 1) * ifelse(K==2, 1, K)
}

## 2b
SHLP <- function(x, hidden, norm=FALSE, Wts=sin(10*1:p(ncol(x)-1, hidden, nlevels(x[,ncol(x)])))) {
  
  ## 2g
  if(norm)
    x[-ncol(x)] <- (scaling <- scale(x[,-ncol(x)])) * (sd <- 1/3) 
  else
    scaling<-sd<-NULL

  net <- nnet(as.formula(paste(names(x)[ncol(x)], "~", ".")), data=x,size=hidden,Wts=Wts,trace = FALSE)
  structure(list(net=net, classes=levels(x[,ncol(x)]),centers=attr(scaling,"scaled:center"),scaling=attr(scaling,"scaled:scale")/sd),class="SHLP")
}

## 2c
is.SHLP <- function(x) inherits(x,"SHLP")
predict.SHLP <- function(o, newdata) {
  stopifnot(is.SHLP(o))
  stopifnot(!is.factor(newdata[[ncol(newdata)]]))
  
  factor(predict(o$net,if (is.null(o$centers)) newdata else scale(newdata,o$centers,o$scaling), type="class"), o$classes)
}

## 2d
heldout <- function(x, newdata=x, method, ...) {
  net <- method(x, ...)
  classified <- predict(net, newdata[,-ncol(newdata)])
  reclassification.rate <- mean(classified != newdata[,ncol(newdata)])*100
  attr(reclassification.rate,"confused") <- table(classified == newdata[,ncol(newdata)],newdata[,ncol(newdata)])
  reclassification.rate
}

## 2f
australia <- load("../australia.rda")
diabetes <- load("../diabetes.rda")
segment <- load("../segment.rda")
vehicle <- load("../vehicle.rda")

error.table <- function(lern, test) {
  H <- c(1, 3, 6, 10, 15, 21, 28)
  lernlern <- sapply(H,function(h) heldout(lern, lern, SHLP, hidden=h, norm=FALSE))
  lerntest <- sapply(H,function(h) heldout(lern, test, SHLP, hidden=h, norm=FALSE))
  error <- rbind(lernlern, lerntest)
  rownames(error) <- c("Lernfehler", "Testfehler")
  colnames(error) <- H
  error
}

error.table(australia.lern, australia.test)
error.table(diabetes.lern, diabetes.test)
error.table(segment.lern, segment.test)
error.table(vehicle.lern, vehicle.test)

# error <- error.table(australia.lern, australia.test)
# error <- rbind(error, error.table(diabetes.lern, diabetes.test))
# error <- rbind(error,error.table(segment.lern, segment.test))
# error <- rbind(error,error.table(vehicle.lern, vehicle.test))
# rownames(error) <- paste(rownames(error), rep(list("australia","diabetes", "segment","vehicle"), each=2))
# write.table(round(error,3), "2f.txt", sep="\t")

## 2h
plot.errorrate <- function(lern, test, main="") {
  H <- c(1, 3, 6, 10, 15, 21, 28)
  lern.test <- sapply(H,function(h) heldout(lern, test, SHLP, hidden=h, norm=FALSE))
  plot(H, lern.test, main=main, ylab="Errorrate [%]", xlab="Numbers of Hiddenneuron", col="red",type = 'b',ylim=c(0, max(lern.test)), log="x")
  lern.lern <- sapply(H, FUN=function(h) heldout(lern, lern, SHLP, hidden=h, norm=FALSE))
  lines(H, lern.lern, lty=2, col="red",type = 'b')
  lern.test.norm <- sapply(H, FUN=function(h) heldout(lern, test, SHLP, hidden=h, norm=TRUE))
  lines(H, lern.test.norm, col="blue",type = 'b')
  lern.lern.norm <- sapply(H, FUN=function(h) heldout(lern, lern, SHLP, hidden=h, norm=TRUE))
  lines(H, lern.lern.norm, lty=2, col="blue",type = 'b')
  legend("bottomleft", legend = c("Testfehler (+nrom)", "Lernfehler (+norm)", "Testfehler(-norm)", "Lernfehler(-norm)"), col=c("blue", "blue", "red", "red"), lty=c(1,2,1,2),bty='n',xpd = NA, cex = 0.65)
  table <- rbind(lern.test, lern.lern)
  table.norm <- rbind(lern.test.norm, lern.lern.norm)
  tables <- rbind(table, table.norm)
  rownames(tables) <- list("Lern/Test", "Lern/Lern", "Lern/Test Norm", "Lern/Lern Norm")
  colnames(tables) <- H
  tables 
}
layout(matrix(1:4, 2, 2, byrow=TRUE))
plot.errorrate(australia.lern, australia.test, "australia")
plot.errorrate(diabetes.lern, diabetes.test, "diabetes")
plot.errorrate(segment.lern, segment.test, "segment")
plot.errorrate(vehicle.lern, vehicle.test, "vehicle")

# error <- plot.errorrate(australia.lern, australia.test)
# error <- rbind(error, diabetes=plot.errorrate(diabetes.lern, diabetes.test))
# error <- rbind(error,plot.errorrate(segment.lern, segment.test))
# error <- rbind(error,plot.errorrate(vehicle.lern, vehicle.test))
# rownames(error) <- paste(rownames(error), rep(list("australia","diabetes", "segment","vehicle"), each=4))
# write.table(round(error,3), "2h.txt", sep="\t")

dev.off()