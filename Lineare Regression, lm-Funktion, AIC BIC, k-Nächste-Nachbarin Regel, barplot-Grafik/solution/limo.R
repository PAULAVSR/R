## Team: Redrum
## Task: 1

graphics.off()

## 1a
dataset.name <- load("../limo.rda")
xydata <- get(dataset.name)

## 1b
polyfun <- function(x,a) {
  sapply(x,function(x) sum(a*x ^ (0:(length(a)-1)))) 
}

## 1c
polyfit <- function(xy,n) {
  stopifnot(n>=0)
  stopifnot(length(xy)>= 2)
  
  left.part <- paste(names(xy)[2],"~")
  right.part <- "1"
  if (n > 0)
    right.part <- paste("I(",names(xy)[1],"^",seq(1,n,length.out = n),")",collapse = "+")
  
  lm(paste(left.part,right.part),xy)
  
}

## 1d
polyfits <- function(xy,order,plot=FALSE) {
  if (length(order) == 1) order <- 1:order
  if (min(order) == 0) offset <- -1 else offset <- 0
  
  list.polyfits <- lapply(order,function(n) polyfit(xy,n))
  
  if (plot) {
    criteria <- sapply(list.polyfits,function(lm) c(AIC(lm),BIC(lm)))
    plot(order,criteria[1,],col = "navyblue",type = "b",pch=16,ylim=range(criteria),xlab="Polynomgrad",ylab = "AIC/BIC Wert",main = paste("AIC&BIC für Ausgleichspolynome der Attribute ",names(xy)[2]," ~ ",names(xy)[1],sep = ""),cex = .75,cex.main=.75)
    lines(order,criteria[2,],col = "tomato",type = "b",pch = 16)
    axis(1,min(order):max(order))
    legend("top",legend = c("AIC","BIC"),pch = 16,col = c("navyblue","tomato"),cex = .75)
    min.AIC = which.min(criteria[1,])+offset
    min.BIC = which.min(criteria[2,])+offset
    if (min.AIC == min.BIC)
      mtext(paste("!!! Gewinner ist Polynomgrad",min.AIC,"!!!"),col = "red",cex=.75)
    else
      mtext(paste("AIC:",min.AIC,"BIC:",min.BIC),cex=.75)
  }
  list.polyfits
}

polyfits(xydata,0:11,plot=TRUE)

## 1e
layout(matrix(1:4,2,2,byrow = T))

plot.polyfits <- function(data,showCriterion=TRUE) {
    
  models <- polyfits(data,0:11,showCriterion)
    
  plot.model <- function(model) {
    ifelse(is.na(model$coeff),model$coeff[which(is.na(model$coeff))] <- 0,model)
    plot(data,pch = 20,cex = .5,col = "red",main = paste("Polynomgrad",length(model$coeff)-1))
    curve(polyfun(x,model$coeff),xlim=c(min(data[,1]),max(data[,1])),n=length(data[,1])*2,add=T,col = "blue")
    mtext(paste(round(signif(model$coeff,2),1),collapse = ",",sep = ""),side=3,cex =.75)
    mtext(paste("BIC=",round(BIC(model),2)),side = 4,cex = .75)
  }

  invisible(mapply(plot.model,models))
}
plot.polyfits(xydata,F)

## 1f
data("cars")
plot.polyfits(cars,T)

## 1g
data("LifeCycleSavings")
attributes <- c("pop15","pop75","dpi")
invisible(apply(cbind(combn(attributes,2),combn(rev(attributes),2)),2,function(i) plot.polyfits(LifeCycleSavings[i])))

#dev.off()