# Team: Redrum
# Task 2
graphics.off()


## 1a
equalize <- function(x,RGB=F) {
  if (RGB == F && is.matrix(x) == T) {
    x[] <- rank(x)/length(x)
    return (as.array(x))
  }
  else {
    stopifnot(is.matrix(x) == F)
    x[] <- apply(x,3,equalize)
    return (x)
  }
  
}

## 1b
data <- load("../equalize.rda")
data.objects <- c("algae","couple","Donald","mri1","soil","turbinate","GUESS")
pics<-list()
for (i in 1:length(data.objects)) {
  if (data.objects[i] %in% data) {
    pics[[i]] <-eval(parse(text = data.objects[i]))
  }
}
names(pics) <- data.objects

## 1c
plot.equalize <- function(x,main = "",K=50) {
  x.equalized <- equalize(x)
  x.lab <- c("Grauwert")
  y.lab <- c("Anzahl Bildpunkte")
  layout(matrix(1:6,3,2))
    
  plot.array(x,main = main)
  
  breaks <- seq(0,1,len = K+1)
  hist.x <- hist(x,breaks,plot = F)
  hist.x.equalized <- hist(x.equalized,breaks,plot = F)
  plot(hist.x,main = "Absolute Häufigkeit",xlab = x.lab,ylab = y.lab)
  
  hist.x.cum <- cumsum(hist.x$counts) / sum(hist.x$counts)
  hist.x.equalized.cum <- cumsum(hist.x.equalized$counts) / sum(hist.x.equalized$counts)
  plot(hist.x.cum,main = "Kumulative Verteilung",xlab = x.lab,ylab = y.lab)

  plot.array(x.equalized,main = "egalisiert")
  plot(hist.x.equalized,main = "Absolute Häufigkeit eq",xlab = x.lab, ylab = y.lab)
  plot(hist.x.equalized.cum,main = "Kumulative Verteilung eq",xlab = x.lab, ylab = y.lab)

}

for (i in 1:length(data.objects)) {
    plot.equalize(pics[[i]],main = paste(names(pics)[i]))
}
## 1d in der Kommentierung

## 1e in der Kommentierung

## 1f
layout(matrix(1:2,1,2))
plot.array(mandrill)
plot.array(equalize(mandrill,RGB = F))
plot.array(mandrill)
plot.array(equalize(mandrill,RGB = T))

dev.off()

