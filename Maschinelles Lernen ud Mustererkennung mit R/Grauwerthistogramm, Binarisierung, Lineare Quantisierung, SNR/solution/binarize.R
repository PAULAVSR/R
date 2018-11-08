# Task 1
graphics.off()

load("../binarize.rda")


pdf("binarize.pdf")

## 1a
par(mfrow = c(2,2))
plot.array(algae)
hist(algae,breaks = seq(min(algae),max(algae),length = 1+25))

## 1b
alga <- algae[210:260,60:110]
plot.array(alga)
hist(alga,breaks = seq(min(alga),max(alga),length.out = 1+25))

## 1c/d
binarize <- function(x,method = 'fixed', threshold=0.5,plot = FALSE) {
  threshold.method <- switch(method,
         fixed = threshold,
         mean = mean(x),
         median = median(x),
         inter = mean(kmeans(matrix(x,ncol = 1),c(min(x),max(x)))$centers)
  )
  binary.select <- ifelse(x<threshold.method,0,1)
  if (plot == T) plot.array(binary.select,main = paste("threshold = ", round(threshold,2)),sub = paste("method = ",method))
  return(x)
}

par(mfrow = c(3,3))

for (i in seq(0.20,0.50,length.out = 9)) {
  binarize(alga,method = 'fixed',threshold = i,plot = T)
}

## 1e
plot.binarized <- function(image = list(...),method = c(...), mfrow = c(1,1)) {
  for (pic in image) {
    par(mfrow = mfrow)
    for (m in method)
      binarize(pic,method = m,plot = T)
  }
}

plot.binarized(image = list(alga,algae,tonga),method = c("fixed","mean","median","inter"),mfrow = c(2,2))

dev.off()
