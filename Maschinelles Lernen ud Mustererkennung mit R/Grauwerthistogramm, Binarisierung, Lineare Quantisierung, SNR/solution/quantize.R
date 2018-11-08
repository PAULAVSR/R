# Team: Redrum
# Task: 2

graphics.off()

## 2a
quantize <- function(x,bits = 8) {
  cells <- 2^bits
  min.value <- (x+0.5)/cells
  max.value <- ((cells-1)+0.5)/cells
  ifelse(x<0,min.value,ifelse(x>=1,max.value,(floor(cells*x)+0.5)/cells))
}

## 2b
test.data <- sin(seq(0,pi,length.out = 1500))
#quantize(test.data)

## 2c
SNR <- function(x,y) {
  10 * log10(mean(x^2)/mean((x-y)^2))
}

SNR.table.generator <- function(table,test.data) {
  table.build <- as.matrix(table)
  new.row <- cbind(SNR(test.data,quantize(test.data)),SNR(test.data,quantize(test.data,12)))
  table.build <- rbind(table.build,new.row)
  return(table.build)
}

SNR.table <- matrix(nrow = 0,ncol = 2)

## 2d
test.data <- sin(seq(0,pi,length.out = 1500))
SNR.table <- SNR.table.generator(SNR.table,test.data)

## 2e
test.data <- runif(1500,0,1)
SNR.table <- SNR.table.generator(SNR.table,test.data)

## 2f
test.data <- rnorm(1500,mean = 0.5,sd = 0.125)
SNR.table <- SNR.table.generator(SNR.table,test.data)

## 2g
test.data <- seq(0,1,length.out = 1500)
SNR.table <- SNR.table.generator(SNR.table,test.data)

## 2h
test.data <- runif(1500,0,1.05)
SNR.table <- SNR.table.generator(SNR.table,test.data)

rownames(SNR.table) <- c("2d","2e","2f","2g","2h")
colnames(SNR.table) <- c("8 Bit", "12 Bit")
SNR.table


graphics.off()
pdf("quantize.pdf")

## 2i
c.vector <- c(seq(1,3.99,length.out = 25),seq(4.1,7,length.out = 25))
results <- matrix(nrow = length(c.vector),ncol = 2)
colnames(results) <- c("C","SNR")

SNR_row <- function(c.vector,results.matrix,method = 'standard',plot = TRUE) {
  for (i in 1:length(c.vector)) {
    if (method == 'set seed') set.seed(21)
    test.data <- rnorm(1500,0.5,1/(2*c.vector[i]))
    results.matrix[i,1] <- c.vector[i]
    results.matrix[i,2] <- SNR(test.data,quantize(test.data))
  }
  
  if (plot) plot(results.matrix, main = "SNR-Werte in Abhängigkeit von C", sub = paste("method = ",method))
  else return(results.matrix)
}

SNR_row(c.vector,results,'set seed',TRUE)
SNR_row(c.vector,results,plot = TRUE)

## Beispiel für vektorisierte SNR_row Funktion
# SNR <- function(x) {
#   y <- quantize(x)
#   10 * log10(mean(x^2)/mean((x-y)^2))
# }
# 
# c.vector <- c(seq(1,3.99,length.out = 25),seq(4.1,7,length.out = 25))
# results <- matrix(nrow = length(c.vector),ncol = 2)
# colnames(results) <- c("C","SNR")
# 
# SNR_row.vectorized <- function(c.vector,results.matrix,method = 'standard',plot = TRUE) {
#   
#   if (method == 'set seed') set.seed(21)
#   test.data <- matrix(nrow = length(c.vector), ncol = 1500)
#   test.data[] <- rnorm(1500,0.5,1/(2*c.vector))
#   results.matrix[,1] <- 1:length(c.vector)
#   results.matrix[,2] <- apply(X = test.data,FUN =SNR,MARGIN = 1)
#   
#   if (plot) plot(results.matrix, main = "SNR-Werte vektorisiert in Abhängigkeit von C", sub = paste("method = ",method))
#   else return(results.matrix)
# }
# 
# SNR_row.vectorized(c.vector,results,'set seed',TRUE)
# SNR_row.vectorized(c.vector,results,plot = TRUE)

dev.off()