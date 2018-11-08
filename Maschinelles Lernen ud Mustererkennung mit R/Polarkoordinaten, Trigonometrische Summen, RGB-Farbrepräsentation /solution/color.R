graphics.off()

load("../color.rda")

#global
pic.arr <- array(dim = c(25,25,3))
intense <- seq(0,1,1/25)
pdf("color.pdf")

#part a
for (i in 1:dim(pic.arr)[1]) {
  for (j in 1:dim(pic.arr)[2]) {
    pic.arr[i,j,1:3] <- intense[i]
  }
}

plot.array(pic.arr,main = "part a")

#part b
pic.arr <- array(dim = c(25,25,3))
for (j in 1:dim(pic.arr)[1]) {
  for (i in 1:dim(pic.arr)[2]) {
    pic.arr[i,j,1] <- 0
    pic.arr[i,j,2] <- intense[j]
    pic.arr[i,j,3] <- 0
  }
}
plot.array(pic.arr, main = "part b")

#part c
for (i in 1:dim(pic.arr)[1]) {
  for (j in 1:dim(pic.arr)[2]) {
    pic.arr[i,j,1] <- intense[i]
    pic.arr[i,j,2] <- 0
    pic.arr[i,j,3] <- intense[j]
  }
}
plot.array(pic.arr,main = "part c")

#part d
part.d <- function (){
  pic.arr <- array(dim = c(25,25,3))
    for (i in 1:dim(pic.arr)[1]) {
      for (j in 1:dim(pic.arr)[2]) {
        pic.arr[i,j,1] <- intense[i]
        pic.arr[i,j,3] <- intense[j]
        
        if (pic.arr[i,j,1]+pic.arr[i,j,3] > 1) {
          pic.arr[i,j,1] <- 0
          pic.arr[i,j,2] <- 0
          pic.arr[i,j,3] <- 0
        } else {
          pic.arr[i,j,2] <- 1- (pic.arr[i,j,1]+pic.arr[i,j,3])
        }
      }
    }
  return(pic.arr)
}

plot.array(part.d(),main = "part d")

#part e
pic.arr <- part.d()
for (i in 1:dim(pic.arr)[1]) {
  for (j in 1:dim(pic.arr)[2]) {
    max.rgb <- max(pic.arr[i,j,])
    if (max.rgb > 0) {
      pic.arr[i,j,] <- pic.arr[i,j,] * 1/max.rgb
    } else {
      pic.arr[i,j,] <- 0
    }
  }
}

plot.array(pic.arr,main = "part e")

#part f
pic.arr <- part.d()
for (i in 1:dim(pic.arr)[1]) {
  for (j in 1:dim(pic.arr)[2]) {
    max.rgb <- max(pic.arr[i,j,])
    if (max.rgb > 0) {
      pic.arr[i,j,] <- pic.arr[i,j,] + 1-max.rgb
    } else {
      pic.arr[i,j,] <- 0
    }
  }
}

plot.array(pic.arr,main = "part f")

dev.off()
