graphics.off()
pdf(file = "sinus.pdf")
layout(matrix(c(1:6),nrow = 2, ncol = 3, byrow = T))

A <- function(n) { 1/n }
sinus.wave <- function(x,A,n,phi) { A(n) * sin(2 * pi * n * x + phi) }


  
#part a
curve(sinus.wave(x,A,1,0), from = 0, to = 1, main = "part a)")  
for (n in 2:5) {
  curve(sinus.wave(x,A,n,0), from = 0, to = 1, col = n, add = T)  
}


#part b
sinus.wave.sum <- function (x,A,n,start,by,phi) {
  sinus.sum <- 0
  for (i in seq(from = start,to = n,by = by)) {
    sinus.sum <- sinus.sum + sinus.wave(x,A,i,phi)
  }
  sinus.sum
  

}

curve(sinus.wave.sum(x,A,5,1,1,0),from = 0, to = 1, main = "part b)", col = 2)

#part c
A2 <- function(n) { (1/2)^n }

curve(sinus.wave.sum(x,A2,5,1,1,0),from = 0, to = 1, main = "part c)", col = 3)

#part d
A3 <- function(n) { 1/sqrt(n) }

curve(sinus.wave.sum(x,A3,5,1,1,0),from = 0, to = 1, main = "part d)", col = 4)

#part e
phi <- pi/3

curve(sinus.wave.sum(x,A,5,1,1,phi),from = 0, to = 1, main = "part e)", col = 5)

#part f

curve(sinus.wave.sum(x,A,10,2,2,0),from = 0, to = 1, main = "part f)", col = 6)


dev.off()
