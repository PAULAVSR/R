# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Bevor das R-Skript ausgef?hrt werden kann, geben Sie eine Simulationsvariante an.! 
# Sie haben zur Auswahl:                                                           !
#                                                                                  !
#     1. (500,1000,5000,10000,15000) Wiederholungen ca. 4 min bei 4 Kernen         !
#                                                                                  !
#     2. (500,1000,10000,50000,100000) Wiederholungen ca. 20 min bei 4 Kernen      !
#                                                                                  !
VARIANTE <- 2 #1                                                                   !
#                                                                                  !
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# OS: Windows only
# Author: Paul Passek

# Package Installer
activate.package <- function(pack) {
  
  if (pack %in% installed.packages())
    require(pack,character.only = T)
  else {
    install.packages(deparse(substitute(pack)))
    require(pack,character.only = T)
  }
}

activate.package(pack="tictoc")
activate.package(pack="doParallel")
activate.package(pack="foreach")
activate.package(pack="ADGofTest")

laeuse.data <- read.csv(text = "../islands.csv")



if (VARIANTE == 1) {
  count = c(500,1000,5000,10000,15000)
  print(paste("Variante 1 gew?hlt. Setting f?r Wiederholungen wurde gesetzt: ", count))
} else if (VARIANTE == 2) {
  count = c(500,1000,10000,50000,100000)
  print(paste("Variante 2 gew?hlt. Setting f?r Wiederholungen wurde gesetzt: ", count))
} else {
  count = c(500,1000,5000,10000,15000)
}

set.seed(203982)
cornumb <- detectCores(logical = T)
cluster <- makeCluster(cornumb)
registerDoParallel(cluster)


n.all<-c(30,50,80,100,120)
result.alln <- matrix(nrow = 512, ncol = 50)
tmper <- matrix(nrow = 512, ncol = 50)
result.all.anb <- matrix(nrow = 50,ncol = 4)
result.ad.test <- matrix(nrow = 26,ncol = 4)
prob.alln <- vector(length = 5)
theta0.alln <- vector(length = 5)
col <- 0
k<-1
tic()
result.alln <- foreach (j = 1:length(n.all),.combine = cbind)%do%{
  print(n.all[j])
  n<-n.all[j]
  tmper [1,col] = n 
  subdata <- laeuse.data[sample(59, n,replace = T), ]
  prob <- sum(subdata$mt.presence)/n
  var <- prob*(1-prob)
  theta0 <- log(prob/(1-prob))
  prob.alln[j] <- prob
  theta0.alln[j] <- theta0
  #Submodel f?r n Stichproben:
  subModelMt <- glm(mt.presence ~ (1+size + no.hab), family = binomial, data = subdata)
  X <- model.matrix(subModelMt)
  
  #Inverse Fisher Info Matrix:
  V <- var * diag(n)
  fisher<- t(X) %*%V%*% X
  #eigen(fisher) # alle Eigenwerte >0 -> positiv definit
  f1 <-solve(fisher)
  
  #Jetzt wiederholung des Experiments m-fach:
  betas <- numeric(0)
  param <- 1 
  result.dens.x <- 0
  result.dens.y <- 0
  result.anb <- matrix(nrow = 10,ncol = 3)
  
  for (c in 1:length(count)){
    #set.seed(c*abs(rnorm(1))*1000)
    col = col + 1
    #betas <- vector(length = count[c])
    counter <- count[c]
    tmp.beta <-0
    row <- 0
    tmp.beta[row] = counter
    betas <- foreach(i = 1:counter,.combine = rbind)%dopar%{
      row = row + 1
      subdata$mt.presence <- rbinom(n,1,prob)
      subModelMt <- glm(mt.presence ~ (1+size + no.hab), family = binomial, data = subdata)
      betaNeu <- coef(subModelMt)[param]
      tmp.beta[row] <- betaNeu
    }
    p.val <- matrix(nrow = 1000,ncol = 1)
    p.val[,] <-foreach (i = 1:1000,.combine = rbind)%dopar%{
      tmp <-(as.numeric(shapiro.test(sample(betas,min(count[c],5000)))[2]))
    }
    #ks.test(betas,"pnorm",mean(betas),sd(betas))
    #rbind(betas,tmp.beta)
    result.ad.test[(col-col/2)+1,1] <- n
    result.ad.test[(col-col/2)+1,2] <- count[c]
    result.ad.test[(col-col/2)+1,3] <- as.numeric(ad.test(betas,distr.fun = pnorm,mean = mean(betas),sd = sd(betas)))[2]
    result.ad.test[(col-col/2)+1,4] <- mean(p.val[,1])
    
    # result muss in ?u?ere Variable speichern
    result.all.anb[0+k,1] <- n
    result.all.anb[0+k,2] <- count[c]
    result.all.anb[0+k,3] <- var(betas)
    result.all.anb[0+k,4] <- mean(betas)
    result.all.anb[1+k,3] <- f1[param,param]
    result.all.anb[1+k,4] <- theta0.alln[j]
    
    k = k+2
    dens.x <- density(betas)$x
    #result.dens.x <-cbind(result.dens.x,dens)
    dens.y <- density(betas)$y
    #result.dens.y <- cbind(result.dens.y,dens)
    
    tmper[,col] <- dens.x
    col = col +1
    tmper[,col] <- dens.y
  }
}
result.alln <- tmper

stopCluster(cluster)
toc()
tmp.names <- c(paste("30:",count[1]),paste("30:",count[1]),paste("30:",count[2]),paste("30:",count[2]),paste("30:",count[3]),paste("30:",count[3]),paste("30:",count[4]),paste("30:",count[4]),paste("30:",count[5]),paste("30:",count[5]),paste("50:",count[1]),paste("50:",count[1]),paste("50:",count[2]),paste("50:",count[2]),paste("50:",count[3]),paste("50:",count[3]),paste("50:",count[4]),paste("50:",count[4]),paste("50:",count[5]),paste("50:",count[5]),paste("80:",count[1]),paste("80:",count[1]),paste("80:",count[2]),paste("80:",count[2]),paste("80:",count[3]),paste("80:",count[3]),paste("80:",count[4]),paste("80:",count[4]),paste("80:",count[5]),paste("80:",count[5]),paste("100:",count[1]),paste("100:",count[1]),paste("100:",count[2]),paste("100:",count[2]),paste("100:",count[3]),paste("100:",count[3]),paste("100:",count[4]),paste("100:",count[4]),paste("100:",count[5]),paste("100:",count[5]),paste("120:",count[1]),paste("120:",count[1]),paste("120:",count[2]),paste("120:",count[2]),paste("120:",count[3]),paste("120:",count[3]),paste("120:",count[4]),paste("120:",count[4]),paste("120:",count[5]),paste("120:",count[5]))
colnames(result.alln) <- tmp.names

x <- seq(-2,2,0.0001)
column <- c(1,11,21,31,41)
#par(mfrow=c(3,2))
col <- c("tomato","orange","darkgreen","skyblue","purple")
for (i in 1:length(column)) {
  y <- dnorm(x, mean=theta0.alln[i], sd=sqrt(f1[param,param]))
  plot(x,y,type = "l",col=1,main = paste("Stichprobe: ",n.all[i]))
  start.col <- column[i]
  points(result.alln[,start.col],result.alln[,start.col+1],type = "l", col = col[1])
  points(result.alln[,start.col+2],result.alln[,start.col+3],type = "l", col = col[2])
  points(result.alln[,start.col+4],result.alln[,start.col+5],type = "l", col = col[3])
  points(result.alln[,start.col+6],result.alln[,start.col+7],type = "l", col = col[4])
  points(result.alln[,start.col+8],result.alln[,start.col+9],type = "l", col = col[5])
  legend("topright",c("500","1.000","10.000","50.000","100.000"),lty = 1,col = col,bty = "n",cex = 0.75)
}

# p-Wert Grafik Anderson-Darling
plot(result.ad.test[,1],result.ad.test[,3],xlab = "Stichprobengr??e", ylab = "p-Wert", main = "p-Werte des Anderson-Darling Tests")
points(result.ad.test[which(result.ad.test[,2] == 500),1],result.ad.test[which(result.ad.test[,2] == 500),3],col = "green",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 1000),1],result.ad.test[which(result.ad.test[,2] == 1000),3],col = "blue",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 10000),1],result.ad.test[which(result.ad.test[,2] == 10000),3],col = "orange",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 50000),1],result.ad.test[which(result.ad.test[,2] == 50000),3],col = "grey",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 100000),1],result.ad.test[which(result.ad.test[,2] == 100000),3],col = "purple",pch = 16)

abline(h = 0.05, col = "red")

legend("center",c("500","1000","10000","50000","100000"),col = c("green","blue","orange","grey","purple"),pch = c(16,16,16,16,16), cex = 0.75)


# p-Wert Shapiro-Wilk
plot(result.ad.test[,1],result.ad.test[,4], xlab = "Stichprobengr??e", ylab = "p-Wert", main = "p-Werte des Shapiro-Wilk Tests")
points(result.ad.test[which(result.ad.test[,2] == 500),1],result.ad.test[which(result.ad.test[,2] == 500),4],col = "green",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 1000),1],result.ad.test[which(result.ad.test[,2] == 1000),4],col = "blue",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 10000),1],result.ad.test[which(result.ad.test[,2] == 10000),4],col = "orange",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 50000),1],result.ad.test[which(result.ad.test[,2] == 50000),4],col = "grey",pch = 16)

points(result.ad.test[which(result.ad.test[,2] == 100000),1],result.ad.test[which(result.ad.test[,2] == 100000),4],col = "purple",pch = 16)

abline(h = 0.05, col = "red")

legend("topleft",c("500","1000","10000","50000","100000"),col = c("green","blue","orange","grey","purple"),pch = c(16,16,16,16,16), cex = 0.75)


# Varianz- und Mittelwertdifferenz
x<-rbind(result.all.anb[which(result.all.anb[,2] == count[5]),3],n.all)
y<-rbind(result.all.anb[which(result.all.anb[,2] == count[5])+1,3],n.all)

a <- rbind(result.all.anb[which(result.all.anb[,2] == count[5]),4],n.all)
b <- rbind(result.all.anb[which(result.all.anb[,2] == count[5])+1,4],n.all)

options(scipen = 10)
plot(x[2,],x[1,]-y[1,],type = "b",xlab = "Stichprobengr??e",ylab = "simulierte Varianz - theoretische Varianz",main= paste("Varianzverlauf bei ",count[5]," Wiederholungen"),pch = 1,col = "darkgreen")

plot(a[2,],a[1,]-b[1,],type = "b", xlab = "Stichprobengr??e", ylab = "simulierter Mittelwert - theoretischer Mittelwert", main = paste("Mittelwertverlauf bei ", count[5]," Wiederholungen"), pch = 1, col = "skyblue")

cbind(x[2,],x[1,]-y[1,])
cbind(a[2,],a[1,]-b[1,])
