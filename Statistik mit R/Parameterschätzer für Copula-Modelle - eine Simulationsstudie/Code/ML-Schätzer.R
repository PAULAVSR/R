#Maximum-Likelihood-Schaetzer
#
#Deklaration
#
# p = rho = Korrelationskoeffizient
# df = degrees of freedom = Freiheitsgrade
# n = Stichprobengröße
# col = Spalte in der Ergebnismatrix
#
#Maximum-Likelihood-Schätzer für:
#           n = 100,250,500,1000
#           Gauß-Copula:  für p = -0.9,-0.7,...,0.7,0.9
#           t-Copula: für p = -0.9,-0.7,...,0.7,0.9 und df = 3,5,10
#           Gumbel-Copula: für p = 1.5,2.5,...,8.5,9.5
#           Clayton-Copula: für p = 0.5,1.5,...8.5,9.5
#
# OS: Windows only!
# Author: Paul Passek

activate.package <- function(pack) {
  
  if (pack %in% installed.packages())
    require(pack,character.only = T)
  else {
    install.packages(deparse(substitute(pack)))
    require(pack,character.only = T)
  }
}

activate.package(pack="copula")
activate.package(pack="foreach")
activate.package(pack="doParallel")
activate.package(pack="tictoc")


cornumb<-makeCluster(2)
registerDoParallel(cornumb)

set.seed(100)
dim = 2
Sim.Length<-c(100,250,500,1000)
p.Normalcop.tcop<-seq(from = -0.9, to = 0.9 , by = 0.2)
p.gumbelcop<-seq(from = 1.5, to = 10 , by = 1)
p.claytoncop<-seq(from = 0.5, to = 10 , by = 1)
df.tcop<-c(3,5,10)


#Gau?-Copula

estimated.overall.Normalcop<- matrix(nrow = 0, ncol = 40) # Ergebnismatrix 1001x114

likelihood.Normalcop = function(param, n.cop_estim){ # Likelihoodfunktion, der Param (=rho) und Stichprobe n.cop_estim ?bergeben wird
  ll<-0
  n.cop<-normalCopula(param = param, dim = dim) # Konstruktion der Referenzcopula mit festem Parameter rho
  ll<--sum(dCopula(n.cop_estim,n.cop,log = TRUE)) # Likelihoodfunktion: negative Summe aus logarithmierten Beobachtungswerten n.cop_estim 
  return(ll)                                      # eingesetzt in die Dichtefunktion der Normal-Copula
}

 # Bef?llung der Ergebnismatrix initialisieren

estimated.overall.Normalcop<-foreach(sim = 1:500,.packages="copula",.combine=rbind)%dopar%{
  col=0
  temp<-matrix(nrow=1,ncol=40)
for (i in 1:length(Sim.Length)){ # Simulation f?r verschiedene Stichprobenumf?nge n
  n<-Sim.Length[i]
  for (j in 1:length(p.Normalcop.tcop)){ # Simulation f?r verschiedene rho p
    p<-p.Normalcop.tcop[j]
    col = col +1 
    
      
      n.cop<-normalCopula(param = p,dim = dim) # Referenzcopula f?r Stichprobe n.cop_estim
      n.cop_estim<-rCopula(n,n.cop) # Simulation der Stichprobencopula f?r Stichprobengr??e n basierend auf Referenzcopula n.cop
      temp[col]<-optim(0, likelihood.Normalcop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, n.cop_estim = n.cop_estim)$par #Optimierung der Likelihoodfunktion
      
    }
  }
  rbind(estimated.overall.Normalcop,temp)
}
stopCluster(cornumb)
write.table(estimated.overall.Normalcop,file = "estimated.overall.Normalcop.txt",sep = " ")

#t-Copula
cornumb<-makeCluster(2)
registerDoParallel(cornumb)

estimated.overall.tcop<- matrix(nrow = 0, ncol = 120) 


likelihood.tcop = function(param, t.cop_estim, df){ 
  ll<-0
  t.cop<-tCopula(param = param, df = df,dim = dim )
  ll<--sum(dCopula(t.cop_estim,t.cop,log = TRUE))
  return(ll)
}

estimated.overall.tcop<-foreach(sim = 1:500,.packages="copula",.combine=rbind)%dopar%{
  temp<-matrix(nrow = 1,ncol = 120)
  col = 0
for (i in 1:length(Sim.Length)){
  n<-Sim.Length[i]  
  for (h in 1:length(df.tcop)){
    v<-df.tcop[h]
    for (j in 1:length(p.Normalcop.tcop)){
      p<-p.Normalcop.tcop[j]
      col = col +1 
      
        
        t.cop<-tCopula(param = p, df = v,dim = dim)
        t.cop_estim<-rCopula(n,t.cop)
        temp[col] <- optim(0, likelihood.tcop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, t.cop_estim = t.cop_estim, df = v)$par
        
      }
    }
  }
  rbind(estimated.overall.tcop,temp)
}
stopCluster(cornumb)
write.table(estimated.overall.tcop,file = "estimated.overall.tcop.txt",sep = " ")

#Gumbel-Copula
cornumb<-makeCluster(2)
registerDoParallel(cornumb)
estimated.overall.gumbelcop<- matrix(nrow = 0, ncol = 36)


likelihood.gumbelcop = function(param, gu.cop_estim){
  ll<-0
  gu.cop<-gumbelCopula(param = param,dim = dim)
  ll<- -sum(dCopula(gu.cop_estim,gu.cop,log = TRUE))
  return(ll)
}


estimated.overall.gumbelcop<-foreach(sim = 1:500,.packages="copula",.combine=rbind)%dopar%{
  col=0
  temp=matrix(nrow=1,ncol=36)
  
  for (i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for (j in 1:length(p.gumbelcop)){
    p<-p.gumbelcop[j]
    col = col +1
    
    gu.cop<-gumbelCopula(param = p,dim = dim)
    gu.cop_estim<-rCopula(n,gu.cop)
    temp[col] = optim(6.5, likelihood.gumbelcop, method = "L-BFGS-B", lower = 1, upper = 100, gu.cop_estim = gu.cop_estim)$par
      
    }
  }
  rbind(estimated.overall.gumbelcop,temp)
}
stopCluster(cornumb)
write.table(estimated.overall.gumbelcop,file = "estimated.overall.gumbelcop.txt",sep = " ")

#Clayton-Copula
cornumb<-makeCluster(2)
registerDoParallel(cornumb)
estimated.overall.claytoncop<- matrix(nrow = 0, ncol = 40)

likelihood.claytoncop = function(param, cl.cop_estim){
  ll<-0
  cl.cop<-claytonCopula(param = param,dim = dim)
  ll<- -sum(dCopula(cl.cop_estim,cl.cop,log = TRUE))
  
  return(ll)
}
tic()
estimated.overall.claytoncop<-foreach(sim = 1:500,.packages="copula",.combine=rbind)%dopar%{
  col = 0
  temp<-matrix(nrow = 1, ncol = 40)
for (i in 1:length(Sim.Length)){
  n<-Sim.Length[i]
  for (j in 1:length(p.claytoncop)){
    p<-p.claytoncop[j]
    col = col + 1
    
    cl.cop<-claytonCopula(param = p,dim = dim)
    cl.cop_estim<-rCopula(n,cl.cop)
    temp[col] <- optim(5.5, likelihood.claytoncop, method = "L-BFGS-B", lower = 0.000001, upper = Inf, cl.cop_estim = cl.cop_estim)$par
      
    }
  }
  rbind(estimated.overall.claytoncop,temp)
}
stopCluster(cornumb)
toc()
write.table(estimated.overall.claytoncop,file = "estimated.overall.claytoncop.txt",sep = " ")
