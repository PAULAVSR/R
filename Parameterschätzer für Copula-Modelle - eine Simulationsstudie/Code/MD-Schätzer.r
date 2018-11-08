# Minimum-Distanz-Sch?tzer
#
# Deklaration
#
# p = rho = Korrelationskoeffizient
# df = degrees of freedom = Freiheitsgrade
# n = Stichprobengr??e
# col = Spalte in der Ergebnismatrix
#
# Minimum-Distanz-Sch?tzer f?r:
#           n = 100,250,500,1000
#           Gau?-Copula:  f?r p = -0.9,0.8,...,0.8,0.9
#           t-Copula: f?r p = -0.9,0.8,...,0.8,0.9 und df = 3,5,10
#           Gumbel-Copula: f?r p = 1.5,2.0,...,9.5,10.0
#           Clayton-Copula: f?r p = 0.5,1.0,...,9.5,10.0
#
# 2 Konzepte:
# 1) empirische Copula (CvM,L1 CvM,KS)
# 2) Kendall?s dependence function  (CvM,L1 CvM,KS)
# 
# OS: Windows only!!!
# Author: Paul Passek
#

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

cornumb<-makeCluster(4)
registerDoParallel(cornumb)

path = ""
setwd(path)
set.seed(4)
dim = 2
Sim.Length<-c(100,250,500,1000)
p.Normalcop.tcop<-seq(from = -0.9, to = 0.9 , by = 0.2)
p.gumbelcop<-seq(from = 1.5, to = 10 , by = 1)
p.claytoncop<-seq(from = 0.5, to = 10 , by = 1)
df.tcop<-c(10)

#CvM Statistik-empirische Copula
#
#Gauss-Copula
estimated.md.Normalcop = matrix(nrow = 0, ncol = 40)


minimizer.Normalcop = function(param,u,U){
  MD = 0
  n.cop = normalCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = sum((ec-pCopula(u,n.cop))^2)  
  return(MD)
}


estimated.md.Normalcop<-foreach(sim = 1:500,.packages="copula",.combine=rbind)%dopar%{
  col=0
  temp<-matrix(nrow=1,ncol=40)
  for(i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.Normalcop.tcop)){
      p<-p.Normalcop.tcop[j]
      col = col+1
      
      n.cop<- normalCopula(param = p, dim = dim)
      U<-rCopula(n,n.cop)
      u<-matrix(runif(n*dim),n,dim)   
      temp[col]<-optim(0,minimizer.Normalcop,method="L-BFGS-B",lower=-0.99999,upper=0.99999,u=u,U=U)$par
    }     
  } 
  rbind(estimated.md.Normalcop,temp)
}
stopCluster(cornumb)


#t-Copula
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.t.cop = matrix(nrow = 0, ncol = 120)


minimizer.t.cop = function(param,u,U,df){
  MD = 0
  t.cop = tCopula(param = param, df = df, dim = dim)
  ec = C.n(u, U=U)
  MD = sum((ec-pCopula(u,t.cop))^2)  
  return(MD)
}

estimated.md.t.cop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind)%dopar%{
  col=0
  temp= matrix(nrow=1,ncol=120)
  for(i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(h in 1:length(df.tcop)){
      v=df.tcop[h]
      for(j in 1:length(p.Normalcop.tcop)){
        p<-p.Normalcop.tcop[j]
        col = col +1 
        
        t.cop<-tCopula(param = p, df = v, dim = dim)
        U<-rCopula(n,t.cop)
        u<-matrix(runif(n*dim),n,dim)
        temp[col] = optim(0,minimizer.t.cop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, df=v, u=u)$par
      }  
    }    
  }  
  rbind(estimated.md.t.cop,temp)
}
stopCluster(cornumb)


#Gumbel-Copula
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.Gumbelcop = matrix(nrow = 0, ncol = 36)


minimizer.Gumbelcop = function(param,u,U){
  MD = 0
  g.cop = gumbelCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = sum((ec-pCopula(u,g.cop))^2)
  
  return(MD)
}



estimated.md.Gumbelcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col=0
  temp= matrix(nrow=1,ncol=36)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.gumbelcop))  {
      col = col +1 
      p<-p.gumbelcop[j]
      
      g.cop<- gumbelCopula(param = p, dim = dim)
      U<-rCopula(n,g.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(5,minimizer.Gumbelcop, method = "L-BFGS-B", lower = 1.000001, upper = Inf, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.Gumbelcop,temp)
}
stopCluster(cornumb)
write.table(estimated.md.Gumbelcop,file = "estimated.md.Gumbelcop.txt",sep = " ")


#Clayton-Copula
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.Claytoncop = matrix(nrow = 0, ncol = 40)


minimizer.Claytoncop = function(param,u,U){
  MD = 0
  c.cop = claytonCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = sum((ec-pCopula(u,c.cop))^2)
  
  return(MD)
}


estimated.md.Claytoncop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col=0
  temp= matrix(nrow=1,ncol=40)
  
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.claytoncop)) { 
      col = col +1 
      p<-p.claytoncop[j]
      
      c.cop<- claytonCopula(param = p, dim = dim)
      U<-rCopula(n,c.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(5,minimizer.Claytoncop, method = "L-BFGS-B", lower = 0.000001, upper = Inf, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.Claytoncop,temp)
}
stopCluster(cornumb)
write.table(estimated.md.Claytoncop,file = "estimated.md.Claytoncop.txt",sep = " ")



# L1 Variante -empirische Copula
#
#Gauss-Copula-L1
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.L1Normalcop = matrix(nrow = 0, ncol = 40)


minimizer.L1Normalcop = function(param,u,U){
  MD = 0
  n.cop = normalCopula(param = param, dim = 2)
  ec = C.n(u, U=U)
  MD = sum(abs(ec-pCopula(u,n.cop)))
  
  return(MD)
}


estimated.md.L1Normalcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col=0
  temp<-matrix(nrow = 1, ncol = 40)
  
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.Normalcop.tcop)) {
      col = col +1 
      p<-p.Normalcop.tcop[j]
      
      n.cop<- normalCopula(param = p, dim = dim)
      U<-rCopula(n,n.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(0,minimizer.L1Normalcop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.L1Normalcop,temp)
}
stopCluster(cornumb)
write.table(estimated.md.L1Normalcop,file = "estimated.md.L1Normalcop.txt",sep = " ")

#t-Copula-L1
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.L1t.cop = matrix(nrow = 0, ncol = 40)


minimizer.L1t.cop = function(param,u,U,df){
  MD = 0
  t.cop = tCopula(param = param, df = df, dim = dim)
  ec = C.n(u, U=U)
  MD = sum(abs(ec-pCopula(u,t.cop)))
  
  return(MD)
}


estimated.md.L1t.cop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(h in 1:length(df.tcop)) {
      v<-df.tcop[h]
      for(j in 1:length(p.Normalcop.tcop)) {
        col = col +1 
        p<-p.Normalcop.tcop[j]
        
        t.cop<-tCopula(param = p, df = v, dim = dim)
        U<-rCopula(n,t.cop)
        u<-matrix(runif(n*dim),n,dim)
        
        temp[col] = optim(0,minimizer.L1t.cop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, df=v, u=u)$par
        
      }
      
    }
    
  }
  rbind(estimated.md.L1t.cop,temp)
}
stopCluster(cornumb)
write.table(estimated.md.L1t.cop,file = "estimated.md.L1t.cop5.txt",sep = " ")

#Gumbel-Copula-L1
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.L1Gumbelcop = matrix(nrow = 0, ncol = 40)



minimizer.L1Gumbelcop = function(param,u,U){
  MD = 0
  g.cop = gumbelCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = sum(abs(ec-pCopula(u,g.cop)))
  
  return(MD)
}

estimated.md.L1Gumbelcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.gumbelcop)) {
      col = col +1 
      p<-p.gumbelcop[j]
      
      g.cop<- gumbelCopula(param = p, dim = dim)
      U<-rCopula(n,g.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(5,minimizer.L1Gumbelcop, method = "L-BFGS-B", lower = 1.000001, upper = Inf, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.L1Gumbelcop,temp)
}
stopCluster(cornumb)



#Clayton-Copula-L1
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.L1Claytoncop = matrix(nrow = 0, ncol = 40)



minimizer.L1Claytoncop = function(param,u,U){
  MD = 0
  c.cop = claytonCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = sum(abs(ec-pCopula(u,c.cop)))
  
  return(MD)
}


estimated.md.L1Claytoncop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.claytoncop)){ 
      col = col +1 
      p<-p.claytoncop[j]
      
      c.cop <-claytonCopula(param = p, dim = dim)
      U<-rCopula(n,c.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(5,minimizer.L1Claytoncop, method = "L-BFGS-B", lower = 0.000001, upper = Inf, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.L1Claytoncop,temp)
}
stopCluster(cornumb)


#Kolmogorov-Smirnoff-empirische Copula
#
#Gauss-Copula-KS
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.KSNormalcop = matrix(nrow = 0, ncol = 40)


minimizer.KSNormalcop = function(param,u,U){
  MD = 0
  n.cop = normalCopula(param = param, dim = 2)
  ec = C.n(u, U=U)
  MD = max(abs(ec-pCopula(u,n.cop)))
  
  return(MD)
}


estimated.md.KSNormalcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow=1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.Normalcop.tcop)){
      col = col +1 
      p<-p.Normalcop.tcop[j]
      
      n.cop<- normalCopula(param = p, dim = dim)
      U<-rCopula(n,n.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(0,minimizer.KSNormalcop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.KSNormalcop,temp)
}
stopCluster(cornumb)


#t-Copula-KS
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.KSt.cop = matrix(nrow = 0, ncol = 120)


minimizer.KSt.cop = function(param,u,U,df){
  MD = 0
  t.cop = tCopula(param = param, df = df, dim = dim)
  ec = C.n(u, U=U)
  MD = max(abs(ec-pCopula(u,t.cop)))
  
  return(MD)
}


estimated.md.KSt.cop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 120)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for (h in 1:length(df.tcop))  {
      v<-df.tcop[h]
      for(j in 1:length(p.Normalcop.tcop)) { 
        col = col +1 
        p<-p.Normalcop.tcop[j]
        
        t.cop<-tCopula(param = p, df = v, dim = dim)
        U<-rCopula(n,t.cop)
        u<-matrix(runif(n*dim),n,dim)
        
        temp[col] = optim(0,minimizer.KSt.cop, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, df=v, u=u)$par
        
      }
      
    }
    
  }
  rbind(estimated.md.KSt.cop,temp)
}
stopCluster(cornumb)


#Gumbel-Copula-KS
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.KSGumbelcop = matrix(nrow = 0, ncol = 36)


minimizer.KSGumbelcop = function(param,u,U){
  MD = 0
  g.cop = gumbelCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = max(abs(ec-pCopula(u,g.cop)))
  
  return(MD)
}


estimated.md.KSGumbelcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 36)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.gumbelcop)) {
      col = col +1 
      p<-p.gumbelcop[j]
      
      g.cop<- gumbelCopula(param = p, dim = dim)
      U<-rCopula(n,g.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(5,minimizer.KSGumbelcop, method = "L-BFGS-B", lower = 1.000001, upper = Inf, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.KSGumbelcop,temp)
}
stopCluster(cornumb)



#Clayton-Copula-KS
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.md.KSClaytoncop = matrix(nrow = 0, ncol = 40)


minimizer.KSClaytoncop = function(param,u,U){
  MD = 0
  c.cop = claytonCopula(param = param, dim = dim)
  ec = C.n(u, U=U)
  MD = max(abs(ec-pCopula(u,c.cop)))
  
  return(MD)
}


estimated.md.KSClaytoncop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.claytoncop))  { 
      col = col +1 
      p<-p.claytoncop[j]
      
      c.cop <-claytonCopula(param = p, dim = dim)
      U<-rCopula(n,c.cop)
      u<-matrix(runif(n*dim),n,dim)
      
      temp[col] = optim(5,minimizer.KSClaytoncop, method = "L-BFGS-B", lower = 0.000001, upper = Inf, u=u, U=U)$par
      
    }
    
  }
  rbind(estimated.md.KSClaytoncop,temp)
}
stopCluster(cornumb)


# Kendall?s dependence function
#
#CvM Statistik-Kendall?s dependence function
#
#Gauss-Copula-Kendall?s dependence function
cornumb<-makeCluster(2)
registerDoParallel(cornumb)
estimated.Kmd.Normalcop = matrix(nrow = 0, ncol = 40)


minimizer.Normalcop.K = function(param,U,Z){
  MD = 0
  n.cop = normalCopula(param = param, dim = dim)
  MD = sum((Z-pCopula(U,n.cop))^2)
  
  return(MD)
}

tic()
estimated.Kmd.Normalcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1, ncol = 40)
  for( i in 1:length(Sim.Length)){#
    n<-Sim.Length[i]
    for(j in 1:length(p.Normalcop.tcop)){#
      col = col +1 
      p<-p.Normalcop.tcop[j]
      
      n.cop<- normalCopula(param = p, dim = dim)
      U<-rCopula(n,n.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(0,minimizer.Normalcop.K,method="L-BFGS-B",lower=-0.99999,upper=0.99999,U=U, Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.Normalcop,temp)
}
stopCluster(cornumb)
write.table(estimated.Kmd.Normalcop,file = "estimated.Kmd.Normalcop.txt",sep = " ")
toc()

#t-Copula-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.t.cop = matrix(nrow = 0, ncol = 40)


minimizer.t.cop.K = function(param,U,df,Z){
  MD = 0
  t.cop = tCopula(param = param, df = df, dim = dim)
  MD = sum((Z-pCopula(U,t.cop))^2)
  
  return(MD)
}


estimated.Kmd.t.cop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(h in 1:length(df.tcop))  {
      v<-df.tcop[h]
      for(j in 1:length(p.Normalcop.tcop)){ 
        col = col +1 
        p<-p.Normalcop.tcop[j]
        
        t.cop<-tCopula(param = p, df = v, dim = dim)
        U<-rCopula(n,t.cop)
        l=rep(1:n,each = n)
        m=rep(1:n,n)
        
        S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
        Z=tapply(S,l,sum)/(n)
        
        temp[col] = optim(0,minimizer.t.cop.K, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, df=v, Z=Z)$par
        
      }
      
    }
    
  }
  rbind(estimated.Kmd.t.cop,temp)
}
stopCluster(cornumb)
write.table(estimated.Kmd.t.cop,file = "estimated.Kmd.t.cop3.txt",sep = " ")

#Gumbel-Copula-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.Gumbelcop = matrix(nrow = 0, ncol = 36)


minimizer.Gumbelcop.K = function(param,U,Z){
  MD = 0
  g.cop = gumbelCopula(param = param, dim = dim)
  MD = sum((Z-pCopula(U,g.cop))^2)
  
  return(MD)
}


estimated.Kmd.Gumbelcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 36)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.gumbelcop)) {
      col = col +1 
      p<-p.gumbelcop[j]
      
      g.cop<- gumbelCopula(param = p, dim = dim)
      U<-rCopula(n,g.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(5,minimizer.Gumbelcop.K, method = "L-BFGS-B", lower = 1.000001, upper = Inf, U=U,Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.Gumbelcop,temp)
}
stopCluster(cornumb)



#Clayton-Copula-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.Claytoncop = matrix(nrow = 0, ncol = 40)


minimizer.Claytoncop.K = function(param,U,Z){
  MD = 0
  c.cop = claytonCopula(param = param, dim = dim)
  MD = sum((Z-pCopula(U,c.cop))^2)
  
  return(MD)
}


estimated.Kmd.Claytoncop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.claytoncop)) { 
      col = col +1 
      p<-p.claytoncop[j]
      
      c.cop<- claytonCopula(param = p, dim = dim)
      U<-rCopula(n,c.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(5,minimizer.Claytoncop.K, method = "L-BFGS-B", lower = 0.000001, upper = Inf, U=U,Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.Claytoncop,temp)  
}
stopCluster(cornumb)


# L1 Variante -Kendall?s dependence function
#
#Gauss-Copula-L1-Kendall?s dependence function
cornumb<-makeCluster(2)
registerDoParallel(cornumb)
estimated.Kmd.L1Normalcop = matrix(nrow = 0, ncol = 40)


minimizer.L1Normalcop.K = function(param,U,Z){
  MD = 0
  n.cop = normalCopula(param = param, dim = dim)
  MD = sum(abs(Z-pCopula(U,n.cop)))
  
  return(MD)
}

tic()
estimated.Kmd.L1Normalcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.Normalcop.tcop)){
      col = col +1 
      p<-p.Normalcop.tcop[j]
      
      n.cop<- normalCopula(param = p, dim = dim)
      U<-rCopula(n,n.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(0,minimizer.L1Normalcop.K, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.L1Normalcop,temp)
}
stopCluster(cornumb)
write.table(estimated.Kmd.L1Normalcop,file = "estimated.Kmd.L1Normalcop.txt",sep = " ")
toc()

#t-Copula-L1-Kendall?s dependence function
cornumb<-makeCluster(2)
registerDoParallel(cornumb)
estimated.Kmd.L1t.cop = matrix(nrow = 0, ncol = 120)



minimizer.L1t.cop.K = function(param,U,df,Z){
  MD = 0
  t.cop = tCopula(param = param, df = df, dim = dim)
  MD = sum(abs(Z-pCopula(U,t.cop)))
  
  return(MD)
}


estimated.Kmd.L1t.cop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow=1,ncol = 120)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(h in 1:length(df.tcop)) {
      v<-df.tcop[h]
      for(j in 1:length(p.Normalcop.tcop)){ 
        col = col +1 
        p<-p.Normalcop.tcop[j]
        
        t.cop<-tCopula(param = p, df = v, dim = dim)
        U<-rCopula(n,t.cop)
        l=rep(1:n,each = n)
        m=rep(1:n,n)
        
        S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
        Z=tapply(S,l,sum)/(n)
        
        temp[col] = optim(0,minimizer.L1t.cop.K, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, df=v, Z=Z)$par
        
      }
      
    }
    
  }
  rbind(estimated.Kmd.L1t.cop,temp)
}
stopCluster(cornumb)
write.table(estimated.Kmd.L1t.cop,file = "estimated.Kmd.L1t.cop.txt", sep = " ")


#Gumbel-Copula-L1-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.L1Gumbelcop = matrix(nrow = 0, ncol = 36)


minimizer.L1Gumbelcop.K = function(param,U,Z){
  MD = 0
  g.cop = gumbelCopula(param = param, dim = dim)
  MD = sum(abs(Z-pCopula(U,g.cop)))
  
  return(MD)
}


estimated.Kmd.L1Gumbelcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 36)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.gumbelcop)){
      col = col +1 
      p<-p.gumbelcop[j]
      
      g.cop<- gumbelCopula(param = p, dim = dim)
      U<-rCopula(n,g.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(5,minimizer.L1Gumbelcop.K, method = "L-BFGS-B", lower = 1.000001, upper = Inf, U=U,Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.L1Gumbelcop,temp)
}
stopCluster(cornumb)



#Clayton-Copula-L1-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.L1Claytoncop = matrix(nrow = 0, ncol = 40)


minimizer.L1Claytoncop.K = function(param,U,Z){
  MD = 0
  c.cop = claytonCopula(param = param, dim = dim)
  MD = sum(abs(Z-pCopula(U,c.cop)))
  
  return(MD)
}


estimated.Kmd.L1Claytoncop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.claytoncop)) { 
      col = col +1 
      p<-p.claytoncop[j]
      
      c.cop <-claytonCopula(param = p, dim = dim)
      U<-rCopula(n,c.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(5,minimizer.L1Claytoncop.K, method = "L-BFGS-B", lower = 0.000001, upper = Inf, U=U, Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.L1Claytoncop,temp)
}
stopCluster(cornumb)


#Kolmogorov-Smirnoff-Kendall?s dependence function
#
#Gauss-Copula-KS-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.KSNormalcop = matrix(nrow = 0, ncol = 40)



minimizer.KSNormalcop.K = function(param,U,Z){
  MD = 0
  n.cop = normalCopula(param = param, dim = 2)
  MD = max(abs(Z-pCopula(U,n.cop)))
  
  return(MD)
}


estimated.Kmd.KSNormalcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1, ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.Normalcop.tcop)){
      col = col +1 
      p<-p.Normalcop.tcop[j]
      
      n.cop<- normalCopula(param = p, dim = dim)
      U<-rCopula(n,n.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(0,minimizer.KSNormalcop.K, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.KSNormalcop,temp)
}
stopCluster(cornumb)

#t-Copula-KS-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.KSt.cop = matrix(nrow = 0, ncol = 40)


minimizer.KSt.cop.K = function(param,U,df,Z){
  MD = 0
  t.cop = tCopula(param = param, df = df, dim = dim)
  MD = max(abs(Z-pCopula(U,t.cop)))
  
  return(MD)
}


estimated.Kmd.KSt.cop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(h in 1:length(df.tcop)){
      v<-df.tcop[h]
      for(j in 1:length(p.Normalcop.tcop)){ 
        col = col +1 
        p<-p.Normalcop.tcop[j]
        
        t.cop<-tCopula(param = p, df = v, dim = dim)
        U<-rCopula(n,t.cop)
        l=rep(1:n,each = n)
        m=rep(1:n,n)
        
        S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
        Z=tapply(S,l,sum)/(n)
        
        temp[col] = optim(0,minimizer.KSt.cop.K, method = "L-BFGS-B", lower = -0.99999, upper = 0.99999, U=U, df=v, Z=Z)$par
        
      }
      
    }
    
  }
  rbind(estimated.Kmd.KSt.cop,temp)
}
stopCluster(cornumb)
write.table(estimated.Kmd.t.cop,file = "estimated.Kmd.KSt.cop3.txt",sep = " ")

#Gumbel-Copula-KS-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.KSGumbelcop = matrix(nrow = 0, ncol = 36)


minimizer.KSGumbelcop.K = function(param,U,Z){
  MD = 0
  g.cop = gumbelCopula(param = param, dim = dim)
  MD = max(abs(Z-pCopula(U,g.cop)))
  
  return(MD)
}


estimated.Kmd.KSGumbelcop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 36)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.gumbelcop)){
      col = col +1 
      p<-p.gumbelcop[j]
      
      g.cop<- gumbelCopula(param = p, dim = dim)
      U<-rCopula(n,g.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(5,minimizer.KSGumbelcop.K, method = "L-BFGS-B", lower = 1.000001, upper = Inf, U=U, Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.KSGumbelcop,temp)
}
stopCluster(cornumb)


#Clayton-Copula-KS-Kendall?s dependence function
cornumb<-makeCluster(4)
registerDoParallel(cornumb)
estimated.Kmd.KSClaytoncop = matrix(nrow = 0, ncol = 40)


minimizer.KSClaytoncop.K = function(param,U,Z){
  MD = 0
  c.cop = claytonCopula(param = param, dim = dim)
  MD = max(abs(Z-pCopula(U,c.cop)))
  
  return(MD)
}


estimated.Kmd.KSClaytoncop<-foreach(sim = 1:500, .packages = "copula",.combine = rbind) %dopar% {
  col = 0
  temp<-matrix(nrow = 1,ncol = 40)
  for( i in 1:length(Sim.Length)){
    n<-Sim.Length[i]
    for(j in 1:length(p.claytoncop)){ 
      col = col +1 
      p<-p.claytoncop[j]
      
      c.cop <-claytonCopula(param = p, dim = dim)
      U<-rCopula(n,c.cop)
      l=rep(1:n,each = n)
      m=rep(1:n,n)
      
      S=((U[l,1]>U[m,1])&(U[l,2]>U[m,2]))
      Z=tapply(S,l,sum)/(n)
      
      temp[col] = optim(5,minimizer.KSClaytoncop.K, method = "L-BFGS-B", lower = 0.000001, upper = Inf, U=U,Z=Z)$par
      
    }
    
  }
  rbind(estimated.Kmd.KSClaytoncop,temp)
}
stopCluster(cornumb)

save.image("MD-Sch?tzer.RData")