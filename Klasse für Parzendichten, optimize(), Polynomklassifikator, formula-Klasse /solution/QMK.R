## Team: Redrum
## Task: 2

graphics.off()

# 2a
cumex <- function(n, dmax) {
	if(dmax == 0) {
		matrix(rep(0, n), ncol=n)
	}
	else {
	if(n == 1) {
		matrix(0:dmax, nrow = dmax+1)
	} else {
		mat <- cumex(n-1, dmax)
		mat <- apply(mat, 2, FUN=function(x) rep(x, each=dmax+1))
		mat <- cbind(mat, 0:dmax)
		matrix(mat[rowSums(mat) <= dmax,], ncol=n)
	}
	}
}

# 2b
polyterms.RHS <- function(n, dmax, vname) {
	paste(apply(cumex(n, dmax), 1, FUN=function(x) str <- paste("I(",paste(vname,"^",x,sep="",collapse = "*"),collapse = "+",")")),collapse = "+")
}

# 2c
QMK <- function(x, dmax) {
	stopifnot(is.data.frame(x))
	stopifnot(is.factor(x[[ncol(x)]]))
	stopifnot(dmax >= 0)
	
	classes <- levels(x[[ncol(x)]])
	label <- names(x)[[ncol(x)]]
	xs <- x[,-ncol(x)]
	vname <- names(xs)
	funs <- polyterms.RHS(ncol(xs), dmax, vname)
	x[,ncol(x)] <- sapply(x[,ncol(x)], FUN=as.numeric) - 1
	model <- lm(paste(label,paste("0+",funs), sep="~"), x)
  
## Since there are some singularities in the fitted lm models: some code to fix warnings
 	# if (length(model$coefficients[!is.na(model$coefficients)][-1]) != 0){
 	#   temp <- paste(names(model$coefficients[!is.na(model$coefficients)])[-1],collapse = "+")
 	#   model <- update(model,formula. = as.formula(paste(label,temp,sep = "~")))
 	# }

	structure(.Data= list(model=model, classes=classes),class = "QMK")
}

# 2d
is.QMK <- function(o) inherits(o,"QMK")

predict.QMK <- function(o, newdata) {
	stopifnot(is.data.frame(newdata))
	stopifnot(is.QMK(o))
	stopifnot(is.factor(newdata[[ncol(newdata)]])==F)
	
	factor(ifelse(predict.lm(o$model, newdata) <= 0.5, o$classes[1], o$classes[2]), o$classes)
} 

# 2e
dim.QMK <- function(o) {
  stopifnot(is.QMK(o))
  
  length(coefficients(o$model))
}

heldout.dim <- function(x, newdata, method=QMK, ...) {
  o <- method(x, ...)
  classes <- predict(o,newdata[,-ncol(newdata)])
  c(dim(o),sum(ifelse(classes == newdata[,ncol(newdata)],0,1))/nrow(newdata))
} 

## 2f
error.table <- function(x,newdata,modelset,...) {
  sol <- sapply(modelset,function(i) heldout.dim(x,newdata,dmax=i),simplify = T,USE.NAMES = T)
  sol <- rbind(sol,sapply(modelset,function(i) heldout.dim(x,x,dmax=i),simplify = T,USE.NAMES = T))[-3,]
  rownames(sol) <- c("modelsize","lern:test","lern:lern")
  colnames(sol) <- paste("Polynomgrad",modelset)
  plot(sol[1,],sol[2,]*100,type = "b",xlim = range(sol[1,]),ylim = range(sol[2:3,])*100,col="tomato",log="x",pch=17,xlab = "Anzahl Polynomterme",ylab="Fehlerrate [%]", ...)
  points(sol[1,],sol[3,]*100,type = "b",col="navyblue",pch=17,lty=2)
  legend("bottomleft",legend = c("Testdaten","Lerndataen"),col = c("tomato","navyblue"),lty = c(1,2))
  signif(sol,3)
}

## 2g
dataset <-c("heart","diabetes","australia","segment","vehicle")
links <- paste("../",dataset,".rda",sep = "")
names(links) <- dataset
data <- lapply(links,function(link) mget(load(link)))

error.table.caller <- function(i,name) {
  sol <- error.table(as.data.frame(i[[1]]), as.data.frame(i[[2]]), 0:3, main=name)
  rownames(sol) <- paste(rownames(sol),name)
  sol
}

solution.raw <- do.call(rbind, mapply(data, dataset, FUN=function(i,n)error.table.caller(i,n), SIMPLIFY =F))
solution.raw


# only for assignment
#solution.dims <- solution.raw[c(1,4,7,10,13),]
#write.table(solution.dims,"dims_2g.txt", sep="\t")
#solution.erros <- signif(solution.raw[-c(1,4,7,10,13),]*100,3)
#rownames(solution.erros) <- paste(rownames(solution.erros),"errors[%]")
#write.table(solution.erros,"errors_2g.txt", sep="\t")

## 2h
x <- error.table(data$diabetes$diabetes.lern,data$diabetes$diabetes.lern,0:5,main="diabetes")
x

#dev.off()
