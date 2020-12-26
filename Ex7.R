library(MASS)
library(plyr)
library(MLmetrics)
library(leaps)
#Exercise 1
baseball <-read.csv("baseball.txt",header=TRUE)
varnames <-names(baseball)[-1]
response <-names(baseball)[1]
unicef.full.model <-as.formula(paste(response,"~",paste(varnames,collapse=" + ")))
fullmod=lm(unicef.full.model,data = baseball)
empty <-lm(baseball$y~1,data=baseball)
backward=step(fullmod,direction="backward", data=unicef)
forward=step(empty,direction="forward",scope=unicef.full.model, data=unicef)

best=lm(y ~ x3 + x7 + x8 + x10 + x11 + x12 + x13 + x14 + x15, data=baseball)
summary(best)

par(mfrow=c(2,2))
plot(best,1)
plot(best,2)
plot(best,3)

box_cox <- boxcox(best, lambda = seq(-1,1, by=.1), main = "Box-Cox Power Transformation")
lamda_bc <- box_cox$x[which.max(box_cox$y)]
roundlamda_bc <-round_any(lamda_bc, .5)  

best_log=lm(log(y) ~ x3 + x7 + x8 + x10 + x11 + x12 + x13 + x14 + x15, data=baseball)
summary(best_log)
plot(best_log,3)

n=337
yhat <- numeric(0)
yhat_log= numeric(0)
for (i in 1:n){
  bs <- lm(y ~ x3 + x7 + x8 + x10 + x11 + x12 + x13 + x14 + x15, data=baseball[-i,])
  yhat[i] <- predict(bs,baseball[i,])
  bs_log=lm(log(y) ~ x3 + x7 + x8 + x10 + x11 + x12 + x13 + x14 + x15, data=baseball[-i,])
  mse=sum(bs_log$residuals^2)/bs_log$df.residual
  yhat_log[i] <- exp(predict(bs_log,baseball[i,]))*exp(mse*0.5)
  
}
sqloss <- sqrt(mean((yhat-baseball$y)^2))
sqloss_log <- sqrt(mean((yhat_log-baseball$y)^2))


#Exercise 2
P=16
bestsub <- regsubsets(unicef.full.model,data=baseball,nvmax=16,method="exhaustive")
sbest <- summary(bestsub)

sqrlossbest_LOO <- numeric(0)
yhat_LOO <- matrix(NA,nrow=P,ncol=n)
modfor <- list()

for (j in 1:P){
  modfor[[j]] <-as.formula(paste(response,"~",paste(varnames[sbest$which[j,2:(P+1)]],
                                                    collapse=" + ")))
}
for (j in 1:P){
  for (i in 1:n){
    fmi <- lm(modfor[[j]], data=baseball[-i,])
    yhat_LOO[j,i] <- predict(fmi,baseball[i,])
  }
  sqrlossbest_LOO[j] <- sqrt(mean((yhat_LOO[j,]-baseball$y)^2))
}

m_LOO=match(min(sqrlossbest_LOO),sqrlossbest_LOO)


###10-fold

#random shuffle data
baseball_shuff<-baseball[sample(nrow(baseball)),]

#Create 10 equally size folds
indexes <- cut(seq(1,nrow(baseball_shuff)),breaks=10,labels=FALSE)
Indexes=list()
for(i in 1:10){
  Indexes[[i]] <- which(indexes==i,arr.ind=TRUE)
  
}



sqrlossbest_10F <- numeric(0)
yhat_10F <- matrix(NA,nrow=P,ncol=n)
for (j in 1:P){
  for (i in 1:10){
    fmi <- lm(modfor[[j]], data=baseball_shuff[-Indexes[[i]],])
    yhat_10F[j,Indexes[[i]]] <- predict(fmi,baseball_shuff[Indexes[[i]],])
  }
  sqrlossbest_10F[j] <- sqrt(mean((yhat_10F[j,]-baseball_shuff$y)^2))
}


#5 Fold cv
baseball_shuff<-baseball[sample(nrow(baseball)),]

#Create 10 equally size folds
indexes <- cut(seq(1,nrow(baseball_shuff)),breaks=5,labels=FALSE)
Indexes=list()
for(i in 1:5){
  Indexes[[i]] <- which(indexes==i,arr.ind=TRUE)
  
}



sqrlossbest_5F <- numeric(0)
yhat_5F <- matrix(NA,nrow=P,ncol=n)
for (j in 1:P){
  for (i in 1:5){
    fmi <- lm(modfor[[j]], data=baseball_shuff[-Indexes[[i]],])
    yhat_5F[j,Indexes[[i]]] <- predict(fmi,baseball_shuff[Indexes[[i]],])
  }
  sqrlossbest_5F[j] <- sqrt(mean((yhat_5F[j,]-baseball_shuff$y)^2))
}


#Exercise 3

inc=function(e1,e2) eval.parent(substitute(e1 <- e1 + e2))

####a)
dataset=list()
dataset2=list()
pb=txtProgressBar(min=0,max=100,style=3)
for (i in 1:1000){
  uniform=list()
  uniform2=list()
  beta=list()
  
  for (u in 1:20){
    uniform[[u]]=runif(50,0,1)
    uniform2[[u]]=runif(50,0,1)
    
  }
  
  
  for (j in 1:21){
    beta[[j]]=rnorm(1,0,0.5)
  }
  
  y_e=rnorm(50,0,1)
  y_e2=rnorm(50,0,1)
  y=c()
  y2=c()
  for (n in 1:50){
    y[n]=0
    y2[n]=0
    for (s in 1:20){
      
      inc(y[n],beta[[s]]*uniform[[s]][n])
      inc(y2[n],beta[[s]]*uniform2[[s]][n])
      
    }
    
    inc(y[n] ,beta[[21]])
    inc(y[n],y_e[n])
    inc(y2[n] ,beta[[21]])
    inc(y2[n],y_e2[n])
    
    
  }
  
  data=c()
  data2=c()
  for(c in 1:20){
    data=cbind(uniform[[c]],data)
    data2=cbind(uniform2[[c]],data2)
  }
  
  dataset[[i]]=cbind(data,y,y_e)
  dataset[[i]]=as.data.frame(dataset[[i]])
  dataset2[[i]]=cbind(data2,y2,y_e2)
  dataset2[[i]]=as.data.frame(dataset2[[i]])
  

  setTxtProgressBar(pb,i)
}

################################################b

#########i)


sqlossfull_b_i=list()
pb=txtProgressBar(min=0,max=100,style=3)
for(i in 1:100){
  
  yhat <- numeric(0)
  for (j in 1:50){
    fmi <- lm(y~1, data=dataset[[i]][-j,])
    yhat[j] <- predict(fmi,dataset[[i]][j,])
  }
  
  sqlossfull_b_i[[i]] <- sqrt(mean((yhat-dataset[[i]][,"y"])^2))
  
  setTxtProgressBar(pb,i)
  
}

#####ii)
sqlossfull_b_ii=list()
pb=txtProgressBar(min=0,max=100,style=3)

for(i in 1:100){
  
  yhat <- numeric(0)
  for (j in 1:50){
    fmi <- lm(y~., data=dataset[[i]][-j,])
    yhat[j] <- predict(fmi,dataset[[i]][j,])
  }
  sqlossfull_b_ii[[i]] <- sqrt(mean((yhat-dataset[[i]][,"y"])^2))
  setTxtProgressBar(pb,i)
  
}

#iii)
sqlossfull_b_iii=list()
pb=txtProgressBar(min=0,max=100,style=3)

for(d in 1:100){
  setTxtProgressBar(pb,i)
  response="y"
  varnames= names(dataset[[d]])[c(-21,-22)]
  full.model <-as.formula(paste(response,"~",paste(varnames,collapse=" + ")))
  P=20
  bestsub <- regsubsets(full.model,data=dataset[[d]],nvmax=20,method="exhaustive")
  sbest <- summary(bestsub)
  
  sqrlossbest_LOO <- numeric(0)
  yhat_LOO <- matrix(NA,nrow=P,ncol=n)
  modfor <- list()
  
  for (j in 1:P){
    modfor[[j]] <-as.formula(paste(response,"~",paste(varnames[sbest$which[j,2:(P+1)]],
                                                      collapse=" + ")))
  }
  for (j in 1:P){
    for (i in 1:n){
      fmi <- lm(modfor[[j]], data=dataset[[d]][-i,])
      yhat_LOO[j,i] <- predict(fmi,dataset[[d]][i,])
    }
    sqrlossbest_LOO[j] <- sqrt(mean((yhat_LOO[j,]-dataset[[d]][,"y"])^2))
  }
  
  sqlossfull_b_iii[[d]]=min(sqrlossbest_LOO)
  setTxtProgressBar(pb,i)
  


}

######################################### DATASET 2###########################

#########i)

pb=txtProgressBar(min=0,max=100,style=3)
sqlossfull_b_i_MOD=list()
for(i in 1:100){
  
  yhat <- numeric(0)
  
    fmi <- lm(y~1, data=dataset[[i]])
    yhat <- predict(fmi,dataset2[[i]])
  
  sqlossfull_b_i_MOD[[i]] <- sqrt(mean((yhat-dataset2[[i]][,"y2"])^2))
}

#####ii)
sqlossfull_b_ii_MOD=list()
for(i in 1:100){
  
  yhat <- numeric(0)
  
    fmi <- lm(y~., data=dataset[[i]])
    yhat <- predict(fmi,dataset2[[i]])
  
  sqlossfull_b_ii_MOD[[i]] <- sqrt(mean((yhat-dataset2[[i]][,"y2"])^2))
}

#iii)
pb=txtProgressBar(min=0,max=100,style=3)
sqlossfull_b_iii_MOD=list()
for(d in 1:100){
  response="y"
  varnames= names(dataset[[d]])[c(-21,-22)]
  full.model <-as.formula(paste(response,"~",paste(varnames,collapse=" + ")))
  P=20
  bestsub <- regsubsets(full.model,data=dataset[[d]],nvmax=20,method="exhaustive")
  sbest <- summary(bestsub)
  
  sqrlossbest_LOO <- numeric(0)
  yhat_LOO <- matrix(NA,nrow=P,ncol=n)
  modfor <- list()
  
  for (j in 1:P){
    modfor[[j]] <-as.formula(paste(response,"~",paste(varnames[sbest$which[j,2:(P+1)]],
                                                      collapse=" + ")))
  }
  for (j in 1:P){
    
      fmi <- lm(modfor[[j]], data=dataset[[d]])
      yhat_LOO[j,] <- predict(fmi,dataset2[[d]])
    
    sqrlossbest_LOO[j] <- sqrt(mean((yhat_LOO[j,]-dataset2[[d]][,"y2"])^2))
  }
  
  sqlossfull_b_iii_MOD[[d]]=min(sqrlossbest_LOO)
  setTxtProgressBar(pb,i)
  
  
}

##############################################a


#########i)


sqlossfull_a_i=list()
for(i in 1:100){
  
  yhat <- numeric(0)
  for (j in 1:50){
    fmi <- lm(y_e~1, data=dataset[[i]][-j,])
    yhat[j] <- predict(fmi,dataset[[i]][j,])
  }
  sqlossfull_a_i[[i]] <- sqrt(mean((yhat-dataset[[i]][,"y_e"])^2))
}

#####ii)
sqlossfull_a_ii=list()
for(i in 1:100){
  
  y_ehat <- numeric(0)
  for (j in 1:50){
    fmi <- lm(y_e~., data=dataset[[i]][-j,])
    y_ehat[j] <- predict(fmi,dataset[[i]][j,])
  }
  sqlossfull_a_ii[[i]] <- sqrt(mean((y_ehat-dataset[[i]][,"y_e"])^2))
}

#iii)
sqlossfull_a_iii=list()
for(d in 1:100){
  response="y_e"
  varnames= names(dataset[[d]])[c(-21,-22)]
  full.model <-as.formula(paste(response,"~",paste(varnames,collapse=" + ")))
  P=20
  bestsub <- regsubsets(full.model,data=dataset[[d]],nvmax=20,method="exhaustive")
  sbest <- summary(bestsub)
  
  sqrlossbest_LOO <- numeric(0)
  y_ehat_LOO <- matrix(NA,nrow=P,ncol=n)
  modfor <- list()
  
  for (j in 1:P){
    modfor[[j]] <-as.formula(paste(response,"~",paste(varnames[sbest$which[j,2:(P+1)]],
                                                      collapse=" + ")))
  }
  for (j in 1:P){
    for (i in 1:n){
      fmi <- lm(modfor[[j]], data=dataset[[d]][-i,])
      y_ehat_LOO[j,i] <- predict(fmi,dataset[[d]][i,])
    }
    sqrlossbest_LOO[j] <- sqrt(mean((y_ehat_LOO[j,]-dataset[[d]][,"y_e"])^2))
  }
  
  sqlossfull_a_iii[[d]]=min(sqrlossbest_LOO)
  
  
}

######################################### DATASET 2###########################

#########i)


sqlossfull_a_i_MOD=c()
for(i in 1:100){
  
  y_ehat <- numeric(0)
  
  fmi <- lm(y_e~1, data=dataset[[i]])
  y_ehat <- predict(fmi,dataset2[[i]])
  
  sqlossfull_a_i_MOD[i] <- sqrt(mean((y_ehat-dataset2[[i]][,"y_e2"])^2))
}

#####ii)
sqlossfull_a_ii_MOD=c()
for(i in 1:100){
  
  y_ehat <- numeric(0)
  
  fmi <- lm(y_e~., data=dataset[[i]])
  y_ehat <- predict(fmi,dataset2[[i]])
  
  sqlossfull_a_ii_MOD[i] <- sqrt(mean((y_ehat-dataset2[[i]][,"y_e2"])^2))
}

#iii)
sqlossfull_a_iii_MOD=c()
for(d in 1:100){
  response="y_e"
  varnames= names(dataset[[d]])[c(-21,-22)]
  full.model <-as.formula(paste(response,"~",paste(varnames,collapse=" + ")))
  P=20
  bestsub <- regsubsets(full.model,data=dataset[[d]],nvmax=20,method="exhaustive")
  sbest <- summary(bestsub)
  
  sqrlossbest_LOO <- numeric(0)
  y_ehat_LOO <- matrix(NA,nrow=P,ncol=n)
  modfor <- list()
  
  for (j in 1:P){
    modfor[[j]] <-as.formula(paste(response,"~",paste(varnames[sbest$which[j,2:(P+1)]],
                                                      collapse=" + ")))
  }
  for (j in 1:P){
    
    fmi <- lm(modfor[[j]], data=dataset[[d]])
    y_ehat_LOO[j,] <- predict(fmi,dataset2[[d]])
    
    sqrlossbest_LOO[j] <- sqrt(mean((y_ehat_LOO[j,]-dataset2[[d]][,"y_e2"])^2))
  }
  
  sqlossfull_a_iii_MOD[d]=min(sqrlossbest_LOO)
  
  
}

sqlossfull_a_i=unlist(sqlossfull_a_i, use.names=FALSE)
sqlossfull_a_ii=unlist(sqlossfull_a_ii, use.names=FALSE)
sqlossfull_a_iii=unlist(sqlossfull_a_iii, use.names=FALSE)
sqlossfull_b_i=unlist(sqlossfull_b_i, use.names=FALSE)
sqlossfull_b_i_MOD=unlist(sqlossfull_b_i_MOD, use.names=FALSE)
sqlossfull_b_ii=unlist(sqlossfull_b_ii, use.names=FALSE)
sqlossfull_b_ii_MOD=unlist(sqlossfull_b_ii_MOD, use.names=FALSE)
sqlossfull_b_iii=unlist(sqlossfull_b_iii)
sqlossfull_b_iii_MOD=unlist(sqlossfull_b_iii_MOD, use.names=FALSE)

sqlossfull_a_i=unlist(sqlossfull_a_i, use.names=FALSE)
sqlossfull_a_ii=unlist(sqlossfull_a_ii, use.names=FALSE)
sqlossfull_a_iii=unlist(sqlossfull_a_iii, use.names=FALSE)
sqlossfull_b_i=unlist(sqlossfull_b_i, use.names=FALSE)
sqlossfull_b_i_MOD=unlist(sqlossfull_b_i_MOD, use.names=FALSE)
sqlossfull_b_ii=unlist(sqlossfull_b_ii, use.names=FALSE)
sqlossfull_b_ii_MOD=unlist(sqlossfull_b_ii_MOD, use.names=FALSE)
sqlossfull_b_iii=unlist(sqlossfull_b_iii)
sqlossfull_b_iii_MOD=unlist(sqlossfull_b_iii_MOD, use.names=FALSE)

matrix=cbind(sqlossfull_a_i,sqlossfull_a_ii,sqlossfull_a_iii,sqlossfull_b_i,
             sqlossfull_b_ii,sqlossfull_b_iii,sqlossfull_b_i_MOD,
             sqlossfull_b_ii_MOD,sqlossfull_b_iii_MOD,sqlossfull_a_i_MOD,
             sqlossfull_a_ii_MOD,sqlossfull_a_iii_MOD)


require(reshape2)
dfcomp<- melt(matrix)

require(ggplot2)
g=ggplot(data = dfcomp, aes(x=Var2, y=value))
g=g +scale_x_discrete(labels=c("sqlossfull_a_i" = "a_i", "sqlossfull_a_ii" = "a_ii",
                               "sqlossfull_a_iii" = "a_iii",
                               "sqlossfull_b_i"="b_i","sqlossfull_b_ii"="b_ii","sqlossfull_b_iii"="b_iii"
                               ,"sqlossfull_b_i_MOD"="b_i_MOD","sqlossfull_b_ii_MOD"="b_ii_MOD",
                               "sqlossfull_b_iii_MOD"="b_iii_MOD","sqlossfull_a_i_MOD"="a_i_MOD",
                               "sqlossfull_a_ii_MOD"="a_ii_MOD","sqlossfull_a_iii_MOD"="a_iii_MOD"))

g=g + theme(axis.text.x = element_text(angle = 90))+ labs(y= "sqlossfunction", x = "")
g=g + geom_boxplot(aes(fill=Var1))
g
