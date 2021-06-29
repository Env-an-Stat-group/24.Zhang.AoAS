library(sp)
library(spdep)
library(Matrix)
library(fda)
library(INLA)
library(lattice)
library(RColorBrewer)
library(plotrix)

############################################################################
############################################################################
######################  Square  ############################################
############################################################################

####################################################
############## Independent    ######################
####################################################

lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(81:85,96:100,111:115,126:130,141:145)
non=c(1:80,86:95,101:110,116:125,131:140,146:225)


beta=rep(0,times=225)
beta[sig]=rnorm(20,mean=10,sd=1)
beta=rep(beta,times=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  coef=c()
  pval=c()
  for(j in 1:225){
    yt=y[loc%in%j]
    xt=cov[loc%in%j]
    mod=lm(yt~xt)
    coef=c(coef,summary(mod)$coefficient[2,1])
    pval=c(pval,summary(mod)$coefficient[2,4])
  }
  betamat[,i]=coef
  indi[,i]=pval
}  


####################################################
##############  Stationary    ######################
####################################################
library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(81:85,96:100,111:115,126:130,141:145)
non=c(1:80,86:95,101:110,116:125,131:140,146:225)



beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,alpha=2)

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}


####################################################
##############  Non-Stationary    ##################
####################################################


library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(81:85,96:100,111:115,126:130,141:145)
non=c(1:80,86:95,101:110,116:125,131:140,146:225)

beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                         theta.prior.mean=c(0, 0, 0),
                         theta.prior.prec=c(0.1, 0.1, 0.1))

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}


############################################################################
############################################################################
###################### zigzag   ############################################
############################################################################


####################################################
############## Independent    ######################
####################################################

lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(65:68,81:84,97:100,113:116,129:132)
non=c(1:64,69:80,85:96,101:112,117:128,133:225)


beta=rep(0,times=225)
beta[sig]=rnorm(20,mean=10,sd=1)
beta=rep(beta,times=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  coef=c()
  pval=c()
  for(j in 1:225){
    yt=y[loc%in%j]
    xt=cov[loc%in%j]
    mod=lm(yt~xt)
    coef=c(coef,summary(mod)$coefficient[2,1])
    pval=c(pval,summary(mod)$coefficient[2,4])
  }
  betamat[,i]=coef
  indi[,i]=pval
}  


####################################################
##############  Stationary    ######################
####################################################
library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(65:68,81:84,97:100,113:116,129:132)
non=c(1:64,69:80,85:96,101:112,117:128,133:225)



beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,alpha=2)

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}


####################################################
##############  Non-Stationary    ##################
####################################################


library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(65:68,81:84,97:100,113:116,129:132)
non=c(1:64,69:80,85:96,101:112,117:128,133:225)

beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                         theta.prior.mean=c(0, 0, 0),
                         theta.prior.prec=c(0.1, 0.1, 0.1))

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}



############################################################################
############################################################################
###################### bar      ############################################
############################################################################


####################################################
############## Independent    ######################
####################################################

lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(80,82,84,86,95,97,99,101,110:116,125,127,129,131,140,142,144,146)
non=c(1:79,81,83,84,87:94,96,98,100,102:109,117:124,126,128,130,132:139,141,143,145,147:225)


beta=rep(0,times=225)
beta[sig]=rnorm(20,mean=10,sd=1)
beta=rep(beta,times=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  coef=c()
  pval=c()
  for(j in 1:225){
    yt=y[loc%in%j]
    xt=cov[loc%in%j]
    mod=lm(yt~xt)
    coef=c(coef,summary(mod)$coefficient[2,1])
    pval=c(pval,summary(mod)$coefficient[2,4])
  }
  betamat[,i]=coef
  indi[,i]=pval
}  


####################################################
##############  Stationary    ######################
####################################################
library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(80,82,84,86,95,97,99,101,110:116,125,127,129,131,140,142,144,146)
non=c(1:79,81,83,84,87:94,96,98,100,102:109,117:124,126,128,130,132:139,141,143,145,147:225)


beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,alpha=2)

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}


####################################################
##############  Non-Stationary    ##################
####################################################


library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(80,82,84,86,95,97,99,101,110:116,125,127,129,131,140,142,144,146)
non=c(1:79,81,83,84,87:94,96,98,100,102:109,117:124,126,128,130,132:139,141,143,145,147:225)

beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                         theta.prior.mean=c(0, 0, 0),
                         theta.prior.prec=c(0.1, 0.1, 0.1))

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}


############################################################################
############################################################################
###################### U shape  ############################################
############################################################################


####################################################
############## Independent    ######################
####################################################

lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(66:71,81:86,100:101,115:116,130:131,141:146,156:161)
non=c(1:65,72:80,87:99,102:114,117:129,132:140,147:155,162:225)


beta=rep(0,times=225)
beta[sig]=rnorm(20,mean=10,sd=1)
beta=rep(beta,times=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  coef=c()
  pval=c()
  for(j in 1:225){
    yt=y[loc%in%j]
    xt=cov[loc%in%j]
    mod=lm(yt~xt)
    coef=c(coef,summary(mod)$coefficient[2,1])
    pval=c(pval,summary(mod)$coefficient[2,4])
  }
  betamat[,i]=coef
  indi[,i]=pval
}  


####################################################
##############  Stationary    ######################
####################################################
library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(66:71,81:86,100:101,115:116,130:131,141:146,156:161)
non=c(1:65,72:80,87:99,102:114,117:129,132:140,147:155,162:225)


beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,alpha=2)

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}


####################################################
##############  Non-Stationary    ##################
####################################################


library(INLA)


lon=seq(1,15,by=1)
lat=seq(1,15,by=1)

grid=expand.grid(lon,lat)

coord=cbind(grid$Var1,grid$Var2)

coordrep=rbind(coord,coord,coord,coord)


cov=rep(c(1,0,1,0),each=225)

repl=rep(c(1,2,3,4),each=225)
loc=rep(c(1:225),times=4)

sig=c(66:71,81:86,100:101,115:116,130:131,141:146,156:161)
non=c(1:65,72:80,87:99,102:114,117:129,132:140,147:155,162:225)

beta=rep(0,times=225)
beta[sig]=rnorm(25,mean=10,sd=1)
beta=rep(beta,times=4)


mesh1<-inla.mesh.create(loc=coordrep)

SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                         theta.prior.mean=c(0, 0, 0),
                         theta.prior.prec=c(0.1, 0.1, 0.1))

A = inla.spde.make.A(mesh = SPDE, loc = coordrep,repl =repl , n.repl=4,weights = cov)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,n.repl=4)

betamat<-matrix(0,ncol=1000,nrow=225)
indi<-matrix(0,ncol=1000,nrow=225)
for(i in 1:1000){
  
  
  y=10+beta*cov+rnorm(900,mean=0,sd=1)
  
  
  stk.b1 = inla.stack(data = list(y = y), A = list(A,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(y)))),tag="est1")
  
  formula = y~-1+b1.intercept+f(b1.field, model = SPDE,replicate = b1.field.repl)
  
  result=inla(formula,family="normal",data=inla.stack.data(stk.b1),
              control.predictor = list(A=inla.stack.A(stk.b1),compute=TRUE),control.compute = list(dic=TRUE),verbose =F)
  
  betares=result$summary.random$b1.field$mean
  
  res<-c()
  for(l in 1:241){
    res[l]<-mean(betares[seq(l,723,by=241)])
  }
  beta025=result$summary.random$b1.field$'0.025quant'
  beta975=result$summary.random$b1.field$'0.975quant'
  b025<-c()
  for(p in 1:241){
    b025[p]<-mean(beta025[seq(p,723,by=241)])
  }
  
  b975<-c()
  for(q in 1:241){
    b975[q]<-mean(beta975[seq(q,723,by=241)])
  }
  
  tt=c()
  for(j in 1:225){
    if(beta025[j]<=0&beta975[j]>0){
      tt[j]<-1
    }
    else{tt[j]<-0}
  }
  
  betamat[,i]=res[1:225]
  indi[,i]=tt
  
}
