library(sp)
library(spdep)
library(Matrix)
library(fda)
library(INLA)


betapbl<-matrix(0,nrow=339,ncol=299)
betares<-matrix(0,nrow=339,ncol=299)
for(i in 1:339){
  for(j in 1:299){
    x<-c(run1[i,j,],run2[i,j,],dat3[i,j,],dat4[i,j,])
    ht<-c(rep(ht1[i,j],times=96),rep(ht2[i,j],times=96))
    dat<-data.frame(x,harmonic,pbl,res,ht)
    mod.h1=lm(x~s1+c1+s2+c2+s3+c3+pbl+res+ht,data=dat) 
    betapbl[i,j]<-summary(mod.h1)$coefficient[8,1]
    betares[i,j]<-summary(mod.h1)$coefficient[9,1]
  }
}


betapbl<-betapbl*mask
betapbl<-ifelse(betapbl==0,NA,betapbl)

cutpts <- seq(-1.4,0.8,by=0.2)
levelplot(betapbl,  at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu"))),main="PBL")

betares<-betares*mask
betares<-ifelse(betares==0,NA,betares)

cutpts <- seq(-1.2,1.2,by=0.3)
levelplot(betares,  at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(9,"Oranges"))),main="Resolution")


#######################################################################
##############  Stationary  ###########################################
#######################################################################

x1<-wind.new$reso
x2<-wind.new$init
ht<-wind.new$ht
Y<-as.numeric(wind.new$wind.rm)

wind.new$rep<-rep(c(1:4),each=71124)

coordrep<-rbind(coord,coord,coord)
coordrep<-rbind(coordrep,coordrep,coordrep,coordrep)

mesh1<-inla.mesh.create(loc=coord)


SPDE = (mesh1,alpha = 2,prior.range = c(3,NA),prior.sigma = c(1,NA))

A1 = inla.spde.make.A(mesh = SPDE, loc = coordrep, repl =wind.new$rep , n.repl=4,weights = x1,group=wind.new$time)
A2 = inla.spde.make.A(mesh = SPDE, loc = coordrep, repl =wind.new$rep , n.repl=4,weights = x2,group=wind.new$time)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,
                                     n.repl=4,n.group=3)
spatial.idx.2 = inla.spde.make.index("b2.field", n.spde=SPDE$n.spde,
                                     n.repl=4,n.group=3)
stk.b1 = inla.stack(data = list(y = Y), A = list(A1,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(Y)))),tag="est1")
stk.b2 = inla.stack(data = list(y = Y), A = list(A2,1), effects = list(spatial.idx.2,list(b2.intercept=rep(1,length(Y)))),tag="est2")
stk.data=inla.stack(stk.b1,stk.b2)
formula = y~-1+b1.intercept+f(b1.field, model = SPDE, replicate = b1.field.repl,control.group = list(model="ar1"))+
  b2.intercept+f(b2.field, model = SPDE, replicate = b2.field.repl,control.group = list(model="ar1"))

result=inla(formula,family="normal",data=inla.stack.data(stk.data),
            control.predictor = list(A=inla.stack.A(stk.data),compute=TRUE),control.compute = list(dic=TRUE),verbose =TRUE)




#######################################################################
##############  Non-Stationary ########################################
#######################################################################
x1<-wind.new$reso
x2<-wind.new$init
ht<-wind.new$ht
Y<-as.numeric(wind.new$wind.rm)

wind.new$rep<-rep(c(1:4),each=71124)

coordrep<-rbind(coord,coord,coord)
coordrep<-rbind(coordrep,coordrep,coordrep,coordrep)

mesh1<-inla.mesh.create(loc=coord)


SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                         theta.prior.mean=c(0, 0, 0),
                         theta.prior.prec=c(0.1, 0.1, 0.1))

A1 = inla.spde.make.A(mesh = SPDE, loc = coordrep, repl =wind.new$rep , n.repl=4,weights = x1,group=wind.new$time)
A2 = inla.spde.make.A(mesh = SPDE, loc = coordrep, repl =wind.new$rep , n.repl=4,weights = x2,group=wind.new$time)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,
                                     n.repl=4,n.group=3)
spatial.idx.2 = inla.spde.make.index("b2.field", n.spde=SPDE$n.spde,
                                     n.repl=4,n.group=3)
stk.b1 = inla.stack(data = list(y = Y), A = list(A1,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(Y)))),tag="est1")
stk.b2 = inla.stack(data = list(y = Y), A = list(A2,1), effects = list(spatial.idx.2,list(b2.intercept=rep(1,length(Y)))),tag="est2")
stk.data=inla.stack(stk.b1,stk.b2)
formula = y~-1+b1.intercept+f(b1.field, model = SPDE, replicate = b1.field.repl,control.group = list(model="ar1"))+
  b2.intercept+f(b2.field, model = SPDE, replicate = b2.field.repl,control.group = list(model="ar1"))

result=inla(formula,family="normal",data=inla.stack.data(stk.data),
            control.predictor = list(A=inla.stack.A(stk.data),compute=TRUE),control.compute = list(dic=TRUE),verbose =TRUE)



##############################################################
##############Bernoulli#######################################
##############################################################
x1<-wind.new$reso
x2<-wind.new$init
ht<-wind.new$ht
Y<-as.numeric(wind.new$wind.bin)

wind.new$rep<-c(rep(c(1:2),each=18162),rep(c(3:4),each=41071))

coordrep<-rbind(coord,coord,coord)
coordrep<-rbind(coordrep,coordrep,coordrep,coordrep)

mesh1<-inla.mesh.create(loc=coord)


SPDE = inla.spde2.matern(mesh1,B.tau=matrix(cbind(0, 1, 0, sin(pi*mesh1$loc[,1])),ncol=4),
                         B.kappa=matrix(c(0, 0, 1, 0),nrow=1,ncol=4),
                         theta.prior.mean=c(0, 0, 0),
                         theta.prior.prec=c(0.1, 0.1, 0.1))

A1 = inla.spde.make.A(mesh = SPDE, loc = coordrep, repl =wind.new$rep , n.repl=4,weights = x1,group=wind.new$time)
A2 = inla.spde.make.A(mesh = SPDE, loc = coordrep, repl =wind.new$rep , n.repl=4,weights = x2,group=wind.new$time)

spatial.idx.1 = inla.spde.make.index("b1.field", n.spde=SPDE$n.spde,
                                     n.repl=4,n.group=3)
spatial.idx.2 = inla.spde.make.index("b2.field", n.spde=SPDE$n.spde,
                                     n.repl=4,n.group=3)
stk.b1 = inla.stack(data = list(y = Y), A = list(A1,1), effects = list(spatial.idx.1,list(b1.intercept=rep(1,length(Y)))),tag="est1")
stk.b2 = inla.stack(data = list(y = Y), A = list(A2,1), effects = list(spatial.idx.2,list(b2.intercept=rep(1,length(Y)))),tag="est2")
stk.data=inla.stack(stk.b1,stk.b2)
formula = y~-1+b1.intercept+f(b1.field, model = SPDE, replicate = b1.field.repl,control.group = list(model="ar1"))+
  b2.intercept+f(b2.field, model = SPDE, replicate = b2.field.repl,control.group = list(model="ar1"))

result=inla(formula,family="binomial",data=inla.stack.data(stk.data),
            control.predictor = list(A=inla.stack.A(stk.data),compute=TRUE),control.compute = list(dic=TRUE),verbose =TRUE)

