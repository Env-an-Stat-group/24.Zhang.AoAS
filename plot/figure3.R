prop.square=c(1:225)
for(i in 1:20){
  prop.square[sig[i]]=sum(indi[sig[i],])/1000
}
for(j in 1:205){
  prop.square[non[j]]=sum(indi[non[j],])/1000
}
(1-mean(prop.square[sig])+mean(prop.square[non]))/2

mat=matrix(prop.square,ncol=15,nrow=15)

cutpts <- c(0,seq(0.9,1,by=0.02))
levelplot(mat,at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu"))),main="Proportion of 95% Credibility Interval Contains 0")


for(i in 1:225){
  indi[i,]<-ifelse(indi[i,]>0.05,1,0)
  
}
prop.square=c(1:225)
for(i in 1:20){
  prop.square[sig[i]]=sum(indi[sig[i],])/1000
}
for(j in 1:205){
  prop.square[non[j]]=sum(indi[non[j],])/1000
}
(1-mean(prop.square[sig])+mean(prop.square[non]))/2
mat=matrix(prop.square,ncol=15,nrow=15)

cutpts <- c(0,seq(0.9,1,by=0.02))
levelplot(mat,at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu"))),main="Proportion of 95% Confidence Interval Contains 0")
