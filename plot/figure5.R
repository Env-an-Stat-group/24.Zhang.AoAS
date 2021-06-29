betares<-result$summary.random$b1.field$mean
betapbl<-result$summary.random$b2.field$mean

res<-betares[1:53362]
pbl<-betapbl[1:53362]

res<-c()
for(i in 1:53362){
  res[i]<-mean(betares[seq(i,106756,by=53362)])
}

pbl<-c()
for(i in 1:53362){
  pbl[i]<-mean(betapbl[seq(i,213512,by=53362)])
}


i=1
tt<-c()
for(j in 1:273951){
  if(mask.row2[j]==0){
    tt[j]<-NA
  }
  else {
    tt[j]<-pbl[i]
    i=i+1
  }
}


mat<-matrix(tt,nrow=549,ncol=499)

cutpts <- seq(-1.4,0.8,by=0.2)
levelplot(mat,  at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu"))),main="PBL")


i=1
tt<-c()
for(j in 1:273951){
  if(mask.row2[j]==0){
    tt[j]<-NA
  }
  else {
    tt[j]<-res[i]
    i=i+1
  }
}
mat<-matrix(tt,nrow=549,ncol=499)

cutpts <- seq(-1.2,1.2,by=0.3)
levelplot(mat,  at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(9,"Oranges"))),main="Resolution")

