#figure 1
mon.mean1<-c()
for(i in 1:48){
  a=ifelse(run1[,,i]==0,NA,run1[,,i])
  mon.mean1[i]=mean(a,na.rm = T)
}

mon.mean2<-c()
for(i in 1:48){
  a=ifelse(run2[,,i]==0,NA,run2[,,i])
  mon.mean2[i]=mean(a,na.rm = T)
}

mon.mean3<-c()
for(i in 1:48){
  a=ifelse(run3[,,i]==0,NA,run3[,,i])
  mon.mean3[i]=mean(a,na.rm = T)
}

mon.mean4<-c()
for(i in 1:48){
  a=ifelse(run4[,,i]==0,NA,run4[,,i])
  mon.mean4[i]=mean(a,na.rm = T)
}

summer=c(6:8,18:20,30:32,42:44)
loc.mean1=matrix(NA,nrow=339,ncol=299)
for(i in 1:339){
  for(j in 1:299){
    loc.mean1[i,j]=mean(run1[i,j,summer])
  }
}

loc.mean1=ifelse(loc.mean1==0,NA,loc.mean1)
cutpts=seq(2,10,by=1)
levelplot(loc.mean1,  at=cutpts,pretty=T, 
          col.regions=(brewer.pal(9,"Reds")),main="PBL")

range(loc.mean1,na.rm = T)


loc.mean2=matrix(NA,nrow=339,ncol=299)
for(i in 1:339){
  for(j in 1:299){
    loc.mean2[i,j]=mean(run2[i,j,summer])
  }
}

loc.mean2=ifelse(loc.mean2==0,NA,loc.mean2)
cutpts=seq(2,10,by=1)
levelplot(loc.mean2,  at=cutpts,pretty=T, 
          col.regions=(brewer.pal(9,"Reds")),main="PBL")


loc.mean3=matrix(NA,nrow=549,ncol=499)
for(i in 1:549){
  for(j in 1:499){
    loc.mean3[i,j]=mean(run3[i,j,summer])
  }
}

loc.mean3=ifelse(loc.mean3==0,NA,loc.mean3)
cutpts=seq(2,10,by=1)
levelplot(loc.mean3,  at=cutpts,pretty=T, 
          col.regions=(brewer.pal(9,"Reds")),main="PBL")


loc.mean4=matrix(NA,nrow=549,ncol=499)
for(i in 1:549){
  for(j in 1:499){
    loc.mean4[i,j]=mean(run4[i,j,summer])
  }
}

loc.mean4=ifelse(loc.mean4==0,NA,loc.mean4)
cutpts=seq(2,10,by=1)
levelplot(loc.mean4,  at=cutpts,pretty=T, 
          col.regions=(brewer.pal(9,"Reds")),main="PBL")
