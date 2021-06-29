diff=c(1:225)
diff[1]=(abs(betamean[1]-betamean[16])+abs(betamean[1]-betamean[2]))/2
diff[15]=(abs(betamean[15]-betamean[14])+abs(betamean[15]-betamean[30]))/2
diff[211]=(abs(betamean[211]-betamean[212])+abs(betamean[211]-betamean[196]))/2
diff[225]=(abs(betamean[225]-betamean[224])+abs(betamean[225]-betamean[210]))/2


for(i in c(2:14)){
  diff[i]=(abs(betamean[i]-betamean[i-1])+abs(betamean[i]-betamean[i+1])+abs(betamean[i]-betamean[i+15]))/3
}

for(i in c(212:224)){
  diff[i]=(abs(betamean[i]-betamean[i-1])+abs(betamean[i]-betamean[i+1])+abs(betamean[i]-betamean[i-15]))/3
}

for(i in seq(16,196,by=15)){
  diff[i]=(abs(betamean[i]-betamean[i-15])+abs(betamean[i]-betamean[i+15])+abs(betamean[i]-betamean[i+1]))/3
}

for(i in seq(30,210,by=15)){
  diff[i]=(abs(betamean[i]-betamean[i-15])+abs(betamean[i]-betamean[i+15])+abs(betamean[i]-betamean[i-1]))/3
}

rest=c(1:225)[-c(1:15,211:225,seq(16,196,by=15),seq(30,210,by=15))]
for(i in 1:181){
  diff[rest[i]]=(abs(betamean[rest[i]]-betamean[rest[i]+1])+abs(betamean[rest[i]]-betamean[rest[i]-1])+abs(betamean[rest[i]]-betamean[rest[i]+15])+abs(betamean[rest[i]]-betamean[rest[i]-15]))/4
}
mat=matrix(diff,ncol=15,nrow=15)
cutpts <- seq(0,5,by=0.5)
levelplot(mat,at=cutpts,pretty=T, 
          col.regions=(rev(brewer.pal(11,"RdBu"))),main="Difference(Independent)")
