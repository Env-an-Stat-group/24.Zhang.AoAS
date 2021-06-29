library(forecast)
library(FNN)
library(gaussDiff)
library(GoFKernel)
library(geoR)
library(ggplot2)


library(abind)

run1 = readRDS("monthly-ws10-d01-run1-2013.RDS")

run1 = abind(run1, readRDS("monthly-ws10-d01-run1-2014.RDS"), along=3)

run1 = abind(run1, readRDS("monthly-ws10-d01-run1-2015.RDS"), along=3)
run1 = abind(run1, readRDS("monthly-ws10-d01-run1-2016.RDS"), along=3)

lat1 = readRDS("LAT_GRID1.RDS")

lon1 = readRDS("LON_GRID1.RDS")

mask = readRDS("saudi_mask-run1.RDS")

#############################################################################################


run2 = readRDS("monthly-ws10-d01-run2-2013.RDS")

run2 = abind(run2, readRDS("monthly-ws10-d01-run2-2014.RDS"), along=3)

run2 = abind(run2, readRDS("monthly-ws10-d01-run2-2015.RDS"), along=3)

run2 = abind(run2, readRDS("monthly-ws10-d01-run2-2016.RDS"), along=3)

#############################################################################################


run3 = readRDS("monthly-ws10-d01-run3-2013.RDS")

run3 = abind(run3, readRDS("monthly-ws10-d01-run3-2014.RDS"), along=3)

run3 = abind(run3, readRDS("monthly-ws10-d01-run3-2015.RDS"), along=3)

run3 = abind(run3, readRDS("monthly-ws10-d01-run3-2016.RDS"), along=3)


lat3 = readRDS("LAT_GRID.RDS")

lon3 = readRDS("LON_GRID.RDS")

#############################################################################################


run4 = readRDS("monthly-ws10-d01-run4-2013.RDS")

run4 = abind(run4, readRDS("monthly-ws10-d01-run4-2014.RDS"), along=3)

run4 = abind(run4, readRDS("monthly-ws10-d01-run4-2015.RDS"), along=3)

run4 = abind(run4, readRDS("monthly-ws10-d01-run4-2016.RDS"), along=3)

#############################################################################################


source("interp_fncs.R") 



dat3 = array(0,dim=c(339,299,48))

for(i in 1:339){
  for(j in 1:299){
    ans = lamb_latlon_to_ij(lat3[1,1], lon3[1,1], 1, 1, 24, 24,47, 6000, lat1[i,j], lon1[i,j], radius= 6370000.0)
    dat3[i,j,] = interp2d_full(run3,ans$i,ans$j,interp_type=2)
  }
}

dat4 = array(0,dim=c(339,299,48))

for(i in 1:339){
  for(j in 1:299){
    ans = lamb_latlon_to_ij(lat3[1,1], lon3[1,1], 1, 1, 24, 24,47, 6000, lat1[i,j], lon1[i,j], radius= 6370000.0)
    dat4[i,j,] = interp2d_full(run4,ans$i,ans$j,interp_type=2)
  }
}

map_mask = function(x){
  for(i in 1:48){
    x[,,i]=x[,,i]*mask 
  }
  return(x)
}


############ make the points outside of Saudi Arabia equals to zero

run1 = map_mask(run1)
run2 = map_mask(run2)
dat3 = map_mask(dat3)
dat4 = map_mask(dat4)

############ 23708 means all the points locate within Saudi Arabia
############ make the dataframe that used in INLA


wind = c(run1[run1!=0],run2[run2!=0],dat3[dat3!=0],dat4[dat4!=0])
time = rep(rep(1:48,each=23708),4)
loc = rep(rep(1:23708,48),4)
reso = rep(c(1,-1),each=2*23708*48)
init = c(rep(c(1,-1),each=23708*48), rev(rep(c(1,-1),each=23708*48)))

############ harmonic term at order 3

h1_sin = sin(2*pi*(0:11)/12)
h1_cos = cos(2*pi*(0:11)/12)

h1_sin_X = rep(rep(h1_sin,each = 23708),16)
h1_cos_X = rep(rep(h1_cos,each = 23708),16)

h2_sin=sin(2*pi*2*(0:11)/12)
h2_cos=cos(2*pi*2*(0:11)/12)

h2_sin_X = rep(rep(h2_sin,each = 23708),16)
h2_cos_X = rep(rep(h2_cos,each = 23708),16)

h3_sin=sin(2*pi*3*(0:11)/12)
h3_cos=cos(2*pi*3*(0:11)/12)

h3_sin_X = rep(rep(h3_sin,each = 23708),16)
h3_cos_X = rep(rep(h3_cos,each = 23708),16)

harmonic = data.frame(s1 = h1_sin_X, c1 = h1_cos_X, s2 = h2_sin_X, c2 = h2_cos_X, s3 = h3_sin_X, c3 = h3_cos_X)

wind_dat = data.frame(wind = wind, location = loc, reso = reso, init = init,time = time)

wind_dat = cbind(wind_dat,harmonic)

rm(run1,run2,run3,run4,dat3,dat4,wind,time,loc,reso,init,harmonic,h1_sin_X,h1_cos_X,h2_sin_X,h2_cos_X,h3_sin_X,h3_cos_X)



##########################################################################
#################    No upscalling   ######################################
##########################################################################

map_mask2 = function(x){
  for(i in 1:48){
    x[,,i]=x[,,i]*mask2
  }
  return(x)
}

run1 = map_mask(run1)
run2 = map_mask(run2)
run3 = map_mask2(run3)
run4 = map_mask2(run4)


wind = c(run1[run1!=0],run2[run2!=0],run3[run3!=0],run4[run4!=0])
time = c(rep(rep(1:48,each=23708),2),rep(rep(1:48,each=53362),2))
loc = c(rep(rep(1:23708,48),2),rep(rep(1:53362,48),2))
reso = c(rep(1,times=2*23708*48),rep(0,times=2*53362*48))
init = c(rep(c(1,0),each=23708*48), rep(c(0,1),each=53362*48))

h1_sin = sin(2*pi*(0:11)/12)
h1_cos = cos(2*pi*(0:11)/12)

h1_sin_X = c(rep(rep(h1_sin,each = 23708),8),rep(rep(h1_sin,each = 53362),8))
h1_cos_X = c(rep(rep(h1_cos,each = 23708),8),rep(rep(h1_cos,each = 53362),8))

h2_sin=sin(2*pi*2*(0:11)/12)
h2_cos=cos(2*pi*2*(0:11)/12)

h2_sin_X = c(rep(rep(h2_sin,each = 23708),8),rep(rep(h1_sin,each = 53362),8))
h2_cos_X = c(rep(rep(h2_cos,each = 23708),8),rep(rep(h1_cos,each = 53362),8))

h3_sin=sin(2*pi*3*(0:11)/12)
h3_cos=cos(2*pi*3*(0:11)/12)

h3_sin_X = c(rep(rep(h3_sin,each = 23708),8),rep(rep(h1_sin,each = 53362),8))
h3_cos_X = c(rep(rep(h3_cos,each = 23708),8),rep(rep(h1_cos,each = 53362),8))

harmonic = data.frame(s1 = h1_sin_X, c1 = h1_cos_X, s2 = h2_sin_X, c2 = h2_cos_X, s3 = h3_sin_X, c3 = h3_cos_X)

wind_dat = data.frame(wind = wind, location = loc, reso = reso, init = init,time = time)

wind_dat = cbind(wind_dat,harmonic)

lat1 = lat1*mask
lon1 = lon1*mask


coord = as.matrix(data.frame(lon = lon1[lon1!=0], lat = lat1[lat1!=0]))

lat3 = lat3*mask2
lon3 = lon3*mask2


coord2 = as.matrix(data.frame(lon = lon3[lon3!=0], lat = lat3[lat3!=0]))
