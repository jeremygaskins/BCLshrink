library(devtools)
source_url("https://raw.githubusercontent.com/jeremygaskins/BCLshrink/main/BCLshrink%20code%20source.R")
library(pgdraw)
library(pscl)
library(coda)
library(GIGrvg)



load('mnlsims1v1v2_df1.RData')
# File is available in the BCLshrink repository on Github
y <- mnlogit_data[,1]
table(y)


X <- as.matrix(mnlogit_data[,-1])
X[1:10,1:10]


## Fitting Horseshoe shared shrinkage (HS SS) model:

fit.HSSS <- BLCshrink.HSSS(y=y, X = X,
	N_it=15000, burn.in=5000, samples=TRUE, progress=FALSE)

round(fit.HSSS$kappa.hat,3)[,1:11]

round(fit.HSSS$kappa.LCI,3)[,1:11]
round(fit.HSSS$kappa.UCI,3)[,1:11]

fit.HSSS$kappa.sig[,1:11]*1

round(apply(fit.HSSS$samples, 2:3, coda::effectiveSize),1)[,1:11]

## Trace plots of the coefficients 
par(mfrow=c(4,2))
for( p in 1:2){
   for( j in 1:4){
      plot( fit.HSSS$samples[,p,j], type='l', xlab='Stored Iteration',
      main=paste0('Trace plot for ',colnames(fit.HSSS$kappa.hat)[j],
	  ' coefficient in class ',p),
      ylab='kappa')
   }
}



## Other shrinkage models

## Normal Gamma shared shrinkage (NG SS) model
fit.NGSS <- BLCshrink.NGSS(y=y, X = X,
	N_it=15000, burn.in=5000, samples=TRUE, progress=FALSE)
round(fit.NGSS$kappa.hat,3)[,1:11]


## Horseshoe unique shrinkage (HS US) model
fit.HSUS <- BLCshrink.HSUS(y=y, X = X,
	N_it=15000, burn.in=5000, samples=TRUE, progress=FALSE) 
round(fit.HSUS$kappa.hat,3)[,1:11]


## Normal Gamma unique shrinkage (NG US) model
fit.NGUS <- BLCshrink.NGUS(y=y, X = X,
	N_it=15000, burn.in=5000, samples=TRUE, progress=FALSE) 
round(fit.NGUS$kappa.hat,3)[,1:11]


## No shrinkage (NS) method
fit.none <- BLCshrink.none(y=y, X = X,
	N_it=15000, burn.in=5000, samples=TRUE, progress=FALSE) 
round(fit.none$kappa.hat,3)[,1:11]
