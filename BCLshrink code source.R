BLCshrink.NGSS <- function(y,X,
	N_it=15000, burn.in=round(.2*N_it), CI.level=0.95, progress=TRUE, 
	c0=0.01, c1=0.01, theta=0.05, samples=FALSE, MCMC.seed=209578){

set.seed(MCMC.seed)
x_mat <- cbind(1,X)
colnames(x_mat)[1] <- "int"
p_mnl <- ncol(x_mat)
m <- nrow(x_mat)
k <- max(y)
Zj <- t(sapply(y, function(s) 1*(s==1:k) ))-.5
prog.count <- seq(0,N_it,length=11)[-1]

kappa <- array(0, dim = c(k,p_mnl)) 
lam_k <- 1/p_mnl 
tausqr <- rep(1,p_mnl)
int.pos <- 1
coef.pos <- c(2:p_mnl)
tausqr[int.pos] <- 100

store.kappa <- array(NA, dim = c(N_it,k,p_mnl))  

   
for( it in 1:N_it){

for ( j in 1:(k-1)){
 	xb=(x_mat)%*%t(kappa[-j,]) 
	max_val <- apply(xb,1,max)
	Dj <- max_val + log(apply(exp(xb-max_val),1,sum))
	Bij <- x_mat%*% kappa[j,]-Dj
	Bij=ifelse(Bij>  1e20,1e20,Bij)
	Bij=ifelse(Bij< -1e20,-1e20,Bij)
	omg <- pgdraw::pgdraw(1,Bij)  
      var_mat <- solve(diag(1/((2*tausqr)/lam_k))+t(x_mat)%*%diag(omg)%*%x_mat)
      mu_vec <- var_mat%*%t(x_mat)%*%( Zj[,j]+diag(omg)%*%Dj)
      kappa[j,] <- MASS::mvrnorm(1,mu_vec,var_mat)        
}

lam_k <- rgamma(1, c0+length(coef.pos)*(k-1)*0.5,  
	c1+0.25* sum(t(t((kappa[,coef.pos]^2))/tausqr[coef.pos])) )
   
for ( p in coef.pos){
	tausqr[p] <- GIGrvg::rgig(1,
		(theta-0.5*(k-1) ),
		0.5*lam_k*sum(kappa[,p]^2),
		2*theta)
}
tausqr <- ifelse(tausqr<1e-10,1e-10, tausqr)
tausqr <- ifelse(tausqr>1e10,1e10, tausqr)

store.kappa[it,,] <- kappa
if( progress ){ if( it %in% prog.count ){
	temp <- which(it==prog.count)
	print(paste0('MCMC progress:  ',temp*10,'%'))
	flush.console()
}}
}


store.kappa <- store.kappa[-(1:burn.in),,]
kappa.hat <- apply(store.kappa,2:3,mean)
kappa.LCI <- apply(store.kappa,2:3,quantile,prob=(1-CI.level)/2)
kappa.UCI <- apply(store.kappa,2:3,quantile,prob=1-(1-CI.level)/2)
kappa.sig <- ( sign(kappa.LCI)*sign(kappa.UCI)==1)
kappa.LCI[k,] <- kappa.UCI[k,] <- kappa.sig[k,] <- NA

colnames(kappa.hat) <- colnames(kappa.LCI) <- colnames(kappa.UCI) <- colnames(kappa.sig) <- colnames(x_mat)
rownames(kappa.hat) <- rownames(kappa.LCI) <- rownames(kappa.UCI) <- rownames(kappa.sig) <- paste0('c',1:k)


if( !samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig)
}
if( samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig, 'samples'=store.kappa)
}

return(result)
}




BLCshrink.NGUS <- function(y,X,
	N_it=15000, burn.in=round(.2*N_it), CI.level=0.95, progress=TRUE, 
	c0=0.01, c1=0.01, theta=0.05, samples=FALSE, MCMC.seed=209578){

set.seed(MCMC.seed)
x_mat <- cbind(1,X)
colnames(x_mat)[1] <- "int"
p_mnl <- ncol(x_mat)
m <- nrow(x_mat)
k <- max(y)
Zj <- t(sapply(y, function(s) 1*(s==1:k) ))-.5
prog.count <- seq(0,N_it,length=11)[-1]

kappa <- array(0, dim = c(k,p_mnl)) 
lam_k <- 1/p_mnl 
tausqr <- array(1,c(k,p_mnl))
int.pos <- 1
coef.pos <- c(2:p_mnl)
tausqr[,int.pos] <- 100

store.kappa <- array(NA, dim = c(N_it,k,p_mnl))  

   
for( it in 1:N_it){

for ( j in 1:(k-1)){
 	xb=(x_mat)%*%t(kappa[-j,])  
	max_val <- apply(xb,1,max)
	Dj <- max_val + log(apply(exp(xb-max_val),1,sum))
	Bij <- x_mat%*% kappa[j,]-Dj
	Bij=ifelse(Bij>  1e20,1e20,Bij)
	Bij=ifelse(Bij< -1e20,-1e20,Bij)
	omg <- pgdraw::pgdraw(1,Bij)  
      var_mat <- solve(diag(1/((2*tausqr[j,])/lam_k))+t(x_mat)%*%diag(omg)%*%x_mat)
      mu_vec <- var_mat%*%t(x_mat)%*%( Zj[,j]+diag(omg)%*%Dj)
      kappa[j,] <- MASS::mvrnorm(1,mu_vec,var_mat)        
}

lam_k <- rgamma(1, c0+length(coef.pos)*(k-1)*0.5,  
	c1+0.25* sum( kappa[,coef.pos]^2 /tausqr[,coef.pos]) )
   
for ( j in 1:(k-1)){ for ( p in coef.pos ){
	tausqr[j,p]=GIGrvg::rgig(1,
		(theta-0.5),
		0.5*lam_k*kappa[j,p]^2,
		2*theta)
}}
tausqr <- ifelse(tausqr<1e-10,1e-10, tausqr)
tausqr <- ifelse(tausqr>1e10,1e10, tausqr)


store.kappa[it,,] <- kappa
if( progress ){ if( it %in% prog.count ){
	temp <- which(it==prog.count)
	print(paste0('MCMC progress:  ',temp*10,'%'))
	flush.console()
}}
}


store.kappa <- store.kappa[-(1:burn.in),,]
kappa.hat <- apply(store.kappa,2:3,mean)
kappa.LCI <- apply(store.kappa,2:3,quantile,prob=(1-CI.level)/2)
kappa.UCI <- apply(store.kappa,2:3,quantile,prob=1-(1-CI.level)/2)
kappa.sig <- ( sign(kappa.LCI)*sign(kappa.UCI)==1)
kappa.LCI[k,] <- kappa.UCI[k,] <- kappa.sig[k,] <- NA

colnames(kappa.hat) <- colnames(kappa.LCI) <- colnames(kappa.UCI) <- colnames(kappa.sig) <- colnames(x_mat)
rownames(kappa.hat) <- rownames(kappa.LCI) <- rownames(kappa.UCI) <- rownames(kappa.sig) <- paste0('c',1:k)


if( !samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig)
}
if( samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig, 'samples'=store.kappa)
}

return(result)
}


BLCshrink.HSSS <- function(y,X,
	N_it=15000, burn.in=round(.2*N_it), CI.level=0.95, progress=TRUE, 
	samples=FALSE, MCMC.seed=209578){


set.seed(MCMC.seed)
x_mat <- cbind(1,X)
colnames(x_mat)[1] <- "int"
p_mnl <- ncol(x_mat)
m <- nrow(x_mat)
k <- max(y)
Zj <- t(sapply(y, function(s) 1*(s==1:k) ))-.5
prog.count <- seq(0,N_it,length=11)[-1]

kappa <- array(0, dim = c(k,p_mnl)) 
xi <- 1
phi_sqr <- 1/p_mnl
delta_sqr <- rep(1,p_mnl)
int.pos <- 1
coef.pos <- c(2:p_mnl)
delta_sqr[int.pos] <- 100
eta <- rep(1,p_mnl)

store.kappa <- array(NA, dim = c(N_it,k,p_mnl))  

   
for( it in 1:N_it){

for ( j in 1:(k-1)){
	delta_phi <- delta_sqr*phi_sqr
	delta_phi[int.pos] <- 100
 	xb <- (x_mat)%*%t(kappa[-j,])  
	max_val <- apply(xb,1,max)
	Dj <- max_val + log(apply(exp(xb-max_val),1,sum))
	Bij <- x_mat%*% kappa[j,]-Dj
	Bij=ifelse(Bij>  1e20,1e20,Bij)
	Bij=ifelse(Bij< -1e20,-1e20,Bij)
	omg <- pgdraw::pgdraw(1,Bij)  
      var_mat <- solve(diag(1/delta_phi)+t(x_mat)%*%diag(omg)%*%x_mat)
      mu_vec <- var_mat%*%t(x_mat)%*%( Zj[,j]+diag(omg)%*%Dj)
      kappa[j,] <- MASS::mvrnorm(1,mu_vec,var_mat)        
}

for ( p in coef.pos){
	delta_sqr[p] <- pscl::rigamma(1,
		k/2,
		(1/eta[p])+((sum(kappa[,p]^2))/(2*phi_sqr) ))
}
delta_sqr <- ifelse( delta_sqr < 1e-10, 1e-10, delta_sqr)
delta_sqr <- ifelse( delta_sqr > 1e10, 1e10, delta_sqr)

for ( p in coef.pos){
	eta[p]=pscl::rigamma(1,1,1+(1/delta_sqr[p]))
}
 
phi_sqr <- pscl::rigamma(1,((k-1)*length(coef.pos) +1)/2,
	(1/xi)+0.5*sum(t(t((kappa[,coef.pos]^2))/delta_sqr[coef.pos])))
phi_sqr <- ifelse( phi_sqr < 1e-10, 1e-10, phi_sqr)
phi_sqr <- ifelse( phi_sqr > 1e10, 1e10, phi_sqr)
xi <- pscl::rigamma(1,1,1+(1/phi_sqr))
  
store.kappa[it,,] <- kappa
if( progress ){ if( it %in% prog.count ){
	temp <- which(it==prog.count)
	print(paste0('MCMC progress:  ',temp*10,'%'))
	flush.console()
}}
}


store.kappa <- store.kappa[-(1:burn.in),,]
kappa.hat <- apply(store.kappa,2:3,mean)
kappa.LCI <- apply(store.kappa,2:3,quantile,prob=(1-CI.level)/2)
kappa.UCI <- apply(store.kappa,2:3,quantile,prob=1-(1-CI.level)/2)
kappa.sig <- ( sign(kappa.LCI)*sign(kappa.UCI)==1)
kappa.LCI[k,] <- kappa.UCI[k,] <- kappa.sig[k,] <- NA

colnames(kappa.hat) <- colnames(kappa.LCI) <- colnames(kappa.UCI) <- colnames(kappa.sig) <- colnames(x_mat)
rownames(kappa.hat) <- rownames(kappa.LCI) <- rownames(kappa.UCI) <- rownames(kappa.sig) <- paste0('c',1:k)


if( !samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig)
}
if( samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig, 'samples'=store.kappa)
}

return(result)
}




BLCshrink.HSUS <- function(y,X,
	N_it=15000, burn.in=round(.2*N_it), CI.level=0.95, progress=TRUE, 
	samples=FALSE, MCMC.seed=209578){


set.seed(MCMC.seed)
x_mat <- cbind(1,X)
colnames(x_mat)[1] <- "int"
p_mnl <- ncol(x_mat)
m <- nrow(x_mat)
k <- max(y)
Zj <- t(sapply(y, function(s) 1*(s==1:k) ))-.5
prog.count <- seq(0,N_it,length=11)[-1]

kappa <- array(0, dim = c(k,p_mnl)) 
xi <- 1
phi_sqr <- 1/p_mnl
delta_sqr <- array(1,c(k,p_mnl))
int.pos <- 1
coef.pos <- c(2:p_mnl)
delta_sqr[,int.pos] <- 100
eta <- array(1,c(k,p_mnl))

store.kappa <- array(NA, dim = c(N_it,k,p_mnl))  

   
for( it in 1:N_it){

for ( j in 1:(k-1)){
	delta_phi <- delta_sqr[j,]*phi_sqr
	delta_phi[int.pos] <- 100
 	xb <- (x_mat)%*%t(kappa[-j,])  
	max_val <- apply(xb,1,max)
	Dj <- max_val + log(apply(exp(xb-max_val),1,sum))
	Bij <- x_mat%*% kappa[j,]-Dj
	Bij=ifelse(Bij>  1e20,1e20,Bij)
	Bij=ifelse(Bij< -1e20,-1e20,Bij)
	omg <- pgdraw::pgdraw(1,Bij)  
      var_mat <- solve(diag(1/delta_phi)+t(x_mat)%*%diag(omg)%*%x_mat)
      mu_vec <- var_mat%*%t(x_mat)%*%( Zj[,j]+diag(omg)%*%Dj)
      kappa[j,] <- MASS::mvrnorm(1,mu_vec,var_mat)        
}

for ( j in 1:(k-1) ){
for ( p in coef.pos){
	delta_sqr[j,p] <- pscl::rigamma(1,1,
		(1/eta[j,p])+((kappa[j,p]^2)/(2*phi_sqr) ))
}}
delta_sqr <- ifelse( delta_sqr < 1e-10, 1e-10, delta_sqr)
delta_sqr <- ifelse( delta_sqr > 1e10, 1e10, delta_sqr)

for ( j in 1:(k-1) ){
for ( p in coef.pos){
	eta[j,p]=pscl::rigamma(1,1,1+(1/delta_sqr[j,p]))
}}
 
phi_sqr <- pscl::rigamma(1,((k-1)*length(coef.pos) +1)/2,
	(1/xi)+0.5*sum( kappa[,coef.pos]^2 / delta_sqr[,coef.pos]))
phi_sqr <- ifelse( phi_sqr < 1e-10, 1e-10, phi_sqr)
phi_sqr <- ifelse( phi_sqr > 1e10, 1e10, phi_sqr)
xi <- pscl::rigamma(1,1,1+(1/phi_sqr))
  
store.kappa[it,,] <- kappa
if( progress ){ if( it %in% prog.count ){
	temp <- which(it==prog.count)
	print(paste0('MCMC progress:  ',temp*10,'%'))
	flush.console()
}}
}


store.kappa <- store.kappa[-(1:burn.in),,]
kappa.hat <- apply(store.kappa,2:3,mean)
kappa.LCI <- apply(store.kappa,2:3,quantile,prob=(1-CI.level)/2)
kappa.UCI <- apply(store.kappa,2:3,quantile,prob=1-(1-CI.level)/2)
kappa.sig <- ( sign(kappa.LCI)*sign(kappa.UCI)==1)
kappa.LCI[k,] <- kappa.UCI[k,] <- kappa.sig[k,] <- NA

colnames(kappa.hat) <- colnames(kappa.LCI) <- colnames(kappa.UCI) <- colnames(kappa.sig) <- colnames(x_mat)
rownames(kappa.hat) <- rownames(kappa.LCI) <- rownames(kappa.UCI) <- rownames(kappa.sig) <- paste0('c',1:k)


if( !samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig)
}
if( samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig, 'samples'=store.kappa)
}

return(result)
}


BLCshrink.none <- function(y,X,
	N_it=15000, burn.in=round(.2*N_it), CI.level=0.95, progress=TRUE, 
	samples=FALSE, MCMC.seed=209578, prior.var=100){


set.seed(MCMC.seed)
x_mat <- cbind(1,X)
colnames(x_mat)[1] <- "int"
p_mnl <- ncol(x_mat)
m <- nrow(x_mat)
k <- max(y)
Zj <- t(sapply(y, function(s) 1*(s==1:k) ))-.5
prog.count <- seq(0,N_it,length=11)[-1]

kappa <- array(0, dim = c(k,p_mnl)) 
delta_sqr <- rep(prior.var,p_mnl)
store.kappa <- array(NA, dim = c(N_it,k,p_mnl))  

   
for( it in 1:N_it){

for ( j in 1:(k-1)){
 	xb <- (x_mat)%*%t(kappa[-j,])  
	max_val <- apply(xb,1,max)
	Dj <- max_val + log(apply(exp(xb-max_val),1,sum))
	Bij <- x_mat%*% kappa[j,]-Dj
	Bij=ifelse(Bij>  1e20,1e20,Bij)
	Bij=ifelse(Bij< -1e20,-1e20,Bij)
	omg <- pgdraw::pgdraw(1,Bij)  
      var_mat <- solve(diag(1/delta_sqr)+t(x_mat)%*%diag(omg)%*%x_mat)
      mu_vec <- var_mat%*%t(x_mat)%*%( Zj[,j]+diag(omg)%*%Dj)
      kappa[j,] <- MASS::mvrnorm(1,mu_vec,var_mat)        
}

  
store.kappa[it,,] <- kappa
if( progress ){ if( it %in% prog.count ){
	temp <- which(it==prog.count)
	print(paste0('MCMC progress:  ',temp*10,'%'))
	flush.console()
}}
}


store.kappa <- store.kappa[-(1:burn.in),,]
kappa.hat <- apply(store.kappa,2:3,mean)
kappa.LCI <- apply(store.kappa,2:3,quantile,prob=(1-CI.level)/2)
kappa.UCI <- apply(store.kappa,2:3,quantile,prob=1-(1-CI.level)/2)
kappa.sig <- ( sign(kappa.LCI)*sign(kappa.UCI)==1)
kappa.LCI[k,] <- kappa.UCI[k,] <- kappa.sig[k,] <- NA

colnames(kappa.hat) <- colnames(kappa.LCI) <- colnames(kappa.UCI) <- colnames(kappa.sig) <- colnames(x_mat)
rownames(kappa.hat) <- rownames(kappa.LCI) <- rownames(kappa.UCI) <- rownames(kappa.sig) <- paste0('c',1:k)


if( !samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig)
}
if( samples ){
	result <- list('kappa.hat'=kappa.hat, 'kappa.LCI'=kappa.LCI, 
		'kappa.UCI'=kappa.UCI, 'kappa.sig'=kappa.sig, 'samples'=store.kappa)
}

return(result)
}




