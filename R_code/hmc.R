data <- read.table("data",header=F)

y <- data[,3]

#constructing A-inverse
ped <- read.table("ped",header=T)
n.ped <- length(ped[,1])
a.mat <- matrix(0,n.ped,n.ped)
a.inv <- matrix(0,n.ped,n.ped)
f <- array(0,n.ped)
tt <- c(1.0,-0.5,-0.5)

for(i in 1:n.ped){
    a.mat[i,i] <- 1
    if(ped[i,2]!=0 && ped[i,3]!=0){
        a.mat[i,i] <- a.mat[i,i]+0.5*a.mat[ped[i,2],ped[i,3]]
        a.mat[1:i-1,i] <- 0.5*(a.mat[1:i-1,ped[i,2]]+a.mat[1:i-1,ped[i,3]])
    }else if(ped[i,2]==0 && ped[i,3]!=0){
        a.mat[1:i-1,i] <- 0.5*(a.mat[1:i-1,ped[i,3]])
    }else if(ped[i,2]!=0 && ped[i,3]==0){
        a.mat[1:i-1,i] <- 0.5*(a.mat[1:i-1,ped[i,2]])
    }
    a.mat[i,1:i-1] <- a.mat[1:i-1,i]
    f[i] <- a.mat[i,i]-1
}

#A-inverse
for(i in 1:n.ped){
    if(ped[i,2]==0 && ped[i,3]==0){
        D <- 1
    }else if(ped[i,2]!=0 && ped[i,3]!=0){
        D <- 0.5-0.25*(f[ped[i,2]]+f[ped[i,3]])
    }else if(ped[i,2]!=0 && ped[i,3]==0){
        D <- 0.75-0.25*(f[ped[i,2]])
    }else if(ped[i,2]==0 && ped[i,3]!=0){
        D <- 0.75-0.25*(f[ped[i,3]])
    }
    for(j in 1:3){
        for(k in 1:3){
            if(ped[i,j]!=0 && ped[i,k]!=0)a.inv[ped[i,j],ped[i,k]] <- a.inv[ped[i,j],ped[i,k]]+tt[j]*tt[k]/D
        }
    }
}

#X-Z-matrix
n <- length(data[,1])
p <- 2
X <- matrix(0,n,p)
Z <- matrix(0,n,n.ped)
for(i in 1:n){
    X[i,data[i,2]] <- 1
    Z[i,data[i,1]] <- 1
}
W <- cbind(X,Z)

#-------------HMC----------
iter <- 10000; burn <- 1000; thin <- 10
sam <- matrix(0,iter,1)
sam[seq(from=(burn+1),to=iter,by=thin)] <- 1

b.sam.hmc <- matrix(0,p,1)
b.sam.hmc2 <- matrix(0,iter,p)

u.sam.hmc <- matrix(0,n.ped,1)
u.sam.hmc2 <- matrix(0,iter,100)

var.e.sam.hmc <- matrix(0,iter,1)
var.u.sam.hmc <- matrix(0,iter,1)

tau <- 10000
v.e <- 10^-6; lambda.e <- 10^-6
v.u <- 10^-6; lambda.u <- 10^-6

accept.var.u <- 0
accept.var.e <- 0
accept.b <- matrix(0,p,1)
accept.u <- matrix(0,n.ped,1)

var.e <- 0.5; var.e.new <- 0
var.u <- 0.5; var.u.new <- 0

b <- matrix(0,p,1)
u <- matrix(0,n.ped,1)

xpx <- matrix(0,p,1)
zpz <- matrix(0,n.ped,1)
for(i in 1:p){
    xpx[i] <- crossprod(X[,i])
}
for(i in 1:n.ped){
    zpz[i] <- crossprod(Z[,i])
}

for(i in 1:iter){
  
# Hamiltonian b
    ee <- y-X%*%b-Z%*%u
    for(jj in 1:p){
        LL.b <- 7
        epsilon.b <- sqrt(1/(xpx[jj]/var.e+1/tau))/(0.1589825*20)

        b.new <- b[jj]
        Xe <- crossprod(ee+X[,jj]*b.new, X[,jj])

        pp <- rnorm(1)
        K1 <- t(pp)%*%pp/2
        U1 <- -((2*Xe*b.new-xpx[jj]*b.new^2)/(2*var.e)-b.new^2/(2*tau))
        H1 <- (U1+K1)

        for(L in 1:LL.b){
            pp <- pp-0.5*epsilon.b*(-((Xe-xpx[jj]*b.new)/var.e-b.new/tau))
            b.new <- b.new+epsilon.b*pp
            pp <- pp-0.5*epsilon.b*(-((Xe-xpx[jj]*b.new)/var.e-b.new/tau))
        }

        K2 <- t(pp)%*%pp/2
        U2 <- -((2*Xe*b.new-xpx[jj]*b.new^2)/(2*var.e)-b.new^2/(2*tau))
        H2 <- (U2+K2)

        if(runif(1) < exp(-H2+H1)){
            accept.b[jj] <- accept.b[jj]+1
        }else{
            b.new <- b[jj]
        }
        
        ee <- ee+X[,jj]*(b[jj]-c(b.new))
        b[jj] <- b.new
    }

# Hamiltonian u
    ee <- y-X%*%b-Z%*%u
    for(jj in 1:n.ped){
        LL.u <- 7
        epsilon.u <- sqrt(1/(zpz[jj]/var.e+a.inv[jj,jj]/var.u))/(0.1589825*20)

        u.new <- u[jj]
        Ze <- crossprod(ee+Z[,jj]*u.new, Z[,jj])
        uG <- crossprod(a.inv[,jj],u)-a.inv[jj,jj]*u.new

        pp <- rnorm(1)
        K1 <- t(pp)%*%pp/2
        U1 <- -((2*Ze*u.new-zpz[jj]*u.new^2)/(2*var.e)-(2*uG*u.new+a.inv[jj,jj]*u.new^2)/(2*var.u))
        H1 <- (U1+K1)

        for(L in 1:LL.u){
            pp <- pp-0.5*epsilon.u*(-((Ze-zpz[jj]*u.new)/var.e-(uG+a.inv[jj,jj]*u.new)/var.u))
            u.new <- u.new+epsilon.u*pp
            pp <- pp-0.5*epsilon.u*(-((Ze-zpz[jj]*u.new)/var.e-(uG+a.inv[jj,jj]*u.new)/var.u))
        }

        K2 <- t(pp)%*%pp/2
        U2 <- -((2*Ze*u.new-zpz[jj]*u.new^2)/(2*var.e)-(2*uG*u.new+a.inv[jj,jj]*u.new^2)/(2*var.u))
        H2 <- (U2+K2)

        if(runif(1) < exp(-H2+H1)){
            accept.u[jj] <- accept.u[jj]+1
        }else{
            u.new <- u[jj]
        }

        ee <- ee+Z[,jj]*(u[jj]-c(u.new))
        u[jj] <- u.new
    }

#leapfrog_for_genetic_variance
    var.u.new <- var.u
    uGu <- t(u)%*%a.inv%*%u

    epsilon.var.u <- sqrt(uGu^2/((n.ped-1)^2*(n.ped-2)))/(0.112485939*20)
    LL.var.u <- 7

    pp <- rnorm(1)
    K1 <- t(pp)%*%pp/2
    U1 <- -(-((n.ped+v.u)/2+1)*log(var.u.new)-(uGu+lambda.u)/(2*var.u.new))
    H1 <- (U1+K1)

    for(L in 1:LL.var.u){
        pp <- pp-0.5*epsilon.var.u*(-(-((n.ped+v.u)/2+1)/var.u.new+0.5*(uGu+lambda.u)/(var.u.new^2)))
        var.u.new <- var.u.new+epsilon.var.u*pp
        pp <- pp-0.5*epsilon.var.u*(-(-((n.ped+v.u)/2+1)/var.u.new+0.5*(uGu+lambda.u)/(var.u.new^2)))
    }

    if(!is.na(var.u.new)&var.u.new > 0){
        K2 <- t(pp)%*%pp/2
        U2 <- -(-((n.ped+v.u)/2+1)*log(var.u.new)-(uGu+lambda.u)/(2*var.u.new))
        H2 <- (U2+K2)

        if(runif(1) < exp(-H2+H1)){
            accept.var.u <- accept.var.u+1
            var.u <- var.u.new
        }
    }

#leapfrog1
    var.e.new <- var.e
    e <- y-X%*%b-Z%*%u
    ee <- t(e)%*%e

    epsilon.var.e <- sqrt(ee^2/((n-1)^2*(n-2)))/(0.112485939*20)
    LL.var.e <- 7

    pp <- rnorm(1)
    K1 <- t(pp)%*%pp/2
    U1 <- -(-((n+v.e)/2+1)*log(var.e.new)-(ee+lambda.e)/(2*var.e.new))
    H1 <- (U1+K1)

    for(L in 1:LL.var.e){
        pp <- pp-0.5*epsilon.var.e*(-(-((n+v.e)/2+1)/var.e.new+0.5*(ee+lambda.e)/(var.e.new^2)))
        var.e.new <- var.e.new+epsilon.var.e*pp
        pp <- pp-0.5*epsilon.var.e*(-(-((n+v.e)/2+1)/var.e.new+0.5*(ee+lambda.e)/(var.e.new^2)))
    }

    if(!is.na(var.e.new)&var.e.new > 0){
        K2 <- t(pp)%*%pp/2
        U2 <- -(-((n+v.e)/2+1)*log(var.e.new)-(ee+lambda.e)/(2*var.e.new))
        H2 <- (U2+K2)

        if(runif(1) < exp(-H2+H1)){
            accept.var.e <- accept.var.e+1
            var.e <- var.e.new
        }
    }

	print(c(as.character(i),var.u,var.e))

	b.sam.hmc <- b.sam.hmc+b*sam[i]
	u.sam.hmc <- u.sam.hmc+u*sam[i]
	var.e.sam.hmc[i] <- var.e
	var.u.sam.hmc[i] <- var.u

    b.sam.hmc2[i,1:2] <- b
    u.sam.hmc2[i,1:100] <- u[501:600]
}

vvv <- cbind(var.e.sam.hmc,var.u.sam.hmc)
colnames(vvv) <- c("R_var","G_var")
write.table(vvv,file="variance_sample_HMC.txt",col.names=T,row.names=F,quote=F)
aaa <- cbind(accept.var.e,accept.var.u)
colnames(aaa) <- c("R_var","G_var")
write.table(aaa,file="variance_sample_accept_HMC.txt",col.names=T,row.names=F,quote=F)

bbb <- b.sam.hmc/sum(sam)
bbb <- cbind(bbb,accept.b)
colnames(bbb) <- c("b_est","Accept")
write.table(bbb,file="b_HMC.txt",col.names=T,row.names=F,quote=F)
write.table(b.sam.hmc2,file="b_sample_HMC.txt",col.names=F,row.names=F,quote=F)

uuu <- u.sam.hmc/sum(sam)
uuu <- cbind(uuu,accept.u)
colnames(uuu) <- c("u_est","Accept")
write.table(uuu,file="u_HMC.txt",col.names=T,row.names=F,quote=F)
write.table(u.sam.hmc2,file="u_sample_HMC.txt",col.names=F,row.names=F,quote=F)

