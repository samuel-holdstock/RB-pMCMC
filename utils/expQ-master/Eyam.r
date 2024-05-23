setwd('Downloads/PhD/MCMC/Inference/Exponential/')
library("Rcpp")
library("RcppArmadillo")
source("Utils/expQ-master/genericQ.r")
## If the library does not work, use sourceCpp
sourceCpp("Utils/expQ-master/rexpQ.cpp")
sourceCpp("Utils/epidemicsRCpp.cpp")
library("expQ")
## Need to install.packages("expm") if package not already installed

## User settings 
ALLDATA=FALSE ## set to FALSE to jump from initial state to end state
## Choose method from uniformisation, scaling and squaring, and expAtv
method="Unif"
##method="SS"  ## beware - this with ALLDATA=FALSE is very slow
##method="expAtv" ## existing package that uses Krylov subspaces
## Other parameters
nits=1  ## number of repeats of the likelihood calculation
prec=1e-15 ## epsilon
thetas=c(0.0196,3.204) ## modal value; i.e. MLE
##thetas=c(0.0178,2.173) ## value used in Ho, Crawford and Suchard (2018)

## Create the data
S=c(254,235,201,153,121,110,97,83)
I=c(7,14,22,29,20,8,8,0)
t=c(0,0.5,1,1.5,2,2.5,3,4)
Eyam=data.frame(t=t,S=S,I=I)

#####################################

## Get the log-likelihood
getll<-function(ds,nuins,Qs,method="Unif",prec=1e-10) {
    ll=0
    for (i in 1:length(ds)) {
        if (method=="Unif") {
            nuout=Unif_v_exp_Q(t(nuins[[i]]),Qs[[i]],prec)
        }
        else if (method=="SS") {
            nuout=SS_v_exp_Q(t(nuins[[i]]),Qs[[i]],prec)
        }
        else if (method=="expAtv") {
            nuout=expm::expAtv(Matrix::t(Qs[[i]]),nuins[[i]])$eAtv
        }
        else {
            stop("Bad method")
        }

        ll=ll+log(nuout[ds[i]])
    }
    return(ll)
}

if (ALLDATA) {
    nobs=length(Eyam$S)-1
    ds=rep(0,nobs)
    Qs<-list(); nuins=list()
    for (i in 1:nobs) {
        Qs[[i]]=SIRgetQslim(as.numeric(Eyam[i,2:3]),as.numeric(Eyam[i+1,2:3]),thetas)*(Eyam$t[i+1]-Eyam$t[i])
        ds[i]=nrow(Qs[[i]])-1
        nuins[[i]]=c(1,rep(0,ds[i]))
    }
} else { ## straight from initial condition to final value
    n=nrow(Eyam)
    x0=as.numeric(Eyam[1,2:3])
    x1=as.numeric(Eyam[n,2:3])
    Qs<-list(); nuins=list()

    tstart=Sys.time()
    Qs[[1]]=SIRgetQslim(x0,x1,thetas)*(Eyam$t[n]-Eyam$t[1])
    tend=Sys.time()
    print("Time for Q")
    print(tend-tstart)

    d=nrow(Qs[[1]])-1; 
    nuin=rep(0,d+1); nuin[1]=1
    nuins[[1]]=nuin
    ds=c(d)
}

print("Starting")
tstart=Sys.time()
for (i in 1:nits) {
    ll=getll(ds,nuins,Qs,method,prec)
}
tend=Sys.time()

options(digits=17)
print((tend-tstart))
print(ll)

get_nu = function(i,lower,upper){
  len = upper - lower + 2
  nu = rep(0,len)
  nu[i] = 1
  return(nu)
}

round_2 = function(x){
  floor(x)+floor(x)%%2
}

thetas = c(1/1000,1/1000,1/10)
width = 10
upper = target + width/2
lower = target - width/2
tau = 1

nuins = get_nu(i,lower,upper)
Qs = BDIGetQ(thetas,c(lower,upper))*tau

microbenchmark(
  a = {
    Unif_v_exp_Q(t(nuins),Qs,prec=1e-10)
  },unit = 'micro'
)
# 100 - 700
# 200 - 1200
microbenchmark(
  a = {
    SS_exp_Q(Qs,prec=1e-10)
  },unit = 'micro'
)
microbenchmark(
  b = { 
    expm(as.matrix(Qs))[i,]
  },unit = 'micro'
)
# 100 - 7000
# 200 - 40000


1400*3*length(taus)*length(widths)/(1000*1000*60)

