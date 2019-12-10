setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
import::here(sSINDy,sSINDy_lm,sSINDy_lm_vectorized,order_coef,
             error_compute,error_compute_lm,.from = 'sSINDy.R')







size_data=100
interval=100
#Creat a data from unifom [-100,100]
x1=runif(size_data,min=-interval,max=interval)
x2=rnorm(size_data,sd=0.1)
x3=runif(size_data,min=-interval,max=interval)
x4=runif(size_data,min=-interval,max=interval)
x5=runif(size_data,min=-interval,max=interval)
x6=runif(size_data,min=-interval,max=interval)
x7=rnorm(size_data)
X=cbind(x1,x2,x3,x4,x5,x6,x7)
# The model is a linear combination of first five variables
# x6 and x7 are data which is not used
y=10*x1+10*x2+0.1*x3+0.1*x4+.1*x5


cv_group=sample(rep(1:10,1000),n_times)
n_fold=10

sSINDy(X,y,10)
lm(y~X[,-7]-1)%>%summary

library(sde)
simulate_BM=function(fun_sim,T,n_times){
  t_diff=T/n_times
  X=sde::BM(x=10, t0=0, T, n_times)
  #X=GBM(x=1, r=10,sigma=1,T=1, N=n_times)
  y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
  X_data=X[1:n_times]
  return(sSINDy_lm_vectorized(cbind(X_data^0,X_data,X_data^2,X_data^3),y_data,10)$coef)
}


simulate_BM(sde::BM,1,100000)
coef_sim=replicate(100,simulate_BM(BM,1,10000))
apply(coef_sim,1,mean)
variation=(X[2:(n_times+1)]-X[1:n_times])*
  (X[2:(n_times+1)]-X[1:n_times])/t_diff

sSINDy(cbind(X_data,X_data^2,X_data^3),variation,10)

plot(BM())

# X=sde::BM(x=10, t0=0, T=1, N=10000)
# library(Sim.DiffProc)
X=Sim.DiffProc::BM(10000,M=1,x0=1,T=1)
# X=read.csv("test.csv")
t_diff=0.0001

#X=GBM(x=1, r=10,sigma=1,T=1, N=n_times)
y_data=(X[2:(n_times+1)]-X[1:n_times])/t_diff
X_data=X[1:n_times]
plot(density(y_data*t_diff))


data_X=X_data
data_y=y_data
cv_group=sample(rep(1:10,1000),n_times)
n_fold=10
