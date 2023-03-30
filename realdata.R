library(MASS)
library(glmnet)
library(matlib)
library(corpcor)
library(MESS)
#=============================================================== read data
data_y=read.csv("/Users/shenwenjun/Desktop/materials/preparation/code/CS_y1.csv",header = FALSE)
data_Exogenous=read.delim("/Users/shenwenjun/Desktop/materials/preparation/code/CS_Exogenous1.txt", sep =",",header = FALSE)
data_Instruments=read.csv("/Users/shenwenjun/Desktop/materials/preparation/code/CS_Instruments1.csv",header = FALSE)
data_Endogenous=read.csv("/Users/shenwenjun/Desktop/materials/preparation/code/CS_Endogenous1.csv",header = FALSE)
#=============================================================== Define variables
## set number of IVs in every submodel and total number of submodels
K1 = 11            
t1 = 7
K1_plus = 11
t1_plus = 7

fn_scale=function(variable)
{
  return(scale(variable,center = TRUE, scale = TRUE))
}
Xo=data_Exogenous[,-40] 
Xo=apply(Xo,2,fn_scale)
Xo=as.matrix(Xo)
XI=data_Instruments[,-c(39,40)]
XI=apply(XI,2,fn_scale)
Xe=data_Endogenous$V1
Xe=scale(Xe,center=TRUE, scale=TRUE)
mean(Xe)   
Y = data_y$V1
Y=scale(Y,center=TRUE, scale=TRUE)

##plot Figure 2 showing the correlation between Xe and the potential instruments
corXIXe = t(cor(Xe, XI))
plot(1:147,corXIXe[1:147,],xlab="Index of Potential Instruments", ylab="Correlation with Xe", main=expression(paste("Correlation between Xe and Potential Instruments in Dataset")),  pch=18, col="blue", xlim=c(1,147), ylim=c(0,0.7))
axis(1, at = c(10,20,30,40,50,60,70,80,90,100,110,120,130,140), labels = c(10,20,30,40,50,60,70,80,90,100,110,120,130,140))
rug(x = c(0, 10,20,30,40,50,60,70,80,90,100,110,120,130,140)+ 5, ticksize = -0.01, side = 1)
###
xIpos = which(abs(cor(Xe, XI))>0.3) 
XI = XI[,xIpos]
dim(XI)
a = as.numeric(cor(Xe,XI))
a=abs(a)

##subsetting the exogenous variables
md1=lm(Y~Xo+Xe+0)
Xo=Xo[,c(which(summary(md1)$coefficients[,4]<=0.1))]

N=dim(data_y)[1]
q=dim(XI)[2]   
p=dim(Xo)[2]   
####==== set the elasticnet mixing parameter in glmnet() for IVs selection by lasso & elastic net
anet1 = 1 
anetel1 = 0.5
anetel2 = 1/3

#==================================== Big function include all methods
fn_big=function(y, xe, xI, xo, n,t,K)
{
  t1=t
  K1=K
  t1_plus=t
  K1_plus=K
  a = as.numeric(cor(xe,xI))
  a=abs(a)
  P1=xo%*%solve(t(xo)%*%xo, tol=1e-30)%*%t(xo)
  In=matrix(0, n, n)
  diag(In)=1
  zI=(In-P1)%*%xI
  M_Ie_tilde=(1/n)*t(zI)%*%xe
  M_II_tilde=(1/n)*t(zI)%*%zI
  omega_hat= (1/n)*t(xI)%*%xI
  v_I = (1/n)*t(xI)%*%xe
  v_o=(1/n)*t(xo)%*%xe
  D_Matrix=(1/n)*t(xI)%*%xo
  Moo=1/n*t(xo)%*%xo
  Moe=(1/n)*t(xo)%*%xe
  a_bias=(1/n)*t(xe)%*%xo%*%solve(t(xo)%*%(xo), tol=1e-30)%*%t(xo)%*%xe
  v=(1/n)*t(xI)%*%xe
  
  #============================ Method 1: Naive
  fit_naive = lm(y~xe + xo + 0)
  beta_M1 = fit_naive$coefficients[1]
  betas1 = fit_naive$coefficients[which(is.na(fit_naive$coefficients)==FALSE)]
  out1 = summary(fit_naive)
  SE_M1 = out1$coefficients[1,2]      
  S1 = '-'
  SXo1 = dim(xo)[2]
  pvalue_M1 = out1$coefficients[1,4]
  fstat1 = '-'
  dfR1 = '-'
  dfE1 = '-'
  pva1 = '-'
  Rsq1 = '-'
  Exe_time1 = '-'
  
  #===== Define function for methods 2: MA on t IVs for K submodels, sampling with equal probability 
  fn_MA = function(y, xe, xI, xo, n, t, K)  {
    t1=t
    K1=K
    p=dim(xo)[2]
    q=dim(xI)[2]
    start_stageone = Sys.time()
    Xej_comb = matrix(0, nrow= n, ncol= K)    
    position = matrix(0, K, t)
    for (J in 1:K) { 
      count = 0
      repeat {
        set.seed(J+count+100)
        count = count + 1
        posJ = sample(1:q, t, replace=FALSE)
        step = 0
        for (h in 1:J) {
          if (sum((sort(posJ)-sort(position[h,]))^2)!=0) {
            step = step + 1
          }
        }
        if (step == J) {
          break
        }
      }
      position[J,] = posJ
      fitJ = lm(xe~xI[,posJ] + xo + 0)
      d.hatJ <- fitJ$fitted.values
      Xej_comb[,J] <- d.hatJ         
    }
    fit_SM = lm(xe~ Xej_comb + 0)         
    Xe_MA = fit_SM$fitted.values
    end_stageone = Sys.time()
    MA_exetime =  difftime(end_stageone, start_stageone, units="secs")
    num_position = as.numeric(position)
    S = length(unique(num_position))
    smy =summary(fit_SM)
    fstatma = smy$fstatistic[1]
    dfRma = smy$fstatistic[2]
    dfEma = smy$fstatistic[3]
    pvama = pf(fstatma, dfRma, dfEma, lower.tail = FALSE)
    Rsq = var(Xe_MA)/var(xe)
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-1) ) 
    psi = matrix(0, K, K)
    v_dk = matrix(0, (t+p)*K, 1)
    d_dk = matrix(0, (t+p)*K, p)
    v_diag = matrix(0, (t+p)*K, K)
    sigma_inv = matrix(0, (t+p)*K, (t+p)*K)
    Ik_vector = matrix(1, K, 1)
    for (I in 1:K) {
      v_dk[(I*(t+p)-(t+p-1)):(I*(t+p)),] = c(v_o,v_I[position[I,],])
      v_diag[(I*(t+p)-(t+p-1)):(I*(t+p)), I] = c(v_o,v_I[position[I,],])
      d_dk[(I*(t+p)-(t+p-1)):(I*(t+p)),] = rbind(Moo,D_Matrix[position[I,],])
      sigma_11=Moo
      sigma_12=t(D_Matrix)[,position[I,]]
      sigma_21=D_Matrix[position[I,],]
      sigma_22=omega_hat[position[I,], position[I,]]
      sigma_II = rbind(cbind(sigma_11,sigma_12),cbind(sigma_21,sigma_22))
      sigma_inv[(I*(t+p)-(t+p-1)):(I*(t+p)), (I*(t+p)-(t+p-1)):(I*(t+p))] = solve(sigma_II, tol=1e-30)
      
      
      for (J in 1:K) {
        vi = c(v_o,v_I[position[I,],])
        vj = c(v_o,v_I[position[J,],])
        sigma_11=Moo
        sigma_12_i=t(D_Matrix)[,position[I,]]
        sigma_21_i=D_Matrix[position[I,],]
        sigma_22_i=omega_hat[position[I,], position[I,]]
        sigma_ii = rbind(cbind(sigma_11,sigma_12_i),cbind(sigma_21_i,sigma_22_i))
        sigma_12_j=t(D_Matrix)[,position[J,]]
        sigma_21_j=D_Matrix[position[J,],]
        sigma_22_j=omega_hat[position[J,], position[J,]]
        sigma_22_ij= omega_hat[position[I,],position[J,]]
        sigma_jj = rbind(cbind(sigma_11,sigma_12_j),cbind(sigma_21_j,sigma_22_j))
        sigma_ii_inv = solve(sigma_ii, tol=1e-30)
        sigma_jj_inv = solve(sigma_jj, tol=1e-30)
        sigma_ij = rbind(cbind(sigma_11,sigma_12_j),cbind(sigma_21_i,sigma_22_ij))
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    d = t(v_diag)%*%sigma_inv%*%d_dk%*%solve(Moo, tol=1e-30)%*%Moe
    bias_beta=a_bias*t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector-1)/(t(u)%*%solve(psi, tol=1e-30)%*%u-a_bias*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector)^2)* beta_MA
    SE_Md=sigma_MA_hat*sqrt(as.numeric(t(u)%*%solve(psi, tol=1e-30)%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%solve(psi, tol=1e-30)%*%u))/(abs(t(u)%*%solve(psi, tol=1e-30)%*%u-a_bias*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector)^2)*sqrt(n))
    #nu=sqrt(3.25)*sqrt(as.numeric(t(u)%*%solve(psi, tol=1e-30)%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%solve(psi, tol=1e-30)%*%u))/(abs(t(u)%*%solve(psi, tol=1e-30)%*%u-a_bias*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector)^2)*sqrt(n))
    beta_nobias=beta_MA-bias_beta
    t_Md = beta_MA/SE_Md
    pvalue_Md = 2*pnorm(abs(t_Md), lower.tail = FALSE)
    result=matrix(c(S, fstatma, dfRma, dfEma, pvama, Rsq, MA_exetime, beta_nobias, SE_Md, pvalue_Md), nrow =1)
    return(result)
  }
  
  #============================ Method 2: MA(t1, K1)
  MA_M2 = fn_MA(y, xe, xI, xo, n, t1, K1)
  S2 = MA_M2[,1]
  fstat2 = MA_M2[,2]
  dfR2 = MA_M2[,3]
  dfE2 = MA_M2[,4]
  pva2 = MA_M2[,5]
  Rsq2 = MA_M2[,6]
  Exe_time2 = MA_M2[,7]
  beta_M2 = MA_M2[,8]
  SE_M2 = MA_M2[,9]
  pvalue_M2 = MA_M2[,10]
  #===== Define function for method 3 : MA on t IVs for K submodels, sampling with unequal probability
  fn_MA_plus = function(y, xe, xI, xo, n, t, K)  {
    a = as.numeric(cor(xe,xI))
    a=abs(a)
    p=dim(xo)[2]
    q=dim(xI)[2]
    start_stageone = Sys.time()
    Xej_comb = matrix(0, nrow= n, ncol= K)    
    position = matrix(0, K, t)
    for (J in 1:K) { 
      count = 0
      repeat {
        set.seed(J+count+100)
        count = count + 1
        posJ = sample(1:q, t, replace=FALSE,prob=a)
        step = 0
        for (h in 1:J) {
          if (sum((sort(posJ)-sort(position[h,]))^2)!=0) {
            step = step + 1
          }
        }
        if (step == J) {
          break
        }
      }
      position[J,] = posJ
      fitJ = lm(xe~xI[,posJ] + xo + 0)
      d.hatJ <- fitJ$fitted.values
      Xej_comb[,J] <- d.hatJ         
    }
    fit_SM = lm(xe~ Xej_comb + 0)         
    Xe_MA = fit_SM$fitted.values
    end_stageone = Sys.time()
    MA_exetime =  difftime(end_stageone, start_stageone, units="secs")
    num_position = as.numeric(position)
    S = length(unique(num_position))
    smy =summary(fit_SM)
    fstatma = smy$fstatistic[1]
    dfRma = smy$fstatistic[2]
    dfEma = smy$fstatistic[3]
    pvama = pf(fstatma, dfRma, dfEma, lower.tail = FALSE)
    Rsq = var(Xe_MA)/var(xe)
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-1) ) 
    psi = matrix(0, K, K)
    v_dk = matrix(0, (t+p)*K, 1)
    d_dk = matrix(0, (t+p)*K, p)
    v_diag = matrix(0, (t+p)*K, K)
    sigma_inv = matrix(0, (t+p)*K, (t+p)*K)
    Ik_vector = matrix(1, K, 1)
    for (I in 1:K) {
      v_dk[(I*(t+p)-(t+p-1)):(I*(t+p)),] = c(v_o,v_I[position[I,],])
      v_diag[(I*(t+p)-(t+p-1)):(I*(t+p)), I] = c(v_o,v_I[position[I,],])
      d_dk[(I*(t+p)-(t+p-1)):(I*(t+p)),] = rbind(Moo,D_Matrix[position[I,],])
      sigma_11=Moo
      sigma_12=t(D_Matrix)[,position[I,]]
      sigma_21=D_Matrix[position[I,],]
      sigma_22=omega_hat[position[I,], position[I,]]
      sigma_II = rbind(cbind(sigma_11,sigma_12),cbind(sigma_21,sigma_22))
      sigma_inv[(I*(t+p)-(t+p-1)):(I*(t+p)), (I*(t+p)-(t+p-1)):(I*(t+p))] = solve(sigma_II, tol=1e-30)
      
      
      for (J in 1:K) {
        vi = c(v_o,v_I[position[I,],])
        vj = c(v_o,v_I[position[J,],])
        sigma_11=Moo
        sigma_12_i=t(D_Matrix)[,position[I,]]
        sigma_21_i=D_Matrix[position[I,],]
        sigma_22_i=omega_hat[position[I,], position[I,]]
        sigma_ii = rbind(cbind(sigma_11,sigma_12_i),cbind(sigma_21_i,sigma_22_i))
        sigma_12_j=t(D_Matrix)[,position[J,]]
        sigma_21_j=D_Matrix[position[J,],]
        sigma_22_j=omega_hat[position[J,], position[J,]]
        sigma_22_ij= omega_hat[position[I,],position[J,]]
        sigma_jj = rbind(cbind(sigma_11,sigma_12_j),cbind(sigma_21_j,sigma_22_j))
        sigma_ii_inv = solve(sigma_ii, tol=1e-30)
        sigma_jj_inv = solve(sigma_jj, tol=1e-30)
        sigma_ij = rbind(cbind(sigma_11,sigma_12_j),cbind(sigma_21_i,sigma_22_ij))
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    d = t(v_diag)%*%sigma_inv%*%d_dk%*%solve(Moo, tol=1e-30)%*%Moe
    bias_beta=a_bias*t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector-1)/(t(u)%*%solve(psi, tol=1e-30)%*%u-a_bias*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector)^2)* beta_MA
    SE_Md=sigma_MA_hat*sqrt(as.numeric(t(u)%*%solve(psi, tol=1e-30)%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%solve(psi, tol=1e-30)%*%u))/(abs(t(u)%*%solve(psi, tol=1e-30)%*%u-a_bias*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector)^2)*sqrt(n))
    #nu=sqrt(3.25)*sqrt(as.numeric(t(u)%*%solve(psi, tol=1e-30)%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%solve(psi, tol=1e-30)%*%u))/(abs(t(u)%*%solve(psi, tol=1e-30)%*%u-a_bias*(t(u)%*%solve(psi, tol=1e-30)%*%Ik_vector)^2)*sqrt(n))
    beta_nobias=beta_MA-bias_beta
    t_Md = beta_MA/SE_Md
    pvalue_Md = 2*pnorm(abs(t_Md), lower.tail = FALSE)
    result=matrix(c(S, fstatma, dfRma, dfEma, pvama, Rsq, MA_exetime, beta_nobias, SE_Md, pvalue_Md), nrow =1)
    return(result)
  }
  
  #============================ Method 3: MA S^+(t1_plus, K1_plus)
  MA_M3 = fn_MA_plus(y, xe, xI, xo, n, t1_plus, K1_plus)
  S3 = MA_M3[,1]
  fstat3 = MA_M3[,2]
  dfR3 = MA_M3[,3]
  dfE3 = MA_M3[,4]
  pva3 = MA_M3[,5]
  Rsq3 = MA_M3[,6] 
  Exe_time3 = MA_M3[,7]
  beta_M3 = MA_M3[,8]
  SE_M3 = MA_M3[,9]
  pvalue_M3 = MA_M3[,10]
  
  #===== Define function for methods 2: MA on t IVs for K submodels, sampling with equal probability 
  fn_MA_ori = function(y, xe, xI, xo, n, t, K)  {
    q=dim(xI)[2]
    start_stageone = Sys.time()
    Xej_comb = matrix(0, nrow= n, ncol= K)    
    position = matrix(0, K, t)
    for (J in 1:K) { 
      count = 0
      repeat {
        set.seed(J+count+100)
        count = count + 1
        posJ = sample(1:q, t, replace=FALSE)     
        step = 0
        for (h in 1:J) {
          if (sum((sort(posJ)-sort(position[h,]))^2)!=0) {
            step = step + 1
          }
        }
        if (step == J) {
          break
        }
      }
      position[J,] = posJ
      fitJ = lm(xe~xI[,posJ] + 0)
      d.hatJ <- fitJ$fitted.values
      Xej_comb[,J] <- d.hatJ         
    }
    fit_SM = lm(xe~ Xej_comb + 0)         
    Xe_MA = fit_SM$fitted.values
    end_stageone = Sys.time()
    MA_exetime =  difftime(end_stageone, start_stageone, units="secs")
    num_position = as.numeric(position)
    S = length(unique(num_position))
    smy =summary(fit_SM)
    fstatma = smy$fstatistic[1]
    dfRma = smy$fstatistic[2]
    dfEma = smy$fstatistic[3]
    pvama = pf(fstatma, dfRma, dfEma, lower.tail = FALSE)
    Rsq = var(Xe_MA)/var(xe)
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-1) ) 
    psi = matrix(0, K, K)
    v_dk = matrix(0, t*K, 1)
    v_diag = matrix(0, t*K, K)
    sigma_inv = matrix(0, t*K, t*K)
    for (I in 1:K) {
      v_dk[((I*t-(t-1)):(I*t)),] = v[position[I,],]
      v_diag[((I*t-(t-1)):(I*t)), I] = v[position[I,],]
      sigma_II = omega_hat[position[I,], position[I,]]
      sigma_inv[((I*t-(t-1)):(I*t)), ((I*t-(t-1)):(I*t))] = Ginv(sigma_II, tol=1e-30)
      
      for (J in 1:K) {
        vi = v[position[I,], ]
        vj = v[position[J,], ]
        sigma_ii_inv = Ginv(omega_hat[position[I,], position[I,]], tol=1e-30)
        sigma_jj_inv = Ginv(omega_hat[position[J,], position[J,]], tol=1e-30)
        sigma_ij = omega_hat[position[I,],position[J,]]
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    SE_Md =  sigma_MA_hat/sqrt(n*t(u)%*%Ginv(psi, tol=1e-30)%*%u) 
    t_Md = beta_MA/SE_Md
    pvalue_Md = 2*pnorm(abs(t_Md), lower.tail = FALSE)
    result=matrix(c(S, fstatma, dfRma, dfEma, pvama, Rsq, MA_exetime, beta_MA, SE_Md, pvalue_Md), nrow =1)
    return(result)
  }
  
  #============================ Method 2: MA S^(t1, K1)
  MA_M4 = fn_MA_ori(y, xe, xI, xo, n, t1, K1)
  S4 = MA_M4[,1]
  fstat4 = MA_M4[,2]
  dfR4 = MA_M4[,3]
  dfE4 = MA_M4[,4]
  pva4 = MA_M4[,5]
  Rsq4 = MA_M4[,6]
  Exe_time4 = MA_M4[,7]
  beta_M4 = MA_M4[,8]
  SE_M4 = MA_M4[,9]
  pvalue_M4 = MA_M4[,10]
  #===== Define function for method 3 : MA on t IVs for K submodels, sampling with unequal probability
  fn_MA_plus_ori = function(y, xe, xI, xo, n, t, K)  {
    a = as.numeric(cor(xe,xI))
    a=abs(a)
    start_stageone = Sys.time()
    Xej_comb = matrix(0, nrow= n, ncol= K)    
    position = matrix(0, K, t)
    for (J in 1:K) { 
      count = 0
      repeat {
        set.seed(J+count+100)
        count = count + 1
        posJ = sample(1:q, t, replace=FALSE, prob=a)    
        step = 0
        for (h in 1:J) {
          if (sum((sort(posJ)-sort(position[h,]))^2)!=0) {
            step = step + 1
          }
        }
        if (step == J) {
          break
        }
      }
      position[J,] = posJ
      fitJ = lm(xe~xI[,posJ] + 0)
      d.hatJ <- fitJ$fitted.values
      Xej_comb[,J] <- d.hatJ         
    }
    fit_SM = lm(xe~ Xej_comb + 0)         
    Xe_MA = fit_SM$fitted.values
    end_stageone = Sys.time()
    MA_exetime =  difftime(end_stageone, start_stageone, units="secs")
    num_position = as.numeric(position)
    S = length(unique(num_position))
    smy =summary(fit_SM)
    fstatma = smy$fstatistic[1]
    dfRma = smy$fstatistic[2]
    dfEma = smy$fstatistic[3]
    pvama = pf(fstatma, dfRma, dfEma, lower.tail = FALSE)
    Rsq = var(Xe_MA)/var(xe)
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-1) ) 
    psi = matrix(0, K, K)
    v_dk = matrix(0, t*K, 1)
    v_diag = matrix(0, t*K, K)
    sigma_inv = matrix(0, t*K, t*K)
    for (I in 1:K) {
      v_dk[((I*t-(t-1)):(I*t)),] = v[position[I,],]
      v_diag[((I*t-(t-1)):(I*t)), I] = v[position[I,],]
      sigma_II = omega_hat[position[I,], position[I,]]
      sigma_inv[((I*t-(t-1)):(I*t)), ((I*t-(t-1)):(I*t))] = Ginv(sigma_II, tol=1e-30)
      
      for (J in 1:K) {
        vi = v[position[I,], ]
        vj = v[position[J,], ]
        sigma_ii_inv = Ginv(omega_hat[position[I,], position[I,]], tol=1e-30)
        sigma_jj_inv = Ginv(omega_hat[position[J,], position[J,]], tol=1e-30)
        sigma_ij = omega_hat[position[I,],position[J,]]
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    SE_Md =  sigma_MA_hat/sqrt(n*t(u)%*%Ginv(psi, tol=1e-30)%*%u) 
    t_Md = beta_MA/SE_Md
    pvalue_Md = 2*pnorm(abs(t_Md), lower.tail = FALSE)
    result=matrix(c(S, fstatma, dfRma, dfEma, pvama, Rsq, MA_exetime, beta_MA, SE_Md, pvalue_Md), nrow =1)
    return(result)
  }
  
  #============================ Method 3: MA S^+(t1_plus, K1_plus)
  MA_M5 = fn_MA_plus_ori(y, xe, xI, xo, n, t1_plus, K1_plus)
  S5 = MA_M5[,1]
  fstat5 = MA_M5[,2]
  dfR5 = MA_M5[,3]
  dfE5 = MA_M5[,4]
  pva5 = MA_M5[,5]
  Rsq5 = MA_M5[,6] 
  Exe_time5 = MA_M5[,7]
  beta_M5 = MA_M5[,8]
  SE_M5 = MA_M5[,9]
  pvalue_M5 = MA_M5[,10]
  
  #===== Define function for methods 4-6: 2SLS using IVs selected by lasso or elastic net in 1st stage
  fn_lasso = function(y, xe, xI, xo, n, anet) {
    start_stageone = Sys.time()
    p=dim(xo)[2]
    q=dim(xI)[2]
    cv.LASSO = cv.glmnet(xI,xe, intercept=FALSE, nfolds = 5, alpha=anet, family="gaussian" )
    LASSO1 = glmnet(xI, xe, lambda =  cv.LASSO$lambda.min, intercept=FALSE, alpha=anet, family="gaussian")  
    alphalasso = which(LASSO1$beta[1:q,]!=0)  
    fitlasso=lm(Xe~xI[,alphalasso]+xo+0)
    LASSOpred = fitlasso$fitted.values
    end_stageone = Sys.time()
    lasso_exetime =  difftime(end_stageone, start_stageone, units="secs")
    S =  length(alphalasso)
    smy = summary(fitlasso)
    fstatL = smy$fstatistic[1]
    dfRL = smy$fstatistic[2]
    dfEL = smy$fstatistic[3]
    pvaL = pf(fstatL, dfRL, dfEL, lower.tail = FALSE)
    Rsq_L = var(LASSOpred)/var(xe)
    fit_lasso = lm(y~LASSOpred + xo + 0)
    beta_lasso = fit_lasso$coefficients[1]  
    sigma_lasso_hat = sqrt( sum( (y-fit_lasso$fitted.values)^2 )/(n-1) )
    xI_LAS = xI[,alphalasso]
    xo_LAS = xo
    P1_LAS=xo_LAS%*%solve(t(xo_LAS)%*%xo_LAS, tol=1e-30)%*%t(xo_LAS)
    In=matrix(0, n, n)
    diag(In)=1
    zI_LAS=(In-P1_LAS)%*%xI_LAS
    M_Ie_tilde_LAS=(1/n)*t(zI_LAS)%*%xe
    M_II_tilde_LAS=(1/n)*t(zI_LAS)%*%zI_LAS
    SE_lasso = sigma_lasso_hat/sqrt(n*t(M_Ie_tilde_LAS)%*%solve(M_II_tilde_LAS, tol=1e-30)%*%M_Ie_tilde_LAS)
    t_ML = beta_lasso/SE_lasso
    pvalue_ML = 2*pnorm(abs(t_ML), lower.tail = FALSE)
    result=matrix(c(S, fstatL, dfRL, dfEL, pvaL, Rsq_L, lasso_exetime, beta_lasso, SE_lasso, pvalue_ML), nrow =1)
    return(result)
  }
  
  #============================ Method 4: pLasso
  lasso_M6 = fn_lasso(y, xe, xI, xo, n, anet1)
  S6 = lasso_M6[,1]
  fstat6 = lasso_M6[,2]
  dfR6 = lasso_M6[,3]
  dfE6 = lasso_M6[,4]
  pva6 = lasso_M6[,5]
  Rsq6 = lasso_M6[,6]
  Exe_time6 = lasso_M6[,7]
  beta_M6 = lasso_M6[,8]
  SE_M6 = lasso_M6[,9]
  pvalue_M6 = lasso_M6[,10]
  #============================ Method 5: pEL_0.5
  lasso_M7 = fn_lasso(y, xe, xI, xo, n, anetel1)
  S7 = lasso_M7[,1]
  fstat7 = lasso_M7[,2]
  dfR7 = lasso_M7[,3]
  dfE7 = lasso_M7[,4]
  pva7 = lasso_M7[,5]
  Rsq7 = lasso_M7[,6]
  Exe_time7 = lasso_M7[,7]
  beta_M7 = lasso_M7[,8]
  SE_M7 = lasso_M7[,9]
  pvalue_M7 = lasso_M7[,10]
  #============================ Method 6: pEL_{1/3}
  lasso_M8 = fn_lasso(y, xe, xI, xo, n, anetel2)
  S8 = lasso_M8[,1]
  fstat8 = lasso_M8[,2]
  dfR8 = lasso_M8[,3]
  dfE8 = lasso_M8[,4]
  pva8 = lasso_M8[,5]
  Rsq8 = lasso_M8[,6]
  Exe_time8 = lasso_M8[,7]
  beta_M8 = lasso_M8[,8]
  SE_M8 = lasso_M8[,9]
  pvalue_M8 = lasso_M8[,10]
  
  #=============================================================== report result
  result=matrix(c(S1, fstat1, dfR1, dfE1, pva1, Rsq1, Exe_time1, beta_M1, SE_M1, pvalue_M1, 
                  S2, fstat2, dfR2, dfE2, pva2, Rsq2, Exe_time2, beta_M2, SE_M2, pvalue_M2,
                  S3, fstat3, dfR3, dfE3, pva3, Rsq3, Exe_time3, beta_M3, SE_M3, pvalue_M3,     
                  S4, fstat4, dfR4, dfE4, pva4, Rsq4, Exe_time4, beta_M4, SE_M4, pvalue_M4, 
                  S5, fstat5, dfR5, dfE5, pva5, Rsq5, Exe_time5, beta_M5, SE_M5, pvalue_M5, 
                  S6, fstat6, dfR6, dfE6, pva6, Rsq6, Exe_time6, beta_M6, SE_M6, pvalue_M6,
                  S7, fstat7, dfR7, dfE7, pva7, Rsq7, Exe_time7, beta_M7, SE_M7, pvalue_M7,
                  S8, fstat8, dfR8, dfE8, pva8, Rsq8, Exe_time8, beta_M8, SE_M8, pvalue_M8
  ),ncol= 10, byrow = TRUE)
  return(result)
}

#========================================== run functions and obtain result
# set case = "ONE" (no adjustment for confounding covariates) or "TWO" (adjust for confounding covariates) in the structural equation
for(t in 5:15){
  for(K in 5:15){
result_full=fn_big(Y, Xe, XI, Xo, N,t,K)      
result_full=as.matrix(result_full)

result_full = rbind(result_full[2:5,], result_full[6:8,],t(as.matrix(result_full[1,])) )

rownames(result_full)=c( "MA", "MA_plus","MA_no_Xo","MA_plus_no_Xo", "pLasso", "pEN_0.5", "pEN_0.33","naive")
colnames(result_full)=c("XI size", "fstat", "dfR", "dfE", "pva", "Rsq", "Exe_time", "beta_e_hat", "SE","p-value" )
## Results in Table 7 
result_7 = result_full[,1:7]
## Results in Table 8
result_8 = result_full[,8:10]
if((result_full[2,10]<0.1)&((result_full[2,8]>0.06))){
  print(t)
  print(K)
  print(result_7[1,1])
  print(result_7[2,1])
print(result_8)}}}
