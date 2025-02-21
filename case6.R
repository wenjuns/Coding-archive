library(MASS)
Sys.setenv(RGL_USE_NULL=TRUE)
library(matlib)
library(corpcor)
library(MESS)
library(foreach)
library(doRNG)
library(doParallel)
library(gcookbook)
library(ggplot2)
library(qqplotr)
library(dplyr)
library(patchwork)
library(ggtext)
cl = makeCluster(detectCores()-1)
registerDoParallel(cl)

#Misspecified Case 1: one of the XIs is (XI)^2, another is sin(Xo) 

nsim=500
rho = 0.2
N=400
q=10                # number of IVs
p = 5               # number of exogeneous variables
sparse = 1*q 
## set number of IVs in every submodel and total number of submodels
K1 = 2           
t1 = 9
K1_plus = 2
t1_plus = 9
####==== coefficient of the endogenous variable Xe
be = -1           
####==== components of the error term e (set as Xu*au+e2) in the reduced form equation
au =  2             
epsilon2 = 1.3     
####==== components of the error term epsilon (set as bu*Xu + e1) in the structural form equation
bu = 1.5           
epsilon1 = 1       
####==== coefficients of XI will be sampled within the interval (0, b)
b = 3.5          
ao_c = 5
####==== set parameters for the covariance matrix of XI
sigmaXIsqfactor = 1
cs = matrix(rho, q, q)
diag(cs)=1
sigma_XI = sigmaXIsqfactor*cs
####==== set parameters for the covariance matrix of Xo
rho_Xo = 0
sigmaXosqfactor = 1
cs_Xo = matrix(rho_Xo, p, p)
diag(cs_Xo)=1
sigma_Xo = sigmaXosqfactor*cs_Xo
####==== Creating matrices to store simulation output
Non_zero_alpha = matrix(0, nrow=nsim, ncol= 8)
allexetime = matrix(0, nrow=nsim, ncol= 8)
B_HAT = matrix(0, nrow=nsim, ncol= 8)
Bias_beta = matrix(0, nrow=nsim, ncol= 8)
Beta_no_bias1 = matrix(0, nrow=nsim, ncol= 8)
Beta_no_bias2 = matrix(0, nrow=nsim, ncol= 8)
Beta_debiased = matrix(0, nrow=nsim, ncol= 8)
SE = matrix(0, nrow=nsim, ncol= 8) 
within_CI_95 = matrix(0, nrow=1, ncol= 8)
#n=500



fn_big=function(y, xe, xI, xo, n)
{ 
  a = as.numeric(cor(xe,xI))
  a=abs(a)
  xo_xo_inv=solve(t(xo)%*%xo,tol=1e-30)
  P1=xo%*%xo_xo_inv%*%t(xo)
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
  a_bias=(1/n)*t(xe)%*%xo%*%xo_xo_inv%*%t(xo)%*%xe
  v=(1/n)*t(xI)%*%xe
  #============================ Method 1: Naive
  start_stageone = Sys.time()
  S1 = 0
  fit_naive = lm(y~xe + xo + 0)
  beta_M1 = fit_naive$coefficients[1]
  bias_beta1=0
  beta_nobias1=beta_M1-bias_beta1
  beta_debiased1=beta_M1
  out1 = summary(fit_naive)
  SE_M1 = out1$coefficients[1,2]      
  within_CI1 =  ifelse( (beta_M1 - 1.96*SE_M1)<be & be< (beta_M1 + 1.96*SE_M1) ,1,0)
  Exe_time1 = 0
  #============================ Method 2: classical 2SLS (same as MA S^(q, 1))
  start_stageone = Sys.time()
  fit = lm(xe~xI+ xo  + 0)
  d.hatss <- fit$fitted.values
  end_stageone = Sys.time()
  Exe_time2 =  difftime(end_stageone, start_stageone, units="secs")
  p = dim(xo)[2]
  S2 = dim(xI)[2]
  fit_2sls = lm(y~d.hatss + xo + 0)
  beta_M2 = fit_2sls$coefficients[1]
  bias_beta2=0
  beta_nobias2=beta_M2-bias_beta2
  beta_debiased2=beta_M2
  sigma_2sls_hat = sqrt( sum( (y-fit_2sls$fitted.values)^2 )/(n-p-1) )
  M_II_tilde_inv=solve(M_II_tilde,tol=1e-30)
  SE_M2 =  sigma_2sls_hat/sqrt(n*t(M_Ie_tilde)%*% M_II_tilde_inv%*%M_Ie_tilde)
  if (is.na(SE_M2))  {
    within_CI2  = 0
  } else {
    within_CI2  = ifelse((beta_M2 - 1.96*SE_M2) <be & be< (beta_M2 + 1.96*SE_M2),1,0)
  }
  
  
  #==================== estimate sigma^2 using LASSO
  fn_lasso_sigma = function(y, xe, xI,xo, n, anet) {
    p=dim(xo)[2]
    start_stageone = Sys.time()
    zI_all= cbind(xo,xI)
    cv.LASSO = cv.glmnet(zI_all,xe, intercept=FALSE, nfolds = 5, alpha=anet, family="gaussian" )
    LASSO1 = glmnet(zI_all, xe, lambda =  cv.LASSO$lambda.min, intercept=FALSE, alpha=anet, family="gaussian")  
    alphalasso = which(LASSO1$beta[1:(p+q),]!=0)
    xo_lasso_ind=which(LASSO1$beta[1:p,]!=0)
    xo_lasso=zI_all[,xo_lasso_ind]
    xI_lasso_ind=which(LASSO1$beta[(1+p):(p+q),]!=0)
    fitlasso=lm(xe~zI_all[,alphalasso]+0)
    LASSOpred = fitlasso$fitted.values
    end_stageone = Sys.time()
    lasso_exetime =  difftime(end_stageone, start_stageone, units="secs")
    S =  length(which(LASSO1$beta[1:(p+q),]!=0))-dim(xo_lasso)[2]
    smy = summary(fitlasso)
    fstatL = smy$fstatistic[1]
    dfRL = smy$fstatistic[2]
    dfEL = smy$fstatistic[3]
    pvaL = pf(fstatL, dfRL, dfEL, lower.tail = FALSE)
    Rsq_L = var(LASSOpred)/var(xe)
    fit_lasso = lm(y~LASSOpred + xo_lasso + 0)
    beta_lasso = fit_lasso$coefficients[1]  
    sigma_lasso_hat = sqrt( sum( (y-fit_lasso$fitted.values)^2 )/(n-1-dim(xo_lasso)[2]) )
    xI_LAS = xI[,xI_lasso_ind]
    xo_LAS = xo_lasso
    P1_LAS=xo_LAS%*%solve(t(xo_LAS)%*%xo_LAS, tol=1e-30)%*%t(xo_LAS)
    In=matrix(0, n, n)
    diag(In)=1
    zI_LAS=(In-P1_LAS)%*%xI_LAS
    M_Ie_tilde_LAS=(1/n)*t(zI_LAS)%*%xe
    M_II_tilde_LAS=(1/n)*t(zI_LAS)%*%zI_LAS
    SE_lasso = sigma_lasso_hat/sqrt(n*t(M_Ie_tilde_LAS)%*%solve(M_II_tilde_LAS, tol=1e-30)%*%M_Ie_tilde_LAS)
    if (is.na(SE_lasso))  {
      within_CI_lasso = 0
    } else {
      within_CI_lasso = ifelse((beta_lasso - 1.96*SE_lasso) <be & be< (beta_lasso + 1.96*SE_lasso),1,0)
    }
    result=matrix(c(S, lasso_exetime, beta_lasso,  SE_lasso, within_CI_lasso,sigma_lasso_hat), nrow =1)
    return(result)
  }
  #==================== Method 5: pLasso
  lasso_M5 = fn_lasso_sigma(y, xe, xI,xo, n, 1)
  S7 = lasso_M5[,1]
  Exe_time7 = lasso_M5[,2]
  beta_M7 = lasso_M5[,3]
  bias_beta7=0
  beta_nobias7=beta_M7
  beta_debiased7=beta_M7
  SE_M7 = lasso_M5[,4]
  within_CI7 = lasso_M5[,5]
  sigma_est=sigma_2sls_hat
  #===== Define function for Method 3: MA on t IVs for K submodels, sampling with equal probability
  fn_MA = function(y, xe, xI, xo, n, t, K)  {
    #t=20
    #K=10
    #estimate of sigma^2 with the 2sls method
    p=dim(xo)[2]
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
      fitJ = lm(xe~xI[,posJ]+xo + 0)
      d.hatJ <- fitJ$fitted.values
      Xej_comb[,J] <- d.hatJ         
    }
    fit_SM = lm(xe~ Xej_comb + 0)         
    Xe_MA = fit_SM$fitted.values
    end_stageone = Sys.time()
    MA_exetime =  difftime(end_stageone, start_stageone, units="secs")
    num_position = as.numeric(position)
    S = length(unique(num_position))
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-p-1) ) 
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
      sigma_inv[(I*(t+p)-(t+p-1)):(I*(t+p)), (I*(t+p)-(t+p-1)):(I*(t+p))] = solve(sigma_II,tol=1e-30)
      
      
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
        sigma_ii_inv = solve(sigma_ii,tol=1e-30)
        sigma_jj_inv = solve(sigma_jj,tol=1e-30)
        sigma_ij = rbind(cbind(sigma_11,sigma_12_j),cbind(sigma_21_i,sigma_22_ij))
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    d = t(v_diag)%*%sigma_inv%*%d_dk%*%solve(Moo, tol=1e-30)%*%Moe
    psi_inv=solve(psi,tol=1e-30)
    c_bias=a_bias*t(u)%*%psi_inv%*%Ik_vector*(t(u)%*%psi_inv%*%Ik_vector-1)/(t(u)%*%psi_inv%*%u-a_bias*(t(u)%*%psi_inv%*%Ik_vector)^2)
    bias_beta1=c_bias* beta_MA
    beta_hat1= beta_MA-bias_beta1
    bias_beta2=c_bias* beta_hat1
    beta_hat2=beta_MA-bias_beta2
    beta_hat_debiased=beta_MA/(1+c_bias)
    nu=sigma_est*sqrt(as.numeric(t(u)%*%psi_inv%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%psi_inv%*%u))/(abs(t(u)%*%psi_inv%*%(u-as.numeric(a_bias)*Ik_vector))*sqrt(n))
    #nu=sqrt(3.25)*sqrt(as.numeric(t(u)%*%ginv(psi)%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%ginv(psi)%*%u))/(abs(t(u)%*%ginv(psi)%*%u-a_bias*(t(u)%*%ginv(psi)%*%Ik_vector)^2)*sqrt(n))
    if (is.na(nu))  {
      within_CI3  = 0
    } else {
      within_CI3  = ifelse((beta_hat2 - 1.96*nu) <be & be< (beta_hat2 + 1.96*nu),1,0)
    }
    #beta_nobias=beta_MA-bias_beta
    #epsilon_hat=(beta_nobias-be)/nu
    #=====================report result
    result=matrix(c(S, MA_exetime, beta_MA, bias_beta2,beta_hat1,beta_hat2,beta_hat_debiased,nu, within_CI3), nrow =1)
    return(result)
  }
  #============================ Method 3: BIC-selected MA S^(t1, K1)
  MA_M3 = fn_MA(y, xe, xI, xo, n, t1, K1)
  S3 = MA_M3[,1]
  Exe_time3 = MA_M3[,2]
  beta_M3 = MA_M3[,3]
  bias_beta3=MA_M3[,4]
  beta_hat1_3=MA_M3[,5]
  beta_hat2_3=MA_M3[,6]
  beta_debiased3=MA_M3[,7]
  SE_M3 = MA_M3[,8]
  within_CI3 = MA_M3[,9]
  #===== Define function for Method 4 : MA on t IVs for K submodels, sampling with unequal probability
  fn_MA_plus = function(y, xe, xI, xo, n, t, K)  {
    p=dim(xo)[2]
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
      fitJ = lm(xe~xI[,posJ]+xo + 0)
      d.hatJ <- fitJ$fitted.values
      Xej_comb[,J] <- d.hatJ         
    }
    fit_SM = lm(xe~ Xej_comb + 0)         
    Xe_MA = fit_SM$fitted.values
    end_stageone = Sys.time()
    MA_exetime =  difftime(end_stageone, start_stageone, units="secs")
    num_position = as.numeric(position)
    S = length(unique(num_position))
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-p-1) ) 
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
      sigma_inv[(I*(t+p)-(t+p-1)):(I*(t+p)), (I*(t+p)-(t+p-1)):(I*(t+p))] = solve(sigma_II,tol=1e-30)
      
      
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
        sigma_ii_inv = solve(sigma_ii,tol=1e-30)
        sigma_jj_inv = solve(sigma_jj,tol=1e-30)
        sigma_ij = rbind(cbind(sigma_11,sigma_12_j),cbind(sigma_21_i,sigma_22_ij))
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    d = t(v_diag)%*%sigma_inv%*%d_dk%*%solve(Moo, tol=1e-30)%*%Moe
    psi_inv=solve(psi,tol=1e-30)
    c_bias=a_bias*t(u)%*%psi_inv%*%Ik_vector*(t(u)%*%psi_inv%*%Ik_vector-1)/(t(u)%*%psi_inv%*%u-a_bias*(t(u)%*%psi_inv%*%Ik_vector)^2)
    bias_beta1=c_bias* beta_MA
    beta_hat1= beta_MA-bias_beta1
    bias_beta2=c_bias* beta_hat1
    beta_hat2=beta_MA-bias_beta2
    beta_hat_debiased=beta_MA/(1+c_bias)
    nu=sigma_est*sqrt(as.numeric(t(u)%*%psi_inv%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%psi_inv%*%u))/(abs(t(u)%*%psi_inv%*%(u-as.numeric(a_bias)*Ik_vector))*sqrt(n))
    #nu=sqrt(3.25)*sqrt(as.numeric(t(u)%*%ginv(psi)%*%(psi-d%*%t(Ik_vector)-Ik_vector%*%t(d)+as.numeric(a_bias)*(Ik_vector%*%t(Ik_vector)))%*%ginv(psi)%*%u))/(abs(t(u)%*%ginv(psi)%*%u-a_bias*(t(u)%*%ginv(psi)%*%Ik_vector)^2)*sqrt(n))
    if (is.na(nu))  {
      within_CI3  = 0
    } else {
      within_CI3  = ifelse((beta_hat2 - 1.96*nu) <be & be< (beta_hat2 + 1.96*nu),1,0)
    }
    #beta_nobias=beta_MA-bias_beta
    #epsilon_hat=(beta_nobias-be)/nu
    #=====================report result
    result=matrix(c(S, MA_exetime, beta_MA, bias_beta2,beta_hat1,beta_hat2,beta_hat_debiased,nu, within_CI3), nrow =1)
    return(result)
  }
  #============================ Method 3: BIC-selected MA S^(t1, K1)
  MA_M4 = fn_MA_plus(y, xe, xI, xo, n, t1, K1)
  S4 = MA_M4[,1]
  Exe_time4 = MA_M4[,2]
  beta_M4 = MA_M4[,3]
  bias_beta4=MA_M4[,4]
  beta_hat1_4=MA_M4[,5]
  beta_hat2_4=MA_M4[,6]
  beta_debiased4=MA_M4[,7]
  SE_M4 = MA_M4[,8]
  within_CI4 = MA_M4[,9]
  
  #===== Define function for Method 5: MA on t IVs for K submodels, sampling with equal probability
  fn_MA5 = function(y, xe, xI, xo, n, t, K)  {
    start_stageone = Sys.time()
    p=dim(xo)[2]
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
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-p-1) ) 
    psi = matrix(0, K, K)
    v_dk = matrix(0, t*K, 1)
    v_diag = matrix(0, t*K, K)
    sigma_inv = matrix(0, t*K, t*K)
    for (I in 1:K) {
      v_dk[((I*t-(t-1)):(I*t)),] = v[position[I,],]
      v_diag[((I*t-(t-1)):(I*t)), I] = v[position[I,],]
      sigma_II = omega_hat[position[I,], position[I,]]
      sigma_inv[((I*t-(t-1)):(I*t)), ((I*t-(t-1)):(I*t))] = solve(sigma_II, tol=1e-30)
      
      for (J in 1:K) {
        vi = v[position[I,], ]
        vj = v[position[J,], ]
        sigma_ii_inv = solve(omega_hat[position[I,], position[I,]], tol=1e-30)
        sigma_jj_inv = solve(omega_hat[position[J,], position[J,]], tol=1e-30)
        sigma_ij = omega_hat[position[I,],position[J,]]
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    SE_Md =  sigma_MA_hat/sqrt(n*t(u)%*%solve(psi, tol=1e-30)%*%u) 
    if (is.na(SE_Md))  {
      within_CI_MA  = 0 
    } else {
      within_CI_MA  = ifelse((beta_MA - 1.96*SE_Md) <be & be< (beta_MA + 1.96*SE_Md),1,0)
    }
    
    #=====================report result
    result=matrix(c(S, MA_exetime, beta_MA,SE_Md,within_CI_MA), nrow =1)
    return(result)
  }
  #============================ Method 5: BIC-selected MA S^(t1, K1)
  MA_M5 = fn_MA5(y, xe, xI, xo, n, t1, K1)
  S5 = MA_M5[,1]
  Exe_time5 = MA_M5[,2]
  beta_M5 = MA_M5[,3]
  bias_beta5=0
  beta_nobias5=beta_M5
  beta_debiased5=beta_M5
  SE_M5 = MA_M5[,4]
  within_CI5 = MA_M5[,5]
  #===== Define function for Method 6 : MA on t IVs for K submodels, sampling with unequal probability
  fn_MA_plus6 = function(y, xe, xI, xo, n, t, K)  {
    start_stageone = Sys.time()
    p=dim(xo)[2]
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
    fit_MA = lm(y~Xe_MA + xo + 0)
    beta_MA = fit_MA$coefficients[1]
    sigma_MA_hat = sqrt( sum( (y-fit_MA$fitted.values)^2 )/(n-p-1) ) 
    psi = matrix(0, K, K)
    v_dk = matrix(0, t*K, 1)
    v_diag = matrix(0, t*K, K)
    sigma_inv = matrix(0, t*K, t*K)
    for (I in 1:K) {
      v_dk[((I*t-(t-1)):(I*t)),] = v[position[I,],]
      v_diag[((I*t-(t-1)):(I*t)), I] = v[position[I,],]
      sigma_II = omega_hat[position[I,], position[I,]]
      sigma_inv[((I*t-(t-1)):(I*t)), ((I*t-(t-1)):(I*t))] = solve(sigma_II, tol=1e-30)
      
      for (J in 1:K) {
        vi = v[position[I,], ]
        vj = v[position[J,], ]
        sigma_ii_inv = solve(omega_hat[position[I,], position[I,]], tol=1e-30)
        sigma_jj_inv = solve(omega_hat[position[J,], position[J,]], tol=1e-30)
        sigma_ij = omega_hat[position[I,],position[J,]]
        psi[I,J] = t(vi)%*%sigma_ii_inv%*%sigma_ij%*%sigma_jj_inv%*%vj
      }  }
    u = t(v_diag)%*%sigma_inv%*%v_dk
    SE_Md =  sigma_MA_hat/sqrt(n*t(u)%*%solve(psi, tol=1e-30)%*%u) 
    if (is.na(SE_Md))  {
      within_CI_MA  = 0 
    } else {
      within_CI_MA  = ifelse((beta_MA - 1.96*SE_Md) <be & be< (beta_MA + 1.96*SE_Md),1,0)
    }
    result=matrix(c(S, MA_exetime, beta_MA,SE_Md,within_CI_MA), nrow =1)
    return(result)
  }
  #============================ Method 6: BIC-selected MA S^+(t1_plus, K1_plus)
  MA_M6 = fn_MA_plus6(y, xe, xI, xo, n, t1_plus, K1_plus)
  S6 = MA_M6[,1]
  Exe_time6 = MA_M6[,2]
  beta_M6 = MA_M6[,3]
  bias_beta6=0
  beta_nobias6=beta_M6
  beta_debiased6=beta_M6
  SE_M6 = MA_M6[,4]
  within_CI6 = MA_M6[,5]
  
  #==================== Method 6: pEL_0.5
  lasso_M6 = fn_lasso_sigma(y, xe, xI,xo, n, 0.5)
  S8 = lasso_M6[,1]
  Exe_time8 = lasso_M6[,2]
  beta_M8 = lasso_M6[,3]
  bias_beta8=0
  beta_nobias8=beta_M8
  beta_debiased8=beta_M8
  SE_M8 = lasso_M6[,4]
  within_CI8 = lasso_M6[,5]
  
  epsilon_hat1=0
  epsilon_hat2=0
  epsilon_hat5=0
  epsilon_hat6=0
  epsilon_hat7=0
  epsilon_hat8=0
  
  #=============================================================== report result
  result=matrix(c(S1, Exe_time1, beta_M1,  bias_beta1, beta_nobias1,beta_nobias1,beta_debiased1,SE_M1,within_CI1,
                  S2, Exe_time2, beta_M2,  bias_beta2, beta_nobias2,beta_nobias2,beta_debiased2,SE_M2,within_CI2,
                  S3, Exe_time3, beta_M3,  bias_beta3, beta_hat1_3,beta_hat2_3,beta_debiased3,SE_M3,within_CI3,
                  S4, Exe_time4, beta_M4,  bias_beta4, beta_hat1_4,beta_hat2_4,beta_debiased4,SE_M4,within_CI4,
                  S5, Exe_time5, beta_M5,  bias_beta5, beta_nobias5,beta_nobias5,beta_debiased5,SE_M5,within_CI5,
                  S6, Exe_time6, beta_M6,  bias_beta6, beta_nobias6,beta_nobias6,beta_debiased6,SE_M6,within_CI6,
                  S7, Exe_time7, beta_M7,  bias_beta7, beta_nobias7,beta_nobias7,beta_debiased7,SE_M7,within_CI7,
                  S8, Exe_time8, beta_M8,  bias_beta8, beta_nobias8,beta_nobias8,beta_debiased8,SE_M8,within_CI8
  ),ncol= 9, byrow = TRUE)
  return(result)
}


fcn_one <- function(n,al_o) {
  XI = mvrnorm(n, mu=replicate(q, 0), sigma_XI)
  Xo = mvrnorm(n, mu=replicate(p, 0), sigma_Xo)
  true_XI=XI
  true_Xo=Xo
  nonlinear_XI5=(XI[,5])^2
  nonlinear_XI4=sin(XI[,4])
  nonlinear_Xo5=exp((Xo[,5]))
  XI[,5]=nonlinear_XI5
  XI[,4]=nonlinear_XI4
  Xo[,5]=nonlinear_Xo5
  Xu = rnorm(n, 0, 1)
  e1=rnorm(n, mean=0, sd = epsilon1)
  e2=rnorm(n, mean=0, sd = epsilon2)
  aI = rep(0, q)
  pos = sample(1:q, sparse, replace=FALSE)
  aI[pos] = sample(seq(0, b, by = 0.1), size=sparse, replace=TRUE)
  bo = sample(c(-5,-4,-3,-2,-1,1,2,3,4,5),size=p,replace=TRUE)  
  ao = runif(p, -ao_c, ao_c)
  #ao=rep(5,p)
  Xe = XI%*%aI  + Xu*au + Xo%*%ao + e2                # Xu*au & e2 make up the error term, e, in the reduced form equation
  Y1 = be*Xe + Xo%*%bo  + bu*Xu + e1        # bu*Xu & e1 make up the error term, epsilon, in the structural equation
  ## Note: sigma^2 in the report equals to sum of bu^2 & epsilon1^2; 
  ## Note: tau^2 in the report equals to sum of au^2 & epsilon2^2; 
  RESULT = fn_big(Y1, Xe, true_XI, true_Xo, n)
  return(RESULT)
}

ite_result <- foreach(f=rep(N, nsim),  .packages=c("MASS", "matlib", "MESS","glmnet"))  %dorng%  fcn_one(n=f)

for (i in 1:nsim) {
  Non_zero_alpha[i,] = t(ite_result[[i]])[1,]
  allexetime[i,] =  t(ite_result[[i]])[2,]
  B_HAT[i,] =  t(ite_result[[i]])[3,]  
  Bias_beta[i,] = t(ite_result[[i]])[4,] 
  Beta_no_bias1[i,] = t(ite_result[[i]])[5,] 
  Beta_no_bias2[i,] = t(ite_result[[i]])[6,]
  Beta_debiased[i,] = t(ite_result[[i]])[7,]
  SE[i,]=t(ite_result[[i]])[8,] 
  within_CI_95[1,] = within_CI_95[1,] + t(ite_result[[i]])[9,]
}


final_sim = rbind(colMeans(Beta_no_bias1)-be,
                  colMeans(Beta_no_bias2)-be,
                  colMeans(Beta_debiased)-be,
                  apply(Beta_no_bias1, 2, sd),
                  apply(Beta_no_bias2, 2, sd),
                  apply(Beta_debiased, 2, sd),
                  colMeans(SE, na.rm = TRUE),
                  within_CI_95/(nsim),
                  colMeans(Non_zero_alpha),
                  colMeans(allexetime))

options(scipen = 200)
final_sim = cbind(final_sim[,2:8],final_sim[,1])
colnames(final_sim)=c("2SLS","MA", "MA_plus","MA_no_Xo","MA_plus_no_Xo","Lasso","elastic_net","Naive")
rownames(final_sim)=c("Bias1","Bias2","Beta_debiased", "SD1", "SD2","SD3","SE", "CP", "XI_size", "T")
round(final_sim,digits = 4)
#qqplot
eps_hat2=Beta_no_bias[,2]
dat_plt2=data.frame(eps_hat2)
g2 = ggplot(dat_plt2,aes(sample=eps_hat2))+
  theme_bw() +
  theme(axis.text = element_text(size = 25),axis.title = element_text(size = 20), plot.caption = element_text(size=25,hjust=0.5))+
  geom_qq_line(size=0.5) +
  geom_qq() +
  labs(caption = "(a) 2SLS", x = "Theoretical Quantiles", y = "Sample Quantiles")
g2
eps_hat3=Beta_no_bias[,3]
dat_plt3=data.frame(eps_hat3)
g3 =ggplot(dat_plt3,aes(sample=eps_hat3))+
  theme_bw() +
  theme(axis.text = element_text(size = 25),axis.title = element_text(size = 20), plot.caption = element_text(size=25,hjust=0.5))+
  geom_qq_line(size=0.5) +
  geom_qq() +
  labs(caption = bquote('(b) MA'^'(t,M)'), x = "Theoretical Quantiles", y = "Sample Quantiles")


eps_hat4=Beta_no_bias[,4]
dat_plt4=data.frame(eps_hat4)
g4 =ggplot(dat_plt4,aes(sample=eps_hat4))+
  theme_bw() +
  theme(axis.text = element_text(size = 25),axis.title = element_text(size = 20), plot.caption = element_text(size=25,hjust=0.5))+
  geom_qq_line(size=0.5) +
  geom_qq() +
  labs(caption = bquote('(c) MA'^'+(t,M)'), x = "Theoretical Quantiles", y = "Sample Quantiles")

eps_hat5=Beta_no_bias[,5]
dat_plt5=data.frame(eps_hat5)
g5 =ggplot(dat_plt5,aes(sample=eps_hat5))+
  theme_bw() +
  theme(axis.text = element_text(size = 25),axis.title = element_text(size = 20), plot.caption = element_text(size=25,hjust=0.5))+
  geom_qq_line(size=0.5) +
  geom_qq() +
  labs(caption = bquote('(d) s'^'(t,M)'), x = "Theoretical Quantiles", y = "Sample Quantiles")

eps_hat6=Beta_no_bias[,6]
dat_plt6=data.frame(eps_hat6)
g6 =ggplot(dat_plt6,aes(sample=eps_hat6))+
  theme_bw() +
  theme(axis.text = element_text(size = 25),axis.title = element_text(size = 20), plot.caption = element_text(size=25,hjust=0.5))+
  geom_qq_line(size=0.5) +
  geom_qq() +
  labs(caption = bquote('(e) s'^'+(t,M)'), x = "Theoretical Quantiles", y = "Sample Quantiles")

eps_hat1=Beta_no_bias[,1]
dat_plt1=data.frame(eps_hat1)
g1 =ggplot(dat_plt1,aes(sample=eps_hat1))+
  theme_bw() +
  theme(axis.text = element_text(size = 25),axis.title = element_text(size = 20), plot.caption = element_text(size=25,hjust=0.5))+
  geom_qq_line(size=0.5) +
  geom_qq() +
  labs(caption = bquote('(f) Naive'), x = "Theoretical Quantiles", y = "Sample Quantiles")

g2+g3+g4+g5+g6+g1

#a <- eps_hat
#t <- (rank(a) -0.5)/length(a)
#q <- qnorm(t)
#plot(q, a)
#abline(mean(a), sd(a), col=2, lwd=2)
