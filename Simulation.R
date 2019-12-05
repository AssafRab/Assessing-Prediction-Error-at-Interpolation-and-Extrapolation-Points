##The following code was sent to a cluster of computers. So for getting the results,
#one should run this code X times.

library(lme4)
library(MASS)
library(reshape)
library(dplyr)

### Parameters

{
  Ve<-c(15,20,25)
  Vb<-c(15,1^2)
  
  beta0<-1
  beta1<-1
  beta2<-1
  beta3<-2
  beta4<-2
  beta5<-2
  beta6<-2
  
  betatime<-0.5
  beta_int<-rbind(beta1,beta2,beta3,beta4,beta5,beta6)

  meanX1<-0.5
  meanX2<-0
  sdX2<-sqrt(1)
  meanX3<-0
  sdX3<-sqrt(1)
  meanX4<-0
  sdX4<-sqrt(1)
  meanX5<-0
  sdX5<-sqrt(1)
  meanX6<-0
  sdX6<-sqrt(1)
  
  Var_cov<-cbind(meanX1,sdX2,sdX3,sdX4,sdX5,sdX6)^2
  }

#Sample size and models setting
{
  Nsubjects<-c(20,100,200)
  
  Nobs_n<-10
  Nobs_m<-2
  Nobs_all<-Nobs_n+Nobs_m
  TestingTimes<-Nobs_n+5*(1:Nobs_m)
  
  length_LLN<-200
  
  models<-list(1:2,1:4,1:6)
  }

### Defining the output matrices
{
  loglike_f<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  loglike_f_est<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  Eloglike_star_f<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  
  C_tAI<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  C_cAI<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  C_AI<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  
  C_tAI_est<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  C_cAI_est<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
  C_AI_est<-array(data=NA,dim=c(length(Nsubjects),length(Ve),length(models)))
}

for (i in 1:length(Nsubjects)) {
  
  {  
    sample_size_n<-Nobs_n*Nsubjects[i]
    sample_size_m<-Nobs_m*Nsubjects[i]
    sample_size_all<-Nobs_all*Nsubjects[i]
    
    Sampling_matrix<-matrix(data=0,nrow=sample_size_all,ncol=Nsubjects[i])
    for(l in 1:Nsubjects[i]){
      rowindex1<-(1+(l-1)*Nobs_all)
      rowindex2<-l*Nobs_all
      
      Sampling_matrix[rowindex1:rowindex2,l]<-rep(1,Nobs_all)
    }
    
    Follow_up_index<-rep(c(rep(0,Nobs_n),rep(1,Nobs_m)),Nsubjects[i])   
  }
  
  # Creating Z 
  {
    Z<-matrix(data=0,nrow=sample_size_all,ncol=2*Nsubjects[i])
    
    for(l in 1:Nsubjects[i])  {
      rowindex1<-(1+(l-1)*Nobs_all)
      rowindex2<-l*Nobs_all
      
      Z[rowindex1:rowindex2,2*l-1]<-rep(1,Nobs_all)
      Z[rowindex1:rowindex2,2*l]<-c(1:Nobs_n,TestingTimes)  
    }
    
    Z_n<-Z[Follow_up_index==0,]
    Z_m<-Z[Follow_up_index==1,]
    
    time<-rep(Z[1:Nobs_all,2],Nsubjects[i])
    
    time_n<-time[Follow_up_index==0]
    time_m<-time[Follow_up_index==1]
    
  }
  
  # Creating X
  {
    X0<-rep(1,Nsubjects[i])
    X1<-rbinom(Nsubjects[i],1,meanX1)
    X2<-rnorm(Nsubjects[i],meanX2,sdX2)
    X3<-rnorm(Nsubjects[i],meanX3,sdX3)
    X4<-rnorm(Nsubjects[i],meanX4,sdX4)
    X5<-rnorm(Nsubjects[i],meanX5,sdX5)
    X6<-rnorm(Nsubjects[i],meanX6,sdX6)
   
    Covariates<-cbind(X1,X2,X3,X4,X5,X6)
  }
  
  
  beta<-rbind(beta0,beta_int,betatime)
  X_n_true<-cbind(Sampling_matrix[Follow_up_index==0,]%*%cbind(X0,Covariates[,1:(length(beta_int))]),time_n)
  X_m_true<-cbind(Sampling_matrix[Follow_up_index==1,]%*%cbind(X0,Covariates[,1:(length(beta_int))]),time_m)
  
  for (k in 1:length(Ve)){
    
    ## Creating y
    {
      e<-rnorm(sample_size_n,mean=0,sd=sqrt(Ve[k])) 
      b <-rnorm(ncol(Z),0,sd=sqrt(Vb))
      
      y<-X_n_true%*%beta+Z_n%*%b+e
    }
    
    
    ## Allocating H and V for the loop below
    {
      R_n<-array(data=NA,dim=c(Nobs_n*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      R_m<-array(data=NA,dim=c(Nobs_m*Nsubjects[i], Nobs_m*Nsubjects[i],length(models)))
      V_n<-array(data=NA,dim=c(Nobs_n*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      V_m<-array(data=NA,dim=c(Nobs_m*Nsubjects[i], Nobs_m*Nsubjects[i],length(models)))
      H_n<-array(data=NA,dim=c(Nobs_n*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      H_m<-array(data=NA,dim=c(Nobs_m*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      
      R_n_est<-array(data=NA,dim=c(Nobs_n*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      R_m_est<-array(data=NA,dim=c(Nobs_m*Nsubjects[i], Nobs_m*Nsubjects[i],length(models)))
      V_n_est<-array(data=NA,dim=c(Nobs_n*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      V_m_est<-array(data=NA,dim=c(Nobs_m*Nsubjects[i], Nobs_m*Nsubjects[i],length(models)))
      H_n_est<-array(data=NA,dim=c(Nobs_n*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
      H_m_est<-array(data=NA,dim=c(Nobs_m*Nsubjects[i], Nobs_n*Nsubjects[i],length(models)))
    }
    

    for(l in 1:length(models))  {
      
      max_beta<-max(unlist(models[l]))
      
      X_n<-cbind(Sampling_matrix[Follow_up_index==0,]%*%cbind(X0,Covariates[,1:max_beta]),time_n)
      X_m<-cbind(Sampling_matrix[Follow_up_index==1,]%*%cbind(X0,Covariates[,1:max_beta]),time_m)
      
      ##Variances
      
      #Mispecification variance
      {
        length_extra<-max_beta-length(beta_int)
        if(length_extra>=0){
          ExtraVar<-0
        } else{
          ExtraVar<-Var_cov[(max_beta+1):length(beta_int)]%*%(beta_int[(max_beta+1):length(beta_int)]^2)
        } 
      }
      
      {
        # The ExtraVar goes to random intercept, this is exactly the idea of random intercept  
        G<-diag(rep(Vb+c(ExtraVar,0),times=Nsubjects[i]))
        
        ZGZ_n<-Z_n%*%G%*%t(Z_n)
        ZGZ_m<-Z_m%*%G%*%t(Z_m)
        
        V_n[,,l]<-(Ve[k]*diag(sample_size_n)+ZGZ_n)
        R_n[,,l] <- Ve[k]*diag(sample_size_n)
        
        V_m[,,l]<-(Ve[k]*diag(sample_size_m)+ZGZ_m)
        R_m[,,l] <-V_m[,,l] - (Z_m%*%G%*%t(Z_n)%*%chol2inv(chol(V_n[,,l]))%*%Z_n%*%G%*%t(Z_m))
      }
      
      ## Hat matrices
      {
        V_n_inv<-chol2inv(chol(V_n[,,l]))
        H_n[,,l]<-(X_n%*%chol2inv(chol(t(X_n)%*%V_n_inv%*%X_n))%*%t(X_n)%*%V_n_inv+
                     Z_n%*%G%*%t(Z_n)%*%V_n_inv%*%(diag(1,sample_size_n)
                                                   -X_n%*%chol2inv(chol(t(X_n)%*%V_n_inv%*%X_n))%*%t(X_n)%*%V_n_inv))
        
        H_m[,,l]<-(X_m%*%chol2inv(chol(t(X_n)%*%V_n_inv%*%X_n))%*%t(X_n)%*%V_n_inv+
                     Z_m%*%G%*%t(Z_n)%*%V_n_inv%*%(diag(1,sample_size_n)
                                                   -X_n%*%chol2inv(chol(t(X_n)%*%V_n_inv%*%X_n))%*%t(X_n)%*%V_n_inv))
      }
      
      # Corrections. (for loss(Opt_t) change C_tAI to w_t, C_cAI to w and delete C_AI)
      {
        
        R_n_inv<-chol2inv(chol(R_n[,,l]))
        R_m_inv<-chol2inv(chol(R_m[,,l]))
        
        C_tAI[i,k,l]<-((1/sample_size_n*sum(diag(R_n_inv%*%H_n[,,l]%*%V_n[,,l]))
                        -1/sample_size_m*sum(diag(R_m_inv%*%H_m[,,l]%*%Z_n%*%G%*%t(Z_m))))
                       +0.5*(1/sample_size_m*sum(diag(R_m_inv%*%V_m[,,l]))
                             -1/sample_size_n*sum(diag(R_n_inv%*%V_n[,,l])))
                       +0.5*(1/sample_size_m*sum(diag(R_m_inv%*%H_m[,,l]%*%V_n[,,l]%*%t(H_m[,,l])))
                             -1/sample_size_n*sum(diag(R_n_inv%*%H_n[,,l]%*%V_n[,,l]%*%t(H_n[,,l]))))
                       +0.5*((1/sample_size_m)*sum(diag(log(R_m[,,l])))-(1/sample_size_n)*sum(diag(log(R_n[,,l])))))
        
        C_cAI[i,k,l]<-(1/sample_size_n)*sum(diag(H_n[,,l]))
        
        C_AI[i,k,l]<-(1/sample_size_n)*ncol(X_n)
        
      }
      
      ##Extracting Estimated variance matrices
        {
        ID<-rep(1:Nsubjects[i], each=10)     
        Prefit <- lmer(y ~X_n-1+(1|ID)+(time_n+-1|ID),REML=T)
        G_est<-diag(rep(c(unclass(VarCorr(Prefit))$'ID'[1],unclass(VarCorr(Prefit))$'ID.1'[1]),times=Nsubjects[i]))
        sig_2<- attr(VarCorr(Prefit),"sc")^2

        ZGZ_n_est<-Z_n%*%G_est%*%t(Z_n)
        ZGZ_m_est<-Z_m%*%G_est%*%t(Z_m)
        
        R_n_est[,,l] <-sig_2*diag(sample_size_n)
        V_n_est[,,l]<-(R_n_est[,,l]+ZGZ_n_est)    
        
        V_m_est[,,l]<-(sig_2*diag(sample_size_m) +ZGZ_m_est)
        R_m_est[,,l] <-V_m_est[,,l] - (Z_m%*%G_est%*%t(Z_n)%*%chol2inv(chol(V_n_est[,,l]))%*%Z_n%*%G_est%*%t(Z_m))
      }  

      ##Creating H_est
      {   
        V_n_est_inv<-chol2inv(chol(V_n_est[,,l]))
        
        H_n_est[,,l]<-(X_n%*%chol2inv(chol(t(X_n)%*%V_n_est_inv%*%X_n))%*%t(X_n)%*%V_n_est_inv+
                         Z_n%*%G_est%*%t(Z_n)%*%V_n_est_inv%*%(diag(1,sample_size_n)
                                                               -X_n%*%chol2inv(chol(t(X_n)%*%V_n_est_inv%*%X_n))%*%t(X_n)%*%V_n_est_inv))
        
        H_m_est[,,l]<-(X_m%*%chol2inv(chol(t(X_n)%*%V_n_est_inv%*%X_n))%*%t(X_n)%*%V_n_est_inv+
                         Z_m%*%G_est%*%t(Z_n)%*%V_n_est_inv%*%(diag(1,sample_size_n)
                                                               -X_n%*%chol2inv(chol(t(X_n)%*%V_n_est_inv%*%X_n))%*%t(X_n)%*%V_n_est_inv))
      }  
      
      ## Correction_est (for loss(Opt_t) change C_tAI_est to w_t_est, C_cAI_est to w_est and delete C_AI_est)
        {
        R_n_est_inv<-chol2inv(chol(R_n_est[,,l]))
        R_m_est_inv<-chol2inv(chol(R_m_est[,,l]))
        
        C_tAI_est[i,k,l]<-((1/sample_size_n*sum(diag(R_n_est_inv%*%H_n_est[,,l]%*%V_n_est[,,l]))
                            -1/sample_size_m*sum(diag(R_m_est_inv%*%H_m_est[,,l]%*%Z_n%*%G_est%*%t(Z_m))))
                           +0.5*(1/sample_size_m*sum(diag(R_m_est_inv%*%V_m_est[,,l]))
                                 -1/sample_size_n*sum(diag(R_n_est_inv%*%V_n_est[,,l])))
                           +0.5*(1/sample_size_m*sum(diag(R_m_est_inv%*%H_m_est[,,l]%*%V_n_est[,,l]%*%t(H_m_est[,,l])))
                                 -1/sample_size_n*sum(diag(R_n_est_inv%*%H_n_est[,,l]%*%V_n_est[,,l]%*%t(H_n_est[,,l]))))
                           +0.5*((1/sample_size_m)*sum(diag(log(R_m_est[,,l])))-(1/sample_size_n)*sum(diag(log(R_n_est[,,l])))))
        
        C_cAI_est[i,k,l]<-(1/sample_size_n)*sum(diag(H_n_est[,,l]))
        
        C_AI_est[i,k,l]<-(1/sample_size_n)*ncol(X_n)
      }
      
      ##loglike_f and loglike_f_est. (for loss(Opt_t) change the likelihood loss to squared errors loss)
      {
        y_p_f_n<-H_n[,,l]%*%y 
        y_p_f_n_est<-H_n_est[,,l]%*%y 
        
        loglike_f[i,k,l]<--1/sample_size_n*(0.5*t(y-y_p_f_n)%*%chol2inv(chol(R_n[,,l]))%*%(y-y_p_f_n))-0.5*(1/sample_size_n)*sum(diag(log(R_n[,,l])))
        loglike_f_est[i,k,l]<--1/sample_size_n*(0.5*t(y-y_p_f_n_est)%*%chol2inv(chol(R_n_est[,,l]))%*%(y-y_p_f_n_est))-0.5*(1/sample_size_n)*sum(diag(log(R_n_est[,,l])))
      }
      
    }
    
    # Allocating loglike_star_f1/2/3
    {
      loglike_star_f1<-matrix(data=NA,nrow=length_LLN,ncol=1)
      loglike_star_f2<-matrix(data=NA,nrow=length_LLN,ncol=1)
      loglike_star_f3<-matrix(data=NA,nrow=length_LLN,ncol=1)
    }
    
    {
      R_m_1_inv<-chol2inv(chol(R_m[,,1]))
      R_m_2_inv<-chol2inv(chol(R_m[,,2]))
      R_m_3_inv<-chol2inv(chol(R_m[,,3]))
    }    
    
    #Creating loglike_star_f1/2/3 (for loss(Opt_t) change the likelihood loss to squared errors loss)
    for (s in 1:length_LLN) {
      e_star<-rnorm(sample_size_m,mean=0,sd=sqrt(Ve[k])) 

      y_star<-X_m_true%*%beta+Z_m%*%b+e_star
      
      y_p_f1_m<-H_m[,,1]%*%y
      y_p_f2_m<-H_m[,,2]%*%y
      y_p_f3_m<-H_m[,,3]%*%y
      
      loglike_star_f1[s]<- -1/sample_size_m*(0.5*t(y_star-y_p_f1_m)%*%R_m_1_inv%*%(y_star-y_p_f1_m))-0.5*(1/sample_size_m)*sum(diag(log(R_m[,,1])))
      loglike_star_f2[s]<- -1/sample_size_m*(0.5*t(y_star-y_p_f2_m)%*%R_m_2_inv%*%(y_star-y_p_f2_m))-0.5*(1/sample_size_m)*sum(diag(log(R_m[,,2])))
      loglike_star_f3[s]<- -1/sample_size_m*(0.5*t(y_star-y_p_f3_m)%*%R_m_3_inv%*%(y_star-y_p_f3_m))-0.5*(1/sample_size_m)*sum(diag(log(R_m[,,3])))
    }
    
    #Averaging to Eloglike_star_f1/2/3    
    {
      Eloglike_star_f[i,k,1]<-1/length_LLN*sum(loglike_star_f1)
      Eloglike_star_f[i,k,2]<-1/length_LLN*sum(loglike_star_f2)
      Eloglike_star_f[i,k,3]<-1/length_LLN*sum(loglike_star_f3)
    }
  }
}  

