library(sp)
library(ggplot2)
library(dplyr)
library(psych)

data(meuse)

### Functions

##Loglikelihood function

{
  Log_like <- function(par,X,x1_dist,x2_dist,Y) {
    a<-par[1]
    b<-par[2]
    c<-par[3]
    d<-par[4]
    
    V<-(a^2)*exp(-1/(2*(b^2))*x1_dist-1/(2*(c^2))*x2_dist)+(d^2)*diag(nrow(x1_dist))
    mean<-X%*%solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%Y
    
    return(log(det(V))+t(Y-mean)%*%solve(V)%*%(Y-mean))
  }
}


###Distance between measurements
{
x1_nor<-meuse$x/max(meuse$x)
x2_nor<-meuse$y/max(meuse$y)

x1_distance<-as.matrix(dist(x1_nor, method = "euclidean",upper=T)^2)
x2_distance<-as.matrix(dist(x2_nor, method = "euclidean",upper=T)^2)
}

X_all<-cbind(model.matrix(~as.factor(meuse$ffreq)*meuse$dist+as.factor(meuse$soil)*meuse$dist))

###Splitting for training and test
{
set.seed(1)
sample <- sample.int(n = nrow(meuse), size = floor(.70*nrow(meuse)), replace = F)

x1_dist_tr<-x1_distance[sample,sample]
x2_dist_tr<-x2_distance[sample,sample]

x1_dist_star<-x1_distance[-sample,-sample]
x2_dist_star<-x2_distance[-sample,-sample]

x1_star_x1_dist<-x1_distance[-sample,sample]
x2_star_x2_dist<-x2_distance[-sample,sample]

X_all_tr<-X_all[sample, ]
X_all_star<-X_all[-sample, ]

y_tr<-log(meuse$lead[sample])
y_star<-log(meuse$lead[-sample])
}  

###Models
model_i<-list(c(1:6),c(1:6,7,8),c(1:6,9,10),c(1:10))
N_models<-length(model_i)

##Defining the output matrices
{
mAI<-matrix(data=NA,nrow=N_models,ncol=1)
cAI<-matrix(data=NA,nrow=N_models,ncol=1)
tAI<-matrix(data=NA,nrow=N_models,ncol=1)
logLike<-matrix(data=NA,nrow=N_models,ncol=1)
logLike_star<-matrix(data=NA,nrow=N_models,ncol=1)
}


for (i in 1:length(model_i)){
{
  X_i<-cbind(X_all_tr[,unlist(model_i[i])])
  X_i_star<-cbind(X_all_star[,unlist(model_i[i])])
}  
  
###Estimating covariance parameters  
{
  result <- optim(par=c(0.2,0.2,0.2,0.2), Log_like, X=X_i,x1_dist=x1_dist_tr,x2_dist=x2_dist_tr,Y=y_tr)
  sig_f<-abs(result$par[1])
  l1<-abs(result$par[2])
  l2<-abs(result$par[3])
  sig_e<-abs(result$par[4])
  } 

###Creating covariance matrices  
  {
  K<-sig_f^2*exp(-1/(2*(l1^2))*x1_dist_tr-1/(2*(l2^2))*x2_dist_tr)
  K_star<-sig_f^2*exp(-1/(2*(l1^2))*x1_dist_star-1/(2*(l2^2))*x2_dist_star)
  k_small_star<-sig_f^2*exp(-1/(2*(l1^2))*x1_star_x1_dist-1/(2*(l2^2))*x2_star_x2_dist)
  R<-sig_e^2*diag(length(y_tr))
  V<-K+R
  V_star<-K_star+sig_e^2*diag(length(y_star))
  R_star<-V_star-k_small_star%*%solve(V)%*%t(k_small_star)
  }  

###Creating covariance matrices  
  {
  H<-(X_i%*%solve(t(X_i)%*%solve(V)%*%X_i)%*%t(X_i)%*%solve(V)+
        K%*%solve(V)%*%(diag(length(y_tr))-X_i%*%solve(t(X_i)%*%solve(V)%*%X_i)%*%t(X_i)%*%solve(V)))
  
  H_star<-(X_i_star%*%solve(t(X_i)%*%solve(V)%*%X_i)%*%t(X_i)%*%solve(V)+
             k_small_star%*%solve(V)%*%(diag(length(y_tr))-X_i%*%solve(t(X_i)%*%solve(V)%*%X_i)%*%t(X_i)%*%solve(V)))
  }

###Results  
  {
  cAI[i]<-(1/length(y_tr))*tr(H)
  tAI[i]<-((1/length(y_tr))*tr(solve(R)%*%H%*%V)
            -(1/length(y_star))*tr(solve(R_star)%*%H_star%*%t(k_small_star))
            +0.5*((1/length(y_star))*tr(solve(R_star)%*%V_star)
                  -(1/length(y_tr))*tr(solve(R)%*%V))
            +0.5*((1/length(y_star))*tr(solve(R_star)%*%H_star%*%V%*%t(H_star))
                  -(1/length(y_tr))*tr(solve(R)%*%H%*%V%*%t(H)))
            +0.5*((1/length(y_star))*sum(log(diag(R_star)))-(1/length(y_tr))*sum(log(diag(R)))))
  mAI[i]<-(1/length(y_tr))*ncol(X_i)
  
  logLike[i]<-(1/length(y_tr))*(-0.5*sum(log(diag(R)))-0.5*t(y_tr-H%*%y_tr)%*%solve(R)%*%(y_tr-H%*%y_tr))
  logLike_star[i]<-(1/length(y_star))*(-0.5*sum(log(diag(R_star)))-0.5*t(y_star-H_star%*%y_tr)%*%solve(R_star)%*%(y_star-H_star%*%y_tr))
  }
}

###Plot
{
Plot_ass<-data.frame(rbind(-logLike+mAI,-logLike+cAI,-logLike+tAI,-logLike_star),c(rep("mAI",length(mAI)),rep("cAI",length(cAI)),rep("tAI",length(tAI)),rep("Oracle",length(logLike_star))),as.factor(rep(1:4)))
names(Plot_ass)<-c("Error","Method","Model")
Plot_ass<-Plot_ass%>%group_by(Method)%>%mutate(best=ifelse(Error==min(Error),"Selected","Not Selected"))
Plot_ass$Method<- factor(Plot_ass$Method,levels=c("tAI", "cAI", "mAI", "Oracle"))
}
 
p_meuse <- ggplot() +
  geom_point(data=Plot_ass, aes(x=as.factor(Model),y=Error,color=Method),size=2,shape=1)+
  labs(x = "Model", y="Prediction Error")+
  scale_colour_manual(values=c("#619CFF","#00BA38", "#F8766D", "#C77CFF"))+
  ggtitle(expression(atop("Meuse Data")))+
  theme(plot.title = element_text(size = 11, face = "bold", colour = "black", vjust = -1,hjust=0.5))+
  theme(legend.text.align=0.5)
p_meuse

