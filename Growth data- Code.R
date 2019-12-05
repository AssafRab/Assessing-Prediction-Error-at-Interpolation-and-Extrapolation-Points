library(lme4)
library(psych)
library(ggplot2)
library(dplyr)


# Output_folder<-"TablePath"
growth.data <- read.table("TablePath",header=T, sep="\t")

###Preparing the data
{
Subject <- rep(growth.data$Subject,4)
Nsubjects <- length(growth.data$Subject)
Sex <- rep(growth.data$Sex,4)
Time <- c(rep(8,Nsubjects),rep(10,Nsubjects),rep(12,Nsubjects),rep(14,Nsubjects))
Excluding_options<-c(8,10,12,14)

y <- c(growth.data$Dist8,growth.data$Dist10,growth.data$Dist12,growth.data$Dist14)
Data<-as.data.frame(cbind(Subject,Sex,Time,y))
}

###Defining the models
model_i<-list(c(1,3),c(1,2,3),c(1,2,3,5))

###Defining the covariates
{
X_all<-data.frame(rep(1,nrow(Data)),Data$Sex,Data$Time,Data$Time^2,
             Data$Sex*Data$Time,Data$Sex*Data$Time^2)
names(X_all)<-c("Intercept","Sex","Time","TimeSq","SexTime","SexTimeSq")
Nobs_subject<-length(unique(Data$Time))
Z_all<-matrix(data=0,nrow=Nobs_subject*Nsubjects,ncol=2*Nsubjects+1)

for(i in 1:Nsubjects)  {
  rowindex1<-(1+(i-1)*Nobs_subject)
  rowindex2<-i*Nobs_subject
  
  Z_all[rowindex1:rowindex2,2*i-1]<-rep(1,Nobs_subject)
  Z_all[rowindex1:rowindex2,2*i]<-c(unique(Data$Time))  
  Z_all[rowindex1:rowindex2,2*Nsubjects+1]<-c(unique(Data$Time))  
}
}
###Defining the output matrices
{
mAI<-matrix(data=NA,nrow=length(Excluding_options),ncol=length(model_i))
cAI<-matrix(data=NA,nrow=length(Excluding_options),ncol=length(model_i))
tAI<-matrix(data=NA,nrow=length(Excluding_options),ncol=length(model_i))
logLike<-matrix(data=NA,nrow=length(Excluding_options),ncol=length(model_i))
logLike_star<-matrix(data=NA,nrow=length(Excluding_options),ncol=length(model_i))
}

for(t in 1:length(Excluding_options)){

{
Exclude<-Excluding_options[t]  
Z_tr<-Z_all[Z_all[,ncol(Z_all)]!=Exclude,1:ncol(Z_all)-1]
Z_star<-Z_all[Z_all[,ncol(Z_all)]==Exclude,1:ncol(Z_all)-1]
}  
  
for (i in 1:length(model_i)){
 
  {
    X_tr_i<-as.matrix(X_all[X_all$Time!=Exclude,unlist(model_i[i])])
    X_star_i<-as.matrix(X_all[X_all$Time==Exclude,unlist(model_i[i])])
    y_tr_i<-Data[Data$Time!=Exclude,4]
    y_star_i<-Data[Data$Time==Exclude,4]
    Data_tr_i<-Data[Data$Time!=Exclude,]
  }

  fit <- lmer(y ~ -1+as.matrix(X_all)+(1+Time|Subject),data=Data,REML=T)

 ###Variance matrices
 {
 G<-diag(Nsubjects)%x%matrix(c(unclass(VarCorr(fit))$'Subject'[1:4]),2,2)
 R_tr<-attr(VarCorr(fit),"sc")^2*diag(length(y_tr_i)) 
 V_tr<-Z_tr%*%G%*%t(Z_tr)+R_tr
 Res_star<-attr(VarCorr(fit),"sc")^2*diag(length(y_star_i))   
 V_star<-Z_star%*%G%*%t(Z_star)+Res_star

 V_tr_inv<-chol2inv(chol(V_tr))

 R_star<-V_star-(Z_star%*%G%*%t(Z_tr)%*%V_tr_inv%*%Z_tr%*%G%*%t(Z_star))  
 R_tr_inv<-chol2inv(chol(R_tr))
 R_star_inv<-chol2inv(chol(R_star))
 }
  
 ###Hat matrices
 {
 H<-(X_tr_i%*%solve(t(X_tr_i)%*%V_tr_inv%*%X_tr_i)%*%t(X_tr_i)%*%V_tr_inv+
       Z_tr%*%G%*%t(Z_tr)%*%V_tr_inv%*%(diag(length(y_tr_i))
                                  -X_tr_i%*%solve(t(X_tr_i)%*%V_tr_inv%*%X_tr_i)%*%t(X_tr_i)%*%V_tr_inv))
 
 H_star<-(X_star_i%*%solve(t(X_tr_i)%*%V_tr_inv%*%X_tr_i)%*%t(X_tr_i)%*%V_tr_inv+
            Z_star%*%G%*%t(Z_tr)%*%V_tr_inv%*%(diag(length(y_tr_i))    
                                  -X_tr_i%*%solve(t(X_tr_i)%*%V_tr_inv%*%X_tr_i)%*%t(X_tr_i)%*%V_tr_inv))
 }
 ###Results
 {
cAI[t,i]<-(1/length(y_tr_i))*tr(H)
tAI[t,i]<-((1/length(y_tr_i))*tr(R_tr_inv%*%H%*%V_tr)
               -(1/length(y_star_i))*tr(R_star_inv%*%H_star%*%Z_tr%*%G%*%t(Z_star))
              +0.5*((1/length(y_star_i))*tr(R_star_inv%*%V_star)
                    -(1/length(y_tr_i))*tr(R_tr_inv%*%V_tr))
              +0.5*((1/length(y_star_i))*tr(R_star_inv%*%H_star%*%V_tr%*%t(H_star))
                    -(1/length(y_tr_i))*tr(R_tr_inv%*%H%*%V_tr%*%t(H)))
              +0.5*((1/length(y_star_i))*log(det(R_star))-(1/length(y_tr_i))*log(det(R_tr))))
mAI[t,i]<-(1/length(y_tr_i))*ncol(X_tr_i)
logLike[t,i]<-(1/length(y_tr_i))*(-0.5*log(det(R_tr))-0.5*t(y_tr_i-H%*%y_tr_i)%*%R_tr_inv%*%(y_tr_i-H%*%y_tr_i))
logLike_star[t,i]<-(1/length(y_star_i))*(-0.5*log(det(R_star))-0.5*t(y_star_i-H_star%*%y_tr_i)%*%R_star_inv%*%(y_star_i-H_star%*%y_tr_i))
 }
 
}

}


####Plots
Data_plot<-data.frame(rbind(cbind(-logLike+tAI,Measure="tAI"),cbind(-logLike+cAI,Measure="cAI"),cbind(-logLike+mAI,Measure="mAI"),cbind(-logLike_star,Measure="Oracle")),rep(c("8","10","12","14")))
names(Data_plot)<-c("1","2","3","Method","Excluded")
Data_plot<-melt(Data_plot,id=c("Method","Excluded"))
names(Data_plot)[3:4]<-c("Model","Error")
Data_plot$Error<-as.numeric(as.matrix(Data_plot[,4]))
Data_plot$Method<- factor(Data_plot$Method,levels=c("tAI", "cAI", "mAI", "Oracle"))
Data_plot$Excluded<- factor(Data_plot$Excluded,levels=c("8", "10", "12", "14"))


p1 <- ggplot() +
  geom_point(data=subset(Data_plot,Excluded==14), aes(x=Model,y=Error,color=Method),size=2,shape=1)+
  scale_colour_manual(values=c("#619CFF","#00BA38", "#F8766D", "#C77CFF"))+
  ggtitle(expression(atop("Growth Data", atop(italic("Holdout Age=14")))))+
  theme(plot.title = element_text(size = 11, face = "bold", colour = "black",hjust = 0.5, vjust = -0.7))+
  labs(y="Prediction Error")+
  theme(legend.text.align=0.5)+
  guides(fill=FALSE)

p1


p2 <- ggplot() +
  geom_point(data=subset(Data_plot,Excluded!=14), aes(x=Model,y=Error,color=Method),size=2,shape=1)+
  facet_wrap(~Excluded)+
  scale_colour_manual(values=c("#619CFF","#00BA38", "#F8766D", "#C77CFF"))+
   ggtitle(expression(atop("Growth Data", atop(italic("Holdout Ages {8, 10, 12}")))))+
  theme(plot.title = element_text(size = 11, face = "bold", colour = "black", vjust = -0.7,hjust = 0.5))+
  theme(legend.text.align=0.5)+
  labs(y="Prediction Error")+
  guides(fill=FALSE)+
  scale_shape_manual(name = "Method",values =c(0,2,1,3))
p2


