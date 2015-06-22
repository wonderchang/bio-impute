NRMSE<-function(M,type,percent)  
{ 
  answer<-NULL
  guess<-NULL
  k<-10      
  m<-2       
  order<-2   
  z_na<-as.matrix(na.omit(M))
  z_com<-as.matrix(z_na)
  
  z_na[sample(1:((nrow(z_na)*ncol(z_na))),(ncol(z_na)*nrow(z_na)*percent/100))]<-NA   #insert na into the completed matrix
  
  answer<-z_com[which(is.na(z_na))]                                 #answer stored the missing numbers
  
  
  if(type==1)
  {    z_zero<-as.matrix(Zero(z_na)) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==2)
  {    z_zero<-as.matrix(RowAverage(z_na)) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==3)
  {    z_zero<-as.matrix(KNN(z_na,15,'EuDist')) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==4)
  {    z_zero<-as.matrix(SKNN(z_na,k,'EuDist')) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==5)
  {    z_zero<-as.matrix(IKNN(z_na,k,'EuDist',m))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==6)
  {    z_zero<-as.matrix(ISKNN(z_na,k,'EuDist',m))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==7)
  {    z_zero<-as.matrix(LS(z_na,k,'EuDist'))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==8)
  {    z_zero<-as.matrix(LLS(z_na,k,'EuDist')) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==9)
  {    z_zero<-as.matrix(SLLS(z_na,k,'EuDist'))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==10)
  {    z_zero<-as.matrix(WLLS(z_na,k,'EuDist',order)) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==11)
  {    z_zero<-as.matrix(ShrLLS(z_na,k,'EuDist')) 
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==12)
  {    z_zero<-as.matrix(ISLLS(z_na,k,'EuDist',m))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==13)
  {    z_zero<-as.matrix(WSLLS(z_na,k,'EuDist',order))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==14)
  {    z_zero<-as.matrix(ShrSLLS(z_na,k,'EuDist'))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==15)
  {    z_zero<-as.matrix(IWLLS(z_na,k,'EuDist',m,order))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==16)
  {    z_zero<-as.matrix(ILLS(z_na,k,'EuDist',m))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==17)
  {    z_zero<-as.matrix(IShrLLS(z_na,k,'EuDist',m))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==18)
  {    z_zero<-as.matrix(ShrWLLS(z_na,k,'EuDist',order))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==19)
  {    z_zero<-as.matrix(IWSLLS(z_na,k,'EuDist',m,order))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==20)
  {    z_zero<-as.matrix(IShrSLLS(z_na,k,'EuDist',m))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==21)
  {    z_zero<-as.matrix(IShrWLLS(z_na,k,'EuDist',m,order))
       guess<-z_zero[which(is.na(z_na))]
  }    
  if(type==22)
  {    z_zero<-as.matrix(ShrWSLLS(z_na,k,'EuDist',order))
       guess<-z_zero[which(is.na(z_na))]
  }
  if(type==23)
  {    z_zero<-as.matrix(IShrWSLLS(z_na,k,'EuDist',m,order))
       guess<-z_zero[which(is.na(z_na))]
  }
  
  nrmse<-(mean((guess-answer)^2,na.rm=1)/var(answer))^0.5                       #calculate the nrmse
  return(nrmse)
}


############################################          #####################################
############################################   PLOT   #######################################
############################################          ######################################

library(ggplot2)
library(MASS)
source('../impute.R')

x<-read.table("../data/BohenSP.txt",header=T,row.names=1,sep="\t");
x.m <- as.matrix(x);

nrmse_mean<-NULL
se<-NULL
ntype<-23                       #from 1 to which type
times<-10                       #calculated for mean
percent<-6.5                    #percent of missing numbers 

type<-matrix(1:ntype,nc=times,nr=ntype,byrow=0)
nrmse<-matrix(1,nc=times,nr=ntype,byrow=0)

for(i in 1:(ncol(type)*nrow(type)))
{
  nrmse[i]<-NRMSE(x.m,type[i],percent)
  print(type[i])
}

nrmse_mean<-as.matrix(rowMeans(nrmse,dims=1))
se<-apply(nrmse,1,function(x) {sqrt(var(x)/length(x))})

df<-data.frame(method=c("Zero","RowAverage","KNN","SKNN","IKNN","ISKNN","LS","LLS","SLLS","WLLS","ShrLLS","ISLLS","WSLLS","ShrSLLS","IWLLS","ILLS","IShrLLS","ShrWLLS","IWSLLS","IShrSLLS","IShrWLLS","ShrWSLLS","IShrWSLLS"),
               nrmse_mean=nrmse_mean,
               se=se)
df$method <- reorder(df$method, -df$nrmse_mean)
ggplot(df, aes(x = method, y = nrmse_mean)) + geom_errorbar(aes(ymin=nrmse_mean-se, ymax=nrmse_mean+se), width=1) + geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

