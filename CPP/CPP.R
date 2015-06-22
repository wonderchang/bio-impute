CPP<-function(M,group,type,percent)
{
  n<-0
  k<-10       
  order<-2    
  m<-2
  z_na<-as.matrix(na.omit(M))
  z_com<-as.matrix(z_na)           #completed matrix: orignal matirx without na
  z_na[sample(1:((nrow(z_na)*ncol(z_na))),(ncol(z_na)*nrow(z_na)*percent/100))]<-NA   
  
  if(type==1)
  {y_com<-as.matrix(Zero(z_na))}
  if(type==2)
  {y_com<-as.matrix(RowAverage(z_na))}
  if(type==3)
  {y_com<-as.matrix(KNN(z_na,15,'EuDist'))}
  if(type==4)
  {y_com<-as.matrix(SKNN(z_na,k,'EuDist'))}
  if(type==5)
  {y_com<-as.matrix(IKNN(z_na,k,'EuDist',m))}
  if(type==6)
  {y_com<-as.matrix(ISKNN(z_na,k,'EuDist',m))}
  if(type==7)
  {y_com<-as.matrix(LS(z_na,k,'EuDist'))}
  if(type==8)
  {y_com<-as.matrix(LLS(z_na,k,'EuDist'))}
  if(type==9)
  {y_com<-as.matrix(SLLS(z_na,k,'EuDist'))}
  if(type==10)
  {y_com<-as.matrix(WLLS(z_na,k,'EuDist',order))}
  if(type==11)
  {y_com<-as.matrix(ShrLLS(z_na,k,'EuDist'))}
  if(type==12)
  {y_com<-as.matrix(ISLLS(z_na,k,'EuDist',m))}
  if(type==13)
  {y_com<-as.matrix(WSLLS(z_na,k,'EuDist',order))}
  if(type==14)
  {y_com<-as.matrix(ShrSLLS(z_na,k,'EuDist'))}
  if(type==15)
  {y_com<-as.matrix(IWLLS(z_na,k,'EuDist',m,order))}
  if(type==16)
  {y_com<-as.matrix(ILLS(z_na,k,'EuDist',m))}
  if(type==17)
  {y_com<-as.matrix(IShrLLS(z_na,k,'EuDist',m))}
  if(type==18)
  {y_com<-as.matrix(ShrWLLS(z_na,k,'EuDist',order))}
  if(type==19)
  {y_com<-as.matrix(IWSLLS(z_na,k,'EuDist',m,order))}
  if(type==20)
  {y_com<-as.matrix(IShrSLLS(z_na,k,'EuDist',m))}
  if(type==21)
  {y_com<-as.matrix(IShrWLLS(z_na,k,'EuDist',m,order))}    
  if(type==22)
  {y_com<-as.matrix(ShrWSLLS(z_na,k,'EuDist',order))}
  if(type==23)
  {y_com<-as.matrix(IShrWSLLS(z_na,k,'EuDist',m,order))}
  
  
  
  y_a<-kmeans(y_com,group)
  y_b<-as.matrix(y_a$cluster) 
  
  z_a<-kmeans(z_com,group) 
  z_b<-as.matrix(z_a$cluster) 
  
  c_m<-matrix(0,nr=group,nc=group)      
  
  for(i in 1:nrow(z_com))
  {c_m[z_b[i],y_b[i]]<-c_m[z_b[i],y_b[i]]+1}
  
  c_max<-apply(c_m,2,function(x) {max(x)}) 
  
  
  return(sum(c_max)/sum(c_m))  
  
}

library(ggplot2)
library(MASS)
source('../impute.R')

x<-read.table("../data/BohenSP.txt",header=T,row.names=1,sep="\t");
x.m <- as.matrix(x);
cpp_mean<-NULL
se<-NULL
ntype<-23
times<-10                       
percent<-6.5                    
group<-10       

type<-matrix(1:ntype,nc=times,nr=ntype,byrow=0)
cpp<-matrix(1,nc=times,nr=ntype,byrow=0)

for(i in 1:(ncol(type)*nrow(type)))
{
  cpp[i]<-CPP(x.m,group,type[i],percent)
  print(type[i])
}

cpp_mean<-as.matrix(rowMeans(cpp,dims=1))
se<-apply(cpp,1,function(x) {sqrt(var(x)/length(x))})

df<-data.frame(method=c("Zero","RowAverage","KNN","SKNN","IKNN","ISKNN","LS","LLS","SLLS","WLLS","ShrLLS","ISLLS","WSLLS","ShrSLLS","IWLLS","ILLS","IShrLLS","ShrWLLS","IWSLLS","IShrSLLS","IShrWLLS","ShrWSLLS","IShrWSLLS"),
               cpp_mean=cpp_mean,
               se=se)
df$method <- reorder(df$method, -df$cpp_mean)
ggplot(df, aes(x = method, y = cpp_mean)) + 
  geom_errorbar(aes(ymin=cpp_mean-se, ymax=cpp_mean+se), width=.1) +
  geom_point() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
