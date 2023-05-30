##模拟数据
#单视图模拟数据 (Fig7.2)
StabilityN=function(N){
n1=N/3;n2=N/3;n3=N/3;set.seed(1);X1<-rnorm(n1,15,1.5);set.seed(2);X2<-rnorm(n2,5,2);set.seed(3);X3<-rnorm(n3,30,2);Xv<-c(X1,X2,X3);data<-matrix(Xv,n1+n2+n3,2);data[1:n1,2]<-1;data[(n1+1):(n1+n2),2]<-2;data[(n1+n2+1):(n1+n2+n3),2]<-3;X<-matrix(data[,1],n1+n2+n3,1);label<-c(data[,2])
return(list(X=X,label=label))}
N=840;Stability=StabilityN(N)
##稳定性分析 (Fig 7.4)
#function of different N
StabilityN=function(N){
n1=N/3;n2=N/3;n3=N/3
set.seed(1);X1<-rnorm(n1,15,1.5);set.seed(2);X2<-rnorm(n2,5,2)
set.seed(3);X3<-rnorm(n3,30,2);Xv<-c(X1,X2,X3);data<-matrix(Xv,n1+n2+n3,2)
data[1:n1,2]<-1;data[(n1+1):(n1+n2),2]<-2;data[(n1+n2+1):(n1+n2+n3),2]<-3
X<-matrix(data[,1],n1+n2+n3,1);label<-c(data[,2])
 return(
list(X=X,label=label)
)
}
###generate
library(ORKM)
StaN<-c(90,210,360,480,600,720,840,960,1080)
StaNMI<-c( );StaP<-c( );StaF<-c( )
StaNMI1<-c( );StaP1<-c( );StaF1<-c( )
for (i in 1:9){
Stability=StabilityN(StaN[i])
set.seed(123)
ptm <- proc.time( )
q1=ORKMeans(X=Stability$X,K=3,V=1,chushi=StaN[i]/2,r=0.5,gamma=0.1,epsilon=1,yita=20,truere=Stability$label,max.iter=1000,method=0)
StaP[i] <- INDEX(q1$result, Stability$label, method = 0)
StaF[i]<- INDEX(q1$result, Stability$label, method = 3)
StaNMI[i]<-q1$NMI

## 敏感性分析 (Fig 7.5)
#fig4 single-view simulation data set
n1=n2=n3=70;set.seed(1);X1<-rnorm(n1,15,1.5);set.seed(2);X2<-rnorm(n2,5,2)
set.seed(3);X3<-rnorm(n3,30,2);Xv<-c(X1,X2,X3);data<-matrix(Xv,n1+n2+n3,2)
data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3
X<-matrix(data[,1],n1+n2+n3,1)
t<-c(50,100,110,120,130,140,150,160,170,180,190,200)
yitasm=c(5,6,7,8,9,10,11,12,13,14,15,20)
smF<-c( );smP<-c( );smNMI<-c( )
esmF<-c( );esmP<-c( );esmNMI<-c( )
e1smF<-c( );e1smP<-c( );e1smNMI<-c( )
for (i in 1:12){
set.seed(123)
q1=ORKMeans(X=X,K=3,V=1,chushi=t[i],r=0.5,gamma=0.1,epsilon=1,yita=20,truere=data[,2],max.iter=1000,method=0)
smP[i] <- INDEX(q1$result, data[,2], method = 0)
smF[i]<- INDEX(q1$result, data[,2], method = 3)
smNMI[i]<-q1$NMI
}
for (i in 1:12){
set.seed(123)
q1=ORKMeans(X=X,K=3,V=1,chushi=150,r=0.5,gamma=0.1,epsilon=1,yita=yitasm[i],truere=data[,2],max.iter=1000,method=0)
esmP[i] <- INDEX(q1$result, data[,2], method = 0)
esmF[i]<- INDEX(q1$result, data[,2], method = 3)
esmNMI[i]<-q1$NMI
}
for (i in 1:12){
set.seed(123)
q1=RKMeans(X=X,K=3,V=1,r=0.5,yita=yitasm[i],truere=data[,2],max.iter=1000,method=0)
e1smP[i] <- INDEX(q1$result, data[,2], method = 0)
e1smF[i]<- INDEX(q1$result, data[,2], method = 3)
e1smNMI[i]<-q1$NMI
}
#多视图模拟数据 (Fig 7.6)
library(ORKM);library(MASS)
n1=n2=n3=70 
set.seed(1) ;X1<-rnorm(n1,20,2); set.seed(2) ;X2<-rnorm(n2,25,1.5) 
set.seed(3) ;X3<-rnorm(n3,30,2) ;Xv<-c(X1,X2,X3) ;data<-matrix(Xv,n1+n2+n3,2) 
data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3 
X<-matrix(data[,1],n1+n2+n3,1) 
lamda1<-0.2;lamda2<-0.8;lamda<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
sol.svd <- svd(lamda);U<-sol.svd$u;D<-sol.svd$d;V<-sol.svd$v
C<-t(U)%*%t(X);Y<-C/D
view<-V%*%Y;view1<-matrix(view[1,]);view2<-matrix(view[2,])
N1<-nrow(view2);J1<-ncol(view2) 
X1<-as.matrix(view1);X2<-as.matrix(view2)
##稳定性分析 (Fig 7.8)
### multi-view simulation data(V=2,K=3,alpha1=0.2,alpha2=0.8)
#function of different N
StabilityN=function(N){
n1=N/3;n2=N/3;n3=N/3;set.seed(1);X1<-rnorm(n1,20,2) 
set.seed(2);X2<-rnorm(n2,25,1.5) ;set.seed(3);X3<-rnorm(n3,30,2) 
Xv<-c(X1,X2,X3) ;data<-matrix(Xv,n1+n2+n3,2) 
data[1:n1,2]<-1;data[(n1+1);(n1+n2),2]<-2.data[(n1+n2+1);(n1+n2+n3),2]<-3
X<-matrix(data[,1],n1+n2+n3,1);lamda1<-0.2;lamda2<-0.8
lamda<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
sol.svd <- svd(lamda);U<-sol.svd$u;D<-sol.svd$d;V<-sol.svd$v
C<-t(U)%*%t(X);Y<-C/D;view<-V%*%Y
view1<-matrix(view[1,]);view2<-matrix(view[2,])
X1<-matrix(view1,n1+n2+n3,1);X2<-matrix(view2,n1+n2+n3,1)
label<-c(data[,2])
 return(
list(X1=X1,X2=X2,label=label))
}
# view1
library(ORKM)
StaN<-c(90,210,360,480,600,720,840,960,1080)
v1StaNMI<-c( );v1StaP<-c( );v1StaF<-c( );v1StaNMI1<-c( );v1StaP1<-c( );v1StaF1<-c( )
for (i in 1:9){
v1Stability=StabilityN(StaN[i])
set.seed(123)
v1q1=ORKMeans(X=v1Stability$X1,K=3,V=2,chushi=StaN[i]/2,r=0.5,gamma=0.01,epsilon=1,yita=20,truere=v1Stability$label,max.iter=1000,method=0)
v1StaP[i] <- INDEX(v1q1$result, v1Stability$label, method = 0)
v1StaF[i]<- INDEX(v1q1$result, v1Stability$label, method = 3)
v1StaNMI[i]<-v1q1$NMI
}
for (i in 1:9){
v1Stability=StabilityN(StaN[i])
set.seed(123)
v1q2 = RKMeans(X = v1Stability$X1, K = 3, V = 2, yita = 20, r = 0.5, max.iter = 1000, truere = v1Stability$label, method = 0)
v1StaNMI1[i] <- v1q2$NMI
v1StaP1[i] <- INDEX(v1q2$result, v1Stability$label, method = 0)
v1StaF1[i] <- INDEX(v1q2$result, v1Stability$label, method = 3)
}
###view2
StaN<-c(90,210,360,480,600,720,840,960,1080)
v2StaNMI<-c( );v2StaP<-c( );v2StaF<-c( );v2StaNMI1<-c( );v2StaP1<-c( );v2StaF1<-c( )
for (i in 8:9){
v2Stability=StabilityN(StaN[i])
set.seed(123)
v2q1=ORKMeans(X=v2Stability$X2,K=3,V=2,chushi=StaN[i]/2,r=0.5,gamma=0.0001,epsilon=1,yita=20,truere=v2Stability$label,max.iter=10,method=0)
v2StaP[i] <- INDEX(v2q1$result, v2Stability$label, method = 0)
v2StaF[i]<- INDEX(v2q1$result, v2Stability$label, method = 3)
v2StaNMI[i]<-v2q1$NMI
}
for (i in 1:9){
v2Stability=StabilityN(StaN[i])
set.seed(123)
v2q2 = RKMeans(X = v2Stability$X2, K = 3, V = 2, yita = 20, r = 0.5, max.iter = 10, truere = v2Stability$label, method = 0)
v2StaNMI1[i] <- v2q2$NMI
v2StaP1[i] <- INDEX(v2q2$result, v2Stability$label, method = 0)
v2StaF1[i] <- INDEX(v2q2$result, v2Stability$label, method = 3)
}

## 敏感性分析 (Fig 7.9)
#fig8 multi-view simulation data set
n1=n2=n3=70 
set.seed(1) ;X1<-rnorm(n1,20,2); set.seed(2) ;X2<-rnorm(n2,25,1.5) 
set.seed(3) ;X3<-rnorm(n3,30,2) ;Xv<-c(X1,X2,X3) ;data<-matrix(Xv,n1+n2+n3,2) 
data[1:70,2]<-1;data[71:140,2]<-2;data[141:210,2]<-3 
X<-matrix(data[,1],n1+n2+n3,1) 
lamda1<-0.2;lamda2<-0.8;lamda<-matrix(c(lamda1,lamda2),nrow=1,ncol=2)
sol.svd <- svd(lamda);U<-sol.svd$u;D<-sol.svd$d;V<-sol.svd$v
C<-t(U)%*%t(X);Y<-C/D
view<-V%*%Y;view1<-matrix(view[1,]);view2<-matrix(view[2,])
N1<-nrow(view2);J1<-ncol(view2) 
X1<-as.matrix(view1);X2<-as.matrix(view2)
t<-c(50,100,110,120,130,140,150,160,170)
yitamvsm=c(1,5,10,30,50,70,90,100,200)
smF<-c( );smP<-c( );smNMI<-c( )
esmF<-c( );esmP<-c( );esmNMI<-c( )
e1smF<-c( );e1smP<-c( );e1smNMI<-c( )
e2smF<-c( );e2smP<-c( );e2smNMI<-c( )
e22smF<-c( );e22smP<-c( );e22smNMI<-c( )
for (i in 1:9){
set.seed(123)
q1=ORKMeans(X=X2,K=3,V=2,chushi=t[i],r=0.5,gamma=0.0001,epsilon=1,yita=20,truere=data[,2],max.iter=1000,method=0)
smP[i] <- INDEX(q1$result, data[,2], method = 0)
smF[i]<- INDEX(q1$result, data[,2], method = 3)
smNMI[i]<-q1$NMI
}
for (i in 1:9){
set.seed(123)
q1=ORKMeans(X=X1,K=3,V=2,chushi=170,r=0.5,gamma=0.0001,epsilon=1,yita=yitamvsm[i],truere=data[,2],max.iter=1000,method=0)
esmP[i] <- INDEX(q1$result, data[,2], method = 0)
esmF[i]<- INDEX(q1$result, data[,2], method = 3)
esmNMI[i]<-q1$NMI
}
for (i in 1:9){
set.seed(123)
q1=RKMeans(X=X1,K=3,V=2,r=0.5,yita=yitamvsm[i],truere=data[,2],max.iter=1000,method=0)
e1smP[i] <- INDEX(q1$result, data[,2], method = 0)
e1smF[i]<- INDEX(q1$result, data[,2], method = 3)
e1smNMI[i]<-q1$NMI
}
for (i in 1:9){
set.seed(123)
q1=ORKMeans(X=X2,K=3,V=2,chushi=170,r=0.5,gamma=0.0001,epsilon=1,yita=yitamvsm[i],truere=data[,2],max.iter=1000,method=0)
e2smP[i] <- INDEX(q1$result, data[,2], method = 0)
e2smF[i]<- INDEX(q1$result, data[,2], method = 3)
e2smNMI[i]<-q1$NMI
}
for (i in 1:9){
set.seed(123)
q1=RKMeans(X=X2,K=3,V=2,r=0.5,yita=yitamvsm[i],truere=data[,2],max.iter=1000,method=0)
e22smP[i] <- INDEX(q1$result, data[,2], method = 0)
e22smF[i]<- INDEX(q1$result, data[,2], method = 3)
e22smNMI[i]<-q1$NMI
}
#敏感性分析Fig 7.10 (只给出V=5的代码) 
StabilityN=function(N){
n1=N/3;n2=N/3;n3=N/3
set.seed(1) ;X1<-rnorm(n1,20,2) ;set.seed(2) ;X2<-rnorm(n2,25,1.5) ;set.seed(3) ;X3<-rnorm(n3,30,2) 
Xv<-c(X1,X2,X3) ;data<-matrix(Xv,n1+n2+n3,2) ;data[1:n1,2]<-1
data[(n1+1):(n1+n2),2]<-2;data[(n1+n2+1):(n1+n2+n3),2]<-3
X<-matrix(data[,1],n1+n2+n3,1) 
lamda1<-0.1;lamda2<-0.1;lamda3<-0.1;lamda4<-0.1;lamda5<-0.1
lamda<-matrix(c(lamda1, lamda2, lamda3, lamda4, lamda5),nrow=1,ncol=1)
sol.svd <- svd(lamda)
U<-sol.svd$u;D<-sol.svd$d;V<-sol.svd$v;C<-t(U)%*%t(X);Y<-C/D
view<-V%*%Y;view1<-matrix(view[1,]);view2<-matrix(view[1,])
view3<-matrix(view[1,]);view4<-matrix(view[1,]);view5<-matrix(view[1,])
X1<-matrix(view1,n1+n2+n3,1);X2<-matrix(view1,n1+n2+n3,1)
X3<-matrix(view1,n1+n2+n3,1);X4<-matrix(view1,n1+n2+n3,1)
X5<-matrix(view1,n1+n2+n3,1);label<-c(data[,2])
 return(
list(main=X2, label=label) ) #the max lamda is the main view
}
V<-c(1,2,3,4,5)
#main view result
z1SenNMI<-c( );z1SenP<-c( );z1SenF<-c( );z2SenNMI<-c( );z2SenP<-c( );z2SenF<-c( )
N=90;Stability=StabilityN(N)
set.seed(123)
z1=ORKMeans(X=Stability$main,K=3,V=3,chushi=N/2,r=0.5,gamma=0.00001,epsilon=1,yita=20,truere=Stability$label,max.iter=10,method=0)
z1SenP[1] <- INDEX(z1$result, Stability$label, method = 0)
z1SenF[1]<- INDEX(z1$result, Stability$label, method = 3)
z1SenNMI[1]<-z1$NMI
set.seed(123)
z2 = RKMeans(X = Stability$main, K = 3, V = 2, yita = 20, r = 0.5, max.iter = 10, truere = Stability$label, method = 0)
z2SenNMI[5] <- z2$NMI
z2SenP[5] <- INDEX(z2$result, Stability$label, method = 0)
z2SenF[5] <- INDEX(z2$result, Stability$label, method = 3)
##真实数据分析
# Seed data set (Fig 7.12(a),Fig 7.13(a))
library(ORKM)
data(seed);data<-seed[,-8];X<-as.matrix(data)
t<-c(100,110,120,130,140,150,160,170,180,190,200);yitaseed=c(0.5,1,5,8,10,15,20)
tseedP<-c( );tseedF<-c( );tseedNMI<-c( )
eseedP<-c( );eseedF<-c( );eseedNMI<-c( )
e1seedP<-c( );e1seedF<-c( );e1seedNMI<-c( )
for (i in 1:11){
set.seed(123)
q1=ORKMeans(X=X,K=3,V=1,chushi=t[i],r=0.5,gamma=0.1,epsilon=1,yita=0.5,truere=seed[,8],max.iter=1000,method=0)
tseedP[i] <- INDEX(q1$result, seed[,8], method = 0)
tseedF[i]<- INDEX(q1$result, seed[,8], method = 3)
tseedNMI[i]<-q1$NMI
}
for (i in 1:7){
set.seed(123)
q1=ORKMeans(X=X,K=3,V=1,chushi=190,r=0.5,gamma=0.1,epsilon=1,yita=yitaseed[i],truere=seed[,8],max.iter=1000,method=0)
eseedP[i] <- INDEX(q1$result, seed[,8], method = 0)
eseedF[i]<- INDEX(q1$result, seed[,8], method = 3)
eseedNMI[i]<-q1$NMI
}
for (i in 1:7){
set.seed(123)
q1=RKMeans(X=X,K=3,V=1,r=0.5,yita=yitaseed[i],truere=seed[,8],max.iter=1000,method=0)
e1seedP[i] <- INDEX(q1$result, seed[,8], method = 0)
e1seedF[i]<- INDEX(q1$result, seed[,8], method = 3)
e1seedNMI[i]<-q1$NMI
}
# Sobar data set  (Fig 7.12(b),Fig 7.13(b))
data(sobar); colnames(sobar) <- NULL;sobar[1:22,20]<-1;sobar[23:72,20]<-2
data<-sobar[,-20];X<-as.matrix(data); t1<-c(45,50,55,60,65);
yitasobar=c(50,80,100,150,200,250,300)
sobarP<-c( );sobarF<-c( );sobarNMI<-c( )
esobarP<-c( );esobarF<-c( );esobarNMI<-c( )
sobarP1<-c( );sobarF1<-c( );sobarNMI1<-c( )
for (i in 1:5){
set.seed(123)
q2=ORKMeans(X=X,K=2,V=1,chushi=t1[i],r=0.5,gamma=0.1,epsilon=1,yita=20,truere=sobar[,20],max.iter=1000,method=0)
sobarP[i] <- INDEX(q2$result, sobar[-73,20], method = 0)
sobarF[i]<- INDEX(q2$result, sobar[-73,20], method = 3)
sobarNMI[i]<-q2$NMI
}
for (i in 1:7){
set.seed(123)
q22=ORKMeans(X=X,K=2,V=1,chushi=45,r=0.5,gamma=0.1,epsilon=1,yita=yitasobar[i],truere=sobar[-73,20],max.iter=1000,method=0)
esobarP[i] <- INDEX(q22$result, sobar[-73,20], method = 0)
esobarF[i]<- INDEX(q22$result, sobar[-73,20], method = 3)
esobarNMI[i]<-q22$NMI
set.seed(123)
q21=RKMeans(X=X,K=2,V=1,r=0.5,yita=yitasobar[i],truere=sobar[-73,20],max.iter=1000,method=0)
sobarP1[i] <- INDEX(q21$result, sobar[-73,20], method = 0)
sobarF1[i]<- INDEX(q21$result, sobar[-73,20], method = 3)
sobarNMI1[i]<-q21$NMI
}
#Cora数据集
# 敏感性分析 (Fig 7.15 Fig 7.16)
t<-c(2000,2100,2200,2300,2400,2500,2600, 2705)
coraF<-c( );coraP<-c( );coraNMI<-c( )
for (i in 1:8){
set.seed(123)
q1=ORKMeans(X=X1,K=7,V=4,chushi=t[i],r=0.5,gamma=0.1,epsilon=1,yita=0.5,truere=label0,max.iter=10,method=0)
coraP[i] <- INDEX(q1$result, label0, method = 0)
coraF[i]<- INDEX(q1$result, label0, method = 3)
coraNMI[i]<-q1$NMI
}
eta<-c(0.5,0.6,0.7,0.8,0.9,1,1.5,2,2.5,3,3.5,4,5)
eNMIcora<-c( );ePcora<-c( );eFcora<-c( );neta<-12
for (o in 1:neta){
q2=ORKMeans(X=X1,K=7,V=4,chushi=2705,r=0.5,gamma=0.1,epsilon=1,yita=eta[o],truere=label0,max.iter=10,method=0)
eNMIcora [o]<-q2$NMI
ePcora [o]<-com_accuracy(q2$result,label0,method=0)
eFcora [o]<-com_accuracy(q2$result,label0,method=3)
}
