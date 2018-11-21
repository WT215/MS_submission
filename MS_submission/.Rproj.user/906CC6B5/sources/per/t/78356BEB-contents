#########LOAD FUNCTIONS############
##########BEGIN OF PREPARE
library("MALDIquant")
library('MALDIquantForeign')
library("readMzXmlData")
library(dendextend)
library(gplots)
readCSV<-function(address){
  setwd(address)
  temp = list.files(path=address,pattern="*.csv")
  a<-list()
  for (i in 1:length(temp)) 
  {
    a[[i]]<-read.csv(temp[i])
    a[[i]]<-na.omit(a[[i]])
  }
  return(a)}
### convert mz int into file read by MALDIquant from list
createSpect<-function(list,cname){
  spectra<-list()
  for(i in 1:length(list)){
    spectra[[i]]<-createMassSpectrum(mass=list[[i]][,1],intensity=list[[i]][,2],metaData=list(name=cname[i]))
  }
  return(spectra)}

changeMZXMLtomzint<-function(templist,address){
  setwd(address)
  setwd(address)
  data<-list()
  for(i in 1:length(templist))
  {
    spec<-readMzXmlFile(templist[i])
    data[[i]]<-cbind(spec$spectrum$mass, spec$spectrum$intensity)
  }
  return(data)}

featureMat<-function(spectra_input,ref){
  
  massList<-lapply(spectra_input,function(x){x@mass})
  intList<-lapply(spectra_input,function(x){x@intensity})
  temp_mat<-extractFeatureMatrix(massList,intList, ref ,length(ref), length(spectra_input))
  rownames(temp_mat)<-ref
  
  return(temp_mat)
}



Rcpp::sourceCpp('E:/MainProject/PROJECT_MRSA_MSSA/C_codes/extractFeatureMatrix.cpp')
########################Begin####################
load("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/LOAD_7_MRSAMSSA_data.RData")

data_name<-c('pc','pa','r1','r2','pa15062016','R1','R2')
All_data_list<-list(pc_list,pa_list,r1_list,r2_list,pa15062016_list,MRMS_050816run1_list,MRMS_050816run2_list)
names(All_data_list)<-data_name
All_data_cname<-list(pc_name,pa_name,r1_name,r2_name,pa15062016_name,MRMS_050816run1_name,MRMS_050816run2_name)




TRUE_LABELS<-c(rep(1,10),rep(3,10))





data_list<-do.call(cbind,(All_data_list))
All_cname<-do.call(c,All_data_cname)
new_spectra<-createSpect(data_list,All_cname)
spectra<-trim(new_spectra)
#spectra<- trim(new_spectra,range=c(1500,1900))
#spectra<- trim(new_spectra,range=c(1600,4000))
spectra3<-smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)


####compare 2 and 1 halfwindowsize
#plot(spectra[[9]],xlim=c(2139,2145))
#lines(spectra3[[9]],xlim=c(2139,2145),col=2)

########test baseline
#baseline <- estimateBaseline(spectra3[[15]], method='TopHat',halfWindowSize=10)
#plot(spectra3[[15]])
#plot(trim(spectra3,range=c(2139,2145))[[15]])
#lines(baseline, col="red", lwd=2)
########end of test baseline
spectra4 <- removeBaseline(spectra3, method='TopHat',halfWindowSize=10)

spectra5<-calibrateIntensity(spectra4, method="TIC")
#spectra6<-alignSpectra(spectra5,SNR=2,noiseMethod='MAD',halfWindowSize = 10,tolerance=0.02,reference=maybe)
#########alignSpectra
SNR_P<-10
peaks<-detectPeaks(spectra5,method="MAD",halfWindowSize=20, SNR=SNR_P)
#peaks<-detectPeaks(spectra5,method="MAD",halfWindowSize=100, SNR=SNR_P)
reff<-referencePeaks(peaks, method=c("strict"), minFrequency=0.5,tolerance=0.02)
fdf<- determineWarpingFunctions(peaks,reference=reff,tolerance=0.02,method="lowess",plot=F,plotInteractive=F)
spectra6<-warpMassSpectra(spectra5, fdf)
length(Select_ind)

peaks <- detectPeaks(spectra6, SNR=2, halfWindowSize=20)
#peaks <- detectPeaks(spectra6, SNR=2, halfWindowSize=100)
peaks <- binPeaks(peaks,method='strict',tolerance=0.02)
peaks<-filterPeaks(peaks,minFrequency=0.75,mergeWhitelists = F)
featureMatrix <- intensityMatrix(peaks, spectra6)
rownames(featureMatrix)<-All_cname
dim(featureMatrix)

TRUE_LABELS_all<-rep(TRUE_LABELS,7)
names(TRUE_LABELS_all)<-All_cname


dim(featureMatrix)

HC_plot(featureMatrix,TRUE_LABELS_all,Method='COR',MAIN='(a) Full features: Routine workflow')


library(glmnet)
train<-sample(seq(1,120),round(120*0.7))
test<-seq(1,120)[-train]


lasso.mod<-glmnet(featureMatrix[train,],TRUE_LABELS_all[train],alpha=1,family="binomial")
set.seed(12300)
cv.out<-cv.glmnet(featureMatrix[train,],TRUE_LABELS_all[train],alpha=1,family="binomial")
plot(cv.out)
bestlam<-cv.out$lambda.min
lasso.pred<-predict(lasso.mod,s=bestlam,newx=featureMatrix[test,],type='class')
table(lasso.pred,TRUE_LABELS_all[test])
out<-glmnet(featureMatrix,TRUE_LABELS_all,alpha=1,family="binomial")
lasso.coef<-predict(out,type='coefficients',s=bestlam)
qq<-as.matrix(lasso.coef)[-1,]
qq<-qq[qq!=0]
qq2<-as.numeric(names(qq))
length(which(qq2<1200))/length(qq2)


dim(featureMatrix)
cv.out<-cv.glmnet(featureMatrix,TRUE_LABELS_all,alpha=1,family="binomial")
bestlam<-cv.out$lambda.min

out<-glmnet(featureMatrix,TRUE_LABELS_all,alpha=1,family="binomial")
lasso.coef<-predict(out,type='coefficients',s=bestlam)
qq<-as.matrix(lasso.coef)[-1,]
qq<-qq[qq!=0]
qq2<-as.numeric(names(qq))
length(which(qq2<1200))/length(qq2)
############gather all data SVM###########
library(e1071)
gather_all<-t(featureMatrix)

gather_all_label<-TRUE_LABELS_all
gather_all2<-rbind(gather_all_label,gather_all)
store<-NULL
set.seed(12300)
for(i in 1:100)
{
  print(i)
  train<-sample(dim(gather_all)[2],dim(gather_all)[2]*0.7)
  test<--train
  
  
  model1<-svm(t(gather_all2[-1,train]), as.factor(gather_all2[1,train]), cost = 10,  scale=F, type="C-classification", kernel="linear" ) 
  p<-predict(model1,t(gather_all2[-1,test]))
  temp_table<-table(p,as.factor(gather_all2[1,test]))
  #print(temp_table)
  acc<-sum(diag(temp_table))/sum(temp_table)
  store<-c(store,acc)
  
}
plot(seq(1:100),store,type='l',xlab='Repeated times',ylab='Accuracy rate',main=paste('Mean accuracy =',round(mean(store),4),', sd =',round(sd(store),4)))


#####binda######
library(binda)
Xall<-featureMatrix
Yall<-TRUE_LABELS_all

thr = optimizeThreshold(Xall, Yall)
Xall.b = dichotomize(Xall, thr)
br<-binda.ranking(Xall.b, Yall)
plot(br, top=100, arrow.col="black", zeroaxis.col="black", ylab="Peaks (m/z)",main="30 Most Differentially Expressed Peaks")






par(mfrow=c(3,1))
data<-Xall.b
#data<-featureMatrix[,Selectind]
DATA=featureMatrix



source("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/FUNCTIONS.R")
HC_plot(featureMatrix,TRUE_LABELS_all,Method='COR',MAIN='(a) Full features: Routine workflow',leg=F)

HC_plot(data,TRUE_LABELS_all,Method='BIN',leg=F,MAIN='(b) Full features: Optimal binary feature matrix')

dim(featureMatrix)
TOP100_mz<-as.numeric(rownames(br[1:100,]))
Selectind<-as.data.frame(br[1:100,])$idx
data<-Xall.b[,Selectind]
HC_plot(data,TRUE_LABELS_all,Method='BIN',leg=F,MAIN='(c) Top 100 features: Optimal binary feature matrix')

plot(br, top=100, arrow.col="black", zeroaxis.col="black", ylab="Peaks (m/z)",main="100 Most Differentially Expressed Peaks")



pr.out=prcomp(data, scale=FALSE)
pr_s<-summary(pr.out)
pr_v<-pr_s$importance[2,]
ff<-as.data.frame(pr.out$x[,c(1,2)])
cool<-TRUE_LABELS_all
ff<-cbind(ff,cool)
plot(ff[,1:2],   pch=19,col=cool,xlab=paste('PC1(',round(pr_v[1],4),'%)',sep=''),ylab=paste('PC2(',round(pr_v[2],4),'%)',sep=''))
legend('bottomleft',legend=c("MRSA","MSSA"),col=c(1,3),pch=16,cex=1)



library(binda)
gather_all<-t(data)
gather_all_label<-TRUE_LABELS_all
gather_all2<-rbind(gather_all_label,gather_all)
store_binda<-NULL
set.seed(12300)
for(i in 1:100)
{
  
  train<-sample(dim(gather_all)[2],dim(gather_all)[2]*0.7)
  test<--train
  
  
  model1<-binda(t(gather_all2[-1,train]), as.factor(gather_all2[1,train]), verbose=FALSE ) 
  p<-predict.binda(model1,t(gather_all2[-1,test]),verbose=F)
  

  
  temp_table<-table(p$class,as.factor(gather_all2[1,test]))
  #print(temp_table)
  acc<-sum(diag(temp_table))/sum(temp_table)
  store_binda<-c(store_binda,acc)
  
}
plot(seq(1:100),store_binda,type='l',xlab='Repeated times',ylab='Accuracy rate',main=paste('Mean accuracy =',round(mean(store_binda),4),', sd =',round(sd(store_binda),4)))









dim(featureMatrix)
#####binda subset######
library(binda)
names(TRUE_LABELS_all)<-All_cname
com_mat<-combn(7, 5)
set.seed(12300)
Reference_ind<-com_mat[,sample(seq(1,dim(com_mat)[2]))]



Train_Samples<-do.call(c,All_data_cname[Reference_ind])
Test_Samples<-setdiff(All_cname,Train_Samples)

Xtrain<-featureMatrix[Train_Samples,]
Ytrain<-TRUE_LABELS_all[Train_Samples]
Xtest<-featureMatrix[Test_Samples,]
Ytest<-TRUE_LABELS_all[Test_Samples]


thr = optimizeThreshold(Xtrain, Ytrain)
Xtrain.b = dichotomize(Xtrain, thr)
br<-binda.ranking(Xtrain.b, Ytrain)
plot(br, top=30, arrow.col="black", zeroaxis.col="black", ylab="Peaks (m/z)",main="30 Most Differentially Expressed Peaks")

Xtest.b<-dichotomize(Xtest, thr)

Selectind<-as.data.frame(br[1:100,])$idx
rownames(br)[1:30]


Inputdat<-Xtrain.b
dmat<-dist(Inputdat, method="binary")
h_temp<-hclust( dmat , method="ward.D2")
memb <- cutree(h_temp, k = 2)
names(memb)
#Ytrain[names(memb[memb==1])]
#Ytrain[names(memb[memb==2])]

leftMSSA<-length(which(Ytrain[names(memb[memb==2])]==3))
leftMRSA<-length(which(Ytrain[names(memb[memb==2])]==1))
rightMSSA<-length(which(Ytrain[names(memb[memb==1])]==3))
rightMRSA<-length(which(Ytrain[names(memb[memb==1])]==1))


hc = as.dendrogram(hclust( dmat , method="ward.D2"))
labels_colors(hc) <- TRUE_LABELS_all[order.dendrogram(hc)]
plot(hc)
legend('top',legend=c(paste("Left:",leftMSSA,"MSSA,",leftMRSA,"MRSA"),paste("Right:",rightMSSA,"MSSA,",rightMRSA,"MRSA")),cex=0.8)


pr.out=prcomp(Inputdat, scale=FALSE)
pr_s<-summary(pr.out)
pr_v<-pr_s$importance[2,]
ff<-as.data.frame(pr.out$x[,c(1,2)])
cool<-Ytrain
ff<-cbind(ff,cool)
plot(ff[,1:2],   pch=19,col=cool,xlab=paste('PC1(',round(pr_v[1],4),'%)',sep=''),ylab=paste('PC2(',round(pr_v[2],4),'%)',sep=''))
legend('bottomleft',legend=c("MRSA","MSSA"),col=c(1,3),pch=16,cex=1)



#test

Inputdat<-Xtest.b
dmat<-dist(Inputdat, method="binary")
h_temp<-hclust( dmat , method="ward.D2")
memb <- cutree(h_temp, k = 2)
names(memb)
#Ytrain[names(memb[memb==1])]
#Ytrain[names(memb[memb==2])]

leftMSSA<-length(which(Ytest[names(memb[memb==2])]==3))
leftMRSA<-length(which(Ytest[names(memb[memb==2])]==1))
rightMSSA<-length(which(Ytest[names(memb[memb==1])]==3))
rightMRSA<-length(which(Ytest[names(memb[memb==1])]==1))


hc = as.dendrogram(hclust( dmat , method="ward.D2"))
labels_colors(hc) <- TRUE_LABELS_all[order.dendrogram(hc)]
plot(hc)
legend('top',legend=c(paste("Left:",leftMSSA,"MSSA,",leftMRSA,"MRSA"),paste("Right:",rightMSSA,"MSSA,",rightMRSA,"MRSA")),cex=0.8)


pr.out=prcomp(Inputdat, scale=FALSE)
pr_s<-summary(pr.out)
pr_v<-pr_s$importance[2,]
ff<-as.data.frame(pr.out$x[,c(1,2)])
cool<-Ytest
ff<-cbind(ff,cool)
plot(ff[,1:2],   pch=19,col=cool,xlab=paste('PC1(',round(pr_v[1],4),'%)',sep=''),ylab=paste('PC2(',round(pr_v[2],4),'%)',sep=''))
legend('bottomleft',legend=c("MRSA","MSSA"),col=c(1,3),pch=16,cex=1)




Inputdat<-Xtest.b
library(binda)
gather_all<-t(Inputdat)
gather_all_label<-Ytest
gather_all2<-rbind(gather_all_label,gather_all)
store_binda<-NULL
set.seed(12300)
for(i in 1:100)
{
  
  train<-sample(dim(gather_all)[2],dim(gather_all)[2]*0.7)
  test<--train
  
  
  model1<-binda(t(gather_all2[-1,train]), as.factor(gather_all2[1,train]), verbose=FALSE ) 
  p<-predict.binda(model1,t(gather_all2[-1,test]),verbose=F)
  
  
  
  temp_table<-table(p$class,as.factor(gather_all2[1,test]))
  #print(temp_table)
  acc<-sum(diag(temp_table))/sum(temp_table)
  store_binda<-c(store_binda,acc)
  
}
plot(seq(1:100),store_binda,type='l',xlab='Repeated times',ylab='Accuracy rate',main=paste('Mean accuracy =',round(mean(store_binda),4),', sd =',round(sd(store_binda),4)))


#####binda subset continue######
library(foreach)
library(binda)
names(TRUE_LABELS_all)<-All_cname

par(mfrow=c(2,3))
Binda_explist<-foreach(j=1:6)%do%{
com_mat<-combn(7, j)
store<-NULL
for(i in 1:dim(com_mat)[2])

{
  print(i)
Reference_ind<-com_mat[,i]
Train_Samples<-do.call(c,All_data_cname[Reference_ind])
Test_Samples<-setdiff(All_cname,Train_Samples)

Xtrain<-featureMatrix[Train_Samples,]
Ytrain<-TRUE_LABELS_all[Train_Samples]
Xtest<-featureMatrix[Test_Samples,]
Ytest<-TRUE_LABELS_all[Test_Samples]


thr = optimizeThreshold(Xtrain, Ytrain)
Xtrain.b = dichotomize(Xtrain, thr)
br<-binda.ranking(Xtrain.b, Ytrain)
#plot(br, top=30, arrow.col="black", zeroaxis.col="black", ylab="Peaks (m/z)",main="30 Most Differentially Expressed Peaks")

Xtest.b<-dichotomize(Xtest, thr)

model1<-binda(Xtrain.b, as.factor(Ytrain), verbose=FALSE ) 
p<-predict.binda(model1,Xtest.b,verbose=F)
temp_table<-table(p$class,as.factor(Ytest))
#print(temp_table)
acc<-sum(diag(temp_table))/sum(temp_table)
store<-c(store,acc)
}


plot(seq(1,dim(com_mat)[2]),store,type='l',xlab='Different subsets of 7 data as training datasets',ylab='Accuracy rates',main=paste('Choose',dim(com_mat)[1],'out of 7 datasets as training datasets'))
legend('bottomright',legend=c(paste('Mean of accuracy rates=',round(mean(store),4)),paste('sd of accuracy rates=',round(sd(store),4))),cex=0.5)

qr<-rbind(com_mat,store)
colnames(qr)<-paste('Combination',seq(1,dim(com_mat)[2]))
#qr[,-dim(qr)[1]]<-as.character(qr[,-dim(qr)[1]])
return(qr)
}


####################Binda_top100#########
library(foreach)
Binda_top100<-foreach(j=1:6)%do%{
  com_mat<-combn(7, j)
  store<-list()
  for(i in 1:dim(com_mat)[2])
    
  {
    print(i)
    Reference_ind<-com_mat[,i]
    Train_Samples<-do.call(c,All_data_cname[Reference_ind])
    Test_Samples<-setdiff(All_cname,Train_Samples)
    
    Xtrain<-featureMatrix[Train_Samples,]
    Ytrain<-TRUE_LABELS_all[Train_Samples]
    Xtest<-featureMatrix[Test_Samples,]
    Ytest<-TRUE_LABELS_all[Test_Samples]
    
    
    thr = optimizeThreshold(Xtrain, Ytrain)
    Xtrain.b = dichotomize(Xtrain, thr)
    br<-binda.ranking(Xtrain.b, Ytrain)
    store[[i]]<-as.numeric(rownames(br)[1:100])
   
  }
  store<-do.call(rbind,store)
   return(store) 
}
    

Binda_top100[[1]]
qq<-do.call(rbind,Binda_top100)
pp<-apply(qq,1,list)

Binary_mat<-foreach(i=1:dim(qq)[1],.combine=rbind)%do%{
  fea<-as.numeric(colnames(featureMatrix))
  Dat<-cbind(fea,rep(0,dim(featureMatrix)[2]))
  Dat[which(fea %in% qq[i,]),2]<-1
  return(Dat[,2])
}
dim(Binary_mat)
colnames(Binary_mat)<-fea

heatmap.2(t(Binary_mat),Rowv=F,dendrogram='none',col = c("black", "white"), trace="column",density.info="none",tracecol='purple')
dim(Binary_mat)

kk<-colSums(Binary_mat)
length(which(kk>70))

Select_ind<-which(kk>1)
Select_ind<-seq(1,dim(featureMatrix)[2])

Xtrain.b.naive = ifelse(is.na(intensityMatrix(peaks)), 0, 1)

length(Select_ind)
data<-Xall.b[,Select_ind]

data<-Xtrain.b.naive

#data<-featureMatrix[,Selectind]
dmat<-dist(data, method="binary")
h_temp<-hclust( dmat , method="ward.D2")
memb <- cutree(h_temp, k = 2)
leftMSSA<-length(which(Ytrain[names(memb[memb==2])]==3))
leftMRSA<-length(which(Ytrain[names(memb[memb==2])]==1))
rightMSSA<-length(which(Ytrain[names(memb[memb==1])]==3))
rightMRSA<-length(which(Ytrain[names(memb[memb==1])]==1))
hc = as.dendrogram(hclust( dmat , method="ward.D2"))
labels_colors(hc)
labels_colors(hc) <- TRUE_LABELS_all[order.dendrogram(hc)]
labels_colors(hc)
plot(hc)
legend('top',legend=c(paste("Left:",leftMSSA,"MSSA,",leftMRSA,"MRSA"),paste("Right:",rightMSSA,"MSSA,",rightMRSA,"MRSA")),cex=0.8)

pr.out=prcomp(data, scale=FALSE)
pr_s<-summary(pr.out)
pr_v<-pr_s$importance[2,]
ff<-as.data.frame(pr.out$x[,c(1,2)])
cool<-TRUE_LABELS_all
ff<-cbind(ff,cool)
plot(ff[,1:2],   pch=19,col=cool,xlab=paste('PC1(',round(pr_v[1],4),'%)',sep=''),ylab=paste('PC2(',round(pr_v[2],4),'%)',sep=''))
legend('bottomleft',legend=c("MRSA","MSSA"),col=c(1,3),pch=16,cex=1)








