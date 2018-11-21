featureMat2<-function(spectra_input,ref,halfwindowsize=1){
  
  massList<-lapply(spectra_input,function(x){x@mass})
  intList<-lapply(spectra_input,function(x){x@intensity})
  temp_mat<-extractFeatureMatrix2(massList,intList, ref ,length(ref), length(spectra_input),halfwindowsize=halfwindowsize)
  rownames(temp_mat)<-ref
  
  return(temp_mat)
}
IntensityMat2<-function(Spectrum_list,RefPeaks,tolerance=0.01){
  featureMatrix<-foreach(i = 1:length(Spectrum_list),.combine=rbind)%do%{

  #RefPeaks<-rep(RefPeaks[1],length(RefPeaks))
  #Spectrum_list[[i]]@mass<-rep(RefPeaks[1],length(Spectrum_list[[i]]@mass))
  match_ind<-match.closest(RefPeaks,Spectrum_list[[i]]@mass,tolerance=tolerance)

  if(sum(duplicated(match_ind[!is.na(match_ind)]))>0){
    match_ind[!is.na(match_ind)][duplicated(match_ind[!is.na( match_ind)])]<-NA
  }

  #match_ind<-match.closest(RefPeaks,Spectrum_list[[i]]@mass)
  kk<-Spectrum_list[[i]]@mass[match_ind]

  #sum(duplicated( match_ind[!is.na( match_ind)]))
  #sum(duplicated( kk[!is.na( kk)]))
  kk_inten<-Spectrum_list[[i]]@intensity[match_ind]

  approxSpectra <- approxfun(Spectrum_list[[i]]@mass,Spectrum_list[[i]]@intensity, yleft=0L, yright=0L)

  Intensity<-kk_inten
  
  approxIntensity <- approxSpectra(RefPeaks[is.na(match_ind)])  
  #if(sum(is.na(approxIntensity))>0)
  #{approxIntensity[is.na(approxIntensity)]<-kk_inten[is.na(approxIntensity)]}
   Intensity[is.na(Intensity)]<-approxIntensity
  return(Intensity)
  
  }
return(featureMatrix)
}


library(dendextend)
HC_plot<-function(DATA,Truelabels,Method='COR',leg=T,MAIN='',input_cex=3){
  if(Method=='COR'){
    dissimarity<-1-abs(cor(t(DATA)))
    dmat<-as.dist(dissimarity)
  } else if(Method=='BIN'){
    dmat<-dist(DATA, method="binary")
  }
  
  h_temp<-hclust( dmat , method="ward.D2")
  memb <- cutree(h_temp, k = 2)
  names(memb)
  leftMSSA<-length(which(Truelabels[names(memb[memb==2])]==3))
  leftMRSA<-length(which(Truelabels[names(memb[memb==2])]==1))
  rightMSSA<-length(which(Truelabels[names(memb[memb==1])]==3))
  rightMRSA<-length(which(Truelabels[names(memb[memb==1])]==1))
  hc = as.dendrogram(hclust( dmat , method="ward.D2"))
  labels_colors(hc) <- Truelabels[order.dendrogram(hc)]
  node_type<-labels_colors(hc)
  #19 solid
  #1 open
  node_type[node_type==3]=19
  node_type[node_type==1]=19
  
  color_type=labels_colors(hc)
  #grey 8
  color_type[color_type==1]=8
  color_type[color_type==3]=1
  
  #set('labels',NULL)%>% 
  hc %>% set("leaves_pch", node_type) %>% set("leaves_cex", 3) %>% set("leaves_col",  color_type )%>%
  #hc %>% set("leaves_pch", color_type) %>% set("leaves_cex", 3) %>% set("leaves_col",  node_type )%>%
    set('clear_branches')%>%
    set('labels_cex',15)%>%

    set("branches_lwd",3)%>% plot(horiz = F,main=MAIN,cex.main=input_cex,adj=0, leaflab = "none", axes=F)
  #title(main=MAIN,adj=0,cex=100)
  
  if(leg){  legend('topright',legend=c("MRSA samples","MSSA samples",paste("Left:",leftMSSA,"MSSA,",leftMRSA,"MRSA"),paste("Right:",rightMSSA,"MSSA,",rightMRSA,"MRSA")),col=c("1","3","0","0"),pch=c(1,19),cex=0.6)}
  
  # ggd1 <- as.ggdend(hc)
  # library(ggplot2)
  # ggplot(ggd1)
}


LASSO_pick<-function(featureMatrix,TRUE_LABELS_all){
  
  cv.out<-cv.glmnet(featureMatrix,TRUE_LABELS_all,alpha=1,family="binomial")
  bestlam<-cv.out$lambda.min
  
  out<-glmnet(featureMatrix,TRUE_LABELS_all,alpha=1,family="binomial")
  lasso.coef<-predict(out,type='coefficients',s=bestlam)
  qq<-as.matrix(lasso.coef)[-1,]
  qq<-qq[qq!=0]
  qq2<-as.numeric(names(qq))
  return(qq2)
}
SVM_pick<-function(featureMatrix,TRUE_LABELS_all){
  
  cv.out<-cv.glmnet(featureMatrix,TRUE_LABELS_all,alpha=1,family="binomial")
  bestlam<-cv.out$lambda.min
  
  out<-svmrfeFeatureRanking(featureMatrix,TRUE_LABELS_all)
  
  return(out)
}



PCA_plot<-function(DATA,Truelabels)
 {
pr.out=prcomp(DATA, scale=FALSE)
pr_v<-summary(pr.out)$importance[2,]

ff<-as.data.frame(pr.out$x[,c(1,2)])
ff<-cbind(ff,Truelabels)
plot(ff[,1:2],   pch=19,main='Technical replicate of the reference data',col=Truelabels,xlab=paste('PC1(',round(pr_v[1],4),'%)',sep=''),ylab=paste('PC2(',round(pr_v[2],4),'%)',sep=''))
legend('bottomleft',legend=c("MRSA","MSSA"),col=c(1,3),pch=16,cex=0.6)
 }


Binda_pick<-function(featureMatrix,TRUE_LABELS_all){
  library(binda)
  Xtrain<-t(featureMatrix)
  Ytrain<-TRUE_LABELS_all
  names(Ytrain)<-Train_cname
  thr = optimizeThreshold(Xtrain, Ytrain)
  Xtrain.b = dichotomize(Xtrain, thr)
  br<-binda.ranking(Xtrain.b, Ytrain,verbose=F)
  #Ref_ind<-as.data.frame(br[1:M1,])$idx
  #TOP_FEA<-as.numeric(rownames(br)[1:M1])
  return(list(br=br,thr=thr))
}





RF_pick<-function(featureMatrix,TRUE_LABELS_all){
  
  library(randomForest)
  library("foreach")
  library("doSNOW")
  cl<-makeCluster(4, type="SOCK")
  registerDoSNOW(cl)
  #Approximately 3 minutes:
  system.time(rf_out <- foreach(ntree = rep(250, 4), .combine = combine, .packages = "randomForest") %dopar%  randomForest(x=t(featureMatrix), y=as.factor(TRUE_LABELS_all), ntree = ntree,mtry=100,importance=T))
  stopCluster(cl)
  
  
  Impotance_out<-rf_out$importance[order(rf_out$importance[,1],decreasing=T),]
  qq<-as.numeric(rownames(Impotance_out))
  return(qq)
}

###routine preprocessing for train and test using MALDIquant

#REF_peaks is for alignment 
FEATUREMATRIX_TRAIN_FUN<-function(Train_list,Train_cname){
  library(MALDIquant)
  new_spectra<-createSpect(Train_list,Train_cname)
  spectra<-trim(new_spectra)
  spectra3<-smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)
  spectra4 <- removeBaseline(spectra3, method='TopHat',halfWindowSize=10)
  spectra5<-calibrateIntensity(spectra4, method="TIC")
  #spectra6<-alignSpectra(spectra5,SNR=2,noiseMethod='MAD',halfWindowSize = 10,tolerance=0.02,reference=maybe)
  #########alignSpectra
  SNR_P<-10
  peaks<-detectPeaks(spectra5,method="MAD",halfWindowSize=20, SNR=SNR_P)
  reff<-referencePeaks(peaks, method=c("strict"), minFrequency=0.75,tolerance=0.02)
  WarpFun1<- determineWarpingFunctions(peaks,reference=reff,tolerance=0.02,method="lowess",plot=F,plotInteractive=F)
  spectra6_t<-warpMassSpectra(spectra5, WarpFun1)
  peaks1 <- detectPeaks(spectra6_t, SNR=2, halfWindowSize=20)
  ################M0
  REF_peaks<-referencePeaks(peaks1,method='strict',minFrequency = 0.75,tolerance = 0.02)
  peaks1 <- binPeaks(peaks1,method='strict',tolerance=0.02)
  peaks1<-filterPeaks(peaks1,minFrequency=0.75,mergeWhitelists = F)
  featureMatrix <- intensityMatrix(peaks1, spectra6_t)
  rownames(featureMatrix)<-Train_cname
  return(list(featureMatrix=featureMatrix,spectra6_t=spectra6_t,REF_peaks=REF_peaks))
}

#For testing data, only return aligned spectrum
FEATUREMATRIX_TEST_FUN<-function(Test_list,Test_cname,REF_peaks){
  new_spectra<-createSpect(Test_list,Test_cname)
  spectra<-trim(new_spectra)
  spectra3<-smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)
  spectra4 <- removeBaseline(spectra3, method='TopHat',halfWindowSize=10)
  spectra5<-calibrateIntensity(spectra4, method="TIC")
  #spectra6<-alignSpectra(spectra5,SNR=2,noiseMethod='MAD',halfWindowSize = 10,tolerance=0.02,reference=maybe)
  #########alignSpectra
  SNR_P<-10
  peaks2<-detectPeaks(spectra5,method="MAD",halfWindowSize=20, SNR=SNR_P)
  WarpFun2<- determineWarpingFunctions(peaks2,reference=REF_peaks,tolerance=0.02,method="lowess",plot=F,plotInteractive=F)
  spectra6<-warpMassSpectra(spectra5, WarpFun2)
  peaks2 <- detectPeaks(spectra6, SNR=2, halfWindowSize=20)
  #featureMatrix2<-featureMat(spectra6,TOP_FEA,Bin_range=Bin_range)
  #colnames(featureMatrix2)<-Test_cname
  #return(featureMatrix2)
  return(list(spectra6=spectra6,peaks2=peaks2))
}

FEATUREMATRIX_TEST_FUN_bench<-function(Test_list,Test_cname,REF_peaks){
  new_spectra<-createSpect(Test_list,Test_cname)
  spectra<-trim(new_spectra)
  spectra3<-smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)
  spectra4 <- removeBaseline(spectra3, method='TopHat',halfWindowSize=10)
  spectra5<-calibrateIntensity(spectra4, method="TIC")
  #spectra6<-alignSpectra(spectra5,SNR=2,noiseMethod='MAD',halfWindowSize = 10,tolerance=0.02,reference=maybe)
  #########alignSpectra
  SNR_P<-10
  peaks2<-detectPeaks(spectra5,method="MAD",halfWindowSize=20, SNR=SNR_P)
  
  #not align to the reference
  WarpFun2<- determineWarpingFunctions(peaks2,tolerance=0.02,method="lowess",plot=F,plotInteractive=F)
  #WarpFun22<- determineWarpingFunctions(peaks2,tolerance=0.02,method="lowess",plot=F,plotInteractive=F,reference=REF_peaks)
  #spectra66<-warpMassSpectra(spectra5, WarpFun22)
  #all.equal(spectra66,spectra6)
  
  spectra6<-warpMassSpectra(spectra5, WarpFun2)
  peaks2 <- detectPeaks(spectra6, SNR=2, halfWindowSize=20)
  #featureMatrix2<-featureMat(spectra6,TOP_FEA,Bin_range=Bin_range)
  #colnames(featureMatrix2)<-Test_cname
  #return(featureMatrix2)
  return(list(spectra6=spectra6,peaks2=peaks2))
}

#Four test FUN#############
##M1=100######
#M1<-100
halfwindowsize<-0.5
#I_RF:
I_RF_FUN<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels,M1=100){
  WF2_result<-FEATUREMATRIX_TRAIN_FUN(Train_list,Train_cname)
  REF_peaks<-WF2_result$REF_peaks
  featureMatrix1<-IntensityMat2( WF2_result$spectra6_t,REF_peaks@mass,tolerance=0.02)
  colnames(featureMatrix1)<-REF_peaks@mass
  rownames(featureMatrix1)<-Train_cname
  FEA_FULL<-as.numeric(colnames(featureMatrix1))
  RF_FEA<-RF_pick(featureMatrix=t(featureMatrix1),TRUE_LABELS_all=Train_labels)
  
  
  #################new data#########
  TEST_RE<-FEATUREMATRIX_TEST_FUN(Test_list,Test_cname,REF_peaks=REF_peaks)
  #directly use aligned spectra to form feature matrix
  featureMatrix2<-IntensityMat2(TEST_RE$spectra6,RefPeaks=REF_peaks@mass,tolerance=0.02)
  colnames(featureMatrix2)<-REF_peaks@mass
  rownames(featureMatrix2)<-Test_cname
  
  if(M1=='FEA_FULL'){M1=length(FEA_FULL)}
  
  Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
  TOP_FEA<-RF_FEA[1:M1]
  RF_model<-randomForest(featureMatrix1[,Ref_ind],as.factor(Train_labels))
  RF_pred<-predict(RF_model,newdata=featureMatrix2[,Ref_ind])
  cm<-table(RF_pred,Test_labels)
  acc<-sum(diag(cm))/sum(cm)
  return(list(acc=acc,TOP_FEA=TOP_FEA))
}
#I_BI:
I_BI_FUN<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels,M1=100){
  WF2_result<-FEATUREMATRIX_TRAIN_FUN(Train_list,Train_cname)
  featureMatrix<-WF2_result$featureMatrix
  REF_peaks<-WF2_result$REF_peaks
  FEA_FULL<-as.numeric(colnames(featureMatrix)) #M0
  featureMatrix<-IntensityMat2( WF2_result$spectra6_t,REF_peaks@mass,tolerance=0.02)
  colnames(featureMatrix)<-REF_peaks@mass
  rownames(featureMatrix)<-Train_cname
  FEA_FULL<-as.numeric(colnames(featureMatrix)) 
  Bind_re<-Binda_pick(t(featureMatrix),Train_labels)
  br<-Bind_re$br
  thr<-Bind_re$thr
  binda_train<-dichotomize(featureMatrix, thr)
  BI_FEA<-as.numeric(rownames(br))
  TEST_RE<-FEATUREMATRIX_TEST_FUN(Test_list,Test_cname,REF_peaks=REF_peaks)
  featureMatrix2<-IntensityMat2(TEST_RE$spectra6,FEA_FULL,tolerance=0.02)
  colnames(featureMatrix2)<-REF_peaks@mass
  rownames(featureMatrix2)<-Test_cname
  binda_test = dichotomize(featureMatrix2, thr)
  
  if(M1=='FEA_FULL'){M1=length(FEA_FULL)}
  
  Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
  TOP_FEA<-BI_FEA[1:M1]
  binda_train_for<-binda_train[,Ref_ind]
  binda_test_for<-binda_test[,Ref_ind]
  binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
  ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
  cm =table(Test_labels, ynew) 
  acc<-sum(diag(cm))/sum(cm)
  return(list(acc=acc,TOP_FEA=TOP_FEA))
  }
#P_RF:
P_RF_FUN<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels,M1=100){
  WF2_result<-FEATUREMATRIX_TRAIN_FUN(Train_list,Train_cname)
  featureMatrix<-WF2_result$featureMatrix
  REF_peaks<-WF2_result$REF_peaks
  FEA_FULL<-as.numeric(colnames(featureMatrix)) #M0
  
  system.time(featureMatrix1<-featureMat2(WF2_result$spectra6_t,FEA_FULL,halfwindowsize=halfwindowsize))
  
  
  
  
  colnames(featureMatrix1)<-Train_cname
  
  RF_FEA<-RF_pick(featureMatrix1,Train_labels)
  TEST_RE<-FEATUREMATRIX_TEST_FUN(Test_list,Test_cname,REF_peaks)
  
  
  featureMatrix2<-featureMat2(TEST_RE$spectra6,FEA_FULL,halfwindowsize=halfwindowsize)
  colnames(featureMatrix2)<-Test_cname
  
  if(M1=='FEA_FULL'){M1=length(FEA_FULL)}
  
  
  Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
  TOP_FEA<-RF_FEA[1:M1]
  RF_model<-randomForest(t(featureMatrix1[Ref_ind,]),as.factor(Train_labels))
  RF_pred<-predict(RF_model,newdata=t(featureMatrix2[Ref_ind,]))
  qq<-table(RF_pred,Test_labels)
  acc<-sum(diag(qq))/sum(qq)
  return(list(acc=acc,TOP_FEA=TOP_FEA))
}
#P_BI:
P_BI_FUN<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels,M1=100){
  WF2_result<-FEATUREMATRIX_TRAIN_FUN(Train_list,Train_cname)
  featureMatrix<-WF2_result$featureMatrix
  REF_peaks<-WF2_result$REF_peaks
  FEA_FULL<-as.numeric(colnames(featureMatrix)) #M0
  system.time(featureMatrix1<-featureMat2(WF2_result$spectra6_t,FEA_FULL,halfwindowsize=halfwindowsize))
  colnames(featureMatrix1)<-Train_cname
  Bind_re<-Binda_pick(featureMatrix1,Train_labels)
  br<-Bind_re$br
  thr<-Bind_re$thr
  BI_FEA<-as.numeric(rownames(br))
  #Test data
  TEST_RE<-FEATUREMATRIX_TEST_FUN(Test_list,Test_cname,REF_peaks)
  featureMatrix2<-featureMat2(TEST_RE$spectra6,FEA_FULL,halfwindowsize=halfwindowsize)
  colnames(featureMatrix2)<-Test_cname
  binda_train<-dichotomize(t(featureMatrix1), thr)
  binda_test<-dichotomize(t(featureMatrix2), thr)
  
  
  if(M1=='FEA_FULL'){M1=length(FEA_FULL)}
  Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
  TOP_FEA<-BI_FEA[1:M1]
  #################new data#########
  binda_test_for<-binda_test[,Ref_ind]
  binda_train_for<-binda_train[,Ref_ind]
  binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
  ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
  qq<-table(ynew,Test_labels)
  acc<-sum(diag(qq))/sum(qq)
  return(list(acc=acc,TOP_FEA=TOP_FEA))
}

#####WORKFLOW1 MALDI########
#MALDI_RF:#############
MALDI_RF_FUN<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels,M1=100){
  new_spectra<-createSpect(c(Train_list,Test_list),c(Train_cname,Test_cname))
  spectra<-trim(new_spectra)
  spectra3<-smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)
  spectra4 <- removeBaseline(spectra3, method='TopHat',halfWindowSize=10)
  spectra5<-calibrateIntensity(spectra4, method="TIC")
  #########alignSpectra
  SNR_P<-10
  peaks<-detectPeaks(spectra5,method="MAD",halfWindowSize=20, SNR=SNR_P)
  reff<-referencePeaks(peaks, method=c("strict"), minFrequency=0.75,tolerance=0.02)
  WarpFun1<- determineWarpingFunctions(peaks,reference=reff,tolerance=0.02,method="lowess",plot=F,plotInteractive=F)
  spectra6_t<-warpMassSpectra(spectra5, WarpFun1)
  peaks1 <- detectPeaks(spectra6_t, SNR=2, halfWindowSize=20)
  ################M0
  REF_peaks<-referencePeaks(peaks1,method='strict',minFrequency = 0.75,tolerance = 0.02)
  peaks1 <- binPeaks(peaks1,method='strict',tolerance=0.02)
  peaks1<-filterPeaks(peaks1,minFrequency=0.75,mergeWhitelists = F)
  featureMatrix <- intensityMatrix(peaks1, spectra6_t)
  rownames(featureMatrix)<-c(Train_cname,Test_cname)
  
  FEA_FULL<-as.numeric(colnames(featureMatrix))
  
  
  if(M1=='FEA_FULL'){M1=length(FEA_FULL)}
  
  RF_FEA<-RF_pick(featureMatrix=t(featureMatrix[Train_cname,]),TRUE_LABELS_all=Train_labels)

  Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
  
  RF_model<-randomForest(featureMatrix[Train_cname,Ref_ind],as.factor(Train_labels))
  RF_pred<-predict(RF_model,newdata=featureMatrix[Test_cname,Ref_ind])
  cm<-table(RF_pred,Test_labels)
  acc<-sum(diag(cm))/sum(cm)
  return(acc)
}



#MALDI_BI:#############
MALDI_BI_FUN<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels,M1=100){
  new_spectra<-createSpect(c(Train_list,Test_list),c(Train_cname,Test_cname))
  spectra<-trim(new_spectra)
  spectra3<-smoothIntensity(spectra, method="MovingAverage",halfWindowSize=2)
  spectra4 <- removeBaseline(spectra3, method='TopHat',halfWindowSize=10)
  spectra5<-calibrateIntensity(spectra4, method="TIC")
  #########alignSpectra
  SNR_P<-10
  peaks<-detectPeaks(spectra5,method="MAD",halfWindowSize=20, SNR=SNR_P)
  reff<-referencePeaks(peaks, method=c("strict"), minFrequency=0.75,tolerance=0.02)
  WarpFun1<- determineWarpingFunctions(peaks,reference=reff,tolerance=0.02,method="lowess",plot=F,plotInteractive=F)
  spectra6_t<-warpMassSpectra(spectra5, WarpFun1)
  peaks1 <- detectPeaks(spectra6_t, SNR=2, halfWindowSize=20)
  ################M0
  REF_peaks<-referencePeaks(peaks1,method='strict',minFrequency = 0.75,tolerance = 0.02)
  peaks1 <- binPeaks(peaks1,method='strict',tolerance=0.02)
  peaks1<-filterPeaks(peaks1,minFrequency=0.75,mergeWhitelists = F)
  featureMatrix <- intensityMatrix(peaks1, spectra6_t)
  rownames(featureMatrix)<-c(Train_cname,Test_cname)
  
  
  Bind_re<-Binda_pick(t(featureMatrix[Train_cname,]),Train_labels)
  br<-Bind_re$br
  BI_FEA<-as.numeric(rownames(br))
  thr<-Bind_re$thr
  
  binda_train<-dichotomize(featureMatrix[Train_cname,], thr)
  binda_test = dichotomize(featureMatrix[Test_cname,], thr)
  FEA_FULL<-as.numeric(colnames(featureMatrix))
  
  
  if(M1=='FEA_FULL'){M1=length(FEA_FULL)}
  
  TOP_FEA<-BI_FEA[1:M1]
  
  Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
  binda_train_for<-binda_train[,Ref_ind]
  binda_test_for<-binda_test[,Ref_ind]
  binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
  ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
  cm =table(Test_labels, ynew) 
  acc<-sum(diag(cm))/sum(cm)
  return(acc)
}
