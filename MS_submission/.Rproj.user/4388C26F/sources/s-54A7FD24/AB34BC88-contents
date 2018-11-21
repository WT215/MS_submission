

#Four test FUN#############
##M1=100######
#M1<-100
halfwindowsize<-0.5
#I_RF:
I_RF_FUN_M1<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  
  TOP_FEA_list<-list()
  
  acc_vec<-foreach(mmm = 1:length(M1_vec),.combine=c)%do%{
    print(mmm)
    
    
    M1=M1_vec[mmm]
    
    if(M1=='FEA_FULL'){M1=length(FEA_FULL)} else{M1=as.numeric(M1)}
    
    Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
    TOP_FEA<-RF_FEA[1:M1]
    TOP_FEA_list[[mmm]]<-TOP_FEA
    
    RF_model<-randomForest(featureMatrix1[,Ref_ind],as.factor(Train_labels))
    RF_pred<-predict(RF_model,newdata=featureMatrix2[,Ref_ind])
    cm<-table(RF_pred,Test_labels)
    acc<-sum(diag(cm))/sum(cm)
    return(acc)
  }
  names(acc_vec)<-M1_vec
  names(TOP_FEA_list)<-M1_vec
  return(list(acc_vec=acc_vec,TOP_FEA_list=TOP_FEA_list))
}
#I_BI:
I_BI_FUN_M1<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  
  
  
  TOP_FEA_list<-list()
  
  acc_vec<-foreach(mmm = 1:length(M1_vec),.combine=c)%do%{
    
    print(mmm)
    
    
    M1=M1_vec[mmm]
    if(M1=='FEA_FULL'){M1=length(FEA_FULL)} else{M1=as.numeric(M1)}
    
    Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
    TOP_FEA<-BI_FEA[1:M1]
    TOP_FEA_list[[mmm]]<-TOP_FEA
    
    binda_train_for<-binda_train[,Ref_ind]
    binda_test_for<-binda_test[,Ref_ind]
    binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
    ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
    cm =table(Test_labels, ynew) 
    acc<-sum(diag(cm))/sum(cm)
    return(acc)
  }
  
  names(acc_vec)<-M1_vec
  names(TOP_FEA_list)<-M1_vec
  return(list(acc_vec=acc_vec,TOP_FEA_list=TOP_FEA_list))
  #return(list(acc=acc,TOP_FEA=TOP_FEA))
}
#P_RF:
P_RF_FUN_M1<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  
  
  TOP_FEA_list<-list()
  
  acc_vec<-foreach(mmm = 1:length(M1_vec),.combine=c)%do%{
    print(mmm)
    
    
    M1=M1_vec[mmm]
    
    
    if(M1=='FEA_FULL'){M1=length(FEA_FULL)} else{M1=as.numeric(M1)}
    
    
    Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
    TOP_FEA<-RF_FEA[1:M1]
    TOP_FEA_list[[mmm]]<-TOP_FEA
    
    RF_model<-randomForest(t(featureMatrix1[Ref_ind,]),as.factor(Train_labels))
    RF_pred<-predict(RF_model,newdata=t(featureMatrix2[Ref_ind,]))
    qq<-table(RF_pred,Test_labels)
    acc<-sum(diag(qq))/sum(qq)
    return(acc)
  }
  
  names(acc_vec)<-M1_vec
  names(TOP_FEA_list)<-M1_vec
  
  return(list(acc_vec=acc_vec,TOP_FEA_list=TOP_FEA_list))
  #return(list(acc=acc,TOP_FEA=TOP_FEA))
}


#P_BI:
P_BI_FUN_M1<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  
  TOP_FEA_list<-list()
  
  acc_vec<-foreach(mmm = 1:length(M1_vec),.combine=c)%do%{
    print(mmm)
    
    
    M1=M1_vec[mmm]
    
    if(M1=='FEA_FULL'){M1=length(FEA_FULL)} else{M1=as.numeric(M1)}
    Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
    TOP_FEA<-BI_FEA[1:M1]
    
    TOP_FEA_list[[mmm]]<-TOP_FEA
    
    #################new data#########
    binda_test_for<-binda_test[,Ref_ind]
    binda_train_for<-binda_train[,Ref_ind]
    binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
    ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
    qq<-table(ynew,Test_labels)
    acc<-sum(diag(qq))/sum(qq)
  }
  
  names(acc_vec)<-M1_vec
  names(TOP_FEA_list)<-M1_vec
  return(list(acc_vec=acc_vec,TOP_FEA_list=TOP_FEA_list))
  
  
  #return(list(acc=acc,TOP_FEA=TOP_FEA))
}

#####WORKFLOW1 MALDI########
#MALDI_RF:#############
MALDI_RF_FUN_M1<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  
  
  acc_vec<-foreach(mmm = 1:length(M1_vec),.combine=c)%do%{
    print(mmm)
    
    
    M1=M1_vec[mmm]
    if(M1=='FEA_FULL'){M1=length(FEA_FULL)} else{M1=as.numeric(M1)}
    
    RF_FEA<-RF_pick(featureMatrix=t(featureMatrix[Train_cname,]),TRUE_LABELS_all=Train_labels)
    
    Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
    
    RF_model<-randomForest(featureMatrix[Train_cname,Ref_ind],as.factor(Train_labels))
    RF_pred<-predict(RF_model,newdata=featureMatrix[Test_cname,Ref_ind])
    cm<-table(RF_pred,Test_labels)
    acc<-sum(diag(cm))/sum(cm)
    return(acc)
  }
  names(acc_vec)<-M1_vec
  return(acc_vec)
  
  
}



#MALDI_BI:#############
MALDI_BI_FUN_M1<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  
  
  
  acc_vec<-foreach(mmm = 1:length(M1_vec),.combine=c)%do%{
    print(mmm)
    
    
    M1=M1_vec[mmm]
    if(M1=='FEA_FULL'){M1=length(FEA_FULL)} else{M1=as.numeric(M1)}
    
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
  names(acc_vec)<-M1_vec
  return(acc_vec)
}
