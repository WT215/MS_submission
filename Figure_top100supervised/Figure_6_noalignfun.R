M1<-100
halfwindowsize<-0.5
#I_RF:
I_RF_FUN2<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
  WF2_result<-FEATUREMATRIX_TRAIN_FUN(Train_list,Train_cname)
  REF_peaks<-WF2_result$REF_peaks
  featureMatrix1<-IntensityMat2( WF2_result$spectra6_t,REF_peaks@mass,tolerance=0.02)
  colnames(featureMatrix1)<-REF_peaks@mass
  rownames(featureMatrix1)<-Train_cname
  FEA_FULL<-as.numeric(colnames(featureMatrix1))
  RF_FEA<-RF_pick(featureMatrix=t(featureMatrix1),TRUE_LABELS_all=Train_labels)
  
  
  #################new data#########
  TEST_RE<-FEATUREMATRIX_TEST_FUN_bench(Test_list,Test_cname,REF_peaks=REF_peaks)
  #directly use aligned spectra to form feature matrix
  featureMatrix2<-IntensityMat2(TEST_RE$spectra6,RefPeaks=REF_peaks@mass,tolerance=0.02)
  colnames(featureMatrix2)<-REF_peaks@mass
  rownames(featureMatrix2)<-Test_cname
  Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
  TOP_FEA<-RF_FEA[1:M1]
  RF_model<-randomForest(featureMatrix1[,Ref_ind],as.factor(Train_labels))
  RF_pred<-predict(RF_model,newdata=featureMatrix2[,Ref_ind])
  cm<-table(RF_pred,Test_labels)
  acc<-sum(diag(cm))/sum(cm)
  return(acc)
}
#I_BI:
I_BI_FUN2<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  TEST_RE<-FEATUREMATRIX_TEST_FUN_bench(Test_list,Test_cname,REF_peaks=REF_peaks)
  featureMatrix2<-IntensityMat2(TEST_RE$spectra6,FEA_FULL,tolerance=0.02)
  colnames(featureMatrix2)<-REF_peaks@mass
  rownames(featureMatrix2)<-Test_cname
  binda_test = dichotomize(featureMatrix2, thr)
  Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
  TOP_FEA<-BI_FEA[1:M1]
  binda_train_for<-binda_train[,Ref_ind]
  binda_test_for<-binda_test[,Ref_ind]
  binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
  ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
  cm =table(Test_labels, ynew) 
  acc<-sum(diag(cm))/sum(cm)
  return(acc)}
#P_RF:
P_RF_FUN2<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
  WF2_result<-FEATUREMATRIX_TRAIN_FUN(Train_list,Train_cname)
  featureMatrix<-WF2_result$featureMatrix
  REF_peaks<-WF2_result$REF_peaks
  FEA_FULL<-as.numeric(colnames(featureMatrix)) #M0
  
  system.time(featureMatrix1<-featureMat2(WF2_result$spectra6_t,FEA_FULL,halfwindowsize=halfwindowsize))
  colnames(featureMatrix1)<-Train_cname
  RF_FEA<-RF_pick(featureMatrix1,Train_labels)
  TEST_RE<-FEATUREMATRIX_TEST_FUN_bench(Test_list,Test_cname,REF_peaks)
  featureMatrix2<-featureMat2(TEST_RE$spectra6,FEA_FULL,halfwindowsize=halfwindowsize)
  colnames(featureMatrix2)<-Test_cname
  Ref_ind<-which(FEA_FULL %in% RF_FEA[1:M1])
  TOP_FEA<-RF_FEA[1:M1]
  RF_model<-randomForest(t(featureMatrix1[Ref_ind,]),as.factor(Train_labels))
  RF_pred<-predict(RF_model,newdata=t(featureMatrix2[Ref_ind,]))
  qq<-table(RF_pred,Test_labels)
  acc<-sum(diag(qq))/sum(qq)
  return(acc)
}
#P_BI:

P_BI_FUN2<-function(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels){
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
  TEST_RE<-FEATUREMATRIX_TEST_FUN_bench(Test_list,Test_cname,REF_peaks)
  featureMatrix2<-featureMat2(TEST_RE$spectra6,FEA_FULL,halfwindowsize=halfwindowsize)
  colnames(featureMatrix2)<-Test_cname
  binda_train<-dichotomize(t(featureMatrix1), thr)
  binda_test<-dichotomize(t(featureMatrix2), thr)
  Ref_ind<-which(FEA_FULL %in% BI_FEA[1:M1])
  TOP_FEA<-BI_FEA[1:M1]
  #################new data#########
  binda_test_for<-binda_test[,Ref_ind]
  binda_train_for<-binda_train[,Ref_ind]
  binda.out = binda(binda_train_for, Train_labels, verbose=FALSE)
  ynew = predict.binda(binda.out, binda_test_for, verbose=FALSE)$class
  qq<-table(ynew,Test_labels)
  acc<-sum(diag(qq))/sum(qq)
  return(acc)
}