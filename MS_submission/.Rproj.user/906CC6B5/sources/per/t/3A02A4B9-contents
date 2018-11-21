source("E:/MainProject/MS_submission/Functions/4BIO_M1_FUNCTIONS.R")

#leave one out MALDIquant first workflow for supervised learning
source("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/MALDI_SUPER_LOOV_funs.R")

#########LOAD FUNCTIONS############
##########BEGIN OF PREPARE
library("MALDIquant")
library('MALDIquantForeign')
library("readMzXmlData")
library(dendextend)
library(gplots)
library(randomForest)
library(foreach)
library("doSNOW")
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

#Rcpp::sourceCpp('E:/MainProject/PROJECT_MRSA_MSSA/C_codes/extractFeatureMatrix.cpp')
library(MALDIrcpp)

load("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/LOAD_7_MRSAMSSA_data.RData")
data_name<-c('pc','pa','r1','r2','pa15062016','R1','R2')
All_data_list<-list(pc_list,pa_list,r1_list,r2_list,pa15062016_list,MRMS_050816run1_list,MRMS_050816run2_list)
names(All_data_list)<-data_name
All_data_cname<-list(pc_name,pa_name,r1_name,r2_name,pa15062016_name,MRMS_050816run1_name,MRMS_050816run2_name)
TRUE_LABELS<-c(rep(1,10),rep(3,10))
source("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/FUNCTIONS.R")
#source("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/LOAD_BASICFUN_DATA.r")


library(foreach)
names(All_data_list)
#1 2
#3 4
#5
#6 7
DATA_IND_LIST<-list(c(1,3,5,6),c(1,4,5,6),c(1,3,5,7),c(1,4,5,7),c(2,3,5,6),c(2,3,4,7),c(2,4,5,6),c(2,4,5,7))





########################Begin analysis####################
names(All_data_list)[3:4]<-c('rr1','rr2')

DATA_IND_LIST




#dir.create("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/split_samples_many/try")




# #generate 10 splits####
# set.seed(12300)
# eightind_traintest_list<-list()
# for(eightind in 1:length(DATA_IND_LIST))
# {
# Bio4<-All_data_list[DATA_IND_LIST[[eightind]]]
# All4bio_cname<-do.call(c,All_data_cname[DATA_IND_LIST[[eightind]]])
# Bio4_samples<-unlist(Bio4,recursive = F)
# 
# TrueLabels_4bio<-rep(c(rep(1,10),rep(3,10)),4)
# names(TrueLabels_4bio)<-names(Bio4_samples)
# 
# 
# trainind_list<-list()
# testind_list<-list()
# Train_labels_list<-list()
# Test_labels_list<-list()
# for(i in 1:10){
#   trainind_list[[i]]<-c(c(sample(seq(1,10),7),sample(seq(11,20),7)),
#                         c(sample(seq(1,10),7)+20,sample(seq(11,20),7)+20),
#                         c(sample(seq(1,10),7)+40,sample(seq(11,20),7)+40),
#                         c(sample(seq(1,10),7)+60,sample(seq(11,20),7)+60)
#   )
#   testind_list[[i]]<-seq(1,length(Bio4_samples))[-trainind_list[[i]]]
#     
#     
#     Train_cname<-All4bio_cname[trainind_list[[i]]]
#   Test_cname<-All4bio_cname[testind_list[[i]]]
#   
#   Train_list<-Bio4_samples[trainind_list[[i]]]
#   Test_list<-Bio4_samples[testind_list[[i]]]
#   Train_labels_list[[i]]<-TrueLabels_4bio[names(Train_list)]
#   Test_labels_list[[i]]<-TrueLabels_4bio[names(Test_list)]
#     
#     
#     eightind_traintest_list[[eightind]]<-list(trainind_list=trainind_list,testind_list=testind_list,Train_labels_list=Train_labels_list,Test_labels_list=Test_labels_list)
# }
# 
# 
# }
# 
# names(eightind_traintest_list)<-unlist(lapply(DATA_IND_LIST,function(x){
#   paste('FourBio_',paste0(data_name[x],collapse = ''),sep='')}))
# 
# 
# eightind_traintest_list$FourBio_pcr1pa15062016R1
# 
# save(eightind_traintest_list,file="E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/split_samples_many/eightind_traintest_list.RData")

load("E:/MainProject/MS_submission/Figure_numfea/eightind_traintest_list.RData")
###begin analysis#########
listt8<-list()

for(i in 1:length(DATA_IND_LIST)){
  print(paste('the',i,'th group'),sep='')
  
  Bio4<-All_data_list[DATA_IND_LIST[[i]]]
  All4bio_cname<-do.call(c,All_data_cname[DATA_IND_LIST[[i]]])
  
  Bio4_samples<-unlist(Bio4,recursive = F)
  
  
  samples_results_list<-list()
  for(sampleind in 1:10){
    
    
    
    Train_ind<-eightind_traintest_list[[i]]$trainind_list[[sampleind]]
    Test_ind<-eightind_traintest_list[[i]]$testind_list[[sampleind]]
    
    Train_cname<-All4bio_cname[Train_ind]
    Test_cname<-All4bio_cname[Test_ind]
    
    Train_list<-Bio4_samples[Train_ind]
    Test_list<-Bio4_samples[Test_ind]
    Train_labels<-eightind_traintest_list[[i]]$Train_labels_list[[sampleind]]
    Test_labels<-eightind_traintest_list[[i]]$Test_labels_list[[sampleind]]
    
    M1_vec<-c(5,10,15,30,50,100,300,500,1000,'FEA_FULL')
  #acc_I_BI<-I_BI_FUN_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)

   acc_I_RF<-I_RF_FUN_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)
  #
  #
  #acc_P_BI<-P_BI_FUN_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)

  #acc_P_RF<-P_RF_FUN_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)
  
  # acc_MALDI_BI<-MALDI_BI_FUN_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)
  # 
  # acc_MALDI_RF<-MALDI_RF_FUN_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)
  
  # acc_MALDI_BI<-MALDI_BI_LOOV_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)
  # 
  # acc_MALDI_RF<-MALDI_RF_LOOV_M1(Train_list,Train_cname,Train_labels,Test_list,Test_cname,Test_labels)
  
  #listt8=list(acc_I_BI=acc_I_BI,acc_I_RF=acc_I_RF,acc_P_BI=acc_P_BI,acc_P_RF=acc_P_RF,acc_MALDI_BI=acc_MALDI_BI,acc_MALDI_RF=acc_MALDI_RF)
   samples_results_list[[sampleind]]<-list(acc_I_RF=acc_I_RF)
  }
  assign(paste('listt8mIRF_',paste0(data_name[DATA_IND_LIST[[i]]],collapse = ''),sep=''), samples_results_list)
  save(list=paste('listt8mIRF_',paste0(data_name[DATA_IND_LIST[[i]]],collapse = ''),sep=''),file=paste('E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/split_samples_many/listtm_IRF/listt8mIRF_',paste0(data_name[DATA_IND_LIST[[i]]],collapse = ''),'.RData',sep=''))
  
}






###FIUGRE4: variation of number of features used in statsitical learning #############


####Figure 4 splitting datasets##########
load("E:/MainProject/PROJECT_MRSA_MSSA/DRAFT/WF2_RFfixbinsize_FIGURES/halfwinsize0.5/Data4_P_RF.RData")
load("E:/MainProject/PROJECT_MRSA_MSSA/DRAFT/WF2_RFfixbinsize_FIGURES/halfwinsize0.5/Data4_P_BI.RData")
load("E:/MainProject/PROJECT_MRSA_MSSA/DRAFT/WF2_PEAKS_FIGURES/Data4_I_BI.RData")
load("E:/MainProject/PROJECT_MRSA_MSSA/DRAFT/WF2_PEAKS_FIGURES/Data4_I_RF.RData")

load("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Figure4/acc_mat32_maldiloov.RData")

#load("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Figure4/ACC_mat_maldi.RData")


ACC_mat_maldi_BI<-t(acc_mat32[,seq(1,64,2)])
ACC_mat_maldi_RF<-t(acc_mat32[,-seq(1,64,2)])

ACC_mat_maldi_BI<-rbind(ACC_mat_maldi_BI,colMeans(ACC_mat_maldi_BI))
ACC_mat_maldi_RF<-rbind(ACC_mat_maldi_RF,colMeans(ACC_mat_maldi_RF))

qq<-Data4_I_BI[1:(dim(Data4_I_BI)[1]-1),]
Data4_I_BI<-rbind(Data4_I_BI,apply(Data4_I_BI[1:(dim(Data4_I_BI)[1]-1),],2,function(x){sd(x)/32}))
Data4_I_RF<-rbind(Data4_I_RF,apply(Data4_I_RF[1:(dim(Data4_I_RF)[1]-1),],2,function(x){sd(x)/32}))
Data4_P_BI<-rbind(Data4_P_BI,apply(Data4_P_BI[1:(dim(Data4_P_BI)[1]-1),],2,function(x){sd(x)/32}))
Data4_P_RF<-rbind(Data4_P_RF,apply(Data4_P_RF[1:(dim(Data4_P_RF)[1]-1),],2,function(x){sd(x)/32}))

rownames(Data4_I_BI)[length(rownames(Data4_I_BI))]<-'se of mean'
rownames(Data4_I_RF)[length(rownames(Data4_I_RF))]<-'se of mean'
rownames(Data4_P_BI)[length(rownames(Data4_P_BI))]<-'se of mean'
rownames(Data4_P_RF)[length(rownames(Data4_P_RF))]<-'se of mean'







library(Rmisc)


#geom_errorbar(aes(ymin=len-se, ymax=len+se),width=.2, position=position_dodge(.9))
MEAN_IND<-which(rownames(Data4_I_BI)=='Average Mean Accuracy rate')
SE_IND<-which(rownames(Data4_I_BI)=='se of mean')

Bar_dat<-cbind(Data4_I_BI[MEAN_IND,],Data4_I_RF[MEAN_IND,],Data4_P_BI[MEAN_IND,],Data4_P_RF[MEAN_IND,],ACC_mat_maldi_BI[33,],ACC_mat_maldi_RF[33,])



SE_dat<-cbind(Data4_I_BI[SE_IND,],Data4_I_RF[SE_IND,],Data4_P_BI[SE_IND,],Data4_P_RF[SE_IND,],apply(ACC_mat_maldi_BI[seq(1,32),],2,sd)/32,apply(ACC_mat_maldi_RF[seq(1,32),],2,sd)/32)
Workflow<-c('I_BI','I_RF','P_BI','P_RF','MALDI_BI','MALDI_RF')



Bar_dat2<-foreach(i=1:dim(Bar_dat)[1],.combine=rbind)%do%{
    
    kk<-Bar_dat[i,]
    kk2<-cbind(kk,SE_dat[i,],Workflow,rep(rownames(Bar_dat)[i],6))
    return(kk2)
}



Bar_dat2<-as.data.frame(Bar_dat2)
colnames(Bar_dat2)<-c('AccuracyRate','SE','Workflow','NumberofTop')
Bar_dat2$AccuracyRate<-as.numeric(as.character(Bar_dat2$AccuracyRate))
Bar_dat2$SE<-as.numeric(as.character(Bar_dat2$SE))
Bar_dat2$AccuracyRate<-round(Bar_dat2$AccuracyRate,3)
Bar_dat2$NumberofTop<-factor(Bar_dat2$NumberofTop,levels=unique(Bar_dat2$NumberofTop))

Bar_dat2$Workflow<-factor(Bar_dat2$Workflow,levels=unique(Bar_dat2$Workflow))


textsize<-15

library(ggplot2)
#Averaged across 32 testing results

dropp<-grep(Bar_dat2$Workflow,pattern='MALDI_')
Bar_dat2<-Bar_dat2[-dropp,]

ori_Fi4A<-ggplot(data=Bar_dat2, aes(x=Bar_dat2$NumberofTop, y=Bar_dat2$AccuracyRate, fill=Bar_dat2$Workflow)) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    scale_y_continuous(limits=c(0,1))+
    geom_errorbar(aes(ymin=Bar_dat2$AccuracyRate-Bar_dat2$SE, ymax=Bar_dat2$AccuracyRate+Bar_dat2$SE),size=0.5,width=.5,position=position_dodge(.9),color='green')+
    geom_text(aes(label=round(Bar_dat2$AccuracyRate,2)), vjust=0.1, color="black", position = position_dodge(0.9), size=4,angle=90,hjust=-0.1)+
    labs(x = "Number of top ranked features used",y='Averaged accuracy rates',fill='Methods')+
    scale_fill_grey()+
    ggtitle('A')+
    theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top",legend.key.size = unit(0.5,"line"),plot.title = element_text(size=30))

ori_Fi4A

####Figure 4 splitting SAMPLES##########
load("E:/MainProject/MS_submission/Figure_numfea/listtm_par1pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_pcr1pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_pcr1pa15062016R2.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_pcr2pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_pcr2pa15062016R2.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_par1r2R2.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_par2pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_numfea/listtm_par2pa15062016R2.RData")

#load correct I_RF results###
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_numfea/listtm_IRF/', pattern="listt8mIRF_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_numfea/listtm_IRF/',x,sep=''),.GlobalEnv)
})



for(i in 1:10){
  listtm_par1pa15062016R1[[i]]$acc_I_RF<-listt8mIRF_par1pa15062016R1[[i]]$acc_I_RF
}
for(i in 1:10){
  listtm_pcr1pa15062016R1[[i]]$acc_I_RF<-listt8mIRF_pcr1pa15062016R1[[i]]$acc_I_RF
}
for(i in 1:10){
  listtm_pcr1pa15062016R2[[i]]$acc_I_RF<-listt8mIRF_pcr1pa15062016R2[[i]]$acc_I_RF
}
for(i in 1:10){
  listtm_pcr2pa15062016R1[[i]]$acc_I_RF<-listt8mIRF_pcr2pa15062016R1[[i]]$acc_I_RF
}
for(i in 1:10){
  listtm_pcr2pa15062016R2[[i]]$acc_I_RF<-listt8mIRF_pcr2pa15062016R2[[i]]$acc_I_RF
}


for(i in 1:10){
  listtm_par1r2R2[[i]]$acc_I_RF<-listt8mIRF_par1r2R2[[i]]$acc_I_RF
}
for(i in 1:10){
  listtm_par2pa15062016R1[[i]]$acc_I_RF<-listt8mIRF_par2pa15062016R1[[i]]$acc_I_RF
}
for(i in 1:10){
  listtm_par2pa15062016R2[[i]]$acc_I_RF<-listt8mIRF_par2pa15062016R2[[i]]$acc_I_RF
}












list_all8<-list(listtm_par1pa15062016R1=listtm_par1pa15062016R1,
                listtm_pcr1pa15062016R1=listtm_pcr1pa15062016R1,
                listtm_pcr1pa15062016R2=listtm_pcr1pa15062016R2,
                listtm_pcr2pa15062016R1=listtm_pcr2pa15062016R1,
                listtm_pcr2pa15062016R2=listtm_pcr2pa15062016R2,
                listtm_par1r2R2=listtm_par1r2R2,
                listtm_par2pa15062016R1=listtm_par2pa15062016R1,
                listtm_par2pa15062016R2=listtm_par2pa15062016R2)


list_all8$listtm_pcr1pa15062016R2[[2]]$acc_MALDI_RF

library(foreach)
acc_I_BI_mat<-foreach(i=1:8,.combine=cbind)%do%{
  
  qq<-do.call(cbind,lapply(list_all8[[i]],function(x){
    x$acc_I_BI$acc_vec
  }))
    return(qq)
    
}

acc_I_RF_mat<-foreach(i=1:8,.combine=cbind)%do%{
  
  qq<-do.call(cbind,lapply(list_all8[[i]],function(x){
    x$acc_I_RF$acc_vec
  }))
  return(qq)
  
}
acc_P_BI_mat<-foreach(i=1:8,.combine=cbind)%do%{
  
  qq<-do.call(cbind,lapply(list_all8[[i]],function(x){
    x$acc_P_BI$acc_vec
  }))
  return(qq)
  
}

acc_P_RF_mat<-foreach(i=1:8,.combine=cbind)%do%{
  
  qq<-do.call(cbind,lapply(list_all8[[i]],function(x){
    x$acc_P_RF$acc_vec
  }))
  return(qq)
  
}


acc_MALDI_BI_mat<-foreach(i=1:8,.combine=cbind)%do%{
  
  qq<-do.call(cbind,lapply(list_all8[[i]],function(x){
    x$acc_MALDI_BI
  }))
  return(qq)
  
}
acc_MALDI_RF_mat<-foreach(i=1:8,.combine=cbind)%do%{
  
  qq<-do.call(cbind,lapply(list_all8[[i]],function(x){
    x$acc_MALDI_RF
  }))
  return(qq)
  
}
M1_vec<-c(5,10,15,30,50,100,300,500,1000,'Full Features')



acc_8mat_list<-list(acc_I_BI_mat=acc_I_BI_mat,
                    acc_I_RF_mat=acc_I_RF_mat,
                    acc_P_BI_mat=acc_P_BI_mat,
                    acc_P_RF_mat=acc_P_RF_mat,
                    acc_MALDI_BI_mat=acc_MALDI_BI_mat,
                    acc_MALDI_RF_mat=acc_MALDI_RF_mat)
Workflow<-c('I_BI','I_RF','P_BI','P_RF','MALDI_BI','MALDI_RF')

BAR_PLOT<-foreach(i=1:length(acc_8mat_list),.combine=rbind)%do%{
  
  qq_mean<-rowMeans(acc_8mat_list[[i]])
  qq_se<-apply(acc_8mat_list[[i]],1,function(x){sd(x)/8})
  qq_out<-cbind(qq_mean,qq_se,M1_vec,rep(Workflow[i],10))
  
}
BAR_PLOT<-as.data.frame(BAR_PLOT)
colnames(BAR_PLOT)<-c('AccuracyRate','SE','NumberofTop','Workflow')


BAR_PLOT$AccuracyRate<-as.numeric(as.character(BAR_PLOT$AccuracyRate))
BAR_PLOT$AccuracyRate<-round(BAR_PLOT$AccuracyRate,3)

BAR_PLOT$SE<-as.numeric(as.character(BAR_PLOT$SE))
BAR_PLOT$NumberofTop<-factor(BAR_PLOT$NumberofTop,levels=unique(BAR_PLOT$NumberofTop))
BAR_PLOT$Workflow<-factor(BAR_PLOT$Workflow,levels=unique(BAR_PLOT$Workflow))





textsize<-15

library(ggplot2)
#Averaged across 32 testing results

droppp2<-grep(BAR_PLOT$Workflow,pattern='MALDI_')
BAR_PLOT<-BAR_PLOT[-droppp2,]
new_Fi4A<-ggplot(data=BAR_PLOT, aes(x=BAR_PLOT$NumberofTop, y=BAR_PLOT$AccuracyRate, fill=BAR_PLOT$Workflow)) +
  geom_bar(stat="identity", position = position_dodge(0.9),width=0.9,show.legend =F)+
  scale_y_continuous(limits=c(0,1))+
  geom_errorbar(aes(ymin=BAR_PLOT$AccuracyRate-BAR_PLOT$SE, ymax=BAR_PLOT$AccuracyRate+BAR_PLOT$SE),size=0.5,width=.5,position=position_dodge(.9),color='green')+
  geom_text(aes(label=round(BAR_PLOT$AccuracyRate,2)), vjust=0.1, color="black", position = position_dodge(0.9), size=4,angle=90,hjust=-0.1)+
  labs(x = "Number of top ranked features used",y='Averaged accuracy rates',fill='Methods')+
  scale_fill_grey()+
  ggtitle('B')+
  theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(size=30))

#jpeg("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Figure4/Figure_4.jpeg", width = 10, height =5.5, units = 'in', res = 500)
new_Fi4A
#dev.off()




library(gridExtra)
jpeg("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/newstart/figures/Figure4.jpeg", width = 10, height =12, units = 'in', res = 300)
grid.arrange(ori_Fi4A,new_Fi4A, ncol=1, nrow=2)
dev.off()








#####jaccard splitting samples#########

M1_vec<-c(5,10,15,30,50,100,300,500,1000)
jaccard_mat_list<-list()
library(foreach)


library(sets)




jaccard_fun<-function(v1,v2){
    qq<-gset_similarity(as.set(v1),as.set(v2), method = "Jaccard")
    return(qq)
}

jaccard_I_BIRF_mat<-foreach(i = 1:8,.combine=cbind)%:%
    foreach(k=1:10,.combine=cbind)%:%
    foreach(j=1:9,.combine=c)%do%{
        kss<-jaccard_fun(list_all8[[i]][[k]]$acc_I_BI$TOP_FEA_list[[j]],list_all8[[i]][[k]]$acc_I_RF$TOP_FEA_list[[j]])
        return(kss)
        
        
    }


jaccard_P_BIRF_mat<-foreach(i = 1:8,.combine=cbind)%:%
    foreach(k=1:10,.combine=cbind)%:%
    foreach(j=1:9,.combine=c)%do%{
        kss<-jaccard_fun(list_all8[[i]][[k]]$acc_P_BI$TOP_FEA_list[[j]],list_all8[[i]][[k]]$acc_P_RF$TOP_FEA_list[[j]])
        return(kss)
    }  



M1_vec<-c(5,10,15,30,50,100,300,500,1000)
jaccard_mat_list<-list(jaccard_I_BIRF_mat=jaccard_I_BIRF_mat,jaccard_P_BIRF_mat=jaccard_P_BIRF_mat)

Comparisons<-c('I_BI*I_RF','P_BI*P_RF')
jaccard_BAR<-foreach(i=1:2,.combine=rbind)%do%{
    qq_mean<-rowMeans(jaccard_mat_list[[i]])
    qq_se<-apply(jaccard_mat_list[[i]],1,function(x){sd(x)/80})
    qq_out<-cbind(qq_mean,qq_se,M1_vec[seq(1,9)],rep(Comparisons[i],9))
}
jaccard_BAR<-as.data.frame(jaccard_BAR)
colnames(jaccard_BAR)<-c('mean','se','sizefeatures','Comparisons')
jaccard_BAR$mean<-as.numeric(as.character(jaccard_BAR$mean))
jaccard_BAR$se<-as.numeric(as.character(jaccard_BAR$se))

jaccard_BAR$sizefeatures<-factor(jaccard_BAR$sizefeatures,levels=unique(jaccard_BAR$sizefeatures))
jaccard_BAR$Comparisons<-factor(jaccard_BAR$Comparisons,levels=unique(jaccard_BAR$Comparisons))


#control#######
qq1<-list_all8[[8]][[10]]$acc_I_BI$TOP_FEA_list$FEA_FULL
qq2<-list_all8[[8]][[10]]$acc_I_RF$TOP_FEA_list$FEA_FULL

pp1<-list_all8[[8]][[10]]$acc_P_BI$TOP_FEA_list$FEA_FULL
pp2<-list_all8[[8]][[10]]$acc_P_RF$TOP_FEA_list$FEA_FULL

#qq2<-list_all8[[8]][[9]]$acc_I_BI$TOP_FEA_list$FEA_FULL

nuuu_vec<-c(5,10,15,30,50,100,300,500,1000)


jaccard_control_I_mat<-foreach(i=1:9,.combine=cbind)%:%
    foreach(j=1:1000,.combine=c)%do%{
        nuuu<-nuuu_vec[i]
        pp<-jaccard_fun(sample(qq1,nuuu),sample(qq2,nuuu))
        return(pp)
    }




jaccard_control_P_mat<-foreach(i=1:9,.combine=cbind)%:%
    foreach(j=1:1000,.combine=c)%do%{
        nuuu<-nuuu_vec[i]
        pp<-jaccard_fun(sample(pp1,nuuu),sample(pp2,nuuu))
        return(pp)
    }

jaccard_control_mean<-c(colMeans(jaccard_control_I_mat),colMeans(jaccard_control_P_mat))
jaccard_control_se<-c(apply(jaccard_control_I_mat,2,function(x){sd(x)/1000}),apply(jaccard_control_P_mat,2,function(x){sd(x)/1000}))

library(ggplot2)
textsize<-15
new_jaccard_plot<-ggplot(data=jaccard_BAR, aes(fill=jaccard_BAR$Comparisons, y=jaccard_BAR$mean, x=jaccard_BAR$sizefeatures)) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9,show.legend =F)+
    #scale_y_continuous(limits=c(0,1))+
    geom_errorbar(aes(ymin=jaccard_BAR$mean-jaccard_BAR$se, ymax=jaccard_BAR$mean+jaccard_BAR$se),size=1,width=.5,position=position_dodge(.9),color='green')+
    
    geom_errorbar(aes(ymin=jaccard_control_mean-jaccard_control_se, ymax=jaccard_control_mean+jaccard_control_se),size=1,width=.5,position=position_dodge(.9),color='blue')+
    
    #geom_text(aes(label=round(jaccard_BAR$mean,2)), vjust=0.1, color="black", position = position_dodge(0.9), size=4,angle=90,hjust=-0.5)+
    labs(fill = "Comparisons",y='Jaccard index',x='# of top ranked features used')+
    scale_fill_grey()+
    ggtitle('B')+
    theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


#####jaccard splitting datasets#########
load("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Figure4/topfea_mat_list_list.RData")



M1_vec<-c(5,10,15,30,50,100,300,500,1000)
ks_mat_list<-list()
library(foreach)

for(i32 in 1:32){
    
    ks_mat<-foreach(i = 1:9,.combine=rbind)%do%{
        
        ks_vec<-NULL
        for(j in 1:3){
            
            for(k in (j+1):4){
                qq<-jaccard_fun(topfea_mat_list_list[[i32]][[i]][,j],topfea_mat_list_list[[i32]][[i]][,k])
                ks_vec<-c(ks_vec,qq)
            }
        }
        names(ks_vec)<-c('I_BI*I_RF','I_BI*P_BI','I_BI*P_RF',
                         'I_RF*P_BI','I_RF*P_RF',
                         'P_BI*P_RF')
        return(ks_vec)
    }
    
    
    
    
    
    
    rownames(ks_mat)<-M1_vec
    ks_mat_list[[i32]]<-ks_mat
    
}
names(ks_mat_list)<-seq(1,32)

I_BI_RF<-foreach(i=1:32,.combine=cbind)%do%{
    qq<-ks_mat_list[[i]]
    return(qq[,1])
}
P_BI_RF<-foreach(i=1:32,.combine=cbind)%do%{
    qq<-ks_mat_list[[i]]
    return(qq[,6])
}

I_BI_RF_rowmean<-rowMeans(I_BI_RF)
P_BI_RF_rowmean<-rowMeans(P_BI_RF)


I_BI_RF_se<-apply(I_BI_RF,1,sd)/sqrt(32)
P_BI_RF_se<-apply(P_BI_RF,1,sd)/sqrt(32)


I_dat<-data.frame(mean=I_BI_RF_rowmean,se=I_BI_RF_se,Comparisons=rep('I_BI*I_RF',9),sizefeatures=as.factor(M1_vec))
P_dat<-data.frame(mean=P_BI_RF_rowmean,se=P_BI_RF_se,Comparisons=rep('P_BI*P_RF',9),sizefeatures=as.factor(M1_vec))
Bar_dat<-rbind(I_dat,P_dat)




#old control#######
load("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/ppqq_OLD.RData")

nuuu_vec<-c(5,10,15,30,50,100,300,500,1000)


old_control_I_mat<-foreach(i=1:9,.combine=cbind)%:%
    foreach(j=1:1000,.combine=c)%do%{
        nuuu<-nuuu_vec[i]
        pp<-jaccard_fun(sample(old_qq1,nuuu),sample(old_qq2,nuuu))
        return(pp)
    }
colMeans(old_control_I_mat)
apply(old_control_I_mat,2,function(x){sd(x)/1000})




old_control_P_mat<-foreach(i=1:9,.combine=cbind)%:%
    foreach(j=1:1000,.combine=c)%do%{
        nuuu<-nuuu_vec[i]
        pp<-jaccard_fun(sample(old_pp1,nuuu),sample(old_pp2,nuuu))
        return(pp)
    }
colMeans(old_control_P_mat)
apply(old_control_P_mat,2,function(x){sd(x)/1000})

old_control_mean<-c(colMeans(old_control_I_mat),colMeans(old_control_P_mat))
old_control_se<-c(apply(old_control_I_mat,2,function(x){sd(x)/1000}),apply(old_control_P_mat,2,function(x){sd(x)/1000}))

library(ggplot2)
textsize<-15
ori_KS_plot<-ggplot(data=Bar_dat, aes(fill=Bar_dat$Comparisons, y=Bar_dat$mean, x=Bar_dat$sizefeatures)) +
    geom_bar(stat="identity", position = position_dodge(0.9),width=0.9)+
    scale_y_continuous(limits=c(0,1))+
    geom_errorbar(aes(ymin=Bar_dat$mean-Bar_dat$se, ymax=Bar_dat$mean+Bar_dat$se),size=1,width=.5,position=position_dodge(.9),color='green')+
    geom_errorbar(aes(ymin=old_control_mean-old_control_se, ymax=old_control_mean+old_control_se),size=1,width=.5,position=position_dodge(.9),color='blue')+
    #geom_text(aes(label=round(Bar_dat$mean,2)), vjust=0.1, color="black", position = position_dodge(0.9), size=4,angle=90,hjust=-0.5)+
    labs(fill = "Comparisons",y='Jaccard index',x='# of top ranked features used')+
    scale_fill_grey()+
    ggtitle('A')+
    theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"),legend.position="top")

ori_KS_plot



g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)}

mylegend<-g_legend(ori_KS_plot)

library(gridExtra)
jpeg("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/newstart/figures/Figure5.jpeg", width = 10, height =5, units = 'in', res = 300)
#grid.arrange(ori_KS_plot+ theme(legend.position="none"),new_jaccard_plot, ncol=2, nrow=1)
p3 <- grid.arrange(arrangeGrob(ori_KS_plot + theme(legend.position="none"),
                               new_jaccard_plot + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))
dev.off()



