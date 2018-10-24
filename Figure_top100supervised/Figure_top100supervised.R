function_matlist<-function(list_all8){
  
  MAT_I_BI<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_I_BI$acc_vec
      return(qq)
    }
  rownames(MAT_I_BI)<-seq(2,14,2)*4
  
  MAT_I_RF<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_I_RF$acc_vec
      return(qq)
    }
  rownames(MAT_I_RF)<-seq(2,14,2)*4
  
  MAT_P_BI<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_P_BI$acc_vec
      return(qq)
    }
  rownames(MAT_P_BI)<-seq(2,14,2)*4
  
  MAT_P_RF<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_P_RF$acc_vec
      return(qq)
    }
  rownames(MAT_P_RF)<-seq(2,14,2)*4
  
  MAT_MALDI_BI<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_MALDI_BI
      return(qq)
    }
  rownames(MAT_MALDI_BI)<-seq(2,14,2)*4
  
  MAT_MALDI_RF<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_MALDI_RF
      return(qq)
    }
  rownames(MAT_MALDI_RF)<-seq(2,14,2)*4
  
  MAT_LIST<-list(MAT_I_BI=MAT_I_BI,
                    MAT_I_RF=MAT_I_RF,
                    MAT_P_BI=MAT_P_BI,
                    MAT_P_RF=MAT_P_RF,
                    MAT_MALDI_BI=MAT_MALDI_BI,
                    MAT_MALDI_RF=MAT_MALDI_RF)
  
  return(MAT_LIST)
}


function_noalignmatlist<-function(list_all8){
  
  MAT_I_BI<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_I_BI
      return(qq)
    }
  rownames(MAT_I_BI)<-seq(2,14,2)*4
  
  MAT_I_RF<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_I_RF
      return(qq)
    }
  rownames(MAT_I_RF)<-seq(2,14,2)*4
  
  MAT_P_BI<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_P_BI
      return(qq)
    }
  rownames(MAT_P_BI)<-seq(2,14,2)*4
  
  MAT_P_RF<-foreach(i=1:length(list_all8),.combine=cbind)%:%
    foreach(j=1:7,.combine=c)%do%{
      qq<-list_all8[[i]][[j]]$acc_P_RF
      return(qq)
    }
  rownames(MAT_P_RF)<-seq(2,14,2)*4

  
  MAT_LIST<-list(MAT_I_BI=MAT_I_BI,
                 MAT_I_RF=MAT_I_RF,
                 MAT_P_BI=MAT_P_BI,
                 MAT_P_RF=MAT_P_RF)
  
  return(MAT_LIST)
}


#aligned v1###########
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_par1pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_pcr1pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_pcr1pa15062016R2.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_pcr2pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_pcr2pa15062016R2.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_par1r2R2.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_par2pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100supervised/aligned/Fig6samples_par2pa15062016R2.RData")


list_all8<-list(
  
                Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
                Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
                Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
                Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
                Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
                Fig6samples_par1r2R2=Fig6samples_par1r2R2,
                Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
                Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
                )



Figure6_traintest_list$FourBio_pcr1pa15062016R1$Train_list_7list[[2]]

Fig6samples_pcr1pa15062016R1[[6]]$acc_P_BI$acc_vec


library(foreach)
MAT_LIST<-function_matlist(list_all8)
lapply(MAT_LIST,function(x){rowMeans(x)})




#####aligned v2#########
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/aligned_v2/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/aligned_v2/',x,sep=''),.GlobalEnv)
})
list_all8_v2<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_LIST_v2<-function_matlist(list_all8_v2)

#####aligned v3#########
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/aligned_v3/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/aligned_v3/',x,sep=''),.GlobalEnv)
})
list_all8_v3<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_LIST_v3<-function_matlist(list_all8_v3)
lapply(MAT_LIST_v3,function(x){rowMeans(x)})




#####aligned v4#########
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/aligned_v4/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/aligned_v4/',x,sep=''),.GlobalEnv)
})
list_all8_v4<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_LIST_v4<-function_matlist(list_all8_v4)
lapply(MAT_LIST_v4,function(x){rowMeans(x)})


#####aligned v5#########
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/aligned_v5/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/aligned_v5/',x,sep=''),.GlobalEnv)
})
list_all8_v5<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_LIST_v5<-function_matlist(list_all8_v5)
lapply(MAT_LIST_v5,function(x){rowMeans(x)})



#plot aligned########

mat_list_list<-list(
  MAT_LIST=MAT_LIST,
                    MAT_LIST_v2=MAT_LIST_v2,
                    MAT_LIST_v3=MAT_LIST_v3,
                    MAT_LIST_v4=MAT_LIST_v4,
                    MAT_LIST_v5=MAT_LIST_v5)



MAT_LIST_ALL<-list()
for(i in 1:6){
  MAT_LIST_ALL[[i]]<-cbind(mat_list_list[[1]][[i]],mat_list_list[[2]][[i]],mat_list_list[[3]][[i]],mat_list_list[[4]][[i]]
                           ,mat_list_list[[5]][[i]]
                           )
}

Workflow<-c('I_BI','I_RF','P_BI','P_RF','MALDI_BI','MALDI_RF')
names(MAT_LIST_ALL)<-Workflow

BAR_PLOT<-foreach(i=1:length(MAT_LIST_ALL),.combine=rbind)%do%{
  seee<-apply(MAT_LIST_ALL[[i]],1,function(x){sd(x)/40})
  qq<-cbind(rowMeans(MAT_LIST_ALL[[i]]),seee,c(8,16,24,32,40,48,56),rep(Workflow[i],7))
  return(qq)
}






BAR_PLOT<-as.data.frame(BAR_PLOT)
colnames(BAR_PLOT)<-c('AccuracyRate','SE','NumberofSam','Workflow')


BAR_PLOT$AccuracyRate<-as.numeric(as.character(BAR_PLOT$AccuracyRate))
BAR_PLOT$AccuracyRate<-round(BAR_PLOT$AccuracyRate,3)

BAR_PLOT$SE<-as.numeric(as.character(BAR_PLOT$SE))
BAR_PLOT$NumberofSam<-factor(BAR_PLOT$NumberofSam,levels=unique(BAR_PLOT$NumberofSam))
BAR_PLOT$Workflow<-factor(BAR_PLOT$Workflow,levels=unique(BAR_PLOT$Workflow))





textsize<-15

library(ggplot2)
#Averaged across 32 testing results

#####Fig6A########
new_Fi6A<-ggplot(data=BAR_PLOT, aes(x=BAR_PLOT$NumberofSam, y=BAR_PLOT$AccuracyRate, fill=BAR_PLOT$Workflow)) +
  geom_bar(stat="identity", position = position_dodge(0.9),width=0.9,show.legend =T)+
  scale_y_continuous(limits=c(0,1))+
  geom_errorbar(aes(ymin=BAR_PLOT$AccuracyRate-BAR_PLOT$SE, ymax=BAR_PLOT$AccuracyRate+BAR_PLOT$SE),size=0.5,width=.5,position=position_dodge(.9),color='green')+
  geom_text(aes(label=round(BAR_PLOT$AccuracyRate,2)), vjust=0.1, color="black", position = position_dodge(0.9), size=4,angle=90,hjust=-0.1)+
  labs(x = "# of samples used for training (80 samples in total)",y='Averaged accuracy rates',fill='Methods')+
  scale_fill_grey()+
  ggtitle('c')+
  theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(size=30))

new_Fi6A




###no align###########

#no align v1
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/noaligned/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/noaligned/',x,sep=''),.GlobalEnv)
})

list_noalign<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_noalignLIST<-function_noalignmatlist(list_noalign)
lapply(MAT_noalignLIST,function(x){rowMeans(x)})



#no align v2
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v2/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v2/',x,sep=''),.GlobalEnv)
})

list_noalign_v2<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_noalignLIST_v2<-function_noalignmatlist(list_noalign_v2)
lapply(MAT_noalignLIST_v2,function(x){rowMeans(x)})



#no align v3
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v3/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v3/',x,sep=''),.GlobalEnv)
})

list_noalign_v3<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_noalignLIST_v3<-function_noalignmatlist(list_noalign_v3)


#no align v4
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v4/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v4/',x,sep=''),.GlobalEnv)
})

list_noalign_v4<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)

MAT_noalignLIST_v4<-function_noalignmatlist(list_noalign_v4)


#no align v5
file_names=as.list(dir(path = 'E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v5/', pattern="Fig6samples_*"))
lapply(file_names,function(x){
  load(paste('E:/MainProject/MS_submission/Figure_top100supervised/noaligned_v5/',x,sep=''),.GlobalEnv)
})

list_noalign_v5<-list(
  
  Fig6samples_par1pa15062016R1=Fig6samples_par1pa15062016R1,
  Fig6samples_pcr1pa15062016R1=Fig6samples_pcr1pa15062016R1,
  Fig6samples_pcr1pa15062016R2=Fig6samples_pcr1pa15062016R2,
  Fig6samples_pcr2pa15062016R1=Fig6samples_pcr2pa15062016R1,
  Fig6samples_pcr2pa15062016R2=Fig6samples_pcr2pa15062016R2,
  Fig6samples_par1r2R2=Fig6samples_par1r2R2,
  Fig6samples_par2pa15062016R1=Fig6samples_par2pa15062016R1,
  Fig6samples_par2pa15062016R2=Fig6samples_par2pa15062016R2
)
MAT_noalignLIST_v5<-function_noalignmatlist(list_noalign_v5)















#plot noaligned########

mat_noalignlist_list<-list(
    MAT_noalignLIST=MAT_noalignLIST,
                    MAT_noalignLIST_v2=MAT_noalignLIST_v2,
                    MAT_noalignLIST_v3=MAT_noalignLIST_v3,
                    MAT_noalignLIST_v4=MAT_noalignLIST_v4,
                    MAT_noalignLIST_v5=MAT_noalignLIST_v5)

MAT_noalignLIST_ALL<-list()
for(i in 1:4){
  MAT_noalignLIST_ALL[[i]]<-cbind(mat_noalignlist_list[[1]][[i]],mat_noalignlist_list[[2]][[i]],mat_noalignlist_list[[3]][[i]],mat_noalignlist_list[[4]][[i]]
                                  ,mat_noalignlist_list[[5]][[i]]
                                  )
}


length(MAT_noalignLIST_ALL)


Workflow<-c('I_BI','I_RF','P_BI','P_RF')
names(MAT_noalignLIST_ALL)<-Workflow

BAR_PLOT2<-foreach(i=1:length(MAT_noalignLIST_ALL),.combine=rbind)%do%{
  seee<-apply(MAT_noalignLIST_ALL[[i]],1,function(x){sd(x)/40})
  qq<-cbind(rowMeans(MAT_noalignLIST_ALL[[i]]),seee,c(8,16,24,32,40,48,56),rep(Workflow[i],7))
  return(qq)
}





BAR_PLOT2<-as.data.frame(BAR_PLOT2)
colnames(BAR_PLOT2)<-c('AccuracyRate','SE','NumberofSam','Workflow')


BAR_PLOT2$AccuracyRate<-as.numeric(as.character(BAR_PLOT2$AccuracyRate))
BAR_PLOT2$AccuracyRate<-round(BAR_PLOT2$AccuracyRate,3)

BAR_PLOT2$SE<-as.numeric(as.character(BAR_PLOT2$SE))
BAR_PLOT2$NumberofSam<-factor(BAR_PLOT2$NumberofSam,levels=unique(BAR_PLOT2$NumberofSam))
BAR_PLOT2$Workflow<-factor(BAR_PLOT2$Workflow,levels=unique(BAR_PLOT2$Workflow))





textsize<-15

library(ggplot2)
#Averaged across 32 testing results

new_Fi6B<-ggplot(data=BAR_PLOT2, aes(x=BAR_PLOT2$NumberofSam, y=BAR_PLOT2$AccuracyRate, fill=BAR_PLOT2$Workflow)) +
  geom_bar(stat="identity", position = position_dodge(0.9),width=0.9,show.legend =T)+
  scale_y_continuous(limits=c(0,1))+
  geom_errorbar(aes(ymin=BAR_PLOT2$AccuracyRate-BAR_PLOT2$SE, ymax=BAR_PLOT2$AccuracyRate+BAR_PLOT2$SE),size=0.5,width=.5,position=position_dodge(.9),color='green')+
  geom_text(aes(label=round(BAR_PLOT2$AccuracyRate,2)), vjust=0.1, color="black", position = position_dodge(0.9), size=4,angle=90,hjust=-0.1)+
  labs(x = "# of samples used for training (80 samples in total)",y='Averaged accuracy rates',fill='Methods')+
  scale_fill_grey()+
  ggtitle('d')+
  theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(size=30))

new_Fi6B



library(ggplot2)
library(gridExtra)
#jpeg("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Figure6_samples/Fig6AB.jpeg", width = 8.5, height =11,units = 'in',res = 500)
grid.arrange(new_Fi6A,new_Fi6B, ncol=1, nrow=2)
#dev.off()

