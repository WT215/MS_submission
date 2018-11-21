
load("E:/MainProject/MS_submission/Figure_top100/listt8_par1pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_par1r2R2.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_par2pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_par2pa15062016R2.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_pcr1pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_pcr1pa15062016R2.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_pcr2pa15062016R1.RData")
load("E:/MainProject/MS_submission/Figure_top100/listt8_pcr2pa15062016R2.RData")
list8all<-list(listt8_par1pa15062016R1=listt8_par1pa15062016R1,
               listt8_par1r2R2=listt8_par1r2R2,
               listt8_par2pa15062016R1=listt8_par2pa15062016R1,
               listt8_par2pa15062016R2=listt8_par2pa15062016R2,
               listt8_pcr1pa15062016R1=listt8_pcr1pa15062016R1,
               listt8_pcr1pa15062016R2=listt8_pcr1pa15062016R2,
               listt8_pcr2pa15062016R1=listt8_pcr2pa15062016R1,
               listt8_pcr2pa15062016R2=listt8_pcr2pa15062016R2)
library(foreach)
Top100_I_BI_list<-foreach(i=1:8,.combine=rbind)%do%{
    qq<-list8all[[i]]$acc_I_BI$TOP_FEA_list$`100`
    return(qq)
}

Top100_I_RF_list<-foreach(i=1:8,.combine=rbind)%do%{
    qq<-list8all[[i]]$acc_I_RF$TOP_FEA_list$`100`
    return(qq)
}

Top100_P_BI_list<-foreach(i=1:8,.combine=rbind)%do%{
    qq<-list8all[[i]]$acc_P_BI$TOP_FEA_list$`100`
    return(qq)
}
Top100_P_RF_list<-foreach(i=1:8,.combine=rbind)%do%{
    qq<-list8all[[i]]$acc_P_RF$TOP_FEA_list$`100`
    return(qq)
}



Top100_4lists<-list(Top100_I_BI_list=Top100_I_BI_list,Top100_I_RF_list=Top100_I_RF_list,Top100_P_BI_list=Top100_P_BI_list,Top100_P_RF_list=Top100_P_RF_list)
###begin analysis#####

Workflow<-c('I_BI','I_RF','P_BI','P_RF')
TOP100<-do.call(rbind,Top100_4lists)
Groupp<-c(rep(Workflow[1],32),rep(Workflow[2],32),rep(Workflow[3],32),rep(Workflow[4],32))
Top100_2<-foreach(i=1:dim(TOP100)[1],.combine=rbind)%do%{
    jj<-cbind(TOP100[i,],rep(paste(i,'sdf'),100),Groupp[i])
    return(jj)
}

Top100_2<-as.data.frame(Top100_2)
colnames(Top100_2)<-c('PeakPos','Run','Group')
Top100_2$PeakPos<-as.numeric(as.character(Top100_2$PeakPos))
class(Top100_2$Run)

levels(Top100_2$Group)



library(ggplot2)
library(reshape2)
#ggplot(Top100_2, aes(x=Top100_2$PeakPos)) + geom_density(aes(group=Top100_2$Run,fill=Top100_2$Group),bw=50,alpha=0.7)+labs(x = "m/z values",y='density',fill='Methods')+theme_minimal()+ ggtitle(paste('Density plot of top 100 ranked features'))+theme(plot.title = element_text(hjust = 0.5))



Top100_2_lists<-foreach(i=1:length(Top100_4lists))%:%
    foreach(j=1:dim(Top100_4lists[[i]])[1],.combine=rbind)%do%{
        kk<-cbind(Top100_4lists[[i]][j,],rep(paste(j),100),rep(Workflow[i],100))

        return(kk)
    }

for(i in 1:4){
    Top100_2_lists[[i]]<-as.data.frame(Top100_2_lists[[i]])
    colnames(Top100_2_lists[[i]])<-c('PeakPos','Run','Method')
    Top100_2_lists[[i]]$PeakPos<-as.numeric(as.character(Top100_2_lists[[i]]$PeakPos))
}
names(Top100_2_lists)<-Workflow



library(gridExtra)
#names(Top100_2_lists[[1]])
textsize<-15
#jpeg("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/TOP100_first7astrain/TOP100.jpeg", width = 25, height = 15, units = 'in', res = 500)
panell<-c('A','B','C','D')
myGrobs <-lapply(seq_along(Top100_2_lists),function(ddd){
    bw=50
    ggplot(Top100_2_lists[[ddd]], aes(x=Top100_2_lists[[ddd]]$PeakPos)) + geom_density(aes(group=Top100_2_lists[[ddd]]$Run),bw=bw)+
        labs(x = "m/z values",y='density',fill='Methods')+
        ggtitle(panell[ddd])+
        theme(axis.text=element_text(size=textsize),axis.title=element_text(size=textsize),legend.text=element_text(size=textsize),legend.title=element_text(size=textsize),plot.title = element_text(size = textsize))+ 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
})
grid.arrange(grobs=myGrobs, ncol=2, nrow=2)

#dev.off()



# rownames(Top100_4lists$Top100_I_BI_list)<-paste('The ',seq(1,8),'th supervised testing',sep='')
# write.csv(t(Top100_4lists$Top100_I_BI_list),file="E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/TOP100_first7astrain/TOP100_I_BI.csv")
# 
# 
# 
# rownames(Top100_4lists$Top100_I_RF_list)<-paste('The ',seq(1,8),'th supervised testing',sep='')
# write.csv(t(Top100_4lists$Top100_I_RF_list),file="E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/TOP100_first7astrain/TOP100_I_RF.csv")
# 
# 
# rownames(Top100_4lists$Top100_P_BI_list)<-paste('The ',seq(1,8),'th supervised testing',sep='')
# write.csv(t(Top100_4lists$Top100_P_BI_list),file="E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/TOP100_first7astrain/TOP100_P_BI.csv")
# 
# 
# 
# rownames(Top100_4lists$Top100_P_RF_list)<-paste('The ',seq(1,8),'th supervised testing',sep='')
# write.csv(t(Top100_4lists$Top100_P_RF_list),file="E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/BMC_rebuttle/Fig_split4BIosamples/TOP100_first7astrain/TOP100_P_RF.csv")




