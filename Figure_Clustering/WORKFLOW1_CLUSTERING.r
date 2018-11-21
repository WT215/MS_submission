
source("D:/MainProject/MS_submission/Functions/FUNCTIONS.R")

#Complete code:E:/MainProject/MS_submission/Figure_Clustering/ANALYSIS_ROUTINE.R
load("D:/MainProject/MS_submission/Figure_Clustering/WORKFLOW1_CLUSTERING.RData")
HC_plot2<-function(DATA,Truelabels,Method='COR',leg=T,MAIN='',input_cex=3){
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
    
    qqq<-All_samples_techrep[order.dendrogram(hc)]
    
    
    
    
    #set('labels',NULL)%>% 
    hc %>% set("leaves_pch", qqq) %>% set("leaves_col", color_type) %>% set("leaves_cex", 1) %>% set('clear_branches')%>% set('labels_cex',12)%>% set("branches_lwd",1.5)%>% plot(horiz = F,main=MAIN,cex.main=input_cex,adj=0, leaflab = "none", axes=F)
    #title(main=MAIN,adj=0,cex=100)
    
    if(leg){  legend('topright',legend=c("MRSA samples","MSSA samples",paste("Left:",leftMSSA,"MSSA,",leftMRSA,"MRSA"),paste("Right:",rightMSSA,"MSSA,",rightMRSA,"MRSA")),col=c("1","3","0","0"),pch=c(1,19),cex=0.6)}
    
    # ggd1 <- as.ggdend(hc)
    # library(ggplot2)
    # ggplot(ggd1)
}


All_samples_techrep<-c(rep(15,40),rep(16,40),rep(17,20),rep(18,40))
names(All_samples_techrep)<-rownames(featureMatrix)

source("E:/MainProject/PROJECT_MRSA_MSSA/R_codes/FUNCTIONS.R")




jpeg("E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/newstart/figures/FIGURE_CLUSTERING_v3_shape.jpeg", width = 15, height = 10, units = 'in', res = 300)
#pdf(file="E:/MainProject/PROJECT_MRSA_MSSA/SUBMISSION/SUBMISSION_02122018/newstart/figures/FIGURE_CLUSTERING_v3_shape_LSA2.pdf", width = 15, height =10)

#pdf(file="E:/MainProject/PROJECT_MRSA_MSSA/FIGURES/FIGURE_CLUSTERING_v3_shape_LSA.pdf", width = 25, height = 20)
#pdf(file="E:/MainProject/PROJECT_MRSA_MSSA/FIGURES/FIGURE_CLUSTERING_v3_shape_LSA2.pdf", width = 9, height =7)

par(mfrow=c(2,1))
data<-Xall.b
#data<-featureMatrix[,Selectind]
DATA=featureMatrix


HC_plot2(DATA=featureMatrix,Truelabels=TRUE_LABELS_all,Method='COR',MAIN='A',leg=F,input_cex=2)

#HC_plot2(data,TRUE_LABELS_all,Method='BIN',leg=F,MAIN='b',input_cex=3)

dim(featureMatrix)
TOP100_mz<-as.numeric(rownames(br[1:100,]))
Selectind<-as.data.frame(br[1:100,])$idx
data<-Xall.b[,Selectind]
HC_plot2(data,TRUE_LABELS_all,Method='BIN',leg=F,MAIN='B',input_cex=2)

dev.off()
