####Script to calculate measures of fit for the phenology model.
##Author:Ian Breckheimer
##Data: E. Theobald
##Last Updated: 1 October 2017

library(ggplot2)
library(reshape2)
library(dplyr)
library(ggmcmc)
library(ROCR)
library(gridExtra)
source("./code/pheno_reassembly_functions.R")

#loads JAGS output
load("./output/th_jagsoutput_allspp_sppplotyr_moist_gdd_perf.Rdata")
obs_all <- Sp_jags_perf$mod$data()$y

ggd <- ggs(Sp_jags_perf$out)
p_log <- grepl("^p",ggd$Parameter)
pred_p <- ggd[p_log,]
p_names <- levels(factor(pred_p$Parameter))
niter <- max(pred_p$Iteration)*max(pred_p$Chain)

p_indices <- as.numeric(regmatches(p_names,gregexpr("[0-9]{1,}",p_names)))
obs <- obs_all[p_indices]

p_preds <- array(dim=c(niter,length(p_names)))
p_obs <- array(dim=c(niter,length(p_names)))

for(i in 1:length(p_names)){
  p_preds[,i] <- pred_p$value[pred_p$Parameter==p_names[i]]
}

#calculates auc
auc <- rep(NA,niter)
tpr <- array(dim=c(niter,length(p_names)+1))
fpr <- array(dim=c(niter,length(p_names)+1))

for (i in 1:niter){
  preds <- p_preds[i,]
  pred <- prediction(preds,factor(obs,levels=c(0,1)))
  perf <- performance(pred,"auc")
  auc[i] <- perf@y.values[[1]]
  perf2 <- performance(pred,measure="tpr",x.measure="fpr")
  fpr[i,1:length(perf2@x.values[[1]])] <- perf2@x.values[[1]]
  tpr[i,1:length(perf2@y.values[[1]])] <- perf2@y.values[[1]]
}

fprm <- melt(fpr,varnames=c("iteration","obs"))
tprm <- melt(tpr,varnames=c("iteration","obs"))

perfdat <- data.frame(Type="Training",
                      iter=fprm$iteration,
                      TPR=tprm$value,
                      FPR=fprm$value)
aucdat <- data.frame(Type="Training",
                     AUC=auc,
                     iter=1:length(auc))

p1 <- ggplot(aucdat)+
  geom_density(aes(x=AUC),fill="grey20")+
  xlim(c(0.5,1))+
  ggtitle("AUC Density")+
  theme_bw()+
  theme(panel.grid=element_blank())


p2 <- ggplot(perfdat)+
        geom_line(aes(x=FPR,y=TPR,group=iter),alpha=0.01)+
        geom_abline(aes(intercept=0,slope=1),linetype="dotted")+
        ggtitle("ROC")+
        theme_bw()+
        theme(panel.grid=element_blank())

#pdf("./figs/roc_auc_1000obs.pdf",width=8,height=4)
grid.arrange(p1,p2,ncol=2)
#dev.off()

####Calculates AUC and ROC curves for testing data without random effects####
tdat <- read.csv("./scratch/final_test_data.csv")

load("./output/th_jagsoutput_allspp_sppplotyr_moist_gdd.Rdata")
out <- Sp_jags_up$out

param.names <- c("height.s",
                 "height_slope.s.sdd",
                 "opt.s",
                 "opt_slope.s.sdd",
                 "opt_slope.s.moist",
                 "opt_slope.s.gdd",
                 "width.s",
                 "width_slope.s.sdd",
                 "width_slope.s.moist",
                 "width_slope.s.gdd")
n.samples=500
spp.indices=1:48

param_array <- create_param_array(out=out,n.sample=n.samples,spp.indices=spp.indices,
                                  param.names=param.names)

test_preds <- matrix(NA,ncol=dim(param_array)[3],nrow=nrow(tdat))

##Predicted flowering probabilities for held out test data for each posterior sample.
for (i in 1:n.samples){
  test_preds[,i] <- jags_fit_fun_moist_gdd(height_int = param_array[tdat$Species_num,1,i],
                                           height_slope_sdd = param_array[tdat$Species_num,2,i],
                                           opt_int = param_array[tdat$Species_num,3,i],
                                           opt_slope_gdd = param_array[tdat$Species_num,4,i],
                                           opt_slope_moist = param_array[tdat$Species_num,5,i],
                                           opt_slope_sdd = param_array[tdat$Species_num,6,i],
                                           width_int = param_array[tdat$Species_num,7,i],
                                           width_slope_gdd = param_array[tdat$Species_num,8,i],
                                           width_slope_moist = param_array[tdat$Species_num,9,i],
                                           width_slope_sdd = param_array[tdat$Species_num,10,i],
                                           x1vec=(tdat$DOY/100 - 1.7),
                                           x2vec=(tdat$SDDDOY/100 - 1.7),
                                           x3vec=tdat$moist_scaled,
                                           x4vec=tdat$gdd_scaled)
}

##Computes ROC and AUC for test data.
niter <- 500
nspp <- 48
nobs <- nrow(tdat)

test_obs <- as.numeric(tdat$PercFlwr > 0)

test_auc <- rep(NA,niter)
test_tpr <- array(dim=c(niter,nobs+1))
test_fpr <- array(dim=c(niter,nobs+1))

for (i in 1:niter){
  test_pred_iter <- test_preds[,i]
  test_pred <- prediction(test_pred_iter,factor(test_obs,levels=c(0,1)))
  test_perf <- performance(test_pred,"auc")
  test_auc[i] <- test_perf@y.values[[1]]
  test_perf2 <- performance(test_pred,measure="tpr",x.measure="fpr")
  test_fpr[i,1:length(test_perf2@x.values[[1]])] <- test_perf2@x.values[[1]]
  test_tpr[i,1:length(test_perf2@y.values[[1]])] <- test_perf2@y.values[[1]]
}

test_fprm <- melt(test_fpr,varnames=c("iteration","obs"))
test_tprm <- melt(test_tpr,varnames=c("iteration","obs"))

test_perfdat <- data.frame(Type="Holdout",
                           iter=test_fprm$iteration,
                           TPR=test_tprm$value,
                           FPR=test_fprm$value)
test_aucdat <- data.frame(Type="Holdout",
                          AUC=test_auc,
                          iter=1:length(test_auc))

####Measures per-species performance####
spp_auc <- array(dim=c(nspp,niter))
spp_tpr <- array(dim=c(nspp,niter,nobs+1))
spp_fpr <- array(dim=c(nspp,niter,nobs+1))

for (j in 1:nspp){
  print(paste("Assessing performance for Species ",j))
  if(nrow(tdat[tdat$Species_num==j,]) > 16){
    for (i in 1:niter){
      spp_pred_iter <- test_preds[tdat$Species_num==j,i]
        spp_pred <- prediction(spp_pred_iter,factor(test_obs[tdat$Species_num==j],levels=c(0,1)))
        spp_perf <- performance(spp_pred,"auc")
        spp_auc[j,i] <- spp_perf@y.values[[1]]
        spp_perf2 <- performance(spp_pred,measure="tpr",x.measure="fpr")
        spp_fpr[j,i,1:length(spp_perf2@x.values[[1]])] <- spp_perf2@x.values[[1]]
        spp_tpr[j,i,1:length(spp_perf2@y.values[[1]])] <- spp_perf2@y.values[[1]]
      }
    }else{
      print(paste("Insufficient observations for species",j))
    }
}

spp_fprm <- melt(spp_fpr,varnames=c("spp","iteration","obs"))
spp_tprm <- melt(spp_tpr,varnames=c("spp","iteration","obs"))

spp_perfdat <- data.frame(Species=spp_fprm$spp,
                           iter=spp_fprm$iteration,
                           TPR=spp_tprm$value,
                           FPR=spp_fprm$value)
spp_perf_grp <- group_by(spp_perfdat,Species,FPR)
spp_perf_sum <- summarise(spp_perf_grp,TPR_mean=mean(TPR,na.rm=TRUE),
                                       TPR_lwr=quantile(TPR,probs=c(0.1),na.rm=TRUE),
                                       TPR_upr=quantile(TPR,probs=c(0.9),na.rm=TRUE))

auc_spp_mean <- apply(spp_auc,FUN=mean,MARGIN=1,na.rm=TRUE)
auc_spp_lwr <- apply(spp_auc,FUN=function(x){quantile(x,probs=0.1,na.rm=TRUE)},MARGIN=1)
auc_spp_upr <- apply(spp_auc,FUN=function(x){quantile(x,probs=0.9,na.rm=TRUE)},MARGIN=1)

auc_spp_all <- data.frame(Species_num=1:48,
                          AUC_mean=auc_spp_mean,
                          AUC_lwr10=auc_spp_lwr,
                          AUC_upr90=auc_spp_upr)
spp_names <- data.frame(Species_num=as.numeric(levels(factor(tdat$Species_num))),
                        Species=levels(tdat$Species))
auc_spp_names <- left_join(auc_spp_all,spp_names,by="Species_num")

####Makes Plots and Summaries####

##Appends performance data
all_perfdat <- rbind(perfdat,test_perfdat)
all_aucdat <- rbind(aucdat,test_aucdat)

##Summarizes performance data for easier plotting.
auc_grp <- group_by(all_aucdat,Type)
auc_sum <- summarise(auc_grp,AUC_mean=mean(AUC,na.rm=TRUE),
                             AUC_lwr=quantile(AUC,probs=0.025,na.rm=TRUE),
                             AUC_upr=quantile(AUC,probs=0.975,na.rm=TRUE))


perf_grp <- group_by(all_perfdat,Type,FPR)
perf_sum <- summarise(perf_grp,TPR_mean=mean(TPR,na.rm=TRUE),
                               TPR_lwr=quantile(TPR,probs = 0.025,na.rm=TRUE),
                               TPR_upr=quantile(TPR,probs = 0.975,na.rm=TRUE))

##Plots training and testing performance together.

p1 <- ggplot()+
  geom_density(aes(x=AUC,fill=Type),color=rgb(0,0,0,0),data=all_aucdat)+
  xlim(c(0.5,1))+
  geom_text(aes(x=(auc_sum$AUC_mean-0.15),y=c(400,200),
                label=paste("Mean:",round(auc_sum$AUC_mean,3)),
                color=Type),data=auc_sum)+
  scale_color_grey(start=0.2,end=0.6)+
  scale_fill_grey(start=0.2,end=0.6)+
  ggtitle("AUC Density")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="none")


p2 <- ggplot(perf_sum)+
  geom_ribbon(aes(x=FPR,ymin=TPR_lwr,ymax=TPR_upr,fill=Type),
              color=rgb(0,0,0,0),alpha=0.4)+
  geom_line(aes(x=FPR,y=TPR_mean,color=Type),linetype="solid",lwd=0.5)+
  geom_abline(aes(intercept=0,slope=1),linetype="dotted")+
  scale_color_grey(start=0.2,end=0.6)+
  scale_fill_grey(start=0.2,end=0.6)+
  ggtitle("ROC")+
  ylab("TPR")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position=c(0.7,0.2))

#pdf("./figs/roc_auc_test_10pct.pdf",width=8,height=4)
grid.arrange(p1,p2,ncol=2)
#dev.off()

## Makes per-species plots of performance.

auc_spp_names$Species_order <- factor(auc_spp_names$Species,
                                      levels=auc_spp_names$Species[order(auc_spp_names$AUC_mean)])

#pdf("./figs/auc_species_cred80.pdf",width=6,height=6)
ggplot(auc_spp_names[complete.cases(auc_spp_names),])+
  geom_point(aes(y=Species_order,x=AUC_mean,color=Species))+
  geom_segment(aes(y=Species_order,yend=Species_order,x=AUC_lwr10,xend=AUC_upr90,
                   color=Species))+
  xlim(c(0.5,1))+
  xlab("AUC")+
  ylab("")+
  theme_bw()+
  theme(legend.position="none")
#dev.off()

