####Script to visualize output of fit phenology model.
####Author: Ian Breckheimer
####Data: Elli Theobald Dissertation
####Date: 23 January 2015

library(ggmcmc)
library(ggplot2)
library(reshape2)
library(dplyr)
library(gridExtra)
library(stringr)
source("./code/pheno_reassembly_functions.R")

#loads JAGS output
load("./output/th_jagsoutput_allspp_sppplotyr_moist_gdd.Rdata")
odata <- read.csv("./scratch/final_train_data.csv")

##Visually checks for convergence.
out <- Sp_jags_up$out
#plot(out)
#summary(out)

##Reads in species traits and merges with prevalence data.
spp <- read.csv("./data/species_number.csv",header=FALSE)
colnames(spp) <- c("Num","Species")
spp$Species <- as.factor(spp$Species)
spp_sum <- read.csv("./data/species_traits_2010_2014.csv")
spp_sum_2015 <- read.csv("./data/species_traits_2015.csv")

##Estimates mean changes in climate drivers that each species experienced in 2015.
spp_env_2015 <- data.frame(Species=spp_sum_2015$Species,
                           Species_num=spp_sum_2015$Species_num,
                           mean_SDD_2015=spp_sum_2015$mean_SDD,
                           mean_gdd_2015=spp_sum_2015$mean_gdd,
                           mean_fdd_2015=spp_sum_2015$mean_fdd,
                           mean_moist_2015=spp_sum_2015$mean_moist)

spp_env <- data.frame(Species=spp_sum$Species,
                      Species_num=spp_sum$Species_num,
                      mean_SDD=spp_sum$mean_SDD,
                      mean_gdd=spp_sum$mean_gdd,
                      mean_fdd=spp_sum$mean_fdd,
                      mean_moist=spp_sum$mean_moist,
                      mean_elev=spp_sum$mean_elev,
                      prevalence=spp_sum$prevalence)
spp_sum_change <- merge(spp_env,spp_env_2015,by=c("Species","Species_num"),all.x=TRUE)

spp_sum_change$SDD_change <- spp_sum_change$mean_SDD_2015 - spp_sum_change$mean_SDD
spp_sum_change$GDD_change <- spp_sum_change$mean_gdd_2015 - spp_sum_change$mean_gdd
spp_sum_change$moist_change <- spp_sum_change$mean_moist_2015 - spp_sum_change$mean_moist

#spp_sum_change <- spp_sum_change[spp_sum_change$prevalence>0.05,]


r1 <- ggplot(spp_sum_change)+
  geom_point(aes(x=mean_elev,y=SDD_change,size=prevalence))+
  geom_smooth(aes(x=mean_elev,y=SDD_change,size=prevalence),
              color="black",method="lm",se=TRUE)+
  ylab("2015 Change in SDD")+
  xlab("Mean Elevation (m)")+
  scale_size_continuous(guide=FALSE)+
  theme_bw()+
  theme(panel.grid=element_blank())
r2 <- ggplot(spp_sum_change)+
  geom_point(aes(x=mean_gdd,y=GDD_change,size=prevalence))+
  geom_smooth(aes(x=mean_gdd,y=GDD_change,size=prevalence),
              color="black",method="lm",se=TRUE)+
  xlim(c(575,670))+
  ylab("2015 Change in post-snow GDD")+
  xlab("Typical (2010-2014) post-snow GDD")+
  scale_size_continuous(guide=FALSE)+
  theme_bw()+
  theme(panel.grid=element_blank())
r3 <- ggplot(spp_sum_change)+
  geom_point(aes(x=mean_moist,y=moist_change,size=prevalence))+
  geom_smooth(aes(x=mean_moist,y=moist_change,size=prevalence),
              color="black",method="lm",se=TRUE)+
  xlab("Typical (2010 - 2014) Soil Moisture Duration (Days)")+
  ylab("2015 Change in Soil Moisture Duration (Days)")+
  scale_size_continuous("Prevalence")+
  theme_bw()+
  theme(panel.grid=element_blank())

#pdf("./figs/spp_exposure_SDD_GDD_moist.pdf",width=12,height=4)
grid.arrange(r1,r2,r3,ncol=3,widths=c(1,1,1.3))
#dev.off()

##Classifies species by their typical environments using kmeans.
spp_kdat <- data.frame(mean_sdd=spp_sum$mean_SDD,
                       mean_elev=spp_sum$mean_elev,
                       mean_gdd=spp_sum$mean_gdd,
                       mean_moist=log(spp_sum$mean_moist))
spp_k <- kmeans(spp_kdat,centers=3,nstart=10)
spp_sum$Cluster <- spp_k$cluster

##Classifies plots by their typical environmental conditions using kmeans.
odata_n2015 <- filter(odata,Year!=2015)
plot_data <- unique(data.frame(Site=odata_n2015$Site,
                                Year=odata_n2015$Year,
                                Year_fact=as.factor(odata_n2015$Year),
                                ElevFact=as.factor(odata_n2015$Elevation),
                                Elevation=odata_n2015$Elevation,
                                SDDDOY=odata_n2015$SDDDOY,
                                gdd50=odata_n2015$gdd50,
                                moist_days=odata_n2015$moist_days,
                                Topo=odata_n2015$Topo))
plot_data_gr <- group_by(plot_data,ElevFact,Topo)
plot_data_sum <- summarise(plot_data_gr,mean_sdd=mean(SDDDOY),
                                        mean_gdd=mean(gdd50),
                                        mean_moist=mean(moist_days))

odata_2015 <- filter(odata,Year==2015)
plot_data_2015 <- unique(data.frame(Site=odata_2015$Site,
                               Year=odata_2015$Year,
                               Year_fact=as.factor(odata_2015$Year),
                               ElevFact=as.factor(odata_2015$Elev),
                               Elevation=odata_2015$Elevation,
                               SDDDOY=odata_2015$SDDDOY,
                               gdd50=odata_2015$gdd50,
                               moist_days=odata_2015$moist_days,
                               Topo=odata_2015$Topo))
plot_data_2015_gr <- group_by(plot_data_2015,ElevFact,Topo)
plot_data_sum_2015 <- summarise(plot_data_2015_gr,mean_sdd_2015=mean(SDDDOY),
                           mean_gdd_2015=mean(gdd50),
                           mean_moist_2015=mean(moist_days))
plot_data_diff <- merge(plot_data_sum,plot_data_sum_2015)
plot_data_diff$ElevTopo <- factor(paste(plot_data_diff$ElevFact,plot_data_diff$Topo,sep="-"))

##Calculates prevalence for each species for each topo band.
odata$ElevFact <- as.factor(odata$Elevation)
odata$ElevTopo <- factor(paste(odata$ElevFact,odata$Topo,sep="-"))
odata_et <- merge(odata,plot_data_diff)
odata_et_gr <- group_by(odata_et,ElevTopo,Species,Species_num)
odata_spp_topo <- summarise(odata_et_gr,prevalence=length(unique(Site[PercFlwr>0]))/5)
spp_topo_prev <- dcast(odata_spp_topo,formula=ElevTopo~Species_num)
spp_topo_prev[is.na(spp_topo_prev)] <- 0
spp_topo_prev_l <- melt(spp_topo_prev,id.vars="ElevTopo")
colnames(spp_topo_prev_l) <- c("ElevTopo","Species_num","Prevalence")
spp_topo_prev_l <- merge(spp_topo_prev_l,plot_data_diff,all.x=TRUE)
spp_topo_prev_l <- merge(spp_topo_prev_l,spp,by.x="Species_num",by.y="Num",all.x=TRUE)
write.csv(spp_topo_prev_l,"./output/spp_topo_prev.csv",row.names=FALSE)

##Gets scaling from original data.
odata <- read.csv("./scratch/final_train_data.csv", header=TRUE)
moist_center <- attributes(scale(odata$moist_days))$'scaled:center'
moist_scale <- attributes(scale(odata$moist_days))$'scaled:scale'
gdd_center <- attributes(scale(odata$gdd50))$'scaled:center'
gdd_scale <- attributes(scale(odata$gdd50))$'scaled:scale'

####Plots typical community curves for each topo band and elevation.
doy <- seq(40,300,by=1)
doy_scaled <- doy / 100 - 1.7

early_sdd <- spp_topo_prev_l$mean_sdd_2015
early_sdd_scaled <- early_sdd / 100 - 1.7
med_sdd <- spp_topo_prev_l$mean_sdd
med_sdd_scaled <- med_sdd / 100 - 1.7

early_gdd <- spp_topo_prev_l$mean_gdd_2015
early_gdd_scaled <- (early_gdd - gdd_center) / gdd_scale
med_gdd <- spp_topo_prev_l$mean_gdd
med_gdd_scaled <- (med_gdd - gdd_center) / gdd_scale

early_sm <- spp_topo_prev_l$mean_moist_2015
early_sm_scaled <- (early_sm - moist_center) / moist_scale
med_sm <- spp_topo_prev_l$mean_moist
med_sm_scaled <- (med_sm - moist_center) / moist_scale

early_plot_preds <- data.frame(Species=spp_topo_prev_l$Species,
                               Spp_num=spp_topo_prev_l$Species_num,
                               ElevTopo=spp_topo_prev_l$ElevTopo,
                               Elev=spp_topo_prev_l$ElevFact,
                               Topo=spp_topo_prev_l$Topo,
                               Prevalence=spp_topo_prev_l$Prevalence,
                               SDD_uns=early_sdd,
                               SDD=early_sdd_scaled,
                               GDD_uns=early_gdd,
                               GDD=early_gdd_scaled,
                               Moist_uns=early_sm,
                               Moist=early_sm_scaled)

early_day_preds <- expand.grid(Melt = "Early",
                               doy = doy_scaled,
                               Species=levels(spp_topo_prev_l$Species),
                               ElevTopo=levels(spp_topo_prev_l$ElevTopo))
early_day_preds$doy_uns <- (early_day_preds$doy + 1.7) * 100
early_jags_preds <- left_join(early_day_preds,early_plot_preds)


med_plot_preds <- data.frame(Species=spp_topo_prev_l$Species,
                               Spp_num=spp_topo_prev_l$Species_num,
                               ElevTopo=spp_topo_prev_l$ElevTopo,
                               Elev=spp_topo_prev_l$ElevFact,
                               Topo=spp_topo_prev_l$Topo,
                               Prevalence=spp_topo_prev_l$Prevalence,
                               SDD_uns=med_sdd,
                               SDD=med_sdd_scaled,
                               GDD_uns=med_gdd,
                               GDD=med_gdd_scaled,
                               Moist_uns=med_sm,
                               Moist=med_sm_scaled)

med_day_preds <- expand.grid(Melt = "Typical",
                               doy = doy_scaled,
                               Species=levels(spp_topo_prev_l$Species),
                               ElevTopo=levels(spp_topo_prev_l$ElevTopo))
med_day_preds$doy_uns <- (med_day_preds$doy + 1.7) * 100
med_jags_preds <- left_join(med_day_preds,med_plot_preds)

jags_pred_spp <- rbind(early_jags_preds,med_jags_preds)

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
spp.indices <- 1:48

param_array <- create_param_array(out=out,n.sample=n.samples,spp.indices=spp.indices,
                                param.names=param.names)
quantfun <- function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))}
param_quants <- apply(param_array,MARGIN=c(1,2),FUN=quantfun)

dimnames(param_quants) <- list(quantile=c("lwr","lwr25","med","upr75","upr"),
                               Spp_num=c(1:48),
                               Parameter=dimnames(param_array)[[2]])
param_meds <- data.frame(param_quants[3,,])
param_meds$Species_num <- rownames(param_meds)

param_quants_m <- melt(param_quants)
param_quants_m <- merge(param_quants_m,spp,by.x="Spp_num",by.y="Num")
write.csv(param_quants_m,"./output/species_param_quantiles.csv")

plot_params_opt <- c("opt_slope.s.sdd","opt_slope.s.moist","opt_slope.s.gdd")
plot_params_width <- c("width_slope.s.sdd","width_slope.s.moist","width_slope.s.gdd")

params_opt_slope_sdd <- filter(param_quants_m, Parameter %in% plot_params_opt[1])
params_opt_slope_sdd <- dcast(params_opt_slope_sdd,formula=Species~quantile)
p1 <- ggplot(params_opt_slope_sdd)+
        geom_point(aes(x=med,y=Species),alpha=0.5)+
        geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=0.5)+
        geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
                     lwd=1.5,alpha=0.5)+
        geom_point(aes(x=med,y=Species),alpha=1,
                   data=filter(params_opt_slope_sdd,lwr > 0 | upr < 0))+
        geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=1,
                     data=filter(params_opt_slope_sdd,lwr > 0 | upr < 0))+
        geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
                     lwd=1.5,alpha=1,data=filter(params_opt_slope_sdd,lwr > 0 | upr < 0))+
        scale_color_discrete(guide=FALSE)+
        scale_x_continuous("Snow Melt")+
        theme_bw()+
        theme(axis.text.y=element_text(face="italic"))

params_opt_slope_moist <- filter(param_quants_m, Parameter %in% plot_params_opt[2])
params_opt_slope_moist <- dcast(params_opt_slope_moist,formula=Species~quantile)
p2 <- ggplot(params_opt_slope_moist)+
  geom_point(aes(x=med,y=Species),alpha=0.2)+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=0.2)+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=0.2)+
  geom_point(aes(x=med,y=Species),alpha=1,
             data=filter(params_opt_slope_moist,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=1,
               data=filter(params_opt_slope_moist,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=1,data=filter(params_opt_slope_moist,lwr > 0 | upr < 0))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  scale_x_continuous("Soil Moist.")+
  scale_color_discrete(guide=FALSE)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())


params_opt_slope_gdd <- filter(param_quants_m, Parameter %in% plot_params_opt[3])
params_opt_slope_gdd <- dcast(params_opt_slope_gdd,formula=Species~quantile)
p3 <- ggplot(params_opt_slope_gdd)+
  geom_point(aes(x=med,y=Species),alpha=0.2)+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=0.2)+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=0.2)+
  geom_point(aes(x=med,y=Species),alpha=1,
             data=filter(params_opt_slope_gdd,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=1,
               data=filter(params_opt_slope_gdd,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=1,data=filter(params_opt_slope_gdd,lwr > 0 | upr < 0))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  scale_x_continuous("GDD")+
  scale_color_discrete(guide=FALSE)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())


params_width_slope_sdd <- filter(param_quants_m, Parameter %in% plot_params_width[1])
params_width_slope_sdd <- dcast(params_width_slope_sdd,formula=Species~quantile)
p4 <- ggplot(params_width_slope_sdd)+
  geom_point(aes(x=med,y=Species),alpha=0.2)+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=0.2)+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=0.2)+
  geom_point(aes(x=med,y=Species),alpha=1,
             data=filter(params_width_slope_sdd,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=1,
               data=filter(params_width_slope_sdd,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=1,data=filter(params_width_slope_sdd,lwr > 0 | upr < 0))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  scale_x_continuous("Snow Melt")+
  scale_color_discrete(guide=FALSE)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())



params_width_slope_moist <- filter(param_quants_m, Parameter %in% plot_params_width[2])
params_width_slope_moist <- dcast(params_width_slope_moist,formula=Species~quantile)
p5 <- ggplot(params_width_slope_moist)+
  geom_point(aes(x=med,y=Species),alpha=0.2)+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=0.2)+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=0.2)+
  geom_point(aes(x=med,y=Species),alpha=1,
             data=filter(params_width_slope_moist,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),alpha=1,
               data=filter(params_width_slope_moist,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=1,data=filter(params_width_slope_moist,lwr > 0 | upr < 0))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  scale_x_continuous("Soil Moist.")+
  scale_color_discrete(guide=FALSE)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

params_width_slope_gdd <- filter(param_quants_m, Parameter %in% plot_params_width[3])
params_width_slope_gdd <- dcast(params_width_slope_gdd,formula=Species~quantile)
p6 <- ggplot(params_width_slope_gdd)+
  geom_point(aes(x=med,y=Species),alpha=0.2)+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),color="grey60")+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=0.2)+
  geom_point(aes(x=med,y=Species),alpha=1,
             data=filter(params_width_slope_gdd,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr,xend=upr,y=Species,yend=Species),
               alpha=1,data=filter(params_width_slope_gdd,lwr > 0 | upr < 0))+
  geom_segment(aes(x=lwr25,xend=upr75,y=Species,yend=Species),
               lwd=1.5,alpha=1,data=filter(params_width_slope_gdd,lwr > 0 | upr < 0))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  scale_x_continuous("GDD")+
  scale_color_discrete(guide=FALSE)+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank())

library(gridExtra)
#pdf("./figs/spp_params_cred_moist_gdd.pdf",width=12,height=6)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=6,widths=c(2.2,1,1,1,1,1))
#dev.off()

spp_nums <- sort(unique(as.numeric(as.character(jags_pred_spp$Spp_num))))

flwr_prob <- matrix(NA,nrow=nrow(jags_pred_spp),ncol=n.samples)
colnames(flwr_prob) <- paste("sample_",1:n.samples,sep="")

for (i in spp_nums){
  print(paste("Predictions for species",i,"of",length(spp_nums)))
  jags_pred_f <- filter(jags_pred_spp,Spp_num==i)
  for (j in 1:dim(param_array)[3]){
    flwr_prob[jags_pred_spp$Spp_num==i,j] <- jags_fit_fun_moist_gdd(height_int=param_array[i,2,j],
                                                                    height_slope_sdd=param_array[i,1,j],
                                                                    opt_int=param_array[i,6,j],
                                                                    opt_slope_sdd=param_array[i,5,j],
                                                                    opt_slope_moist=param_array[i,4,j],
                                                                    opt_slope_gdd=param_array[i,3,j],
                                                                    width_int=param_array[i,10,j],
                                                                    width_slope_sdd=param_array[i,9,j],
                                                                    width_slope_moist=param_array[i,8,j],
                                                                    width_slope_gdd=param_array[i,7,j],
                                                                    x1vec=jags_pred_f$doy,
                                                                    x2vec=jags_pred_f$SDD,
                                                                    x3vec=jags_pred_f$Moist,
                                                                    x4vec=jags_pred_f$GDD)
  }
}

jags_pred_attr <- data.frame(jags_pred_spp,flwr_prob)
save(jags_pred_attr,file="./output/spp_ElevTopo_early_late_samples_attrib.Rdata",compress=TRUE)
#load("./output/spp_ElevTopo_early_late_samples_attrib.Rdata")

##Calculates quantiles for each species-site combo.
quantfun <- function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)}
flwr_quant <- t(apply(flwr_prob,MARGIN=1,FUN=quantfun))
colnames(flwr_quant) <- c("flwr_lwr","flwr_lwr25","flwr_med","flwr_upr75","flwr_upr")
jags_pred_quant <- data.frame(jags_pred_spp,flwr_quant)

jags_pred_spp$flwr_prob <- flwr_quant[,3]
jags_pred_spp$Elev <- factor(jags_pred_spp$Elev,levels=c("1901","1791","1680","1570","1490"))
jags_pred_spp$Topo <- factor(jags_pred_spp$Topo,levels=c("r","s","c"),labels=c("Ridge","Slope","Cove"))
jags_pred_spp$DSS <- jags_pred_spp$doy_uns - jags_pred_spp$SDD_uns
jags_pred_spp$flwr_prob_pres <- jags_pred_spp$flwr_prob * jags_pred_spp$Prevalence
jags_pred_spp$prob_sq <- NA
jags_pred_spp$prob_sq[jags_pred_spp$Melt=="Typical"] <- jags_pred_spp$flwr_prob_pres[jags_pred_spp$Melt=="Typical"] / 2 + 0.5
jags_pred_spp$prob_sq[jags_pred_spp$Melt=="Early"] <- jags_pred_spp$flwr_prob_pres[jags_pred_spp$Melt=="Early"] / 2
jags_pred_pres <- filter(jags_pred_spp,Prevalence>0)
jags_pred_gr <- group_by(jags_pred_pres,Melt,Elev,Topo,Spp_num,ElevTopo,Species,Prevalence)
jags_pred_max <- summarise(jags_pred_gr,max_prob=max(prob_sq),
                           max_DSS=DSS[which.max(prob_sq)],
                           start_DSS=obs_intervals(data.frame(dss=DSS,pred=flwr_prob),threshold=0.1)[1],
                           end_DSS=obs_intervals(data.frame(dss=DSS,pred=flwr_prob),threshold=0.1)[2],
                           max_DOY=doy_uns[which.max(prob_sq)],
                           start_DOY=obs_intervals(data.frame(dss=doy_uns,pred=flwr_prob),threshold=0.1)[1],
                           end_DOY=obs_intervals(data.frame(dss=doy_uns,pred=flwr_prob),threshold=0.1)[2])
jags_pred_change <- data.frame(jags_pred_max[jags_pred_max$Melt=="Early",],
                          jags_pred_max[jags_pred_max$Melt=="Typical",8:14])
colnames(jags_pred_change) <- c("Melt","Elev","Topo","Spp_num","ElevTopo","Species","Prevalence","max_prob_early","max_DSS_early","start_DSS_early",
                                "end_DSS_early","max_DOY_early","start_DOY_early","end_DOY_early",
                                "max_prob_late","max_DSS_late","start_DSS_late",
                                "end_DSS_late","max_DOY_late","start_DOY_late","end_DOY_late")
jags_pred_change$diff_start <- jags_pred_change$start_DSS_late - jags_pred_change$start_DSS_early
jags_pred_change$diff_end <- jags_pred_change$end_DSS_late - jags_pred_change$end_DSS_early
jags_pred_change$diff_max <- jags_pred_change$max_DSS_late - jags_pred_change$max_DSS_early
jags_pred_change$diff_maxprob <- jags_pred_change$max_prob_late - jags_pred_change$max_prob_early
jags_pred_change$diff_length <- (jags_pred_change$end_DSS_late - jags_pred_change$start_DSS_late)-
  (jags_pred_change$end_DSS_early - jags_pred_change$start_DSS_early) 
jags_pred_change_pres <- filter(jags_pred_change,Prevalence>0)
jags_pred_change_gr <- group_by(jags_pred_change_pres,Species,Spp_num)
jags_pred_change_spp <- summarise(jags_pred_change_gr,
                                 opt_DSS_early=mean(max_DSS_early),
                                 start_DSS_early=mean(start_DSS_early),
                                 end_DSS_early=mean(end_DSS_early),
                                 opt_DSS_med=mean(max_DSS_late),
                                 start_DSS_med=mean(start_DSS_late),
                                 end_DSS_med=mean(end_DSS_late),
                                 opt_DOY_early=mean(max_DOY_early),
                                 start_DOY_early=mean(start_DOY_early),
                                 end_DOY_early=mean(end_DOY_early),
                                 opt_DOY_med=mean(max_DOY_late),
                                 start_DOY_med=mean(start_DOY_late),
                                 end_DOY_med=mean(end_DOY_late))
jags_pred_change_spp$diff_start <- jags_pred_change_spp$start_DSS_med - jags_pred_change_spp$start_DSS_early
jags_pred_change_spp$diff_end <- jags_pred_change_spp$end_DSS_med - jags_pred_change_spp$end_DSS_early
jags_pred_change_spp$diff_max <- jags_pred_change_spp$opt_DSS_med - jags_pred_change_spp$opt_DSS_early
jags_pred_change_spp$diff_length <- (jags_pred_change_spp$end_DSS_med - jags_pred_change_spp$start_DSS_med)-
  (jags_pred_change_spp$end_DSS_early - jags_pred_change_spp$start_DSS_early) 
jags_pred_change_spp$length_med <-  (jags_pred_change_spp$end_DSS_med - jags_pred_change_spp$start_DSS_med)
jags_pred_change_spp$length_early <-  (jags_pred_change_spp$end_DSS_early - jags_pred_change_spp$start_DSS_early)
jags_pred_change_spp$start_diff_DOY <-  (jags_pred_change_spp$start_DOY_early - jags_pred_change_spp$start_DOY_med)
jags_pred_change_spp$end_diff_DOY <-  (jags_pred_change_spp$end_DOY_early - jags_pred_change_spp$end_DOY_med)
jags_pred_change_spp$opt_diff_DOY <-  (jags_pred_change_spp$opt_DOY_early - jags_pred_change_spp$opt_DOY_med)


Spp_order <- levels(jags_pred_change_spp$Species)[order(jags_pred_change_spp$opt_DSS_med)]
jags_pred_change$Spp_ord <- factor(jags_pred_change$Species,levels=Spp_order)
jags_pred_change_pres$Spp_ord <- factor(jags_pred_change_pres$Species,levels=Spp_order)
jags_pred_spp$Spp_ord <- factor(jags_pred_spp$Species,levels=Spp_order)
jags_pred_pres$Spp_ord <- factor(jags_pred_pres$Species,levels=Spp_order)
jags_pred_change_spp$Spp_ord <- factor(jags_pred_change_spp$Species,levels=Spp_order)


##Merges species attribute data
jags_pred_change_attr <- merge(jags_pred_change_spp,spp_sum[,-2],by.x="Spp_num",
                               by.y="Species_num",all.x=TRUE)
jags_pred_change_attr <- merge(jags_pred_change_attr,param_meds,
                               by.x="Spp_num",by.y="Species_num")
jags_pred_change_attr$Spp_ord <- factor(jags_pred_change_attr$Species,
                                        levels=Spp_order)
write.csv(jags_pred_change_attr,"./output/theobald_spp_early_late_topo_params.csv")
jags_pred_change_attr <- read.csv("./output/theobald_spp_early_late_topo_params.csv")

####Community statistics for manuscript.
min_shift_opt <- min(jags_pred_change_attr$opt_diff_DOY)
max_shift_opt <- max(jags_pred_change_attr$opt_diff_DOY)
mean_shift_opt <- mean(jags_pred_change_attr$opt_diff_DOY)
sd_shift_opt <- sd(jags_pred_change_attr$opt_diff_DOY)

melt_sens_min <- min(jags_pred_change_attr$opt_slope.s.sdd)
melt_sens_max <- max(jags_pred_change_attr$opt_slope.s.sdd)
melt_sens_mean <- mean(jags_pred_change_attr$opt_slope.s.sdd)
melt_sens_sd <- sd(jags_pred_change_attr$opt_slope.s.sdd)

moist_sens_min <- min(jags_pred_change_attr$opt_slope.s.moist)
moist_sens_max <- max(jags_pred_change_attr$opt_slope.s.moist)
moist_sens_mean <- mean(jags_pred_change_attr$opt_slope.s.moist)
moist_sens_sd <- sd(jags_pred_change_attr$opt_slope.s.moist)

gdd_sens_min <- min(jags_pred_change_attr$opt_slope.s.gdd)
gdd_sens_max <- max(jags_pred_change_attr$opt_slope.s.gdd)
gdd_sens_mean <- mean(jags_pred_change_attr$opt_slope.s.gdd)
gdd_sens_sd <- sd(jags_pred_change_attr$opt_slope.s.gdd)

diff_length_min <- min(jags_pred_change_attr$diff_length)
diff_length_max <- max(jags_pred_change_attr$diff_length)
diff_length_mean <- mean(jags_pred_change_attr$diff_length)
diff_length_sd <- sd(jags_pred_change_attr$diff_length)
jags_pred_change_attr$lengthen <- as.numeric(jags_pred_change_attr$diff_length < 0)
nspp_lengthen <- sum(jags_pred_change_attr$lengthen)
nspp_shorten <- 48 - nspp_lengthen
pct_lengthen <- (nspp_lengthen / 48) * 100

melt_sens_min_len <- min(jags_pred_change_attr$width_slope.s.sdd)
melt_sens_max_len <- max(jags_pred_change_attr$width_slope.s.sdd)
melt_sens_mean_len <- mean(jags_pred_change_attr$width_slope.s.sdd)
melt_sens_sd_len <- sd(jags_pred_change_attr$width_slope.s.sdd)

moist_sens_min_len <- min(jags_pred_change_attr$width_slope.s.moist)
moist_sens_max_len <- max(jags_pred_change_attr$width_slope.s.moist)
moist_sens_mean_len <- mean(jags_pred_change_attr$width_slope.s.moist)
moist_sens_sd_len <- sd(jags_pred_change_attr$width_slope.s.moist)

gdd_sens_min_len <- min(jags_pred_change_attr$width_slope.s.gdd)
gdd_sens_max_len <- max(jags_pred_change_attr$width_slope.s.gdd)
gdd_sens_mean_len <- mean(jags_pred_change_attr$width_slope.s.gdd)
gdd_sens_sd_len <- sd(jags_pred_change_attr$width_slope.s.gdd)


##Calculates flowering species richness for each topo band.
jags_pred_topo <- group_by(jags_pred_spp,Melt,Elev,Topo,ElevTopo,DSS)
jags_pred_topo_rich <- summarise(jags_pred_topo,nspp=length(unique(Species[Prevalence>0])),
                                                flwr_rich=sum(flwr_prob_pres))

#pdf("./figs/flwr_richness_topo_elev.pdf",height=8,width=5)
ggplot(jags_pred_topo_rich)+
  geom_line(aes(x=DSS,y=flwr_rich,linetype=Melt))+
  facet_grid(facets=Elev~Topo)+
  xlim(c(-10,90))+
  scale_y_continuous("Number of Species in Flower",breaks=seq(0,8,by=2))+
  xlab("Days Since Snow Melt")+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()


g1 <- ggplot(jags_pred_change_attr)+
  geom_point(aes(x=opt_DSS_med,y=opt_DSS_early,size=prevalence),
             position=position_jitter(width=0.2))+
  geom_smooth(aes(x=opt_DSS_med,y=opt_DSS_early),color="black",method="lm",se=TRUE)+
  geom_abline(aes(slope=1,intercept=0),linetype="dotted")+
  scale_size_continuous("Prevalence",limits=c(0,1),range=c(0.1,4),guide=FALSE)+
  xlim(c(0,80))+
  ylim(c(0,80))+
  ylab("Peak Flower DSS (2015)")+
  xlab("Peak Flower DSS (2010-2014)")+
  theme_bw()+
  theme(panel.grid=element_blank())

g2 <- ggplot(jags_pred_change_attr)+
  geom_point(aes(x=start_DSS_med,y=start_DSS_early,size=prevalence),
             position=position_jitter(width=0.2))+
  geom_smooth(aes(x=start_DSS_med,y=start_DSS_early),color="black",method="lm",se=TRUE)+
  geom_abline(aes(slope=1,intercept=0),linetype="dotted")+
  scale_size_continuous("Prevalence",limits=c(0,1),range=c(0.1,4),guide=FALSE)+
  #xlim(c(10,40))+
  #ylim(c(10,40))+
  ylab("First Flower DSS (2015)")+
  xlab("First Flower DSS (2010-2014)")+
  theme_bw()+
  theme(panel.grid=element_blank())

g3 <- ggplot(jags_pred_change_attr)+
  geom_point(aes(x=end_DSS_med,y=end_DSS_early,size=prevalence),
             position=position_jitter(width=0.2))+
  geom_smooth(aes(x=end_DSS_med,y=end_DSS_early),color="black",method="lm",se=TRUE)+
  geom_abline(aes(slope=1,intercept=0),linetype="dotted")+
  scale_size_continuous("Prevalence",limits=c(0,1),range=c(0.1,4))+
  #xlim(c(10,40))+
  #ylim(c(10,40))+
  ylab("Last Flower DSS (2015)")+
  xlab("Last Flower DSS (2010-2014)")+
  theme_bw()+
  theme(panel.grid=element_blank())


library(gridExtra)
#pdf("./figs/SDD_opt_width_early_late.pdf",width=12,height=4)
grid.arrange(g1,g2,g3,nrow=1,widths=c(1,1,1.3))
#dev.off()

##Figure showing flowering periods for erythronium, valerian, and gentian.
ggplot(filter(jags_pred_change_attr,Species %in% c("Erythronium montanum",
                                     "Valeriana sitchensis",
                                     "Gentiana calycosa")))+
    geom_point(aes(x=opt_DSS_early,xend=opt_DSS_early,
                     y=as.numeric(factor(Species))+0.1,yend=as.numeric(factor(Species))+0.1),
                   shape=21)+
    geom_segment(aes(x=start_DSS_early,xend=end_DSS_early,
                   y=as.numeric(factor(Species))+0.1,yend=as.numeric(factor(Species))+0.1),
               linetype="solid",show.legend=TRUE)+
    geom_point(aes(x=opt_DSS_med,xend=opt_DSS_med,
                   y=as.numeric(factor(Species))-0.1,yend=as.numeric(factor(Species))-0.1),
               shape=21)+
    geom_segment(aes(x=start_DSS_med,xend=end_DSS_med,
                     y=as.numeric(factor(Species))-0.1,yend=as.numeric(factor(Species))-0.1),
                 linetype="dotted",show.legend=TRUE)+
    scale_y_continuous("",breaks=c(1:3),labels=c("Erythronium montanum",
                                                 "Gentiana calycosa",
                                                 "Valeriana sitchensis"))+
    xlab("Days Since Snow")+
    scale_linetype_manual("Year",values=c("solid","dotted"),
                          labels=c("2015","2010-2014"),
                          guide=guide_legend(position="right"))+
  theme_bw()
  

p1 <- ggplot(jags_pred_spp)+
  geom_hline(aes(yintercept=0.5),linetype=1)+
  geom_vline(aes(xintercept=0),linetype=2)+
  geom_line(aes(x=DSS,y=prob_sq,color=Spp_ord,linetype=Melt),
            data=jags_pred_spp)+
  scale_x_continuous("Days Since Snow Melt",limits=c(-10,100))+
  scale_y_continuous("Flowering Probability",breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0,0.25,0,0.25,1))+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(legend.position="right")
#pdf("./figs/spp_curves_elev_topo.pdf",width=15,height=10)
print(p1)
#dev.off()


month_breaks <- c(1,32,60,91,121,152,182,213,244,274,305,335)
month_labels <- c("Jan.","Feb.","Mar.","Apr.","May.","Jun.",
                  "Jul.","Aug.","Sep.","Oct.","Nov.","Dec.")

p2 <-  ggplot(jags_pred_pres)+
  geom_line(aes(x=doy_uns,y=prob_sq,color=Spp_ord,linetype=Melt),
            data=jags_pred_pres)+
  geom_hline(aes(yintercept=0.5),linetype=1)+
  geom_segment(aes(x=SDD_uns,xend=SDD_uns,y=0,yend=0.5),linetype=2,
               data=filter(jags_pred_pres,Melt=="Early"))+
  geom_segment(aes(x=SDD_uns,xend=SDD_uns,y=0.5,yend=1),linetype=2,
               data=filter(jags_pred_pres,Melt=="Typical"))+
  scale_x_continuous("Day of Year",limits=c(90,290),
                     breaks=month_breaks,labels=month_labels)+
  scale_y_continuous("Flowering Probability",breaks=c(0,0.25,0.5,0.75,1),
                     labels=c(0,0.25,0,0.25,1))+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(legend.position="bottom")

#pdf("./figs/spp_curves_elev_topo_DOY.pdf",width=15,height=10)
print(p2)
#dev.off()


p3 <- ggplot(jags_pred_change_attr)+
  geom_vline(aes(xintercept=0),linetype=2)+
  geom_segment(aes(x=start_DSS_early,xend=end_DSS_early,
                   y=as.numeric(Spp_ord)+0.2,yend=as.numeric(Spp_ord)+0.2,color=Spp_ord),
               linetype=1,data=jags_pred_change_attr)+
  geom_point(aes(x=opt_DSS_early,y=as.numeric(Spp_ord)+0.2,color=Spp_ord),
             data=jags_pred_change_attr)+
  geom_segment(aes(x=start_DSS_med,xend=end_DSS_med,
                   y=as.numeric(Spp_ord)-0.2,yend=as.numeric(Spp_ord)-0.2,color=Species),
               linetype=2,alpha=0.6,data=jags_pred_change_attr)+
  geom_point(aes(x=opt_DSS_med,y=as.numeric(Spp_ord)-0.2,color=Spp_ord),
             alpha=0.5,data=jags_pred_change_attr)+
  scale_y_continuous("Species",breaks=1:48,labels=levels(jags_pred_change_attr$Spp_ord))+
  scale_x_continuous("Days Since Snow",limits=c(-10,100))+
  scale_color_discrete(guide=FALSE)+
  #facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(axis.text.y=element_text(size=11,face="italic"))

#pdf("./figs/spp_lines_all.pdf",width=6,height=8)
print(p3)
#dev.off()
