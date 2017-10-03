####Script to measure reassembly using fit phenology model.
####Author: Ian Breckheimer
####Data: Elli Theobald
####Last Updated: 1 October 2017

library(dplyr)
library(reshape2)
library(stringr)
library(ggplot2)

source("./code/pheno_reassembly_functions.R")
load("./output/spp_ElevTopo_early_late_samples_attrib.Rdata")

##Subsets by elevation to avoid memory problems
jags_pred_1490 <- tbl_df(jags_pred_attr[jags_pred_attr$ElevTopo %in% 
                                          c("1490-r","1490-s","1490-c"),])
jags_pred_1570 <- tbl_df(jags_pred_attr[jags_pred_attr$ElevTopo %in%
                                          c("1570-r","1570-s","1570-c"),])
jags_pred_1680 <- tbl_df(jags_pred_attr[jags_pred_attr$ElevTopo %in%
                                          c("1680-r","1680-s","1680-c"),])
jags_pred_1791 <- tbl_df(jags_pred_attr[jags_pred_attr$ElevTopo %in%
                                          c("1791-r","1791-s","1791-c"),])
jags_pred_1901 <- tbl_df(jags_pred_attr[jags_pred_attr$ElevTopo %in%
                                          c("1901-r","1901-s","1901-c"),])
rm(jags_pred_attr)

####Calculates posterior quantiles for flowering species richness####

quants_rich_1490 <- jags_pred_1490 %>%
                      select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
                             contains("sample")) %>%
                      melt(id.vars=c("Melt","doy_uns","ElevTopo",
                                     "SDD_uns","Spp_num","Prevalence")) %>%
                      mutate(flwr_prob_pres=value*Prevalence,
                             DSS=doy_uns - SDD_uns) %>%
                      group_by(Melt,ElevTopo,DSS,Spp_num,Prevalence) %>%
                      summarise(flwr_med=quantile(flwr_prob_pres,0.5),
                                flwr_lwr=quantile(flwr_prob_pres,0.1),
                                flwr_upr=quantile(flwr_prob_pres,0.9)) %>%
                      group_by(Melt,ElevTopo,DSS) %>%
                      summarise(rich_lwr=sum(flwr_lwr),
                                rich_med=sum(flwr_med),
                                rich_upr=sum(flwr_upr))
quants_rich_1570 <- jags_pred_1570 %>%
                    select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
                           contains("sample")) %>%
                    melt(id.vars=c("Melt","doy_uns","ElevTopo",
                                   "SDD_uns","Spp_num","Prevalence")) %>%
                    mutate(flwr_prob_pres=value*Prevalence,
                           DSS=doy_uns - SDD_uns) %>%
                    group_by(Melt,ElevTopo,DSS,Spp_num,Prevalence) %>%
                    summarise(flwr_med=quantile(flwr_prob_pres,0.5),
                              flwr_lwr=quantile(flwr_prob_pres,0.1),
                              flwr_upr=quantile(flwr_prob_pres,0.9)) %>%
                    group_by(Melt,ElevTopo,DSS) %>%
                    summarise(rich_lwr=sum(flwr_lwr),
                              rich_med=sum(flwr_med),
                              rich_upr=sum(flwr_upr))
quants_rich_1680 <- jags_pred_1680 %>%
          select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
                 contains("sample")) %>%
          melt(id.vars=c("Melt","doy_uns","ElevTopo",
                         "SDD_uns","Spp_num","Prevalence")) %>%
          mutate(flwr_prob_pres=value*Prevalence,
                 DSS=doy_uns - SDD_uns) %>%
          group_by(Melt,ElevTopo,DSS,Spp_num,Prevalence) %>%
          summarise(flwr_med=quantile(flwr_prob_pres,0.5),
                    flwr_lwr=quantile(flwr_prob_pres,0.1),
                    flwr_upr=quantile(flwr_prob_pres,0.9)) %>%
          group_by(Melt,ElevTopo,DSS) %>%
          summarise(rich_lwr=sum(flwr_lwr),
                    rich_med=sum(flwr_med),
                    rich_upr=sum(flwr_upr))
quants_rich_1791 <- jags_pred_1791 %>%
          select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
                 contains("sample")) %>%
          melt(id.vars=c("Melt","doy_uns","ElevTopo",
                         "SDD_uns","Spp_num","Prevalence")) %>%
          mutate(flwr_prob_pres=value*Prevalence,
                 DSS=doy_uns - SDD_uns) %>%
          group_by(Melt,ElevTopo,DSS,Spp_num,Prevalence) %>%
          summarise(flwr_med=quantile(flwr_prob_pres,0.5),
                    flwr_lwr=quantile(flwr_prob_pres,0.1),
                    flwr_upr=quantile(flwr_prob_pres,0.9)) %>%
          group_by(Melt,ElevTopo,DSS) %>%
          summarise(rich_lwr=sum(flwr_lwr),
                    rich_med=sum(flwr_med),
                    rich_upr=sum(flwr_upr))
quants_rich_1901 <- jags_pred_1901 %>%
        select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
               contains("sample")) %>%
        melt(id.vars=c("Melt","doy_uns","ElevTopo",
                       "SDD_uns","Spp_num","Prevalence")) %>%
        mutate(flwr_prob_pres=value*Prevalence,
               DSS=doy_uns - SDD_uns) %>%
        group_by(Melt,ElevTopo,DSS,Spp_num,Prevalence) %>%
        summarise(flwr_med=quantile(flwr_prob_pres,0.5),
                  flwr_lwr=quantile(flwr_prob_pres,0.1),
                  flwr_upr=quantile(flwr_prob_pres,0.9)) %>%
        group_by(Melt,ElevTopo,DSS) %>%
        summarise(rich_lwr=sum(flwr_lwr),
                  rich_med=sum(flwr_med),
                  rich_upr=sum(flwr_upr))

quants_rich <- rbind(quants_rich_1490,
                     quants_rich_1570,
                     quants_rich_1680,
                     quants_rich_1791,
                     quants_rich_1901)

quants_rich$Elev <- str_split_fixed(quants_rich$ElevTopo,pattern="-",n=2)[,1]
quants_rich$Elev <- factor(quants_rich$Elev,levels=c("1901","1791","1680","1570","1490"),
                           labels=c("1901m","1791m","1680m","1570m","1490m"))

quants_rich$Topo <- str_split_fixed(quants_rich$ElevTopo,pattern="-",n=2)[,2]
quants_rich$Topo <- factor(quants_rich$Topo,levels=c("r","s","c"),
                           labels=c("Ridge","Slope","Cove"))

##Calculates lag and season length.
quants_rich$season <- as.factor(as.numeric(quants_rich$rich_med >= 0.5))
quants_seas <- filter(quants_rich,season=="1")
quants_seas_grp <- group_by(quants_seas,Elev,Topo,Melt)
quants_seas_sum <- summarise(quants_seas_grp,start=min(DSS),
                                             peak=DSS[which.max(rich_med)],
                                             end=max(DSS))
quants_seas_sum$y_val <- 11
quants_seas_sum$y_val[quants_seas_sum$Melt=="Early"] <- 12

#pdf("./figs/fl_richness_elev_topo_cred80.pdf",width=6,height=5)
ggplot(quants_rich) +
  geom_ribbon(aes(x=DSS,ymin=rich_lwr,ymax=rich_upr,fill=Melt),alpha=0.3)+
  # geom_point(aes(x=start,y=y_val,color=Melt),shape=23,
  #            show.legend=FALSE,data=quants_seas_sum)+
  # geom_point(aes(x=peak,y=y_val,color=Melt),shape=20,
  #            show.legend=FALSE,data=quants_seas_sum)+
  # geom_point(aes(x=end,y=y_val,color=Melt),shape=23,
  #            show.legend=FALSE,data=quants_seas_sum)+
  # geom_segment(aes(x=start,xend=end,y=y_val,yend=y_val,color=Melt),linetype="dotted",
  #            show.legend=FALSE,data=quants_seas_sum)+
  geom_line(aes(x=DSS,y=rich_lwr,color=Melt),alpha=0.3)+
  geom_line(aes(x=DSS,y=rich_upr,color=Melt),alpha=0.3)+
  geom_line(aes(x=DSS,y=rich_med,color=Melt))+
  geom_vline(aes(xintercept=0),linetype="dotted")+
  xlim(c(-5,90))+
  xlab("Days Since Snow Melt")+
  scale_y_continuous("Number of Species in Flower",
                     breaks=c(0,5,10),limits=c(0,11))+
  scale_color_brewer("Melt",type="qual",palette=2,direction=-1)+
  scale_fill_brewer("Melt",type="qual",palette=2,direction=-1)+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()

####Calculates bray-curtis distances between early and late-melt conditions for each site.
dist_1490 <- site_bc_dist(jags_pred_1490)
dist_1570 <- site_bc_dist(jags_pred_1570)
dist_1680 <- site_bc_dist(jags_pred_1680)
dist_1791 <- site_bc_dist(jags_pred_1791)
dist_1901 <- site_bc_dist(jags_pred_1901)

##Combines estimates for plotting
dist_1490_dss <- dist_1490$quants_DSS
dist_1490_dss$Elev <- "1490"
dist_1570_dss <- dist_1570$quants_DSS
dist_1570_dss$Elev <- "1570"
dist_1680_dss <- dist_1680$quants_DSS
dist_1680_dss$Elev <- "1680"
dist_1791_dss <- dist_1791$quants_DSS
dist_1791_dss$Elev <- "1791"
dist_1901_dss <- dist_1901$quants_DSS
dist_1901_dss$Elev <- "1901"
dist <- rbind(dist_1490_dss,dist_1570_dss,dist_1680_dss,dist_1791_dss,dist_1901_dss)

dist$Topo <- factor(dist$Topo,levels=c("r","s","c"),
                    labels=c("Ridge","Slope","Cove"))
dist$Elev <- factor(dist$Elev,levels=c("1901","1791","1680","1570","1490"),
                    labels=c("1901m","1791m","1680m","1570m","1490m"))

##Extracts seasonal averages.
dist_1490_seas <- data.frame(dist_1490$quants_seas)
dist_1490_seas$quant <- factor(rownames(dist_1490_seas))
dist_1490_seas$Elev <- "1490"
dist_1570_seas <- data.frame(dist_1570$quants_seas)
dist_1570_seas$quant <- factor(rownames(dist_1570_seas))
dist_1570_seas$Elev <- "1570"
dist_1680_seas <- data.frame(dist_1680$quants_seas)
dist_1680_seas$quant <- factor(rownames(dist_1680_seas))
dist_1680_seas$Elev <- "1680"
dist_1791_seas <- data.frame(dist_1791$quants_seas)
dist_1791_seas$quant <- factor(rownames(dist_1791_seas))
dist_1791_seas$Elev <- "1791"
dist_1901_seas <- data.frame(dist_1901$quants_seas)
dist_1901_seas$quant <- factor(rownames(dist_1901_seas))
dist_1901_seas$Elev <- "1901"
dist_seas <- rbind(dist_1490_seas,dist_1570_seas,dist_1680_seas,dist_1791_seas,dist_1901_seas)
dist_seas <- melt(dist_seas,id.vars=c("quant","Elev"))
colnames(dist_seas) <- c("quant","Elev","Topo","value")
dist_seas_c <- dcast(dist_seas,formula=Elev+Topo~quant)

dist_seas_c$Topo <- factor(dist_seas_c$Topo,levels=c("r","s","c"),
                           labels=c("Ridge","Slope","Cove"))

####Calculates community novelty between early and late-melt conditions for each site.
nov_1490 <- site_bc_nov(jags_pred_1490)
nov_1570 <- site_bc_nov(jags_pred_1570)
nov_1680 <- site_bc_nov(jags_pred_1680)
nov_1791 <- site_bc_nov(jags_pred_1791)
nov_1901 <- site_bc_nov(jags_pred_1901)

nov_1490_dss <- nov_1490$quants_DSS
nov_1490_dss$Elev <- "1490"
nov_1570_dss <- nov_1570$quants_DSS
nov_1570_dss$Elev <- "1570"
nov_1680_dss <- nov_1680$quants_DSS
nov_1680_dss$Elev <- "1680"
nov_1791_dss <- nov_1791$quants_DSS
nov_1791_dss$Elev <- "1791"
nov_1901_dss <- nov_1901$quants_DSS
nov_1901_dss$Elev <- "1901"
nov <- rbind(nov_1490_dss,nov_1570_dss,nov_1680_dss,nov_1791_dss,nov_1901_dss)

nov$Topo <- factor(nov$Topo,levels=c("r","s","c"),
                    labels=c("Ridge","Slope","Cove"))
nov$Elev <- factor(nov$Elev,levels=c("1901","1791","1680","1570","1490"),
                    labels=c("1901m","1791m","1680m","1570m","1490m"))

##Extracts seasonal averages.
nov_1490_seas <- data.frame(nov_1490$quants_seas)
nov_1490_seas$quant <- factor(rownames(nov_1490_seas))
nov_1490_seas$Elev <- "1490"
nov_1570_seas <- data.frame(nov_1570$quants_seas)
nov_1570_seas$quant <- factor(rownames(nov_1570_seas))
nov_1570_seas$Elev <- "1570"
nov_1680_seas <- data.frame(nov_1680$quants_seas)
nov_1680_seas$quant <- factor(rownames(nov_1680_seas))
nov_1680_seas$Elev <- "1680"
nov_1791_seas <- data.frame(nov_1791$quants_seas)
nov_1791_seas$quant <- factor(rownames(nov_1791_seas))
nov_1791_seas$Elev <- "1791"
nov_1901_seas <- data.frame(nov_1901$quants_seas)
nov_1901_seas$quant <- factor(rownames(nov_1901_seas))
nov_1901_seas$Elev <- "1901"
nov_seas <- rbind(nov_1490_seas,nov_1570_seas,nov_1680_seas,nov_1791_seas,nov_1901_seas)
nov_seas <- melt(nov_seas,id.vars=c("quant","Elev"))
colnames(nov_seas) <- c("quant","Elev","Topo","value")
nov_seas_c <- dcast(nov_seas,formula=Elev+Topo~quant)

nov_seas_c$Topo <- factor(nov_seas_c$Topo,levels=c("r","s","c"),
                           labels=c("Ridge","Slope","Cove"))

##Extracts study-wide median novelty for each day.
all_median_samps <- array(NA,dim=c(5,91,500))
all_median_samps[1,,] <- nov_1490$nov_median_samples
all_median_samps[2,,] <- nov_1570$nov_median_samples
all_median_samps[3,,] <- nov_1680$nov_median_samples
all_median_samps[4,,] <- nov_1791$nov_median_samples
all_median_samps[5,,] <- nov_1901$nov_median_samples

quantfun <- function(x){ quantile(x,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
all_median_avg <- apply(all_median_samps,MARGIN=c(2,3),FUN=mean)
all_median_dss <- data.frame(t(apply(all_median_avg,MARGIN=1,FUN=quantfun)))
all_median_dss$DSS <- 0:90

####Computes posterior quantiles for lag and season length####
quants_lag_1490 <- jags_pred_1490 %>%
  select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
         contains("sample")) %>%
  melt(id.vars=c("Melt","doy_uns","ElevTopo",
                 "SDD_uns","Spp_num","Prevalence")) %>%
  mutate(flwr_prob_pres=value*Prevalence,
         DSS=doy_uns - SDD_uns) %>%
  group_by(Melt,ElevTopo,DSS,variable) %>%
  summarise(rich=sum(flwr_prob_pres,na.rm=TRUE)) %>%
  group_by(Melt,ElevTopo,variable) %>%
  summarise(lag=min(DSS[rich > 0.5]),
            end=max(DSS[rich > 0.5]),
            length=end - lag) %>%
  group_by(Melt,ElevTopo) %>%
  summarise(lag_med=quantile(lag,0.5,na.rm=TRUE),
            lag_lwr=quantile(lag,0.025,na.rm=TRUE),
            lag_upr=quantile(lag,0.975,na.rm=TRUE),
            lag_upr75=quantile(lag,0.75,na.rm=TRUE),
            lag_lwr25=quantile(lag,0.25,na.rm=TRUE),
            end_med=quantile(end,0.5,na.rm=TRUE),
            end_lwr=quantile(end,0.025,na.rm=TRUE),
            end_upr=quantile(end,0.975,na.rm=TRUE),
            end_lwr25=quantile(end,0.25,na.rm=TRUE),
            end_upr75=quantile(end,0.75,na.rm=TRUE),
            length_med=quantile(end,0.5,na.rm=TRUE),
            length_lwr=quantile(end,0.025,na.rm=TRUE),
            length_upr=quantile(end,0.975,na.rm=TRUE),
            length_lwr25=quantile(end,0.25,na.rm=TRUE),
            length_upr75=quantile(end,0.75,na.rm=TRUE))
quants_lag_1570 <- jags_pred_1570 %>%
  select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
         contains("sample")) %>%
  melt(id.vars=c("Melt","doy_uns","ElevTopo",
                 "SDD_uns","Spp_num","Prevalence")) %>%
  mutate(flwr_prob_pres=value*Prevalence,
         DSS=doy_uns - SDD_uns) %>%
  group_by(Melt,ElevTopo,DSS,variable) %>%
  summarise(rich=sum(flwr_prob_pres,na.rm=TRUE)) %>%
  group_by(Melt,ElevTopo,variable) %>%
  summarise(lag=min(DSS[rich > 0.5]),
            end=max(DSS[rich > 0.5]),
            length=end - lag) %>%
  group_by(Melt,ElevTopo) %>%
  summarise(lag_med=quantile(lag,0.5,na.rm=TRUE),
            lag_lwr=quantile(lag,0.025,na.rm=TRUE),
            lag_upr=quantile(lag,0.975,na.rm=TRUE),
            lag_upr75=quantile(lag,0.75,na.rm=TRUE),
            lag_lwr25=quantile(lag,0.25,na.rm=TRUE),
            end_med=quantile(end,0.5,na.rm=TRUE),
            end_lwr=quantile(end,0.025,na.rm=TRUE),
            end_upr=quantile(end,0.975,na.rm=TRUE),
            end_lwr25=quantile(end,0.25,na.rm=TRUE),
            end_upr75=quantile(end,0.75,na.rm=TRUE),
            length_med=quantile(end,0.5,na.rm=TRUE),
            length_lwr=quantile(end,0.025,na.rm=TRUE),
            length_upr=quantile(end,0.975,na.rm=TRUE),
            length_lwr25=quantile(end,0.25,na.rm=TRUE),
            length_upr75=quantile(end,0.75,na.rm=TRUE))
quants_lag_1680 <- jags_pred_1680 %>%
  select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
         contains("sample")) %>%
  melt(id.vars=c("Melt","doy_uns","ElevTopo",
                 "SDD_uns","Spp_num","Prevalence")) %>%
  mutate(flwr_prob_pres=value*Prevalence,
         DSS=doy_uns - SDD_uns) %>%
  group_by(Melt,ElevTopo,DSS,variable) %>%
  summarise(rich=sum(flwr_prob_pres,na.rm=TRUE)) %>%
  group_by(Melt,ElevTopo,variable) %>%
  summarise(lag=min(DSS[rich > 0.5]),
            end=max(DSS[rich > 0.5]),
            length=end - lag) %>%
  group_by(Melt,ElevTopo) %>%
  summarise(lag_med=quantile(lag,0.5,na.rm=TRUE),
            lag_lwr=quantile(lag,0.025,na.rm=TRUE),
            lag_upr=quantile(lag,0.975,na.rm=TRUE),
            lag_upr75=quantile(lag,0.75,na.rm=TRUE),
            lag_lwr25=quantile(lag,0.25,na.rm=TRUE),
            end_med=quantile(end,0.5,na.rm=TRUE),
            end_lwr=quantile(end,0.025,na.rm=TRUE),
            end_upr=quantile(end,0.975,na.rm=TRUE),
            end_lwr25=quantile(end,0.25,na.rm=TRUE),
            end_upr75=quantile(end,0.75,na.rm=TRUE),
            length_med=quantile(end,0.5,na.rm=TRUE),
            length_lwr=quantile(end,0.025,na.rm=TRUE),
            length_upr=quantile(end,0.975,na.rm=TRUE),
            length_lwr25=quantile(end,0.25,na.rm=TRUE),
            length_upr75=quantile(end,0.75,na.rm=TRUE))
quants_lag_1791 <- jags_pred_1791 %>%
  select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
         contains("sample")) %>%
  melt(id.vars=c("Melt","doy_uns","ElevTopo",
                 "SDD_uns","Spp_num","Prevalence")) %>%
  mutate(flwr_prob_pres=value*Prevalence,
         DSS=doy_uns - SDD_uns) %>%
  group_by(Melt,ElevTopo,DSS,variable) %>%
  summarise(rich=sum(flwr_prob_pres,na.rm=TRUE)) %>%
  group_by(Melt,ElevTopo,variable) %>%
  summarise(lag=min(DSS[rich > 0.5]),
            end=max(DSS[rich > 0.5]),
            length=end - lag) %>%
  group_by(Melt,ElevTopo) %>%
  summarise(lag_med=quantile(lag,0.5,na.rm=TRUE),
            lag_lwr=quantile(lag,0.025,na.rm=TRUE),
            lag_upr=quantile(lag,0.975,na.rm=TRUE),
            lag_upr75=quantile(lag,0.75,na.rm=TRUE),
            lag_lwr25=quantile(lag,0.25,na.rm=TRUE),
            end_med=quantile(end,0.5,na.rm=TRUE),
            end_lwr=quantile(end,0.025,na.rm=TRUE),
            end_upr=quantile(end,0.975,na.rm=TRUE),
            end_lwr25=quantile(end,0.25,na.rm=TRUE),
            end_upr75=quantile(end,0.75,na.rm=TRUE),
            length_med=quantile(end,0.5,na.rm=TRUE),
            length_lwr=quantile(end,0.025,na.rm=TRUE),
            length_upr=quantile(end,0.975,na.rm=TRUE),
            length_lwr25=quantile(end,0.25,na.rm=TRUE),
            length_upr75=quantile(end,0.75,na.rm=TRUE))
quants_lag_1901 <- jags_pred_1901 %>%
  select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
         contains("sample")) %>%
  melt(id.vars=c("Melt","doy_uns","ElevTopo",
                 "SDD_uns","Spp_num","Prevalence")) %>%
  mutate(flwr_prob_pres=value*Prevalence,
         DSS=doy_uns - SDD_uns) %>%
  group_by(Melt,ElevTopo,DSS,variable) %>%
  summarise(rich=sum(flwr_prob_pres,na.rm=TRUE)) %>%
  group_by(Melt,ElevTopo,variable) %>%
  summarise(lag=min(DSS[rich > 0.5]),
            end=max(DSS[rich > 0.5]),
            length=end - lag) %>%
  group_by(Melt,ElevTopo) %>%
  summarise(lag_med=quantile(lag,0.5,na.rm=TRUE),
            lag_lwr=quantile(lag,0.025,na.rm=TRUE),
            lag_upr=quantile(lag,0.975,na.rm=TRUE),
            lag_upr75=quantile(lag,0.75,na.rm=TRUE),
            lag_lwr25=quantile(lag,0.25,na.rm=TRUE),
            end_med=quantile(end,0.5,na.rm=TRUE),
            end_lwr=quantile(end,0.025,na.rm=TRUE),
            end_upr=quantile(end,0.975,na.rm=TRUE),
            end_lwr25=quantile(end,0.25,na.rm=TRUE),
            end_upr75=quantile(end,0.75,na.rm=TRUE),
            length_med=quantile(end,0.5,na.rm=TRUE),
            length_lwr=quantile(end,0.025,na.rm=TRUE),
            length_upr=quantile(end,0.975,na.rm=TRUE),
            length_lwr25=quantile(end,0.25,na.rm=TRUE),
            length_upr75=quantile(end,0.75,na.rm=TRUE))

quants_lag <- rbind(quants_lag_1490,
                     quants_lag_1570,
                     quants_lag_1680,
                     quants_lag_1791,
                     quants_lag_1901)

quants_lag$Elev <- str_split_fixed(quants_lag$ElevTopo,pattern="-",n=2)[,1]
quants_lag$Elev <- factor(quants_lag$Elev,levels=c("1901","1791","1680","1570","1490"),
                           labels=c("1901m","1791m","1680m","1570m","1490m"))

quants_lag$Topo <- str_split_fixed(quants_lag$ElevTopo,pattern="-",n=2)[,2]
quants_lag$Topo <- factor(quants_lag$Topo,levels=c("r","s","c"),
                           labels=c("Ridge","Slope","Cove"))

####Computes summary statistics.####

##Overall changes in lag and length.
load("./output/spp_ElevTopo_early_late_samples_attrib.Rdata")
jags_pred_pres <- filter(jags_pred_attr,Prevalence > 0)
DSS_uns <- jags_pred_pres$doy_uns - jags_pred_pres$SDD_uns
jags_pred_dss <- data.frame(DSS_uns,jags_pred_pres)
jags_pred_dss <- filter(jags_pred_dss,DSS_uns >= -10 & DSS_uns <= 90)

jags_attr <- jags_pred_dss[,1:16]
jags_samps <- jags_pred_dss[,17:516]

length_samps <- matrix(NA,ncol=3,nrow=500)
colnames(length_samps) <- c("Early","Typical","Diff")
lag_samps <- matrix(NA,ncol=3,nrow=500)
colnames(lag_samps) <- c("Early","Typical","Diff")

for(i in 1:dim(jags_samps)[2]){
  jags_attr_prob <- data.frame(jags_attr,prob=jags_samps[,i])
  jags_attr_prob$prob_pres <- jags_attr_prob$prob * jags_attr_prob$Prevalence
  jags_attr_grp <- group_by(jags_attr_prob,DSS_uns,Melt,ElevTopo)
  jags_attr_rich <- summarise(jags_attr_grp,rich=sum(prob_pres))
  jags_attr_grp2 <- group_by(jags_attr_rich,Melt,ElevTopo)
  jags_attr_ElevTopo <- summarise(jags_attr_grp2,start=min(DSS_uns[rich > 0.5]),
                                  end=max(DSS_uns[rich > 0.5]),
                                  length=end - start)
  jags_attr_grp3 <- group_by(jags_attr_ElevTopo,Melt)
  jags_attr_melt <- summarise(jags_attr_grp3,avg_lag=mean(start),
                              avg_length=mean(end))
  jags_attr_diff <- data.frame(Melt="Diff",avg_lag=jags_attr_melt$avg_lag[1] - jags_attr_melt$avg_lag[2],
                               avg_length=jags_attr_melt$avg_length[1] - jags_attr_melt$avg_length[2])
  jags_attr_change <- rbind(jags_attr_melt,jags_attr_diff)
  lag_samps[i,] <- jags_attr_change$avg_lag
  length_samps[i,] <- jags_attr_change$avg_length
}

quantfun <- function(x) {quantile(x,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
lag_quants <- apply(lag_samps,FUN=quantfun,MARGIN = 2)
length_quants <- apply(length_samps,FUN=quantfun,MARGIN = 2)

##CIs for low-elevation sites
length_samps_low <- matrix(NA,ncol=3,nrow=500)
colnames(length_samps_low) <- c("Early","Typical","Diff")
lag_samps_low <- matrix(NA,ncol=3,nrow=500)
colnames(lag_samps_low) <- c("Early","Typical","Diff")

for(i in 1:dim(jags_samps)[2]){
  jags_attr_prob <- data.frame(jags_attr,prob=jags_samps[,i])
  jags_attr_prob <- filter(jags_attr_prob,Elev %in% c("1490","1570","1680"))
  jags_attr_prob$prob_pres <- jags_attr_prob$prob * jags_attr_prob$Prevalence
  jags_attr_grp <- group_by(jags_attr_prob,DSS_uns,Melt,ElevTopo)
  jags_attr_rich <- summarise(jags_attr_grp,rich=sum(prob_pres))
  jags_attr_grp2 <- group_by(jags_attr_rich,Melt,ElevTopo)
  jags_attr_ElevTopo <- summarise(jags_attr_grp2,start=min(DSS_uns[rich > 0.5]),
                                  end=max(DSS_uns[rich > 0.5]),
                                  length=end - start)
  jags_attr_grp3 <- group_by(jags_attr_ElevTopo,Melt)
  jags_attr_melt <- summarise(jags_attr_grp3,avg_lag=mean(start),
                              avg_length=mean(end))
  jags_attr_diff <- data.frame(Melt="Diff",avg_lag=jags_attr_melt$avg_lag[1] - jags_attr_melt$avg_lag[2],
                               avg_length=jags_attr_melt$avg_length[1] - jags_attr_melt$avg_length[2])
  jags_attr_change <- rbind(jags_attr_melt,jags_attr_diff)
  lag_samps_low[i,] <- jags_attr_change$avg_lag
  length_samps_low[i,] <- jags_attr_change$avg_length
}

quantfun <- function(x) {quantile(x,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
lag_quants_low <- apply(lag_samps_low,FUN=quantfun,MARGIN = 2)
length_quants_low <- apply(length_samps_low,FUN=quantfun,MARGIN = 2)

##CIs for all topo bands and elevations
lag_samps_topo <- matrix(NA,nrow=500,ncol=15)
colnames(lag_samps_topo) <- levels(jags_attr$ElevTopo)
length_samps_topo <- matrix(NA,nrow=500,ncol=15)
colnames(length_samps_topo) <- levels(jags_attr$ElevTopo)

for(i in 1:dim(jags_samps)[2]){
  jags_attr_prob <- data.frame(jags_attr,prob=jags_samps[,i])
  jags_attr_prob$prob_pres <- jags_attr_prob$prob * jags_attr_prob$Prevalence
  jags_attr_grp <- group_by(jags_attr_prob,DSS_uns,Melt,Elev,Topo,ElevTopo)
  jags_attr_rich <- summarise(jags_attr_grp,rich=sum(prob_pres,na.rm=TRUE))
  jags_attr_grp2 <- group_by(jags_attr_rich,Melt,Elev,Topo,ElevTopo)
  jags_attr_ElevTopo <- summarise(jags_attr_grp2,start=min(DSS_uns[rich > 0.5],na.rm=TRUE),
                                  end=max(DSS_uns[rich > 0.5],na.rm=TRUE),
                                  length=end - start)
  jags_attr_early <- filter(jags_attr_ElevTopo,Melt=="Early")
  jags_attr_late <- filter(jags_attr_ElevTopo,Melt=="Typical")
  lag_samps_topo[i,] <- jags_attr_early$start - jags_attr_late$start
  length_samps_topo[i,] <- jags_attr_early$length - jags_attr_late$length
}

quantfun <- function(x) {quantile(x,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
lag_quants_topo <- apply(lag_samps_topo,FUN=quantfun,MARGIN = 2)
length_quants_topo <- apply(length_samps_topo,FUN=quantfun,MARGIN = 2)
lag_quants_melt <- melt(lag_quants_topo)
colnames(lag_quants_melt) <- c("Quantile","ElevTopo","Value")
lag_quants_melt$Elev <- as.factor(substr(lag_quants_melt$ElevTopo,1,4))
lag_quants_melt$Topo <- as.factor(substr(lag_quants_melt$ElevTopo,6,6))
lag_quants_melt$Topo <- factor(lag_quants_melt$Topo,levels=c("r","s","c"),
                          labels=c("Ridge","Slope","Cove"))
lag_quants_cast <- dcast(lag_quants_melt,formula=Elev+Topo~Quantile,value.var="Value")

length_quants_melt <- melt(length_quants_topo)
colnames(length_quants_melt) <- c("Quantile","ElevTopo","Value")
length_quants_melt$Elev <- as.factor(substr(length_quants_melt$ElevTopo,1,4))
length_quants_melt$Topo <- as.factor(substr(length_quants_melt$ElevTopo,6,6))
length_quants_melt$Topo <- factor(length_quants_melt$Topo,levels=c("r","s","c"),
                               labels=c("Ridge","Slope","Cove"))
length_quants_cast <- dcast(length_quants_melt,formula=Elev+Topo~Quantile,value.var="Value")

##Computes credible intervals for average novelty across the flowering season.
nov_med_samps <- rep(NA,500)


for(i in 1:dim(jags_samps)[2]){
  jags_attr_prob <- data.frame(jags_attr,prob=jags_samps[,i])
  jags_attr_prob$prob_pres <- jags_attr_prob$prob * jags_attr_prob$Prevalence
  jags_attr_early <- 
  jags_attr_grp <- group_by(jags_attr_prob,DSS_uns,Melt,ElevTopo)
  jags_attr_rich <- summarise(jags_attr_grp,nov=bc_nov())
  jags_attr_grp2 <- group_by(jags_attr_rich,Melt,ElevTopo)
  jags_attr_ElevTopo <- summarise(jags_attr_grp2,start=min(DSS_uns[rich > 0.5]),
                                  end=max(DSS_uns[rich > 0.5]),
                                  length=end - start)
  jags_attr_grp3 <- group_by(jags_attr_ElevTopo,Melt)
  jags_attr_melt <- summarise(jags_attr_grp3,avg_lag=mean(start),
                              avg_length=mean(end))
  jags_attr_diff <- data.frame(Melt="Diff",avg_lag=jags_attr_melt$avg_lag[1] - jags_attr_melt$avg_lag[2],
                               avg_length=jags_attr_melt$avg_length[1] - jags_attr_melt$avg_length[2])
  jags_attr_change <- rbind(jags_attr_melt,jags_attr_diff)
  lag_samps[i,] <- jags_attr_change$avg_lag
  length_samps[i,] <- jags_attr_change$avg_length
}

quantfun <- function(x) {quantile(x,probs=c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
lag_quants <- apply(lag_samps,FUN=quantfun,MARGIN = 2)
length_quants <- apply(length_samps,FUN=quantfun,MARGIN = 2)


####Makes plots####

#pdf("./figs/reassembly_bcdist_elev_topo.pdf",width=6,height=8)
ggplot(dist)+
  geom_line(aes(x=DSS,y=value),
            data=filter(dist,quant=="med"))+
  # stat_summary(aes(x=DSS,y=value),geom="text",fun.y=mean,label="Mean",
  #              data=filter(dist,quant=="med"))+
  geom_line(aes(x=DSS,y=value),
            data=filter(dist,quant=="upr"),linetype="dotted")+
  geom_line(aes(x=DSS,y=value),
            data=filter(dist,quant=="lwr"),linetype="dotted")+
  scale_y_continuous("Reassembly Index",breaks=c(0,0.5,1))+
  scale_x_continuous("Days Since Snow Melt",limits=c(0,75))+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()

#pdf("./figs/reassembly_bcdistmin_elev_topo.pdf",width=6,height=8)
ggplot(dist)+
  geom_line(aes(x=DSS,y=value),
            data=filter(dist,quant=="lwr"))+
  scale_y_continuous("Min. Reassembly Index",breaks=c(0,0.5,1))+
  scale_x_continuous("Days Since Snow Melt",limits=c(0,75))+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()

#pdf("./figs/novelty_bcdist_elev_topo.pdf",width=6,height=5)
ggplot(nov)+
  geom_line(aes(x=DSS,y=value),
            data=filter(nov,quant=="med"))+
  # stat_summary(aes(x=DSS,y=value),geom="text",fun.y=mean,label="Mean",
  #              data=filter(dist,quant=="med"))+
  geom_line(aes(x=DSS,y=value),
            data=filter(nov,quant=="upr"),linetype="dotted")+
  geom_line(aes(x=DSS,y=value),
            data=filter(nov,quant=="lwr"),linetype="dotted")+
  scale_y_continuous("Novelty Index",breaks=c(0,0.5,1))+
  scale_x_continuous("Days Since Snow Melt",limits=c(0,75))+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()

#pdf("./figs/seas_median_bcdist_elev_topo.pdf",width=5,height=4)
ggplot(dist_seas_c)+
  geom_linerange(aes(x=Elev,ymin=lwr,ymax=upr,color=Topo),
                 position=position_dodge(width=0.5),lwd=0.8)+
  geom_linerange(aes(x=Elev,ymin=lwr25,ymax=upr75,color=Topo),
                 position=position_dodge(width=0.5),lwd=1.5)+
  geom_point(aes(x=Elev,y=med,color=Topo),shape=21,fill="white",
             position=position_dodge(width=0.5))+
  scale_color_grey("Topo.\nPosition")+
  scale_x_discrete("Elevation (m)")+
  scale_y_continuous("Median Reassembly Index",limits=c(0,1))+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()

#pdf("./figs/seas_median_bcnov_elev_topo.pdf",width=5,height=4)
ggplot(nov_seas_c)+
  geom_linerange(aes(x=Elev,ymin=lwr,ymax=upr,color=Topo),
                 position=position_dodge(width=0.5),lwd=0.8)+
  geom_linerange(aes(x=Elev,ymin=lwr25,ymax=upr75,color=Topo),
                 position=position_dodge(width=0.5),lwd=1.5)+
  geom_point(aes(x=Elev,y=med,color=Topo),shape=21,fill="white",
             position=position_dodge(width=0.5))+
  scale_color_grey("Topo.\nPosition")+
  scale_x_discrete("Elevation (m)")+
  scale_y_continuous("Median Novelty Index",limits=c(0,0.3))+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()

##Median novelty across elevations and topo positions.
#pdf("./figs/bcnov_all_elev.pdf",width=4,height=4)
ggplot(all_median_dss)+
  geom_ribbon(aes(x=DSS,ymin=X25.,ymax=X75.),fill="grey60")+
  geom_line(aes(x=DSS,y=X2.5.),linetype="dotted")+
  geom_line(aes(x=DSS,y=X97.5.),linetype="dotted")+
  geom_line(aes(x=DSS,y=X50.),linetype="solid")+
  scale_y_continuous("Novelty Index",breaks=c(0,0.5,1),limits=c(0,1.01))+
  scale_x_continuous("Days Since Snow Melt",limits=c(0,90))+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()


#pdf("./figs/lag_elev_topo.pdf",width=5,height=10)
ggplot(quants_lag)+
  geom_linerange(aes(x=Melt,ymin=lag_lwr,ymax=lag_upr,color=Melt),
                 position=position_dodge(width=0.5),lwd=0.8)+
  geom_linerange(aes(x=Melt,ymin=lag_lwr25,ymax=lag_upr75,color=Melt),
                 position=position_dodge(width=0.5),lwd=1.5)+
  geom_point(aes(x=Melt,y=lag_med,color=Melt),shape=21,fill="white",
             position=position_dodge(width=0.5))+
  scale_color_grey("Topo.\nPosition")+
  scale_x_discrete("Elevation (m)")+
  scale_y_continuous("Lag (days)")+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()


#pdf("./figs/length_elev_topo.pdf",width=5,height=10)
ggplot(quants_lag)+
  geom_linerange(aes(x=Melt,ymin=length_lwr,ymax=length_upr,color=Melt),
                 position=position_dodge(width=0.5),lwd=0.8)+
  geom_linerange(aes(x=Melt,ymin=length_lwr25,ymax=length_upr75,color=Melt),
                 position=position_dodge(width=0.5),lwd=1.5)+
  geom_point(aes(x=Melt,y=length_med,color=Melt),shape=21,fill="white",
             position=position_dodge(width=0.5))+
  scale_color_grey("Topo.\nPosition")+
  scale_x_discrete("Elevation (m)")+
  scale_y_continuous("Flowering Season Length (days)")+
  facet_grid(facets=Elev~Topo)+
  theme_bw()+
  theme(panel.grid=element_blank())
#dev.off()


p1 <- ggplot(lag_quants_cast)+
  geom_linerange(aes(x=Elev,ymin=lag_quants_cast$'2.5%',ymax=lag_quants_cast$'97.5%',color=Topo),
                 position=position_dodge(width=0.5),lwd=0.8)+
  geom_linerange(aes(x=Elev,ymin=lag_quants_cast$'25%',ymax=lag_quants_cast$'75%',color=Topo),
                 position=position_dodge(width=0.5),lwd=1.5)+
  geom_point(aes(x=Elev,y=lag_quants_cast$'50%',color=Topo),shape=21,fill="white",
             position=position_dodge(width=0.5))+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  scale_color_grey("Topo.\nPosition")+
  scale_x_discrete("Elevation (m)")+
  scale_y_continuous("Change in Lag (Days)",limits=c(-11,28))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        legend.position="none")

p2 <- ggplot(length_quants_cast)+
  geom_linerange(aes(x=Elev,ymin=length_quants_cast$'2.5%',ymax=length_quants_cast$'97.5%',color=Topo),
                 position=position_dodge(width=0.5),lwd=0.8)+
  geom_linerange(aes(x=Elev,ymin=length_quants_cast$'25%',ymax=length_quants_cast$'75%',color=Topo),
                 position=position_dodge(width=0.5),lwd=1.5)+
  geom_point(aes(x=Elev,y=length_quants_cast$'50%',color=Topo),shape=21,fill="white",
             position=position_dodge(width=0.5))+
  geom_hline(aes(yintercept=0),linetype="dotted")+
  scale_color_grey("Topo.\nPosition")+
  scale_x_discrete("Elevation (m)")+
  scale_y_continuous("Change in Length (Days)",limits=c(-11,28))+
  theme_bw()+
  theme(panel.grid=element_blank())

library(gridExtra)
#pdf("./figs/lag_length_elev_topo.pdf",width=8,height=4)
grid.arrange(p1,p2,ncol=2,widths=c(0.9,1.1))
#dev.off()


