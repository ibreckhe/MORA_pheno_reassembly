####Script to fit a Bayesian nonlinear mixed-effects model to phenology data ####

##Analysis supporting "Theobald, Breckheimer, and HilleRisLambers. 2017. Climate drives phenological reassembly 
##of a mountain wildflower meadow community. Ecology."

##Author: Ian Breckheimer
##Data: Elli Theobald
####Updated: 1 October 2017

#### Load data and packages ####
library(rjags)
library(ggmcmc)

data <- read.csv("./data/Phenology_Env_2010to2015_gapfilled.csv", header=TRUE)
table(data$Species)

##Sources phenology functions.
source("./code/pheno_reassembly_functions.R")

#### Prep data ####

# separates data into training and testing subsets.
set.seed(42)
f <- length(data[,1])
random <- runif(f,0,1)
train_data <- data[random<=0.9,]
write.csv(train_data,"./scratch/final_train_data.csv")
test_data <- data[random>0.9,]
head(train_data)

length(train_data$Year)

##Creates a column for days since snow for each observation.
train_data$dss <- train_data$DOY-train_data$SDDDOY

##Unique species in train_dataset.
UniSpecies <- unique(train_data$Species)

## Prepares vectors of scaled and centered train_data for JAGS.
SpFlwr <- train_data$PercFlwr
SpFlwr_jags <- as.numeric(train_data$PercFlwr > 0)

##Scaled and centered days since snow.
Spdss <- scale(train_data$dss)
Spdss_scale <- attr(Spdss,"scaled:scale")
Spdss_center <- attr(Spdss,"scaled:center")
Spdss <- as.numeric(Spdss)

##Scaled and centered snow dissapearance day.
sdd <- (train_data$SDDDOY/100)-1.7

##Scaled and centered DOY # not scaled and centered yet
doy <- (train_data$DOY/100)-1.7

##Scaled and centered soil moisture window.
moist_center <- attributes(scale(train_data$moist_days))$'scaled:center'
moist_scale <- attributes(scale(train_data$moist_days))$'scaled:scale'
moist <- as.numeric(scale(train_data$moist_days))

##Scaled and centered growing-degree days.
gdd_center <- attributes(scale(train_data$gdd50))$'scaled:center'
gdd_scale <- attributes(scale(train_data$gdd50))$'scaled:scale'
gdd <- as.numeric(scale(train_data$gdd50))

##Categorical train_data for possible random effects.
Year <- train_data$Year
length(Year)
SiteYear <- paste(train_data$Site,train_data$Year,sep="-")
SppSiteYear <- paste(train_data$Species,train_data$Site,train_data$Year,sep="-")
SppSite <- paste(train_data$Site,train_data$Species,sep="-")
Site <- train_data$Site
Species <- train_data$Species

##Applies training scaling to test data.
test_data$Spp_num_jags <- as.numeric(factor(test_data$Species))
test_data$gdd_scaled <- (test_data$gdd50 - gdd_center) / gdd_scale
test_data$moist_scaled <- (test_data$moist_days - moist_center) / moist_scale
write.csv(test_data,"./scratch/final_test_data.csv")

##Plots train_data to gut-check.
jagsdat <- data.frame(sdd=train_data$SDDDOY,doy=train_data$DOY,moist=train_data$moist_days,
                      year=train_data$Year,site=train_data$Site,gdd=train_data$gdd50,
                      spp=train_data$Species,flwr=train_data$PercFlwr)
pdf("./figs/theobald_allspp_train_data_gdd_ez.pdf",width=15,height=10)
ggplot(data=jagsdat)+
  geom_point(aes(x=sdd,y=doy,alpha=flwr,color=moist),
             position=position_jitter(width=0.2,height=0),
             shape=19,size=0.1)+
  geom_abline(slope=1,intercept=0,linetype="dotted")+
  scale_color_gradientn("Soil Moisture Window (days)",colors=rainbow(20))+
  scale_alpha_continuous("% Flowering",range=c(0.05,1))+
  xlab("Snow Melt Day of Year")+
  ylab("Day of Year")+
  facet_wrap(facets=~spp)+
  theme_bw()
dev.off()

#### Fits a nonlinear logistic mixed-model to presence-absence data in JAGS #####

##This part of the script uses functions from the file pheno_reassembly_functions.R##

##Fits the model and produces posterior samples (takes a few hours)
Sp_jags <- fit.jags.mixed.moist.gdd(x1=doy,x2=sdd,x3=moist,x4=gdd,y=SpFlwr_jags,species=Species,
                                 groups=SppSiteYear,
                                 nsamples=20000)
Sp_jags_up <- update.jags.mixed.moist.gdd(Sp_jags,n.update=1000,n.iter=10000,thin=10,
                                params=c("height.s","height.s.mu","height.s.sigma",
                                         "opt.g[5]","opt.g.sigma",
                                         "opt.s","opt.s.mu","opt.s.sigma",
                                         "width.s","width.s.mu","width.s.sigma",
                                         "height_slope.s.sdd","height_slope.s.sdd.mu","height_slope.s.sdd.sigma",
                                         "opt_slope.s.sdd","opt_slope.s.sdd.mu","opt_slope.s.sdd.sigma",
                                         "width_slope.s.sdd","width_slope.s.sdd.mu","width_slope.s.sdd.sigma",
                                         "opt_slope.s.moist","opt_slope.s.moist.mu","opt_slope.s.moist.sigma",
                                         "width_slope.s.moist","width_slope.s.moist.mu","width_slope.s.moist.sigma",
                                         "opt_slope.s.gdd","opt_slope.s.gdd.mu","opt_slope.s.gdd.sigma",
                                         "width_slope.s.gdd","width_slope.s.gdd.mu","width_slope.s.gdd.sigma"))

##Collects posterior samples of predicted flower probability for a random set of 1000 observations.
set.seed(50)
perf_obs <- base::sample(1:(nrow(jagsdat)-5000),size=5000,replace=FALSE)
obs_labels <- paste("p[",perf_obs,"]",sep="")

Sp_jags_perf <- update.jags.mixed.moist.gdd(Sp_jags,n.update=1000,n.iter=5000,thin=10,
                                          params=obs_labels)

##Saves model output to disk.
save(Sp_jags_up,file="./output/th_jagsoutput_allspp_sppplotyr_moist_gdd.Rdata")
save(Sp_jags_perf,file="./output/th_jagsoutput_allspp_sppplotyr_moist_gdd_perf.Rdata")
