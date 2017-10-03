####Functions for the phenology analysis####
####Author: Ian Breckheimer
####October 1st 2017.

##Logit and antilogit function.####
antilogit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}

####Function to fit a Bayesian mixed-effects model in JAGS with species, plot, and year as random effects.####
fit.jags.mixed.moist.gdd <- function(x1,x2,x3,x4,y,groups,species,
                                     nsamples=10000){
  require(rjags)
  
  ##Removes negative and NA values from x, y and group vectors.
  data <- data.frame(x1,x2,x3,x4,y,groups,species) 
  data_complete <- data[complete.cases(data),]
  
  ##Prepares input.
  group <- as.numeric(factor(data_complete$groups))
  ngroups <- length(unique(group))
  species <- as.numeric(factor(data_complete$species))
  nspp <- length(unique(species))
  x1 <- data_complete$x1
  x2 <- data_complete$x2 # want this to be sdd
  x3 <- data_complete$x3
  x4 <- data_complete$x4
  y <- as.numeric(data_complete$y > 0)
  n <- length(x1)
  
  ##Fits model.
  cat(
    "
    model{
    # priors
    #height_int ~ dunif(-20,20)
    #height_slope ~ dunif(-5,5)
    #opt_int ~ dunif(-20,20)
    #opt_slope ~ dunif(-5,5)
    #width_int ~ dunif(-20,20)
    #width_slope ~ dnorm(0,0.1) 
    
    # priors .g=group, .s=species, .y=years
    #height.g.sigma ~ dunif(0.001,3.5)
    #height.g.tau <- pow(height.g.sigma,-2)
    width.g.sigma ~ dunif(0.01,15)
    width.g.tau <- pow(width.g.sigma,-2)
    #opt.g.sigma ~ dunif(0.01,5)
    #opt.g.tau <- pow(opt.g.sigma,-2)
    
    # priors .g=group, .s=species, .y=years
    #height.y.sigma ~ dunif(0,3.5)
    #height.y.tau <- pow(height.y.sigma,-2)
    #width.g.sigma ~ dunif(0.001,15)
    #width.g.tau <- pow(width.g.sigma,-2)
    opt.g.sigma ~ dunif(0.001,5)
    opt.g.tau <- pow(opt.g.sigma,-2)
    
    ##Community mean parameters for snow
    height_slope.s.sdd.mu ~ dnorm(0,0.1)
    height_slope.s.sdd.sigma ~ dunif(0.001,10)
    height_slope.s.sdd.tau <- pow(height_slope.s.sdd.sigma,-2)
    opt_slope.s.sdd.mu ~ dnorm(1,0.1)
    opt_slope.s.sdd.sigma ~ dunif(0.001,5)
    opt_slope.s.sdd.tau <- pow(opt_slope.s.sdd.sigma,-2)
    width_slope.s.sdd.mu ~ dnorm(0,0.01)
    width_slope.s.sdd.sigma ~ dunif(0.001,10)
    width_slope.s.sdd.tau <- pow(width_slope.s.sdd.sigma,-2)
    
    ##Community mean parameters for soil moisture.
    opt_slope.s.moist.mu ~ dnorm(0,0.1)
    opt_slope.s.moist.sigma ~ dunif(0.001,5)
    opt_slope.s.moist.tau <- pow(opt_slope.s.moist.sigma,-2)
    width_slope.s.moist.mu ~ dnorm(0,0.1)
    width_slope.s.moist.sigma ~ dunif(0.001,5)
    width_slope.s.moist.tau <- pow(width_slope.s.moist.sigma,-2)
    
    ##Community mean parameters for growing degree-days.
    opt_slope.s.gdd.mu ~ dnorm(0,0.1)
    opt_slope.s.gdd.sigma ~ dunif(0.001,5)
    opt_slope.s.gdd.tau <- pow(opt_slope.s.gdd.sigma,-2)
    width_slope.s.gdd.mu ~ dnorm(0,0.1)
    width_slope.s.gdd.sigma ~ dunif(0.001,5)
    width_slope.s.gdd.tau <- pow(width_slope.s.gdd.sigma,-2)
    
    ##Species mean intercepts
    height.s.mu ~ dnorm(0,0.01)
    height.s.sigma ~ dunif(0.001,10)
    height.s.tau <- pow(height.s.sigma,-2)
    opt.s.mu ~ dnorm(0,0.01)
    opt.s.sigma ~ dunif(0.001,10)
    opt.s.tau <- pow(opt.s.sigma,-2)
    width.s.mu ~ dnorm(100,0.001)T(0.001,1000)
    width.s.mu.log <- log(width.s.mu)
    width.s.sigma ~ dunif(0.001,400)
    width.s.tau <- pow(width.s.sigma,-2)
    
    ##Group random effects (plot).
    for (j in 1:ngroups){
    opt.g[j] ~ dnorm(0, opt.g.tau)
    width.g[j] ~ dnorm(0,width.g.tau)
    }
    
    ##Year random effects (yr).
    #for (j in 1:nyears){
    #height.y[j] ~ dnorm(0,height.y.tau)
    #opt.g[j] ~ dnorm(0, opt.g.tau)
    #width.g[j] ~ dnorm(0,width.g.tau)
    #}
    
    ##Species random effects on intercept
    for (l in 1:nspp){
    height.s[l] ~ dnorm(height.s.mu,height.s.tau)
    opt.s[l] ~ dnorm(opt.s.mu, opt.s.tau)
    width.s[l] ~ dnorm(width.s.mu,width.s.tau)T(0.001,1000)
    width.s.log[l] <- log(width.s[l])
    }
    
    ##Species random effects on slope
    for (l in 1:nspp){
    height_slope.s.sdd[l] ~ dnorm(height_slope.s.sdd.mu, height_slope.s.sdd.tau)
    opt_slope.s.sdd[l] ~ dnorm(opt_slope.s.sdd.mu, opt_slope.s.sdd.tau)
    width_slope.s.sdd[l] ~ dnorm(width_slope.s.sdd.mu, width_slope.s.sdd.tau)
    opt_slope.s.moist[l] ~ dnorm(opt_slope.s.moist.mu, opt_slope.s.moist.tau)
    width_slope.s.moist[l] ~ dnorm(width_slope.s.moist.mu, width_slope.s.moist.tau)
    opt_slope.s.gdd[l] ~ dnorm(opt_slope.s.gdd.mu, opt_slope.s.gdd.tau)
    width_slope.s.gdd[l] ~ dnorm(width_slope.s.gdd.mu, width_slope.s.gdd.tau)
    }
    
    ##Likelihood
    for (i in 1:n){
    y[i] ~ dbinom(p[i],1)
    
    height[i] <-  height_slope.s.sdd[spp[i]] * x2[i] + height.s[spp[i]]
    opt[i] <- opt.s[spp[i]] + opt_slope.s.sdd[spp[i]] * x2[i] + opt_slope.s.moist[spp[i]] * x3[i] + 
    opt_slope.s.gdd[spp[i]] * x4[i] + opt.g[group[i]]
    width[i] <-  exp(width.s.log[spp[i]] + width_slope.s.sdd[spp[i]] * x2[i] + width_slope.s.moist[spp[i]] * x3[i] +
    width_slope.s.gdd[spp[i]] * x4[i]) * -1
    logit(p[i]) <- width[i]*(x1[i] - opt[i])^2 + height[i]
    }
    }  
    ", file="jagsmodel_log_vertexform_allspp_sppsiteyear_group_moist_gdd.txt"
  )
  
  jd <- list(x1=x1, x2=x2,x3=x3,x4=x4,y=y,n=n,spp=species,nspp=nspp,group=group,ngroups=ngroups)
  mod <- jags.model("jagsmodel_log_vertexform_allspp_sppsiteyear_group_moist_gdd.txt", data= jd, n.chains=3, n.adapt=1000)
  update(mod,n.iter=nsamples)
  out <- coda.samples(mod, c("height.s","height.s.mu","height.s.sigma",
                             "opt.g[5]","opt.g.sigma",
                             "opt.s","opt.s.mu","opt.s.sigma",
                             "width.s","width.s.mu","width.s.sigma",
                             "height_slope.s.sdd","height_slope.s.sdd.mu","height_slope.s.sdd.sigma",
                             "opt_slope.s.sdd","opt_slope.s.sdd.mu","opt_slope.s.sdd.sigma",
                             "width_slope.s.sdd","width_slope.s.sdd.mu","width_slope.s.sdd.sigma",
                             "opt_slope.s.moist","opt_slope.s.moist.mu","opt_slope.s.moist.sigma",
                             "width_slope.s.moist","width_slope.s.moist.mu","width_slope.s.moist.sigma",
                             "opt_slope.s.gdd","opt_slope.s.gdd.mu","opt_slope.s.gdd.sigma",
                             "width_slope.s.gdd","width_slope.s.gdd.mu","width_slope.s.gdd.sigma"),
                      n.iter=nsamples,thin=10)
  #diags <- gelman.diag(out)
  return(list(mod=mod,out=out))
  }

####Function to update a JAGS model.####
update.jags.mixed.moist.gdd <- function(jagsmodel,n.update,n.iter,thin,
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
                                                 "width_slope.s.gdd","width_slope.s.gdd.mu","width_slope.s.gdd.sigma")){
  mod <- jagsmodel$mod
  update(mod,n.iter=n.update)
  out <- coda.samples(mod, params,
                      n.iter=n.iter,thin=thin)
  return(list(mod=mod,out=out))
}



##### Function to find the values that define a given percentage of the area under the phenology curve.####

##Defaults to +- one standard deviation.

obs_intervals <- function(preds,threshold=0.1586553){
  
  # Temporary variables.
  preds <- preds[order(preds$dss),]
  dss <- preds$dss
  pred <- preds$pred
  
  # Gets the bin width from the first pred interval.
  bin_width <- dss[2] - dss[1]
  
  # Assume all NA predictions are zero
  pred[is.na(pred)] <- 0
  
  # Total area under the curve.
  total_area <- sum(pred*bin_width,na.rm=TRUE)
  
  # Computes cumulative proportions in both directions
  cumprop_up <- cumsum(pred)*bin_width/total_area
  cumprop_down <- rev(cumsum(rev(pred))*bin_width)/total_area
  
  # Finds the indices of the first and last values greater than 0.158
  lwr_index <- min(which(cumprop_up >= threshold))
  upr_index <- max(which(cumprop_down >= threshold))
  
  # Finds the corresponding values of dss.
  lwr_bound <- dss[lwr_index]
  upr_bound <- dss[upr_index]
  bounds <- c(lwr_bound,upr_bound)
  names(bounds) <- c("lwr_bound","upr_bound")
  
  # Output
  return(bounds)
}

##Function to parse JAGS output into a matrix of functions representing posterior samples of response curves
##for each species.
make.fun.matrix <- function(jags.out,par.names=c("width.s","opt.s","height.s"),
                            n.samples=100){
  gg.out <- ggmcmc::ggs(jags.out)
  ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  samples <- sample(1:dim(jags.out[[1]])[1],size=n.samples,replace=FALSE)
  chains <- sample(1:3,size=n.samples,replace=TRUE)
  gg.sample <- gg.out[gg.out$Iteration %in% samples & gg.out$Chain == chains,]
  fun.list <- replicate(ngroups*n.samples,function(x){x})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    pname1 <- paste(par.names[1],"[",j,"]",sep="")
    pname2 <- paste(par.names[2],"[",j,"]",sep="")
    pname3 <- paste(par.names[3],"[",j,"]",sep="")
    p1.gr <- gg.sample[gg.sample$Parameter == pname1,4]
    p2.gr <- gg.sample[gg.sample$Parameter == pname2,4]
    p3.gr <- gg.sample[gg.sample$Parameter == pname3,4]
    np <- length(p3.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x - ",np),p2.gr,rep(")^2 + ",np),
                      p3.gr,rep(")",np),sep="")
    for (i in 1:n.samples){
      body(fun.matrix[i,j][[1]]) <- parse(text=fun.body[i])
    }
  }
  return(fun.matrix)
}

##Function to parse JAGS output into a matrix of functions representing posterior samples of response curves
##for each species.
make.med.fun.matrix <- function(jags.out,par.names=c("width.s","opt.s","height.s"),
                                n.samples=100,spp=TRUE){
  gg.out <- ggmcmc::ggs(jags.out)
  if(spp==FALSE){
    ngroups <- 1
  }else{
    ngroups <- sum(grepl(paste(par.names[1],"[",sep=""),unique(gg.out$Parameter),fixed=TRUE))
  }
  fun.list <- replicate(ngroups*n.samples,function(x){x})
  fun.matrix <- matrix(fun.list,ncol=ngroups)
  print("Making function matrix for species ")
  for(j in 1:ngroups){
    print(paste(j,"of ",ngroups))
    if(spp==FALSE){
      pname1 <- par.names[1]
      pname2 <- par.names[2]
      pname3 <- par.names[3]
    }else{
      pname1 <- paste(par.names[1],"[",j,"]",sep="")
      pname2 <- paste(par.names[2],"[",j,"]",sep="")
      pname3 <- paste(par.names[3],"[",j,"]",sep="")
    }
    p1.gr <- quantile(gg.out[gg.out$Parameter == pname1,4],probs=c(0.5))
    p2.gr <- quantile(gg.out[gg.out$Parameter == pname2,4],probs=c(0.5))
    p3.gr <- quantile(gg.out[gg.out$Parameter == pname3,4],probs=c(0.5))
    np <- length(p3.gr)
    fun.body <- paste(rep("antilogit(",np),p1.gr,rep(" * (x - ",np),p2.gr,rep(")^2 + ",np),
                      p3.gr,rep(")",np),sep="")
    body(fun.matrix[1,j][[1]]) <- parse(text=fun.body)
  }
  return(fun.matrix)
}

make.community.fun.matrix <- function(jags.out,par.names=c("width.s.mu","opt.s.mu","height.s.mu"),
                                      n.samples=1000){
  gg.out <- ggmcmc::ggs(jags.out)
  samples <- sample(1:dim(jags.out[[1]])[1],size=n.samples,replace=FALSE)
  chains <- sample(1:3,size=n.samples,replace=TRUE)
  gg.sample <- gg.out[gg.out$Iteration %in% samples & gg.out$Chain == chains,]
  fun.list <- replicate(n.samples,function(x){x})
  fun.matrix <- matrix(fun.list,ncol=1)
  print("Making function matrix for iteration.")
  p1 <- gg.sample[gg.sample$Parameter == par.names[1],4]
  p2 <- gg.sample[gg.sample$Parameter == par.names[2],4]
  p3 <- gg.sample[gg.sample$Parameter == par.names[3],4]
  np <- length(p3)
  fun.body <- paste(rep("antilogit(",np),p1,rep(" * (x - ",np),p2,rep(")^2 + ",np),
                    p3,rep(")",np),sep="")
  for (i in 1:n.samples){
    print(paste(i,"of ",n.samples))
    body(fun.matrix[i,1][[1]]) <- parse(text=fun.body[i])
  }
  return(fun.matrix)
}

##Inverse logit function
inv.logit <- function(x){exp(x)/(1+exp(x))}


##Function to predict the mean value of the response.
jags_fit_fun <- function(height_int,height_slope,
                         opt_int,opt_slope,
                         width_int,width_slope,
                         x1vec,x2vec){
  #Checks inputs.
  stopifnot(length(x1vec)==length(x2vec))
  
  #Calculates the value of the response variable (probability of flowering)
  height <- height_int + height_slope * x2vec
  opt <- opt_int + opt_slope * x2vec
  width<- exp(width_int + width_slope * x2vec) * -1
  mu <- inv.logit(width * (x1vec - opt)^2 + height)
  return(mu)
}

##Function to predict the mean value of the response.
jags_fit_fun_sq <- function(height_int,height_slope,
                            opt_int,opt_slope,opt_slope_sq,
                            width_int,width_slope,width_slope_sq,
                            x1vec,x2vec){
  #Checks inputs.
  stopifnot(length(x1vec)==length(x2vec))
  
  #Calculates the value of the response variable (probability of flowering)
  height <- height_int + height_slope * x2vec
  opt <- opt_int + opt_slope * x2vec + opt_slope_sq * x2vec * x2vec
  width<- exp(width_int + width_slope * x2vec + width_slope_sq * x2vec * x2vec) * -1
  mu <- inv.logit(width * (x1vec - opt)^2 + height)
  return(mu)
}

jags_fit_fun_moist_gdd <- function(height_int,height_slope_sdd,
                                    opt_int,opt_slope_sdd,
                                    opt_slope_moist,opt_slope_gdd,
                                    width_int,width_slope_sdd,
                                    width_slope_moist,width_slope_gdd,
                                    x1vec,x2vec,x3vec,x4vec){
  #Checks inputs.
  stopifnot(length(x1vec)==length(x2vec))
  stopifnot(length(x3vec)==length(x4vec))
  stopifnot(length(x2vec)==length(x4vec))
  
  #Calculates the value of the response variable (probability of flowering)
  height <- height_int + height_slope_sdd * x2vec
  opt <- opt_int + opt_slope_sdd * x2vec + opt_slope_moist * x3vec + 
    opt_slope_gdd * x4vec
  width<- exp(log(width_int) + width_slope_sdd * x2vec + width_slope_moist * x3vec + 
                + width_slope_gdd * x4vec) * -1
  mu <- inv.logit(width * (x1vec - opt)^2 + height)
  return(mu)
}


##Calculates pairwise distances between corresponding rows in two matrices.
bc_dist <- function(m1,m2){
  stopifnot(dim(m1)==dim(m2))
  bc_row <- function(r1,r2){(sum(abs(r1-r2))/sum(r1+r2))}
  dist_vec <- rep(NA,dim(m1)[1])
  for(i in 1:dim(m1)[1]){
    dist_vec[i] <- bc_row(m1[i,],m2[i,])
  }
  return(dist_vec)
}

##Calculates minimum pairwise distances between a given row in m1 and any row in m2.
bc_nov <- function(m1,m2){
  stopifnot(dim(m1)==dim(m2))
  bc_row <- function(r1,r2){(sum(abs(r1-r2))/sum(r1+r2))}
  dist_mat <- matrix(NA,nrow=dim(m1)[1],ncol=dim(m1)[1])
  for(i in 1:dim(m1)[1]){
    for(j in 1:dim(m1)[1]){
      dist_mat[i,j] <- bc_row(m1[i,],m2[j,])
    }
  }
  dist_vec <- apply(dist_mat,MARGIN=1,FUN=min)  
  return(dist_vec)
}


euc_dist <- function(m1,m2){
  stopifnot(dim(m1)==dim(m2))
  euc_row <- function(r1,r2){sqrt(sum(r1-r2)^2)}
  dist_vec <- rep(NA,dim(m1)[1])
  for(i in 1:dim(m1)[1]){
    dist_vec[i] <- euc_row(m1[i,],m2[i,])
  }
  return(dist_vec)
}

create_param_array <- function(out,n.samples,spp.indices,param.names){
  require(ggmcmc)
  require(stringr)
  require(reshape2)
  
  nspp <- length(spp.indices)
  nparams <- length(param.names)
  jags.out <- ggmcmc::ggs(out)
  samples <- sample(unique(jags.out$Iteration),size=n.samples,replace=FALSE)
  chains <- sample(unique(jags.out$Chain),size=n.samples,replace=TRUE)
  gg.sample <- jags.out[jags.out$Iteration %in% samples & jags.out$Chain == chains,]
  gg.sample$Species <- as.numeric(stringr::str_extract(gg.sample$Parameter, "[0-9]+"))
  gg.sample <- gg.sample[!is.na(gg.sample$Species),]
  gg.sample$Param <- stringr::str_split_fixed(gg.sample$Parameter,"\\[[0-9]+\\]+",n=2)[,1]
  gg.sample <- filter(gg.sample,Param %in% param.names)
  gg.sample <- filter(gg.sample,Species %in% spp.indices)
  gg.sample$Param <- factor(gg.sample$Param)
  gg.sample <- gg.sample[,c(1,4,5,6)]
  gg.array <- reshape2::acast(gg.sample,Species~Param~Iteration,value.var="value")
  
  ##Returns a named array
  return(gg.array)
}

####Calculates reassembly index for each site-melt combination.
site_bc_dist <- function(jags_pred){
  require(dplyr)
  require(reshape2)
  quants_spp_early <-  jags_pred %>%
    select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
           contains("sample")) %>%
    filter(Melt=="Early") %>%
    melt(id.vars=c("Melt","doy_uns","ElevTopo",
                   "SDD_uns","Spp_num","Prevalence")) %>%
    mutate(flwr_prob_pres=value*Prevalence,
           DSS=(round(doy_uns) - round(SDD_uns))) %>%
    filter(DSS >= 0 & DSS <= 90) %>%
    select(ElevTopo,Spp_num,DSS,variable,flwr_prob_pres) %>%
    acast(formula=ElevTopo~Spp_num~DSS~variable)
  quants_spp_med <-  jags_pred %>%
    select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
           contains("sample")) %>%
    filter(Melt=="Typical") %>%
    melt(id.vars=c("Melt","doy_uns","ElevTopo",
                   "SDD_uns","Spp_num","Prevalence")) %>%
    mutate(flwr_prob_pres=value*Prevalence,
           DSS=(round(doy_uns) - round(SDD_uns))) %>%
    filter(DSS >= 0 & DSS <= 90) %>%
    select(ElevTopo,Spp_num,DSS,variable,flwr_prob_pres) %>%
    acast(formula=ElevTopo~Spp_num~DSS~variable)
  
  ##Loops through each array and calculates dissimilarity
  dist_mat <- array(NA,dim=c(dim(quants_spp_med)[1],
                             dim(quants_spp_med)[3],
                             dim(quants_spp_med)[4]))
  
  ##Converts all species with less than 1% prob to 0 to avoid artifacts early and late
  quants_spp_early[quants_spp_early<0.001] <- 0
  quants_spp_med[quants_spp_med<0.001] <- 0
  
  for(i in 1:dim(quants_spp_med)[1]){
    print(paste("Processing",i,"of 3 sites"))
    for(j in 1:dim(quants_spp_med)[4]){
      early_r_mat <- t(quants_spp_early[i,,,j])
      med_r_mat <- t(quants_spp_med[i,,,j])
      dist_mat[i,,j] <- bc_dist(early_r_mat,med_r_mat)
    }
  }
  dist_mat[is.na(dist_mat)] <- 0
  dist_mat[is.nan(dist_mat)] <- 0
  quantfun <- function(x){quantile(x,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
  quants <- apply(dist_mat,MARGIN=c(1,2),FUN=quantfun)
  dimnames(quants) <- list(quant=c("lwr","lwr10","lwr25","med","upr75","upr10","upr"),
                           Topo=c("c","r","s"),
                           DSS=seq(0,90,by=1))
  quants_m <- melt(quants)
  quants_seas_avg <- apply(dist_mat,MARGIN=c(1,3),FUN=quantfun)
  quants_seas_quant <- apply(quants_seas_avg,MARGIN=c(1,2),FUN=quantfun)
  quants_seas_med <- quants_seas_quant[,4,]
  quants_seas_dss <- apply(dist_mat,MARGIN=c(2,3),FUN=quantfun)
  dimnames(quants_seas_med) <- list(quant=c("lwr","lwr10","lwr25","med","upr75","upr10","upr"),
                                    Topo=c("c","r","s"))
  out <- list(quants_DSS=quants_m,
              quants_seas=quants_seas_med,
              nov_median_samples=quants_seas_dss)
  return(out)
}


site_bc_nov <- function(jags_pred){
  require(dplyr)
  require(reshape2)
  quants_spp_early <-  jags_pred %>%
    select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
           contains("sample")) %>%
    filter(Melt=="Early") %>%
    melt(id.vars=c("Melt","doy_uns","ElevTopo",
                   "SDD_uns","Spp_num","Prevalence")) %>%
    mutate(flwr_prob_pres=value*Prevalence,
           DSS=(round(doy_uns) - round(SDD_uns))) %>%
    filter(DSS >= 0 & DSS <= 90) %>%
    select(ElevTopo,Spp_num,DSS,variable,flwr_prob_pres) %>%
    acast(formula=ElevTopo~Spp_num~DSS~variable)
  quants_spp_med <-  jags_pred %>%
    select(Melt,doy_uns,ElevTopo,SDD_uns,Spp_num,Prevalence,
           contains("sample")) %>%
    filter(Melt=="Typical") %>%
    melt(id.vars=c("Melt","doy_uns","ElevTopo",
                   "SDD_uns","Spp_num","Prevalence")) %>%
    mutate(flwr_prob_pres=value*Prevalence,
           DSS=(round(doy_uns) - round(SDD_uns))) %>%
    filter(DSS >= 0 & DSS <= 90) %>%
    select(ElevTopo,Spp_num,DSS,variable,flwr_prob_pres) %>%
    acast(formula=ElevTopo~Spp_num~DSS~variable)
  
  ##Loops through each array and calculates dissimilarity
  dist_mat <- array(NA,dim=c(dim(quants_spp_med)[1],
                             dim(quants_spp_med)[3],
                             dim(quants_spp_med)[4]))
  
  ##Converts all species with less than 1% prob to 0 to avoid artifacts early and late
  quants_spp_early[quants_spp_early<0.001] <- 0
  quants_spp_med[quants_spp_med<0.001] <- 0
  
  for(i in 1:dim(quants_spp_med)[1]){
    print(paste("Processing",i,"of 3 sites"))
    for(j in 1:dim(quants_spp_med)[4]){
      early_r_mat <- t(quants_spp_early[i,,,j])
      med_r_mat <- t(quants_spp_med[i,,,j])
      dist_mat[i,,j] <- bc_nov(early_r_mat,med_r_mat)
    }
  }
  dist_mat[is.na(dist_mat)] <- 0
  dist_mat[is.nan(dist_mat)] <- 0
  quantfun <- function(x){quantile(x,c(0.025,0.1,0.25,0.5,0.75,0.9,0.975))}
  quants <- apply(dist_mat,MARGIN=c(1,2),FUN=quantfun)
  dimnames(quants) <- list(quant=c("lwr","lwr10","lwr25","med","upr75","upr10","upr"),
                           Topo=c("c","r","s"),
                           DSS=seq(0,90,by=1))
  quants_m <- melt(quants)
  quants_seas_avg <- apply(dist_mat,MARGIN=c(1,3),FUN=quantfun)
  quants_seas_quant <- apply(quants_seas_avg,MARGIN=c(1,2),FUN=quantfun)
  quants_seas_med <- quants_seas_quant[,4,]
  dimnames(quants_seas_med) <- list(quant=c("lwr","lwr10","lwr25","med","upr75","upr10","upr"),
                                    Topo=c("c","r","s"))
  quant_seas_dss <- apply(dist_mat,MARGIN=c(2,3),FUN=quantfun)[4,,]
  out <- list(quants_DSS=quants_m,
              quants_seas=quants_seas_med,
              nov_median_samples=quant_seas_dss)
  return(out)
}


