# This script was run after the completion of the 2023 EFH project. As a time and storage 
# saving measure for the project, it was useful to discard failed models and not retain their full results.
# However, for the paper, it was necessary to go back and reproduce those results. This script goes through and checks 
# if models were successful in the past, and if so, retrieves the necessay information for the paper. If they were discarded
# in the past, it reruns only those parts that are necessary to the paper. 

rm(list = ls())

library(EFHSDM)
library(raster)
library(maxnet)

EFH.path <- "/Users/kisei.tanaka/Desktop/10211730"

masterplan<-read.csv("/Users/kisei.tanaka/Desktop/10211730/Masterplan.csv")

out.table <- data.frame(
  Region = vector(), 
  Species = vector(), 
  Lifestage = vector(), 
  Model = vector(), 
  N = vector(), 
  Weight = vector(), 
  Scale = vector(), 
  RMSE = vector(), 
  cvRMSE = vector(), 
  Rho = vector(), 
  cvRho = vector(), 
  AUC = vector(), 
  PDE = vector(), 
  cvPDE = vector(), 
  Prob_Area = vector(), 
  Abund_Area = vector(), 
  C_Area = vector()
)

out.table <- read.csv(paste0(EFH.path,"/Trawl_Models2/Metrics_all_models.csv"))

for(r in 1:3){
  
  r = 1
  
  region <- c("AI","EBS","GOA")[r]
  
  # First the setup
  bathy <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Bathy"))
  slope <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Slope"))
  tmax <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Tmax"))
  btemp <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Btemp"))
  btemp <- raster::crop(x = btemp, y = bathy)
  BPI <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/BPI"))
  BPI <- raster::crop(x = BPI, y = bathy)
  Curve <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Curve_Mean"))
  AspectE <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Aspect_East"))
  AspectN <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Aspect_North"))
  
  lat <- raster::init(bathy, v = 'y')
  lat <- raster::mask(lat, bathy, overwrite = FALSE)
  lon <- raster::init(bathy, v = 'x')
  lon <- raster::mask(lon, bathy, overwrite = FALSE)
  coral <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Coralfactor"))
  sponge <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Spongefactor"))
  whips <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/Whipsfactor"))
  
  east <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/ROMSbcurrentEastings"))
  north <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/ROMSbcurrentNorthings"))
  eastSD <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/ROMSbEastingsSD"))
  northSD <- raster::raster(paste0(EFH.path, "/Variables/Variables_", region, "_1km/ROMSbNorthingsSD"))
  
  if (region == "EBS") {
    
    phi <- raster::raster(paste0(EFH.path, "/Variables/Variables_EBS_1km/phi"))
    
    raster.stack <- raster::stack(
      lon,
      lat,
      bathy,
      slope,
      # AspectE,
      # AspectN,
      # Curve,
      btemp,
      # east,
      # north,
      # eastSD,
      # northSD,
      tmax,
      phi,
      BPI,
      sponge,
      coral,
      whips
    )
    
    names(raster.stack) <- c(
      "lon",
      "lat",
      "bdepth",
      "slope",
      # "aspectE",
      # "aspectN",
      # "curve",
      "btemp",
      # "bcurrentU",
      # "bcurrentV",
      # "bcurrentUSD",
      # "bcurrentVSD",
      "tmax",
      "phi",
      "BPI",
      "sponge",
      "coral",
      "pen"
    )
    
  }
  
  
  # GOA and AI don't have sediment grabs to calculate phi, so there is a "rockiness" variable instead
  if(region%in%c("AI","GOA")){
    
    rocky <- raster::raster(paste0(EFH.path,"/Variables/Variables_",region,"_1km/rocky"))
    raster.stack <- raster::stack(lon,lat,bathy,slope,AspectE,AspectN,Curve,btemp,east,north,eastSD,northSD,tmax,rocky,BPI, sponge, coral, whips)
    
    raster.stack <- raster::stack(lon,lat,bathy,slope,btemp,tmax,BPI, sponge, coral, whips)
    
    names(raster.stack) <- c("lon","lat","bdepth","slope","btemp","tmax","BPI","sponge","coral","pen")
    
  }
  
  ak.raster <- raster::raster(paste0(EFH.path,"/Variables/", region, "_Alaska_raster"))
  
  if(region=="GOA"){
    
    GOA.mask<-raster::raster(paste0(EFH.path,"/Variables/GOA_mask_raster"))
    raster.stack<-raster::mask(raster.stack,GOA.mask)
    ak.raster <- NULL
    
  }
  
  # Load the data
  # region.data<-utils::read.csv(paste0(EFH.path,"/Trawl_Models2/",region,"/all_",region,"_data_2021.csv"))
  region.data <- subset(region_data_all, year >= 2012)
  region.data$sponge <- as.integer(region.data$sponge > 0)
  region.data$logarea <- log(region.data$area)
  
  # small adjustment for the REBS complex
  region.data$j_rebs <- region.data$j_rebs + region.data$j_rough + region.data$j_bspot
  region.data$a_rebs <- region.data$a_rebs + region.data$a_rough + region.data$a_bspot
  
  # For now, we will leave the weights of the SFIs in the dataset, but we are using it as a binary variable for the analysis
  region.data$sponge <- as.factor(as.integer(region.data$sponge > 0))
  region.data$coral <- as.factor(as.integer(region.data$coral > 0))
  region.data$pen <- as.factor(as.integer(region.data$pen > 0))
  
  region.data$logarea <- log(region.data$area)
  
  off.raster <- raster::raster(bathy)
  off.raster <- raster::setValues(off.raster, values = mean(region.data$logarea))
  names(off.raster) <- "logarea"
  raster.stack2 <- stack(raster.stack, off.raster)
  
  r.vals <- getValues(raster.stack2)
  r.vals2 <- as.data.frame(r.vals)
  r.vals2$sponge <- as.factor(r.vals2$sponge)
  r.vals2$coral <- as.factor(r.vals2$coral)
  r.vals2$pen <- as.factor(r.vals2$pen)
  
  rm(lat,lon,bathy,slope,tmax,btemp,BPI,AspectN,AspectE,sponge,coral,whips,Curve,phi,rocky,east,eastSD,north,northSD,r.vals)
  
  covars2d <- list(c("lon","lat"),
                   c("bcurrentU","bcurrentV"),
                   c("bcurrentUSD","bcurrentVSD"))
  
  if (region == "EBS") {
    
    covars <- c("bdepth", "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "phi", "BPI")
    maxnet.covars <- c("bcurrentU", "bcurrentV", "bcurrentUSD", "bcurrentVSD", "bdepth",
                       "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "phi", "BPI")
    
  } else {
    
    covars <- c("bdepth", "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "rocky", "BPI")
    maxnet.covars <- c("bcurrentU", "bcurrentV", "bcurrentUSD", "bcurrentVSD", "bdepth",
                       "slope", "aspectE", "aspectN", "curve", "btemp", "tmax", "rocky", "BPI")
    
  }
  
  cofactors <- c("sponge", "coral", "pen")
  
  maxnet2d <- list(c("bcurrentU", "bcurrentV"), 
                   c("bcurrentUSD", "bcurrentVSD"))
  
  # this looks cumbersome, but automates the GAM formulas
  basic.gam.table <- data.frame(
    type = c(rep("smooth", length(covars2d) + length(covars)), rep("factor", length(cofactors)), "offset"),
    dims = c(rep(2, length(covars2d)), rep(1, length(covars) + length(cofactors) + 1)),
    term = c(unlist(lapply(covars2d, FUN = function(x) { return(x[1]) })), covars, cofactors, "logarea"),
    term2 = c(unlist(lapply(covars2d, FUN = function(x) { return(x[2]) })), rep(NA, length(covars) + length(cofactors) + 1)),
    bs = c(rep("ds", length(covars2d)), rep("tp", length(covars)), rep(NA, length(cofactors) + 1)),
    k = c(rep(10, length(covars2d)), rep(4, length(covars)), rep(NA, length(cofactors) + 1)),
    m = c(rep(1, length(covars2d) + length(covars)), rep(NA, length(cofactors) + 1)),
    m2 = c(rep(0.5, length(covars2d)), rep(NA, length(covars) + length(cofactors) + 1))
  )
  
  alt.gam.table <- basic.gam.table
  alt.gam.table$bs[alt.gam.table$bs == "tp"] <- "cr"
  
  #####################################################
  # Now look at each species
  
  region.plan <- subset(masterplan, Region == region)
  
  # region.species <- list.files(paste0("Y:/RACE_EFH_variables/Trawl_Models2/", region))
  # 
  # valid.species <- region.species[region.species %in% region.plan$File_name2]
  # 
  # valid.species <- valid.species[!valid.species %in% c("adult_rougheyerockfish", 
  #                                                      "subadult_rougheyerockfish",
  #                                                      "adult_blackspottedrockfish",
  #                                                      "subadult_blackspottedrockfish")]
  
  valid.species <- c("adult_rougheyerockfish", 
                     "subadult_rougheyerockfish",
                     "adult_blackspottedrockfish",
                     "subadult_blackspottedrockfish")
  
  for(s in 1:length(valid.species)){
    
    s = 1
    
    sp.match <- which(region.plan$File_name2 == valid.species[s])
    
    spec <- region.plan$Species[sp.match]
    lstage <- region.plan$Lifestage[sp.match]
    styear <- region.plan$Start_Year[sp.match]
    abbr <- region.plan$Abbreviation[sp.match]
    
    spec.path <- paste0(EFH.path, "/Trawl_Models2/", region, "/", valid.species[s])
    
    intable <- any(out.table$Region == region & out.table$Lifestage == lstage & out.table$Species == spec)
    
    if(intable == F){
      
      # spec.errors <- read.csv(paste0(spec.path,"/cv_error_data.csv"))
      # e.table <- XML::htmlParse(paste0("file:///Y:/RACE_EFH_variables/Trawl_Models2/",region,"/",valid.species[s],"/ensemble_table.html"))
      # e.table2 <- XML::readHTMLTable(e.table,as.data.frame = F)[[1]]
      # ensemble.table <- as.data.frame(e.table2,stringsAsFactors = F)
      # names(ensemble.table) <- names(e.table2)
      
      ensemble.errors <- subset(spec.errors,Model == "ensemble")
      already.done <- unique(spec.errors$Model)
      
      species.data <- subset(region.data, year >= styear)
      
      # haul.match <- match(species.data$hauljoin,ensemble.errors$hauljoin)
      # species.data$Folds <- ensemble.errors$Folds[haul.match]
      
      n.pres <- sum(species.data[,abbr]>0)
      
      basic.gam.formula <- AssembleGAMFormula(yvar = abbr, gam.table = basic.gam.table)
      alt.gam.formula <- AssembleGAMFormula(yvar = abbr, gam.table = alt.gam.table)
      basic.hgam.formula <- AssembleGAMFormula(yvar = abbr, gam.table = basic.gam.table, hgam = T)
      alt.hgam.formula <- AssembleGAMFormula(yvar = abbr, gam.table = alt.gam.table, hgam = T)
      
      #############################################
      # Maxnet
      
      print(paste("Running MaxEnt model for",region,lstage,spec))
      
      if(any(already.done=="maxnet")){
        
        maxnet.model<-readRDS(paste0(spec.path,"/maxnet_model.rds"))
        maxnet.abund<-raster::raster(paste0(spec.path,"/maxnet_abundance"))
        maxnet.errors<-subset(spec.errors,Model=="maxnet")
        maxnet.scale<-ensemble.table$`Scale Factor`[ensemble.table$Model=="maxnet"]
        maxnet.converge<-T
        
      }else{
        
        # Initialize variables for maxnet modeling
        maxnet_converge <- rep(FALSE, 6)
        maxnet_abund_check <- rep(FALSE, 6)
        maxnet_model_list <- list()
        maxnet_abund_list <- list()
        maxnet_cv_model_list <- list()
        maxnet_error_list <- list()
        maxnet_scale_vec <- rep(NA, 6)
        maxnet_rmse <- rep(NA, 6)
        
        # loop through and check different multiplication constants
        for(r in 1:6){
          
          r = 1
          
          r.mult <- c(.5,1,1.5,2,2.5,3)[r]
          
          # check out a prospective maxnet model
          try(maxnet.model0 <- FitMaxnet(data = species.data,
                                         species = abbr,
                                         vars = maxnet.covars,
                                         facs = cofactors,
                                         regmult = r.mult,
                                         reduce = T))
          
          if(exists("maxnet.model0")){
            
            # maxnet.converge0[r] <- T
            
            maxnet.scale.vec[r] <- mean(species.data[,abbr])/
              mean(exp(stats::predict(maxnet.model0, 
                                      newdata = species.data,
                                      type = "link") + maxnet.model0$ent))
            
            maxnet.abund.list[[r]] <- MakeMaxEntAbundance(model = maxnet.model0,
                                                          maxent.stack = raster.stack,
                                                          scale.fac = maxnet.scale.vec[r],
                                                          type = "cloglog",
                                                          land = ak.raster,
                                                          filename = "")
            
            maxnet.abund.check0[r] <- raster::cellStats(maxnet.abund.list[[r]], max) < (max(species.data[,abbr]) * 10)
            
          }
          
          # if it converges, do the checks and crossvalidation
          if(maxnet.converge0[r]){
            
            maxnet.model.list[[r]]<-maxnet.model0
            
            maxnet.cv <- CrossValidateModel(model = maxnet.model0,
                                            data = species.data,species = abbr,
                                            group = "Folds",
                                            model.type = "maxnet",
                                            key="hauljoin",
                                            scale.preds = T,
                                            regmult = r.mult)
            
            maxnet.error.list[[r]] <- maxnet.cv[[1]]
            
            maxnet.rmse0[r] <- max(c(RMSE(pred = maxnet.error.list[[r]]$pred,obs = maxnet.error.list[[r]]$abund),
                                     RMSE(pred = maxnet.error.list[[r]]$cvpred,obs = maxnet.error.list[[r]]$abund)))
            
            # will also discard a model any of the cv folds fails, as it becomes impossible to calculate RMSE
            if(any(is.na(maxnet.cv[[1]]$cvpred))){
              
              maxnet.model.list[[r]]<-NA
              maxnet.cv.model.list[[r]]<-NA
              maxnet.error.list[[r]]<-NA
              maxnet.rmse0[r]<-NA
              
            }
            
          }else{
            
            maxnet.model.list[[r]]<-NA
            maxnet.cv.model.list[[r]]<-NA
            maxnet.error.list[[r]]<-NA
            maxnet.rmse0[r]<-NA
          }
          
          rm(maxnet.model0)
        }
        
        if(any(maxnet.abund.check0)){
          
          passed.abund<-which(maxnet.abund.check0)
          which.best<-passed.abund[which.min(maxnet.rmse0[passed.abund])]
          
        }else{
          
          which.best<-which.min(maxnet.rmse0)
          
        }
        
        maxnet.model <- maxnet.model.list[[which.best]]
        maxnet.abund <- maxnet.abund.list[[which.best]]
        maxnet.errors <- maxnet.error.list[[which.best]]
        maxnet.scale <- maxnet.scale.vec[which.best]
        maxnet.converge <- T
        
        saveRDS(maxnet.model,file=paste0(spec.path,"/maxnet_model.rds"))
      }
      
      if(maxnet.converge){
        
        # Generate probability raster using the MaxEnt model
        prob.raster <- MakeMaxEntAbundance(model = maxnet.model,
                                           maxent.stack = raster.stack,
                                           land = ak.raster,
                                           type = "maxnet")
        
        # Create a data frame for AUC calculation
        auc.dat <- data.frame(Index = 1:nrow(maxnet.errors),
                              Abundance = maxnet.errors$abund,
                              Probability = maxnet.errors$prob,
                              CV_Probability = maxnet.errors$cvprob)
        
        # Find the breakpoint for EFH using cumulative method
        cbreak <- FindEFHbreaks(abund.raster = maxnet.abund,
                                method = "cumulative",
                                quantiles = 0.05)[2]
        
        # Create a data frame for output data
        out.dat <- data.frame(Region = region,
                              Species = spec,
                              Lifestage = lstage,
                              Model = "maxnet",
                              N = n.pres,
                              Weight = ensemble.table$weight[ensemble.table$Model == "maxnet"],
                              Scale = maxnet.scale,
                              RMSE = RMSE(pred = maxnet.errors$pred, obs = maxnet.errors$abund),
                              cvRMSE = RMSE(pred = maxnet.errors$cvpred, obs = maxnet.errors$abund),
                              Rho = cor(round(maxnet.errors$pred, 2), maxnet.errors$abund, method = "spearman"),
                              cvRho = cor(round(maxnet.errors$cvpred, 2), maxnet.errors$abund, method = "spearman"),
                              AUC = PresenceAbsence::auc(auc.dat, which.model = 1)[1, 1],
                              cvAUC = PresenceAbsence::auc(auc.dat, which.model = 2)[1, 1],
                              PDE = PDE(pred = maxnet.errors$pred, obs = maxnet.errors$abund),
                              cvPDE = PDE(pred = maxnet.errors$cvpred, obs = maxnet.errors$abund),
                              Prob_Area = round(sum(getValues(prob.raster) >= 0.05, na.rm = TRUE), -2),
                              Abund_Area = round(sum(getValues(maxnet.abund) >= 0.0513, na.rm = TRUE), -2),
                              C_Area = round(sum(getValues(maxnet.abund) >= cbreak, na.rm = TRUE), -2))
        
      }else{
        
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="maxnet",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="maxnet"],
                            Scale=NA,RMSE=NA,cvRMSE=NA,Rho=NA,cvRho=NA,AUC=NA,cvAUC=NA,
                            PDE=NA,cvPDE=NA,Prob_Area=NA,Abund_Area=NA,C_Area=NA)
        
      }
      
      out.table <- rbind(out.table, out.dat)
      
      ##########################################################################################
      # Cloglog
      
      print(paste("Running Cloglog model for",region,lstage,spec))
      
      if(any(already.done=="cloglog")){
        
        cloglog.model<-readRDS(paste0(spec.path,"/cloglog_model.rds"))
        cloglog.abund<-raster::raster(paste0(spec.path,"/cloglog_abundance"))
        cloglog.errors<-subset(spec.errors,Model=="cloglog")
        cloglog.scale<-ensemble.table$`Scale Factor`[ensemble.table$Model=="cloglog"]
        cloglog.converge<-T
        
      }else{
        
        try(cloglog.model <- FitGAM(gam.formula = basic.gam.formula,
                                    data = species.data,
                                    family.gam = "binomial",
                                    link.fx = "cloglog",
                                    reduce = T,
                                    select = T,
                                    verbose = F))
        
        cloglog.converge <- exists("cloglog.model") & 
          any(is.infinite(predict(cloglog.model,type = "response"))) == F &
          any(is.na(predict(cloglog.model,type="response"))) == F
        
        if(cloglog.converge){
          
          cloglog.scale <- mean(species.data[,abbr])/mean(exp(predict(cloglog.model, type = "link")))
          
          cloglog.abund <- MakeGAMAbundance(model = cloglog.model,
                                            r.stack = raster.stack,
                                            scale.factor = cloglog.scale,
                                            # land = ak.raster,
                                            filename = "")
          
          cloglog.cv<-CrossValidateModel(model = cloglog.model,model.type = "cloglog",
                                         data = species.data,scale.preds = T,
                                         species = abbr,group = "Folds",key="hauljoin")
          
          cloglog.errors <- cloglog.cv[[1]]
          
          saveRDS(cloglog.model, file = paste0(spec.path, "/cloglog_model.rds"))
          
        }
      }
      if(cloglog.converge){
        prob.raster<-raster::predict(raster.stack2, cloglog.model,factors = list(coral=c(0,1),sponge=c(0,1),pen=c(0,1)),
                                     type = "response", newdata.guaranteed = TRUE)
        
        auc.dat<-data.frame(1:nrow(cloglog.errors),cloglog.errors$abund,cloglog.errors$prob,cloglog.errors$cvprob)
        
        cbreak<-FindEFHbreaks(abund.raster = cloglog.abund,method = "cumulative",quantiles = .05)[2]
        
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="cloglog",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="cloglog"],
                            Scale=cloglog.scale,
                            RMSE=RMSE(pred = cloglog.errors$pred,obs = cloglog.errors$abund),
                            cvRMSE=RMSE(pred = cloglog.errors$cvpred,obs = cloglog.errors$abund),
                            Rho=cor(round(cloglog.errors$pred,2),cloglog.errors$abund,method="spearman"),
                            cvRho=cor(round(cloglog.errors$cvpred,2),cloglog.errors$abund,method="spearman"),
                            AUC=PresenceAbsence::auc(auc.dat,which.model=1)[1,1],
                            cvAUC=PresenceAbsence::auc(auc.dat,which.model=2)[1,1],
                            PDE=PDE(pred = cloglog.errors$pred,obs = cloglog.errors$abund),
                            cvPDE=PDE(pred = cloglog.errors$cvpred,obs = cloglog.errors$abund),
                            Prob_Area=round(sum(getValues(prob.raster)>=.05,na.rm=T),-2),
                            Abund_Area=round(sum(getValues(cloglog.abund)>=.0513,na.rm=T),-2),
                            C_Area=round(sum(getValues(cloglog.abund)>=cbreak,na.rm=T),-2))
      }else{
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="maxnet",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="maxnet"],
                            Scale=NA,RMSE=NA,cvRMSE=NA,Rho=NA,cvRho=NA,AUC=NA,cvAUC=NA,
                            PDE=NA,cvPDE=NA,Prob_Area=NA,Abund_Area=NA,C_Area=NA)
      }
      out.table<-rbind(out.table,out.dat)
      
      ##############################################################################################
      #Hpoisson
      
      print(paste("Running HGAM model for",region,lstage,spec))
      
      # special handling for some species that caused problems
      doit<-T
      if(region=="AI" & lstage=="adult" & spec=="English_sole"){doit<-F}
      if(region=="AI" & lstage=="earlyjuvenile" & spec=="flathead_sole"){doit<-F}
      if(region=="EBS" & lstage=="earlyjuvenile" & spec=="Pacific_cod"){doit<-F}
      if(region=="EBS" & lstage=="subadult" & spec=="whiteblotched_skate"){doit<-F}
      
      if(doit){  
        if(any(already.done=="hpoisson")){
          hpoisson.model<-readRDS(paste0(spec.path,"/hpoisson_model.rds"))
          hpoisson.abund<-raster::raster(paste0(spec.path,"/hpoisson_abundance"))
          hpoisson.errors<-subset(spec.errors,Model=="hpoisson")
          hpoisson.scale<-ensemble.table$`Scale Factor`[ensemble.table$Model=="hpoisson"]
          hpoisson.converge<-T
        }else{
          try(hpoisson.model<-FitHurdleGAM(density.formula = basic.hgam.formula[[1]],prob.formula = basic.hgam.formula[[2]],
                                           data = species.data,verbose = F,select = T,reduce = T))
          
          hpoisson.converge<-exists("hpoisson.model") & any(is.infinite(predict(hpoisson.model,type="response")))==F &
            any(is.na(predict(hpoisson.model,type="response")))==F
          
          if(hpoisson.converge){
            hpoisson.scale<-mean(species.data[,abbr])/mean(predict(hpoisson.model,type="response"))
            hpoisson.abund<-MakeGAMAbundance(model = hpoisson.model,r.stack = raster.stack,scale.factor = hpoisson.scale,
                                             land = ak.raster,filename = "")
            
            hpoisson.cv<-CrossValidateModel(model = hpoisson.model,model.type = "hgam",data = species.data,
                                            species = abbr,group = "Folds",key="hauljoin",scale.preds = T)
            
            hpoisson.errors<-hpoisson.cv[[1]]
            saveRDS(hpoisson.model,file=paste0(spec.path,"/hpoisson_model.rds"))
          }
        }
      }else{
        hpoisson.converge<-F
      }
      if(hpoisson.converge){
        
        hpoisson.prob.vals<-as.data.frame(predict(hpoisson.model,newdata=as.data.frame(r.vals2)))
        
        auc.dat<-data.frame(1:nrow(hpoisson.errors),hpoisson.errors$abund,hpoisson.errors$prob,hpoisson.errors$cvprob)
        
        cbreak<-FindEFHbreaks(abund.raster = hpoisson.abund,method = "cumulative",quantiles = .05)[2]
        
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="hpoisson",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="hpoisson"],
                            Scale=hpoisson.scale,
                            RMSE=RMSE(pred = hpoisson.errors$pred,obs = hpoisson.errors$abund),
                            cvRMSE=RMSE(pred = hpoisson.errors$cvpred,obs = hpoisson.errors$abund),
                            Rho=cor(round(hpoisson.errors$pred,2),hpoisson.errors$abund,method="spearman"),
                            cvRho=cor(round(hpoisson.errors$cvpred,2),hpoisson.errors$abund,method="spearman"),
                            AUC=PresenceAbsence::auc(auc.dat,which.model=1)[1,1],
                            cvAUC=PresenceAbsence::auc(auc.dat,which.model=2)[1,1],
                            PDE=PDE(pred = hpoisson.errors$pred,obs = hpoisson.errors$abund),
                            cvPDE=PDE(pred = hpoisson.errors$cvpred,obs = hpoisson.errors$abund),
                            Prob_Area=round(sum(hpoisson.prob.vals>=.05,na.rm=T),-2),
                            Abund_Area=round(sum(getValues(hpoisson.abund)>=.0513,na.rm=T),-2),
                            C_Area=round(sum(getValues(hpoisson.abund)>=cbreak,na.rm=T),-2))
      }else{
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="hpoisson",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="hpoisson"],
                            Scale=NA,RMSE=NA,cvRMSE=NA,Rho=NA,cvRho=NA,AUC=NA,cvAUC=NA,
                            PDE=NA,cvPDE=NA,Prob_Area=NA,Abund_Area=NA,C_Area=NA)
      }
      out.table<-rbind(out.table,out.dat)
      
      #################################################################################################
      #Poisson
      print(paste("Running Poisson model for",region,lstage,spec))
      if(any(already.done=="poisson")){
        poisson.model<-readRDS(paste0(spec.path,"/poisson_model.rds"))
        poisson.abund<-raster::raster(paste0(spec.path,"/poisson_abundance"))
        poisson.errors<-subset(spec.errors,Model=="poisson")
        poisson.scale<-ensemble.table$`Scale Factor`[ensemble.table$Model=="poisson"]
        poisson.converge<-T
      }else{
        try(poisson.model<-FitGAM(gam.formula = basic.gam.formula,data = species.data,verbose = F,
                                  reduce = T,select = T,family.gam = "poisson"))
        
        poisson.converge<-exists("poisson.model") & any(is.infinite(predict(poisson.model,type="response")))==F &
          any(is.na(predict(poisson.model,type="response")))==F
        
        if(poisson.converge){
          poisson.scale<-mean(species.data[,abbr])/mean(predict(poisson.model,type="response"))
          poisson.abund<-MakeGAMAbundance(model = poisson.model,r.stack = raster.stack,scale.factor = poisson.scale,
                                          land = ak.raster,filename = "")
          
          poisson.cv<-CrossValidateModel(model = poisson.model,model.type = "gam",species = abbr,data = species.data,
                                         group = "Folds",key="hauljoin",scale.preds = T)
          poisson.errors<-poisson.cv[[1]]
          
          saveRDS(poisson.model,file=paste0(spec.path,"/poisson_model.rds"))
        } 
      }
      if(poisson.converge){
        poisson.prob.vals<-1-dpois(0,lambda = getValues(poisson.abund))
        
        auc.dat<-data.frame(1:nrow(poisson.errors),poisson.errors$abund,poisson.errors$prob,poisson.errors$cvprob)
        
        cbreak<-FindEFHbreaks(abund.raster = poisson.abund,method = "cumulative",quantiles = .05)[2]
        
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="poisson",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="poisson"],
                            Scale=poisson.scale,
                            RMSE=RMSE(pred = poisson.errors$pred,obs = poisson.errors$abund),
                            cvRMSE=RMSE(pred = poisson.errors$cvpred,obs = poisson.errors$abund),
                            Rho=cor(round(poisson.errors$pred,2),poisson.errors$abund,method="spearman"),
                            cvRho=cor(round(poisson.errors$cvpred,2),poisson.errors$abund,method="spearman"),
                            AUC=PresenceAbsence::auc(auc.dat,which.model=1)[1,1],
                            cvAUC=PresenceAbsence::auc(auc.dat,which.model=2)[1,1],
                            PDE=PDE(pred = poisson.errors$pred,obs = poisson.errors$abund),
                            cvPDE=PDE(pred = poisson.errors$cvpred,obs = poisson.errors$abund),
                            Prob_Area=round(sum(poisson.prob.vals>=.05,na.rm=T),-2),
                            Abund_Area=round(sum(getValues(poisson.abund)>=.0513,na.rm=T),-2),
                            C_Area=round(sum(getValues(poisson.abund)>=cbreak,na.rm=T),-2))
      }else{
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="poisson",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="poisson"],
                            Scale=NA,RMSE=NA,cvRMSE=NA,Rho=NA,cvRho=NA,AUC=NA,cvAUC=NA,
                            PDE=NA,cvPDE=NA,Prob_Area=NA,Abund_Area=NA,C_Area=NA)
      }
      out.table<-rbind(out.table,out.dat)
      
      #####################################################################
      # Negative Binomial
      print(paste("Running Negative Binomial model for",region,lstage,spec))
      
      if(any(already.done=="negbin")){
        
        negbin.model<-readRDS(paste0(spec.path,"/negbin_model.rds"))
        negbin.abund<-raster::raster(paste0(spec.path,"/negbin_abundance"))
        negbin.errors<-subset(spec.errors,Model=="negbin")
        negbin.scale<-ensemble.table$`Scale Factor`[ensemble.table$Model=="negbin"]
        negbin.converge<-T
        
      }else{
        
        try(negbin.model<-FitGAM(gam.formula = basic.gam.formula,
                                 data = species.data,verbose = F,
                                 reduce = T,select = T,family.gam = "nb"))
        
        negbin.converge <- exists("negbin.model") & 
          any(is.infinite(predict(negbin.model,type="response"))) == F &
          any(is.na(predict(negbin.model,type="response"))) == F
        
        if(negbin.converge){
          
          negbin.scale<-mean(species.data[,abbr])/mean(predict(negbin.model,type="response"))
          negbin.abund<-MakeGAMAbundance(model = negbin.model,r.stack = raster.stack,scale.factor = negbin.scale,
                                         land = ak.raster,filename = "")
          
          negbin.cv<-CrossValidateModel(model = negbin.model,model.type = "gam",species = abbr,data = species.data,
                                        group = "Folds",key="hauljoin",scale.preds = T)
          negbin.errors<-negbin.cv[[1]]
          
          saveRDS(negbin.model,file=paste0(spec.path,"/negbin_model.rds"))
        } 
      }
      if(negbin.converge){
        theta<-as.numeric(strsplit(negbin.model$family[[1]],split="[()]")[[1]][2])
        negbin.prob.vals<-1-dnbinom(x = 0,mu = getValues(negbin.abund),size=theta)
        
        auc.dat<-data.frame(1:nrow(negbin.errors),negbin.errors$abund,negbin.errors$prob,negbin.errors$cvprob)
        
        cbreak<-FindEFHbreaks(abund.raster = negbin.abund,method = "cumulative",quantiles = .05)[2]
        
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="negbin",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="negbin"],
                            Scale=negbin.scale,
                            RMSE=RMSE(pred = negbin.errors$pred,obs = negbin.errors$abund),
                            cvRMSE=RMSE(pred = negbin.errors$cvpred,obs = negbin.errors$abund),
                            Rho=cor(round(negbin.errors$pred,2),negbin.errors$abund,method="spearman"),
                            cvRho=cor(round(negbin.errors$cvpred,2),negbin.errors$abund,method="spearman"),
                            AUC=PresenceAbsence::auc(auc.dat,which.model=1)[1,1],
                            cvAUC=PresenceAbsence::auc(auc.dat,which.model=2)[1,1],
                            PDE=PDE(pred = negbin.errors$pred,obs = negbin.errors$abund),
                            cvPDE=PDE(pred = negbin.errors$cvpred,obs = negbin.errors$abund),
                            Prob_Area=round(sum(negbin.prob.vals>=.05,na.rm=T),-2),
                            Abund_Area=round(sum(getValues(negbin.abund)>=.0513,na.rm=T),-2),
                            C_Area=round(sum(getValues(negbin.abund)>=cbreak,na.rm=T),-2))
      }else{
        out.dat<-data.frame(Region=region,Species=spec,Lifestage=lstage,Model="negbin",N=n.pres,
                            Weight=ensemble.table$weight[ensemble.table$Model=="negbin"],
                            Scale=NA,RMSE=NA,cvRMSE=NA,Rho=NA,cvRho=NA,AUC=NA,cvAUC=NA,
                            PDE=NA,cvPDE=NA,Prob_Area=NA,Abund_Area=NA,C_Area=NA)
      }
      out.table<-rbind(out.table,out.dat)
      
      ##############################################################
      # Ensemble
      
      ensemble.abund <- raster::raster(paste0(spec.path, "/ensemble_abundance"))
      ensemble.prob.vals <- 1 - dpois(0, lambda = getValues(ensemble.abund))
      
      auc.dat <- data.frame(1:nrow(ensemble.errors), ensemble.errors$abund, ensemble.errors$prob, ensemble.errors$cvprob)
      
      cbreak <- FindEFHbreaks(abund.raster = ensemble.abund, method = "cumulative", quantiles = 0.05)[2]
      
      out.dat <- data.frame(Region = region, Species = spec, Lifestage = lstage, Model = "ensemble", N = n.pres,
                            Weight = 1, Scale = 1,
                            RMSE = RMSE(pred = ensemble.errors$pred, obs = ensemble.errors$abund),
                            cvRMSE = RMSE(pred = ensemble.errors$cvpred, obs = ensemble.errors$abund),
                            Rho = cor(round(ensemble.errors$pred, 2), ensemble.errors$abund, method = "spearman"),
                            cvRho = cor(round(ensemble.errors$cvpred, 2), ensemble.errors$abund, method = "spearman"),
                            AUC = PresenceAbsence::auc(auc.dat, which.model = 1)[1, 1],
                            cvAUC = PresenceAbsence::auc(auc.dat, which.model = 2)[1, 1],
                            PDE = PDE(pred = ensemble.errors$pred, obs = ensemble.errors$abund),
                            cvPDE = PDE(pred = ensemble.errors$cvpred, obs = ensemble.errors$abund),
                            Prob_Area = round(sum(ensemble.prob.vals >= 0.05, na.rm = T), -2),
                            Abund_Area = round(sum(getValues(ensemble.abund) >= 0.0513, na.rm = T), -2),
                            C_Area = round(sum(getValues(ensemble.abund) >= cbreak, na.rm = T), -2))
      out.table <- rbind(out.table, out.dat)
      write.csv(out.table, file = paste0(EFH.path, "/Trawl_Models2/Metrics_all_models.csv"), row.names = F)
      print(out.table[nrow(out.table),])
      
      # clear out previous run
      rm(maxnet.model,maxnet.abund,maxnet.errrors,maxnet.scale,maxnet.converge,
         cloglog.model,cloglog.abund,cloglog.errrors,cloglog.scale,cloglog.converge,
         hpoisson.model,hpoisson.abund,hpoisson.errrors,hpoisson.scale,hpoisson.converge,
         poisson.model,poisson.abund,poisson.errrors,poisson.scale,poisson.converge,
         negbin.model,negbin.abund,negbin.errrors,negbin.scale,negbin.converge)
      
    }
  }
}



