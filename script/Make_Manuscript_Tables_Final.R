# In this script, we are going to make the various tables for the manuscript,
library(ggplot2)

masterplan<-read.csv("Y:/RACE_EFH_variables/Masterplan.csv")

ai.data<-read.csv("Y:/RACE_EFH_variables/Trawl_Models2/AI/all_AI_data_2021.csv")
ebs.data<-read.csv("E:/Documents/EFH/all_EBS_data_2021.csv")
goa.data<-read.csv("Y:/RACE_EFH_variables/Trawl_Models2/GOA/all_GOA_data_2021.csv")


metrics.table<-read.csv("E:/Documents/EFH/Metrics_all_models.csv")
metrics.table<-subset(metrics.table,Species%in%c("dark_rockfish","yellow_Irish_lord")==F)

bad.rows<-which(metrics.table$RMSE>5000)
metrics.table[bad.rows,7:18]<-NA

maxnet.table<-subset(metrics.table,Model=="maxnet")
cloglog.table<-subset(metrics.table,Model=="cloglog")
hpoisson.table<-subset(metrics.table,Model=="hpoisson")
poisson.table<-subset(metrics.table,Model=="poisson")
negbin.table<-subset(metrics.table,Model=="negbin")
ensemble.table<-subset(metrics.table,Model=="ensemble")


inverts<-c("Tanner_crab","snow_crab","red_king_crab","golden_king_crab","blue_king_crab","giant_octopus")
roundfish<-c("Pacific_cod","walleye_pollock","Atka_mackerel","sablefish")
skates<-c(unique(metrics.table$Species[unlist(lapply(strsplit(metrics.table$Species,split="_"),
                                                     FUN=function(x){any(x=="skate")}))]),"spiny_dogfish")
rockfish<-unique(metrics.table$Species[unlist(lapply(strsplit(metrics.table$Species,split="_"),
                                                     FUN=function(x){any(x=="rockfish")}))])
rockfish<-c(rockfish,"Pacific_ocean_perch","rougheye_blackspotted_complex","shortspine_thornyhead")
flatfish<-unique(metrics.table$Species)[unique(metrics.table$Species)%in%c(inverts,roundfish,skates,rockfish)==F]

################################################################################
# This section makes a table of the metrics by taxonomic group and then by region
# Table 2 in paper

rmse.table<-data.frame(Region=ensemble.table$Region,Species=ensemble.table$Species,Lifestage=ensemble.table$Lifestage,Type=NA,
                       maxnet=ifelse(maxnet.table$cvRMSE>maxnet.table$RMSE,maxnet.table$cvRMSE,maxnet.table$RMSE),
                       cloglog=ifelse(cloglog.table$cvRMSE>cloglog.table$RMSE,cloglog.table$cvRMSE,cloglog.table$RMSE),
                       hpoisson=ifelse(hpoisson.table$cvRMSE>hpoisson.table$RMSE,hpoisson.table$cvRMSE,hpoisson.table$RMSE),
                       poisson=ifelse(poisson.table$cvRMSE>poisson.table$RMSE,poisson.table$cvRMSE,poisson.table$RMSE),
                       negbin=ifelse(negbin.table$cvRMSE>negbin.table$RMSE,negbin.table$cvRMSE,negbin.table$RMSE),
                       ensemble=ensemble.table$cvRMSE)


rmse.table$Type[rmse.table$Species%in%flatfish]<-"flatfish"
rmse.table$Type[rmse.table$Species%in%roundfish]<-"roundfish"
rmse.table$Type[rmse.table$Species%in%rockfish]<-"rockfish"
rmse.table$Type[rmse.table$Species%in%skates]<-"skates"
rmse.table$Type[rmse.table$Species%in%inverts]<-"inverts"

rmse.table2<-na.omit(rmse.table)

type.table2<-data.frame(Type=vector(),Model=vector(),N=vector(),RMSE=vector())

for(i in 1:5){
  species<-list(flatfish,roundfish,rockfish,skates,inverts)[[i]]
  
  species.table<-subset(rmse.table2,Species%in%species)
  
  sp.N<-apply(species.table[,5:10],MARGIN = 2,FUN = function(x){sum(!is.na(x))})
  sp.rmse<-apply(species.table[,5:10],MARGIN = 2,FUN = median,na.rm=T)
  
  type.dat<-data.frame(Type=c("flatfish","roundfish","rockfish","skates","inverts")[i],
                       Model=c("maxnet","cloglog","hpoisson","poisson","negbin","ensemble"),
                       N=sp.N,
                       RMSE=round(sp.rmse,2))
  type.table2<-rbind(type.table2,type.dat)
}
write.csv(type.table2,file="Y:/RACE_EFH_variables/Trawl_Models2/metrics_by_type2.csv",row.names = F)

###############################################################################
# RMSE histogram 
# Figure 3 in paper

rmse.diff1<-data.frame(x=log(ifelse(maxnet.table$RMSE>maxnet.table$cvRMSE,maxnet.table$RMSE,maxnet.table$cvRMSE))-
                         log(ensemble.table$cvRMSE),group="A")
rmse.diff2<-data.frame(x=log(ifelse(cloglog.table$RMSE>cloglog.table$cvRMSE,cloglog.table$RMSE,cloglog.table$cvRMSE))-
                         log(ensemble.table$cvRMSE),group="B")
rmse.diff3<-data.frame(x=log(ifelse(hpoisson.table$RMSE>hpoisson.table$cvRMSE,hpoisson.table$RMSE,hpoisson.table$cvRMSE))-
                         log(ensemble.table$cvRMSE),group="C")
rmse.diff4<-data.frame(x=log(ifelse(poisson.table$RMSE>poisson.table$cvRMSE,poisson.table$RMSE,poisson.table$cvRMSE))-
                         log(ensemble.table$cvRMSE),group="D")
rmse.diff5<-data.frame(x=log(ifelse(negbin.table$RMSE>negbin.table$cvRMSE,negbin.table$RMSE,negbin.table$cvRMSE))-
                         log(ensemble.table$cvRMSE),group="E")

gg.rmse<-rbind(rmse.diff1,rmse.diff2,rmse.diff3,rmse.diff4,rmse.diff5)
gg.rmse$group<-factor(gg.rmse$group,labels = c("MaxEnt","paGAM","hGAM","GAM[P]","GAM[nb]"))

meds<-aggregate(gg.rmse$x,FUN=median,by=list(gg.rmse$group),na.rm=T)[,2]
meds<-data.frame(meds=meds,group=LETTERS[1:5])
meds$group<-factor(meds$group,labels = c("MaxEnt","paGAM","hGAM","GAM[P]","GAM[nb]"))

labs<-data.frame(lab=paste0("median = ",round((meds$meds),2)),group=LETTERS[1:5],x=-1.5,y=110)
labs$group<-factor(labs$group,labels = c("MaxEnt","paGAM","hGAM","GAM[P]","GAM[nb]"))

png("E:/Documents/EFH/RMSE_histogram3.png",width=6.5,height=6,units="in",res=300)
ggplot()+
  geom_histogram(data=gg.rmse,aes(x=x),binwidth = .1)+geom_vline(xintercept=0,linetype=3)+ylab("Frequency")+
  geom_vline(data=meds,aes(xintercept=meds),col=2,linetype=2)+
  geom_label(data=labs,aes(x=x,y=y,label=lab),size=3,hjust=.07)+xlim(c(-1.5,1.5))+
  facet_wrap(~group,labeller=label_parsed)+xlab("log(RMSE of Constituent/RMSE of Ensemble)")+
  theme_bw()
dev.off()

###############################################################################
# this section is just to figure out which species/lifestages are going to get highlighted
ebs.data<-metrics.table[metrics.table$Region=="EBS",]

ebs.areas<-data.frame(Region=vector(),Species=vector(),Lifestage=vector(),maxnet=vector(),cloglog=vector(),
                      hpoisson=vector(),poisson=vector(),negbin=vector(),ensemble=vector())
ebs.sp<-unique(paste(ebs.data$Species,ebs.data$Lifestage))
for(i in 1:length(ebs.sp)){
  spec<-strsplit(ebs.sp[i],split=" ")[[1]][1]
  ls<-strsplit(ebs.sp[i],split=" ")[[1]][2]
  
  species.data<-subset(ebs.data,Species==spec & Lifestage==ls)
  
  species.dat<-data.frame(Region="EBS",Species=spec,Lifestage=ls,
                          maxnet=species.data$Abund_Area[1],cloglog=species.data$Abund_Area[2],
                          hpoisson=species.data$Abund_Area[3],poisson=species.data$Abund_Area[4],
                          negbin=species.data$Abund_Area[5],ensemble=species.data$Abund_Area[6])
  ebs.areas<-rbind(ebs.areas,species.dat)
}

maxnet.diffs<-log(ebs.areas$ensemble/ebs.areas$maxnet)
hist(maxnet.diffs)
ebs.areas[order(maxnet.diffs)[1:5],]    # maybe rex
ebs.areas[order(maxnet.diffs,decreasing = T)[1:5],] # subadult atf looks good

cloglog.diffs<-log(ebs.areas$ensemble/ebs.areas$cloglog)
hist(cloglog.diffs)
ebs.areas[order(cloglog.diffs)[1:5],]    # maybe rex
ebs.areas[order(cloglog.diffs,decreasing = T)[1:5],] # goodmaybe adult ssth

hpoisson.diffs<-log(ebs.areas$ensemble/ebs.areas$hpoisson)
hist(hpoisson.diffs)
ebs.areas[order(hpoisson.diffs)[1:5],]    # maybe adult dover
ebs.areas[order(hpoisson.diffs,decreasing = T)[1:10],] # I think adult pop

poisson.diffs<-log(ebs.areas$ensemble/ebs.areas$poisson)
hist(poisson.diffs)
ebs.areas[order(poisson.diffs)[1:5],]    # maybe adult ssth
ebs.areas[order(poisson.diffs,decreasing = T)[1:5],] # I think adult pop

negbin.diffs<-log(ebs.areas$ensemble/ebs.areas$negbin)
hist(negbin.diffs)
ebs.areas[order(negbin.diffs)[1:5],]    # maybe adult ssth
ebs.areas[order(negbin.diffs,decreasing = T)[1:5],] # I think adult pop

# going to try adult rex, subadult pop, and tanner crab
#################################################################################
# Some of these species/lifestages had individual SDMs that didn't converge, so we'll need to make the maps here
# calling in the necessary rasters
library(raster)
library(EFHSDM)
library(magrittr)
library(akgfmaps)
library(maxnet)
library(patchwork)

source("E:/Documents/EFH/Functions_Maxent.R")
source("E:/Documents/EFH/Functions_GAMModel.R")
source("E:/Documents/EFH/Functions_LoadMap.R")
source("E:/Documents/EFH/Functions_akgfmaps.R")



region.data<-read.csv("E:/Documents/EFH/all_EBS_data_2021.csv")
region.data$logarea<-log(region.data$area)

# (network location is \\akc0ss-n086/SEA_Programs/RACE_EFH_variables)
bathy <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Bathy")
slope <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Slope")
tmax <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Tmax")
btemp <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Btemp")
btemp<-raster::crop(x = btemp,y=bathy)
BPI <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/BPI")
BPI<-raster::crop(x = BPI,y=bathy)
Curve <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Curve_Mean")
AspectE <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Aspect_East")
AspectN <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Aspect_North")

lat <- raster::init(bathy, v ='y')
lat <- raster::mask(lat, mask = bathy,overwrite = F)
lon <- raster::init(bathy, v ='x')
lon <- raster::mask(lon, mask = bathy,overwrite = F)
coral <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Coralfactor")
sponge <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Spongefactor")
whips <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/Whipsfactor")

east<-raster::raster("E:/Documents/EFH/Variables_EBS_1km/ROMSbcurrentEastings")
north<-raster::raster("E:/Documents/EFH/Variables_EBS_1km/ROMSbcurrentNorthings")
eastSD<-raster::raster("E:/Documents/EFH/Variables_EBS_1km/ROMSbEastingsSD")
northSD<-raster::raster("E:/Documents/EFH/Variables_EBS_1km/ROMSbNorthingsSD")
phi <- raster::raster("E:/Documents/EFH/Variables_EBS_1km/phi")

raster.stack <- raster::stack(lon,lat,bathy,slope,AspectE,AspectN,Curve,btemp,east,north,eastSD,northSD,tmax,phi,BPI, sponge, coral, whips)
names(raster.stack) <- c("lon","lat","bdepth","slope","aspectE","aspectN","curve","btemp","bcurrentU","bcurrentV",
                         "bcurrentUSD","bcurrentVSD", "tmax","phi","BPI","sponge","coral","pen")

ak.raster<-raster::raster("Y:/RACE_EFH_variables/Variables/EBS_Alaska_raster")

off.raster<-raster::raster(bathy)
off.raster<-raster::setValues(off.raster,values = mean(region.data$logarea))
names(off.raster)<-"logarea"
raster.stack2<-stack(raster.stack,off.raster)
rm(bathy,slope,tmax,btemp,BPI,Curve,AspectE,AspectN,lat,lon,coral,sponge,whips,east,north,eastSD,northSD,phi,off.raster)
################################################################################
#######################################################
# Figure 5: Area comparison histogram

area.diffs<-log(ensemble.table$C_Area)-log(ensemble.table$Abund_Area)
area.diffs2<-ensemble.table$C_Area/ensemble.table$Abund_Area
hist(area.diffs,breaks=40)
hist(area.diffs,breaks=40,plot=F)
abline(v=median(area.diffs,na.rm=T),col=2)
median(log(ensemble.table$Abund_Area/ensemble.table$C_Area),na.rm=T)

ebs.area.diffs<-area.diffs[ensemble.table$Region=="EBS"]
ebs.area.diffs2<-cbind(ebs.area.diffs,ensemble.table[ensemble.table$Region=="EBS",3:4])

ssth.x<-ebs.area.diffs2[ebs.area.diffs2$Species=="shortspine_thornyhead" & ebs.area.diffs2$Lifestage=="adult",1]
dover.x<-ebs.area.diffs2[ebs.area.diffs2$Species=="Dover_sole" & ebs.area.diffs2$Lifestage=="subadult",1]
poll.x<-ebs.area.diffs2[ebs.area.diffs2$Species=="walleye_pollock" & ebs.area.diffs2$Lifestage=="adult",1]

sum(area.diffs<(-.1))/length(area.diffs)


area.comp.hist<-ggplot()+
  geom_histogram(data=data.frame(x=data.frame(x=area.diffs)),aes(x=x),binwidth=.2)+
  geom_vline(xintercept=0,linetype=2,alpha=.2)+ylab("Frequency")+
  xlab("log(Area Cumulative Method/Area Probability Method)")+theme_bw()+xlim(c(-2,2))+
  geom_vline(xintercept = median(area.diffs,na.rm=T),col=2,linetype=2,linewidth=1.2)+
  geom_segment(data=data.frame(x=c(ssth.x,dover.x,poll.x-.02),xend=c(ssth.x,dover.x,poll.x-.02),
                               yend=c(10,2,32),y=c(30,18,40)),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(angle=20,length = unit(.1,"inches")),linewidth=1)+
  geom_label(data=data.frame(x=-2,y=42,label="median = -0.55"),aes(x=x,y=y,label=label),hjust=.23,vjust=.06)+
  annotate(geom="text",x=ssth.x-.1,y=32,label="adult ssth",size=6)+
  annotate(geom="text",x=dover.x,y=20,label="subadult dover sole",size=6)+
  annotate(geom="text",x=poll.x+.8,y=42,label="adult walleye pollock",size=6)

png(filename = "E:/Documents/EFH/Area_comp_histogram3.png",width = 6.5,height=6,units="in",res=300)
area.comp.hist
dev.off()

###################################################################################################################
# Figure 6
# let's try adult ssth, subadult dover, and adult poll
#ssth
abund.raster<-raster("Y:/RACE_EFH_variables/Trawl_Models2/EBS/adult_shortspinethornyhead/ensemble_abundance")

ssth.per.ao<-abund.raster>.05
ssth.cum.ao<-abund.raster>FindEFHbreaks(abund.raster,method = "cumulative",quantiles = .05)[2]


AO.plot<-function(ao,fig.name){
  
  poly0<-stars::st_as_stars(ao)
  poly<-sf::st_as_sf(poly0,merge = TRUE)
  poly<-poly[poly$layer == 1, ]
  
  ao.plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = survey.sf, fill = "grey95") +
    ggplot2::geom_sf(data = EBS$akland, fill = "grey40") +
    ggplot2::geom_sf(data = EBS$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = EBS$bathymetry, color = "grey60")+
    ggplot2::geom_sf(data = poly,fill="steelblue",alpha=.5,col="black") +
    ggplot2::coord_sf(xlim = EBS$plot.boundary$x, ylim = EBS$plot.boundary$y) +
    ggplot2::scale_x_continuous(name = "Longitude", breaks = EBS$lon.breaks) +
    ggplot2::scale_y_continuous(name = "Latitude", breaks = EBS$lat.breaks) +
    ggplot2::theme_bw() +
    ggplot2:: theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = c(0.2, 0.12),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA))+
    ggplot2::geom_label(data = data.frame(x = -550000, y = 800000, label = fig.name),
                        ggplot2::aes(x = x, y = y, label = label, hjust = 0, vjust = 1), size = 5)
  
  return(ao.plot)
}


ssth.per.map<-AO.plot(ao = ssth.per.ao,fig.name = "Adult Shortspine Thornyhead\nProbability Method")
ssth.cum.map<-AO.plot(ao = ssth.cum.ao,fig.name = "Adult Shortspine Thornyhead\nCumulative Method")

#subadult dover
abund.raster<-raster("Y:/RACE_EFH_variables/Trawl_Models2/EBS/subadult_doversole/ensemble_abundance")

dover.per.ao<-abund.raster>.05
dover.cum.ao<-abund.raster>FindEFHbreaks(abund.raster,method = "cumulative",quantiles = .05)[2]


dover.per.map<-AO.plot(ao = dover.per.ao,fig.name = "Subadult Dover Sole\nProbability Method")
dover.cum.map<-AO.plot(ao = dover.cum.ao,fig.name = "Subadult Dover Sole \nCumulative Method")

#subadult dover
abund.raster<-raster("Y:/RACE_EFH_variables/Trawl_Models2/EBS/adult_walleyepollock/ensemble_abundance")

poll.per.ao<-abund.raster>.05
poll.cum.ao<-abund.raster>FindEFHbreaks(abund.raster,method = "cumulative",quantiles = .05)[2]


poll.per.map<-AO.plot(ao = poll.per.ao,fig.name = "Adult Walleye Pollock\nProbability Method")
poll.cum.map<-AO.plot(ao = poll.cum.ao,fig.name = "Adult Walleye Pollock\nCumulative Method")

png(filename = "Y:/RACE_EFH_variables/Trawl_Models2/Area_comp_panel.png",width = 12,height=18,units="in",res=300)
(ssth.per.map|ssth.cum.map)/
  (dover.per.map|dover.cum.map)/
  (poll.per.map|poll.cum.map)
dev.off()


############################################################################################################
############################################################################################################
############################################################################################################
############################################################################################################
# Figure 4

# adult rex sole

maxnet.model.rex<-readRDS("E:/Documents/EFH/adult_rexsole/maxnet_model.rds")
maxnet.rex<-MakeMaxEntAbundance(model = maxnet.model.rex,maxent.stack = raster.stack,scale.fac = 8,
                                land = NULL,type = "cloglog",filename = "")
cloglog.rex<-raster("E:/Documents/EFH/adult_rexsole/cloglog_abundance")
hpoisson.rex<-raster("E:/Documents/EFH/adult_rexsole/hpoisson_abundance")
poisson.rex<-raster("E:/Documents/EFH/adult_rexsole/poisson_abundance")
negbin.model.rex<-readRDS("E:/Documents/EFH/adult_rexsole/negbin_model.rds")
negbin.rex<-MakeGAMAbundance(model=negbin.model.rex,r.stack=raster.stack,scale.factor = .5547)
ensemble.rex<-raster("E:/Documents/EFH/adult_rexsole/ensemble_abundance")

# subadult pop
maxnet.pop<-raster("E:/Documents/EFH/subadult_pacificoceanperch/maxnet_abundance")
cloglog.pop<-raster("E:/Documents/EFH/subadult_pacificoceanperch/cloglog_abundance")
hpoisson.pop<-raster("E:/Documents/EFH/subadult_pacificoceanperch/hpoisson_abundance")
poisson.model.pop<-readRDS("E:/Documents/EFH/subadult_pacificoceanperch/poisson_model.rds")
poisson.pop<-MakeGAMAbundance(model=poisson.model.pop,r.stack=raster.stack,scale.factor = 1)
negbin.pop<-raster("E:/Documents/EFH/subadult_pacificoceanperch/negbin_abundance")
ensemble.pop<-raster("E:/Documents/EFH/subadult_pacificoceanperch/ensemble_abundance")

# Tanner crab
maxnet.tanner<-raster("E:/Documents/EFH/all_tannercrab/maxnet_abundance")
cloglog.tanner<-raster("E:/Documents/EFH/all_tannercrab/cloglog_abundance")
hpoisson.model.tanner<-readRDS("E:/Documents/EFH/all_tannercrab/hpoisson_model.rds")
hpoisson.tanner<-MakeGAMAbundance(model=hpoisson.model.tanner,r.stack=raster.stack,scale.factor = .0004)
poisson.tanner<-raster("E:/Documents/EFH/all_tannercrab/poisson_abundance")
negbin.tanner<-raster("E:/Documents/EFH/all_tannercrab/negbin_abundance")
ensemble.tanner<-raster("E:/Documents/EFH/all_tannercrab/ensemble_abundance")

map.list<-list()

abund.list<-list(maxnet.rex,cloglog.rex,hpoisson.rex,poisson.rex,negbin.rex,maxnet.pop,cloglog.pop,hpoisson.pop,
                 poisson.pop,negbin.pop,maxnet.tanner,cloglog.tanner,hpoisson.tanner,poisson.tanner,negbin.tanner)
ensemble.abund.list<-list(ensemble.rex,ensemble.pop,ensemble.tanner)


EBS<-akgfmaps::get_base_layers(select.region = "ebs",set.crs = "auto")

survey.sf1<-stars::st_as_stars(is.na(ensemble.rex))
survey.sf2<-sf::st_as_sf(survey.sf1,merge = TRUE)
survey.sf3<-sf::st_cast(survey.sf2,"POLYGON") # cast the polygons to polylines

survey.sf <- sf::st_transform(survey.sf3, sf::st_crs(EBS$akland))[1:(nrow(survey.sf3) - 1), ]
rm(survey.sf1,survey.sf2,survey.sf3)

fig.names<-c("atop('Adult Rex Sole','MaxEnt')","atop('Adult Rex Sole','paGAM')","atop('Adult Rex Sole','hGAM')",
             "atop('Adult Rex Sole','GAM'[P])","atop('Adult Rex Sole','GAM'[nb])",
             "atop('Subadult POP','MaxEnt')","atop('Subadult POP','paGAM')","atop('Subadult POP','hGAM')",
             "atop('Subadult POP','GAM'[P])","atop('Subadult POP','GAM'[nb])",
             "atop('Tanner Crab','MaxEnt')","atop('Tanner Crab','paGAM')","atop('Tanner Crab','hGAM')",
             "atop('Tanner Crab','GAM'[P])","atop('Tanner Crab','GAM'[nb])")

lon.labels<-c("180\u00B0E","175\u00B0E","170\u00B0E","165\u00B0E","160\u00B0E","155\u00B0E")
lat.labels<-c("54\u00B0N","56\u00B0N","58\u00B0N","60\u00B0N","62\u00B0N","64\u00B0N","66\u00B0N")

EBS$graticule$degree_label<-c(lon.labels,lat.labels)

for(i in 1:15){
  e.abund<-ensemble.abund.list[[ceiling(i/5)]]
  e.abund.vals<-getValues(e.abund)
  e.prob.vals<-1-dpois(x = 0,lambda = e.abund.vals)
  e.prob<-setValues(e.abund,values = e.prob.vals)
  
  abund.vals<-getValues(abund.list[[i]])
  prob.vals<-1-dpois(x = 0,lambda = abund.vals)
  prob<-setValues(e.abund,values = prob.vals)
  
  spec<-c("Adult Rex Sole","Subadult POP","Tanner Crab")[ceiling(i/5)]
  mod.name<-rep(c("MaxEnt","paGAM","hGAM","GAM[P]","GAM[nb]"),3)[i]
  fig.name<-fig.names[i]
  
  model.ao<-prob>.05
  ensemble.ao<-e.prob>.05
  
  model.poly0<-stars::st_as_stars(model.ao)
  model.poly<-sf::st_as_sf(model.poly0,merge = TRUE)
  model.poly<-model.poly[model.poly$layer == 1, ]
  
  ensemble.poly0<-stars::st_as_stars(ensemble.ao)
  ensemble.poly<-sf::st_as_sf(ensemble.poly0,merge = TRUE)
  ensemble.poly<-ensemble.poly[ensemble.poly$layer == 1, ]
  
  ao.plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(data = survey.sf, fill = "grey95") +
    ggplot2::geom_sf(data = EBS$akland, fill = "grey40") +
    ggplot2::geom_sf(data = ensemble.poly,col="red",fill=rgb(1,0,0,.2))+
    ggplot2::geom_sf(data = model.poly,fill="white") +
    ggplot2::geom_sf(data = EBS$graticule, color = "grey70", alpha = 0.5) +
    ggplot2::geom_sf(data = EBS$bathymetry, color = "grey60")+
    ggplot2::geom_sf(data = model.poly,fill="steelblue",alpha=.5,col="black") +
    ggplot2::geom_sf(data=ensemble.poly,fill=NA,col="red")+
    ggplot2::scale_x_continuous( breaks = EBS$lon.breaks,labels = function(x) paste0(abs(x), '\u00B0', "E")) +
    ggplot2::scale_y_continuous(breaks=EBS$lat.breaks,labels = function(x) paste0(x, '\u00B0', "N")) +
    ggplot2::coord_sf(xlim = EBS$plot.boundary$x, ylim = EBS$plot.boundary$y-c(40000,0))+
    ggplot2::theme_bw() +
    ggplot2:: theme(
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      panel.background = ggplot2::element_rect(fill = NA, color = "black"),
      legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
      legend.position = c(0.2, 0.12),
      axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
      plot.background = ggplot2::element_rect(fill = NA, color = NA))
  
  if(i%in%c(1,6,11)){
    
    ao.plot <- ggplot2::ggplot() +
      ggplot2::geom_sf(data = survey.sf, fill = "grey95") +
      ggplot2::geom_sf(data = EBS$akland, fill = "grey40") +
      ggplot2::geom_sf(data = ensemble.poly,col="red",fill=rgb(1,0,0,.2))+
      ggplot2::geom_sf(data = model.poly,fill="white") +
      ggplot2::geom_sf(data = EBS$graticule, color = "grey70", alpha = 0.5) +
      ggplot2::geom_sf(data = EBS$bathymetry, color = "grey60")+
      ggplot2::geom_sf(data = model.poly,fill="steelblue",alpha=.5,col="black") +
      ggplot2::geom_sf(data=ensemble.poly,fill=NA,col="red")+
      ggplot2::scale_x_continuous( breaks = EBS$lon.breaks,labels = function(x) paste0(abs(x), '\u00B0', "E")) +
      ggplot2::scale_y_continuous(breaks=EBS$lat.breaks,labels = function(x) paste0(x, '\u00B0', "N")) +
      ggplot2::coord_sf(xlim = EBS$plot.boundary$x, ylim = EBS$plot.boundary$y-c(40000,0))+
      ggplot2::theme_bw() + ggtitle(c("Adult Rex Sole","Subadult POP","Tanner Crab")[ceiling(i/5)])+
      ggplot2:: theme(
        panel.border = ggplot2::element_rect(color = "black", fill = NA),
        panel.background = ggplot2::element_rect(fill = NA, color = "black"),
        legend.key = ggplot2::element_rect(fill = NA, color = "grey30"),
        legend.position = c(0.2, 0.12),
        axis.title = ggplot2::element_blank(), axis.text = ggplot2::element_text(size = 12),
        legend.text = ggplot2::element_text(size = 12), legend.title = ggplot2::element_text(size = 12),
        plot.background = ggplot2::element_rect(fill = NA, color = NA),plot.title = element_text(hjust = 0.5,size=18))
    
  }
  map.list[[i]]<-ao.plot
}
## now make the histograms

maxnet.metrics<-subset(metrics.table,Model=="maxnet")
cloglog.metrics<-subset(metrics.table,Model=="cloglog")
hpoisson.metrics<-subset(metrics.table,Model=="hpoisson")
poisson.metrics<-subset(metrics.table,Model=="poisson")
negbin.metrics<-subset(metrics.table,Model=="negbin")
ensemble.metrics<-subset(metrics.table,Model=="ensemble")

median(maxnet.metrics$Abund_Area/ensemble.metrics$Abund_Area,na.rm=T)
median(cloglog.metrics$Abund_Area/ensemble.metrics$Abund_Area,na.rm=T)
median(hpoisson.metrics$Abund_Area/ensemble.metrics$Abund_Area,na.rm=T)
median(poisson.metrics$Abund_Area/ensemble.metrics$Abund_Area,na.rm=T)
median(log(negbin.metrics$Abund_Area/ensemble.metrics$Abund_Area),na.rm=T)



rex.spot<-which(ensemble.metrics$Region=="EBS" & ensemble.metrics$Species=="rex_sole" & ensemble.metrics$Lifestage=="adult")
pop.spot<-which(ensemble.metrics$Region=="EBS" & ensemble.metrics$Species=="Pacific_ocean_perch" & ensemble.metrics$Lifestage=="subadult")
tan.spot<-which(ensemble.metrics$Region=="EBS" & ensemble.metrics$Species=="Tanner_crab" & ensemble.metrics$Lifestage=="all")

# maxnet
rex.x<-log(maxnet.metrics$Abund_Area/ensemble.metrics$Abund_Area)[rex.spot]
pop.x<-log(maxnet.metrics$Abund_Area/ensemble.metrics$Abund_Area)[pop.spot]
tan.x<-log(maxnet.metrics$Abund_Area/ensemble.metrics$Abund_Area)[tan.spot]

maxnet.hist<-ggplot()+
  geom_histogram(data=data.frame(x=log(maxnet.metrics$Abund_Area/ensemble.metrics$Abund_Area)),aes(x=x),binwidth=.1)+
  geom_vline(xintercept=0,linetype=2)+ylab("Frequency")+
  xlab("Area relative to ensemble (log ratio)")+xlim(-2,2)+ylim(0,120)+theme_bw()+
  ylab("MaxEnt")+theme(axis.title.y=element_text(size=18))+
  geom_segment(data=data.frame(x=c(rex.x+.03,pop.x+.05,tan.x+.05),xend=c(rex.x+.03,pop.x+.05,tan.x+.05),
                               yend=c(4,14,33),y=c(30,55,85)),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(angle=20,length = unit(.1,"inches")),size=1)+
  annotate(geom="text",x=rex.x+.02,y=36,label="rex",size=7)+
  annotate(geom="text",x=pop.x+.05,y=60,label="pop",size=7,hjust=0)+
  annotate(geom="text",x=tan.x+.05,y=90,label="tan",size=7,hjust=0)+
  annotate(geom="text",x=-1.4,y=75,label="Ensemble\nlarger",size=6)+
  annotate(geom="text",x=1.4,y=75,label="Constituent\nlarger",size=6)


# cloglog

rex.x<-log(cloglog.metrics$Abund_Area/ensemble.metrics$Abund_Area)[rex.spot]
pop.x<-log(cloglog.metrics$Abund_Area/ensemble.metrics$Abund_Area)[pop.spot]
tan.x<-log(cloglog.metrics$Abund_Area/ensemble.metrics$Abund_Area)[tan.spot]

cloglog.hist<-ggplot()+
  geom_histogram(data=data.frame(x=log(cloglog.metrics$Abund_Area/ensemble.metrics$Abund_Area)),aes(x=x),binwidth=.1)+
  geom_vline(xintercept=0,linetype=2)+ylab("Frequency")+
  xlab("Area relative to ensemble (log ratio)")+xlim(-2,2)+ylim(0,120)+theme_bw()+
  ylab("paGAM")+theme(axis.title.y=element_text(size=18))+
  geom_segment(data=data.frame(x=c(rex.x+.05,pop.x+.04,tan.x-.6),xend=c(rex.x+.05,pop.x+.04,tan.x-.04),
                               yend=c(7,49,110),y=c(35,85,110)),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(angle=20,length = unit(.1,"inches")),size=1)+
  annotate(geom="text",x=rex.x+.05,y=40,label="rex",size=7,hjust=0)+
  annotate(geom="text",x=pop.x+.05,y=90,label="pop",size=7,hjust=0)+
  annotate(geom="text",x=tan.x-.55,y=116,label="tan",size=7,hjust=.5)

# hpoisson
rex.x<-log(hpoisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)[rex.spot]
pop.x<-log(hpoisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)[pop.spot]
tan.x<-log(hpoisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)[tan.spot]

hpoisson.hist<-ggplot()+
  geom_histogram(data=data.frame(x=log(hpoisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)),aes(x=x),binwidth=.1)+
  geom_vline(xintercept=0,linetype=2)+ylab("Frequency")+
  xlab("Area relative to ensemble (log ratio)")+xlim(-4,2)+ylim(0,120)+theme_bw()+
  ylab("hGAM")+theme(axis.title.y=element_text(size=18))+
  geom_segment(data=data.frame(x=c(rex.x+.01,pop.x,tan.x),xend=c(rex.x+.01,pop.x,tan.x),
                               yend=c(40,4,3),y=c(60,28,45)),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(angle=20,length = unit(.1,"inches")),size=1)+
  annotate(geom="text",x=rex.x-.02,y=66,label="rex",size=7,hjust=.9)+
  annotate(geom="text",x=pop.x,y=35,label="pop",size=7,hjust=.1)+
  annotate(geom="text",x=tan.x,y=51,label="tan",size=7,hjust=.1)

# poisson
rex.x<-log(poisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)[rex.spot]
pop.x<-log(poisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)[pop.spot]
tan.x<-log(poisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)[tan.spot]

poisson.hist<-ggplot()+
  geom_histogram(data=data.frame(x=log(poisson.metrics$Abund_Area/ensemble.metrics$Abund_Area)),aes(x=x),binwidth=.1)+
  geom_vline(xintercept=0,linetype=2)+ylab("Frequency")+
  xlab("Area relative to ensemble (log ratio)")+xlim(-2,2)+ylim(0,120)+theme_bw()+
  ylab(expression("GAM"["P"]))+theme(axis.title.y=element_text(size=18))+
  geom_segment(data=data.frame(x=c(rex.x-.04,pop.x+.03),xend=c(rex.x-.04,pop.x+.03),
                               yend=c(44,6),y=c(75,35)),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(angle=20,length = unit(.1,"inches")),size=1)+
  annotate(geom="text",x=rex.x-.15,y=78,label="rex\ntan",size=7,hjust=.9)+
  annotate(geom="text",x=pop.x-.02,y=41,label="pop",size=7,hjust=.25)

# negbin
rex.x<-log(negbin.metrics$Abund_Area/ensemble.metrics$Abund_Area)[rex.spot]
pop.x<-log(negbin.metrics$Abund_Area/ensemble.metrics$Abund_Area)[pop.spot]
tan.x<-log(negbin.metrics$Abund_Area/ensemble.metrics$Abund_Area)[tan.spot]

negbin.hist<-ggplot()+
  geom_histogram(data=data.frame(x=log(negbin.metrics$Abund_Area/ensemble.metrics$Abund_Area)),aes(x=x),binwidth=.1)+
  geom_vline(xintercept=0,linetype=2)+ylab("Frequency")+
  xlab("Area relative to ensemble (log ratio)")+xlim(-2,2)+ylim(0,120)+theme_bw()+
  ylab(expression("GAM"["nb"]))+theme(axis.title.y=element_text(size=18))+
  geom_segment(data=data.frame(x=c(rex.x,pop.x),xend=c(rex.x,pop.x),
                               yend=c(34,4),y=c(65,35)),
               aes(x=x,xend=xend,y=y,yend=yend),arrow=arrow(angle=20,length = unit(.1,"inches")),size=1)+
  annotate(geom="text",x=rex.x-.05,y=78,label="rex\ntan",size=7,hjust=.5)+
  annotate(geom="text",x=pop.x,y=41,label="pop",size=7)


png(filename = "E:/Documents/EFH/Hist_maps_panel4.png",width = 16,height=20,units="in",res=600)
(maxnet.hist|map.list[[1]]|map.list[[6]]|map.list[[11]])/
  (cloglog.hist|map.list[[2]]|map.list[[7]]|map.list[[12]])/
  (hpoisson.hist|map.list[[3]]|map.list[[8]]|map.list[[13]])/
  (poisson.hist|map.list[[4]]|map.list[[9]]|map.list[[14]])/
  (negbin.hist|map.list[[5]]|map.list[[10]]|map.list[[15]])
dev.off()

median(log(maxnet.metrics$Abund_Area/ensemble.metrics$Abund_Area),na.rm=T)
median(log(cloglog.metrics$Abund_Area/ensemble.metrics$Abund_Area),na.rm=T)
median(log(hpoisson.metrics$Abund_Area/ensemble.metrics$Abund_Area),na.rm=T)
median(log(poisson.metrics$Abund_Area/ensemble.metrics$Abund_Area),na.rm=T)
median(log(negbin.metrics$Abund_Area/ensemble.metrics$Abund_Area),na.rm=T)

######################################################################################


