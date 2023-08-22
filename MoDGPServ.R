################Lathyrus modelling #########
setwd("C:\\Users\\desalecg\\OneDrive - Universitetet i Oslo\\Documents\\MoDGPI")
# install.packages
install.packages(c("sdm","raster","maptools","rgdal","rgbif","dplyr"))
devtools::install_github("babaknaimi\\sdm")#if you want the recent version
library(sp)
library(sdm)
library(rgdal)
library(raster)
installAll() ##done only the first time - no need to run it every time
library(rgbif)
library(dplyr)
#########################
getmethodNames() ####to check the lists of available sdm algorithms
##############################
name_backbone(name = "Lathyrus", rank = "genus", kingdom = "plante") # to obtain taxon key
occ_count(taxonKey=10356062) # number of occurrences at genus level #1.542.715
occ_search(scientificName = "Lathyrus") #to obtain metadata
Lathyrus$meta #
####################### to sign into gbif ## for tracking DOI for dowlaoded data
options(gbif_user=rstudioapi::askForPassword("my gbif username"))#desalegn
options(gbif_email=rstudioapi::askForPassword("my registred gbif e-mail"))#desdchala@gmail.com
options(gbif_pwd=rstudioapi::askForPassword("my gbif password"))
######################## downloading ---gives guide next steps  --only those with coordinates
lathyrus_occ <- occ_download(pred("taxonKey", 10356062), 
                             pred("hasCoordinate", TRUE),
                             format = "SIMPLE_CSV")
############################################ to download 
##########this provides a download code of of 0303638-220831081235567 ###########
LatOcc2 = occ_download_get("0303638-220831081235567")%>%
  occ_download_import()
##########exporting the data and local saving
write.csv(LatOcc2,"C:\\Users\\desalecg\\Dropbox (UiO)\\BioDT\\2023\\MoDGPI\\occurrence\\0303638220831081235567.csv")
head(LatOcc2)
nrow(LatOcc2)
### to make summary by group
install.packages("janitor")
library(janitor)
###########################
freqsp = tabyl(LatOcc2, species)
dim(freqsp)###162 species
#########################avoid duplicates
Uniqloc = LatOcc2[ ,c("species","decimalLatitude", "decimalLongitude")]
dim(Uniqloc)
Uniqloc2 = unique(Uniqloc)
dim(Uniqloc2)
##############
# export as uniLoc2
################ unique locations - after duplicates are removed
uniLoc2 = read.csv("uniLoc2.csv")
############ select 10,000 points where other species occur 
###but the model target is lacking
set.seed(42)
Abs_ls <- lapply(unique(uniLoc2$species), FUN = function(x){
  NotSpecies_df <- uniLoc2[uniLoc2$species != x, ]
  SampledAbs <- sample(1:nrow(NotSpecies_df), size = 1e5, replace = FALSE)
  Report_df <- NotSpecies_df[SampledAbs,c("decimalLongitude", "decimalLatitude")]
})
names(Abs_ls) <- unique(uniLoc2$species)
head(Abs_ls,3)
#######
####################
###################### aabsence  ##########
Abs_ls_presence=mapply(cbind, Abs_ls, "presence"=0, SIMPLIFY=F)
dim(Abs_ls_presence[[1]])
head(Abs_ls_presence)
class(Abs_ls_presence)
names(Abs_ls_presence)
LsativusA = which(names(Abs_ls_presence)=="Lathyrus sativus")
dim(LsativusA)
head(LsativusA,3)
LsativusA=Abs_ls_presence[[27]]
head(LsativusA,3)
dim(LsativusA)
LsativusA=cbind(species="Lathyrus sativus",LsativusA)
names(LsativusA)=c("species","lon","lat","presence")
################### presence ##############
head(uniLoc2)
dim(uniLoc2)
class(uniLoc2)
uniLoc2_presence=cbind(uniLoc2,"presence"=1)
head(uniLoc2_presence)
ncol(uniLoc2_presence)
uniLoc2_presence=uniLoc2_presence[ ,-1]
names(uniLoc2_presence)=c("species","decimalLongitude","decimalLatitude","presence")
LsativusP=uniLoc2_presence[uniLoc2_presence$species=="Lathyrus sativus", ]
head(LsativusP,3)
##############
LsativusPA = rbind(LsativusP,LsativusA)
############ convering to list
uniLoc2_presencel = uniLoc2_presence %>%
  group_by(species)  %>%
  list()
head(uniLoc2_presence[[1]])
head(uniLoc2_presence)
class(uniLoc2_presence)
################################
########### downloading bioclim #################
#bio = raster:: getData('worldclim', var='bio', res=2.5, lon=12, lat=38)
#bb=getData("SRTM", lon >= -180 & lon <=180)###didn't work
#bio2.5 = getData('worldclim', var='bio', res=2.5)
bio2.5=list.files(path="C:\\Users\\desalecg\\OneDrive - Universitetet i Oslo\\Documents\\wc2-5\\bio_2-5m_bil", pattern="bil", full.names=TRUE)
bio2.5=raster::stack(bio2.5)
plot(bio2.5[[1]])
points(uniLoc2[uniLoc2$species=="Lathyrus sativus", c(3,4)])
for (i in row.names(spp_to_model)){
  points(uniLoc2[uniLoc2$species==i, c(3,4)])
} ##########it works
#####In the case of res=0.5, you must also provide a lon and lat argument#####
###for a tile; for the lower resolutions global data will be downloaded.######
head(uniLoc2)
uniLoc2=uniLoc2[ ,-1]
names(uniLoc2)=c("species", "lon","lat")
head(uniLoc2,3)
presence=1
uniLoc2=cbind(uniLoc2,presence)
head(uniLoc2,3)
#################
########### multicollinearity ###########
##########to check vif from the usdam package#######
install.packages("usdm")
library(usdm)
vif(bio2.5)
#############################
###########vifcor removes the one with the highest with vif (0.9 threshold default)######
#vifcor(bio)
###############but better to extract it at each of the presence points###
LsativusP=LsativusP[ ,c(2:4)] ####presence points ###1165
names(LsativusP)=c("lon","lat","species")
coordinates(LsativusP)=c("decimalLongitude","decimalLatitude")
coordinates(LsativusA)=c("decimalLongitude","decimalLatitude")
bg=LsativusA
ex = raster::extract(bio2.5,LsativusPA[, c(1,2)])
head(ex)
v = vifcor(bio2.5, th=0.7) #variable inflation factor
############ No need to do this ##############
corTree = function(dataframe){
  corMat = cor(dataframe, use="pairwise.complete.obs", method="pearson")
  disVals = abs(as.dist(corMat))
  cluVals = hclust(1 - disVals)
  plot(cluVals)
  abline(h=0.3, col="red")
}

corTree(ex)
####################
######now exclude those with high cor and vif #######
biomod = exclude(bio2.5,v)
biomod ##############
############ compare it with the tree ########
names(biomod)
######################
pred=ex[ ,names(biomod)]
head(pred,3)
pred=cbind(LsativusPA,pred)
preda=raster::extract(biomod,LsativusA[ ,c(2,3)])
preda=cbind(LsativusA,preda)
##########################
bg=LsativusA
plot(bio2.5[[1]])
points(bg, col="red")
##########################
Abs_ls_presence[[27]]$species=names(Abs_ls_presence[27])
head(Abs_ls_presence[[27]],3)
Abs_ls_presence[[27]] = Abs_ls_presence[[27]][, -4]
##################
length(Abs_ls_presence) #87
names(Abs_ls_presence)
###########################
for (i in names(Abs_ls_presence)){
  Abs_ls_presence[[i]]$species=names(Abs_ls_presence[i])
}
################################
head(uniLoc2_presence,3)
class(uniLoc2_presence)
nrow(uniLoc2_presence)
################converting Abs_ls_absence to data frame
Abs_ls_absence=rbindlist(rbind,Abs_ls_presence)
do.call(rbind, Abs_ls_presence)
df=as.data.frame(do.call(cbind, Abs_ls_presence))
table(df$`Lathyrus linifolius.decimalLongitude`)
table(df$presence)
table(names(df$presence))
########################
uniLoc2_presence=uniLoc2_presence[ ,c(2,3,4,1)]
################list
uniLoc2_presenceL= split.data.frame(uniLoc2_presence, uniLoc2_presence$species)
names(uniLoc2_presenceL[[27]])
df2=as.data.frame(do.call(cbind, uniLoc2_presenceL))
# pre-allocate a list and fill it with a loop
xy.list <- vector("list", nrow(xy.df))
for (i in 1:nrow(xy.df)) {
  xy.list[[i]] <- xy.df[i,]
}
#########################

######################## it considers presence as species
data=rbind(pred,preda)
d <- sdmData(~.+coords(decimalLongitude+decimalLatitude),
             train=pred)  
d
d = sdmData( presence ~ .,
             LsativusPPr, 
            #predictors = biomod,
            bg = LsativusAPr) ##need to be substituted by the points we have generated
            #bg = list (n = 100000, method = "gRandom"))
d
d@species.names #"presence"
d@groups
########## it considers as presence only
dd= sdmData(species ~ . + coords(lon +lat),
             LsativusPA, 
             predictors = biomod)
dd
dd@species.names #"Lathyrus sativus"
dd@factors
##################
LsativusPA = as.data.frame(LsativusPA)
#############
for(i in nrow(LsativusPA)){
        mm=sample(i, 5000, replace = FALSE)
  }
nn=LsativusPA[mm, ]
summary(nn$presence)
sum(nn$presence==1)
sum(nn$presence==0)
nn%in%LsativusPA
for(i in nrow(LsativusPA)){
  for(j in ncol(LsativusPA)){
    if (nn%in%LsativusPA)
      rownames(LsativusPA)="trial"
    }
  }
######################
present=uniLoc2_presence
absent=Abs_ls_presence[-59]
head(present)
head(absent)
###############
#absentdf <-  data.frame(matrix(unlist(absent), nrow=length(absent), byrow=TRUE))
presentdf <- present
library (plyr)
absentdf <- ldply (absent, data.frame)
##################################
head(absentdf,3)
absentdf=absentdf[ ,-1]
head(absentdf)
head(present)
Data=rbind(present,absentdf)
head(Data)
tail(Data)
Data <- na.omit(Data)


Trial <- reshape(data = Data,
                 idvar = c("decimalLongitude", "decimalLatitude"),
                 timevar = "species",
                 direction = "wide"
)
View(head(Trial))



summary(Indometh) # data in long format

## long to wide (direction = "wide") requires idvar and timevar at a minimum
reshape(Indometh, direction = "wide", idvar = "Subject", timevar = "time")
##################
######## use getmethodsnames to see the methods######
install.packages("parallel")
library(parallel)
getmethodNames() ###it mentions the methods that are available
m = sdm(~., d, methods = c ("maxent","gbm","GAM","RF"), 
        replications = c("sub", "boot"), test.p = 25, n = 3,
        parallelSetting = list(ncore = 4, method = "parallel"))
#########################################
m@models$species$brt$`6` #######to check different 
gui(m)
###############
save(d, file=dlsat.RData)
save(m,file=mlsat.RData)
###### prediction ########filename is to identify the output####
####### you can run the model on smaller scale and then project on larger one###
#p = predict(m, biomod, filename = "Lsat.img", overwrite = TRUE)
p1 <- predict(m, biomod)
#set extent and projection
extent(p1) <- extent(biomod)
projection(p1) <- projection(biomod)
#write the predicted raster to a file
writeRaster(p1, filename="Lsat.tif", format="GTiff", overwrite=TRUE)

#################################################
en = ensemble(m,bio,filename = "enpr.img", setting = list(method = "weighted",stat="tss", opt =2), overwrite = TRUE)
ev = getEvaluation(m, stat = c("TSS","threshold"))
#########################
############saving #######
save(p,file="p.Rdata")
save(ev,file="ev.Rdata")
####################################################
################### Soil ############################
###############SOIL:https://soilgrids.org/
install.packages('XML')###
install.packages('rgdal')####
install.packages('gdalUtils')
install.packages('sf')
install.packages('dplyr')
install_github("envirometrix/landmap")
install.packages("leaflet")
install.packages("mapview")
library(devtools)
devtools::install_github("https://github.com/be-marc/leaflet.opacity", dependencies=TRUE)
#################libraries ################
library(XML)
library(rgdal)
library(gdalUtils)
library(raster)
library(sf)
library(dplyr)
library(RColorBrewer)
library(leaflet.opacity)
library(leaflet)
library(leaflet.opacity)
library(mapview)
########################################
pH = raster("https://files.isric.org/soilgrids/latest/data/phh2o/phh2o_15-30cm_mean.vrt")
pH
plot(pH)
#####################
nitro = raster("https://files.isric.org/soilgrids/latest/data/nitrogen/nitrogen_15-30cm_mean.vrt")
nitro
plot(nitro)
############
projection(pH) <- projection(bio2.5)
extent(pH) <- extent(bio2.5)
resolution(pH)=resolution(bio2.5)#didn't work
