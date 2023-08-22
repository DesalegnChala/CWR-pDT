################Lathyrus modelling #########
setwd("C:\\Users\\desalecg\\OneDrive - Universitetet i Oslo\\Documents\\MoDGPI")
# install.packages(c("sdm","gbif"))
devtools::install_github("babaknaimi\\sdm")#if you want the recent version
library(sp)
library(sdm)
library(rgdal)
library(raster)
installAll() ##done only the first time
library(rgbif)
library(dplyr)
name_backbone(name = "Lathyrus", rank = "genus", kingdom = "plante")
options(gbif_user=rstudioapi::askForPassword("my gbif username"))#desalegn
options(gbif_email=rstudioapi::askForPassword("my registred gbif e-mail"))#desdchala@gmail.com
options(gbif_pwd=rstudioapi::askForPassword("my gbif password"))##Password
lathyrus_occ <- occ_download(pred("taxonKey", 10356062), ##taxonkey=10356062
                             pred("hasCoordinate", TRUE),###only these which have coordinates
                             format = "SIMPLE_CSV")##format
####After it finishes, use
LatOcc2 = occ_download_get("0303638-220831081235567")%>% #provided during downlaoding
  occ_download_import()
write.csv(LatOcc2,"C:\\Users\\desalecg\\OneDrive - Universitetet i Oslo\\Documents\\MoDGPI\\occurrence_0303638220831081235567.csv")
summary(LatOcc2$year)#1600    1990    2006    1997    2015    2023 #dim=1269112, 50
LatOcc2 = LatOcc2[LatOcc2$year>1973, ] #filter only the one matches current climate data
dim(LatOcc2) ### dim= 1140027, 50
#######
install.packages("janitor") ### to function "tbyl" #frequency of occurrence points
library(janitor)
freqsp = tabyl(LatOcc2, species)
dim(freqsp)
uniLoc = unique(LatOcc2[ ,c('species', 'decimalLongitude', 'decimalLatitude')]) #only returns three columns
uniLoc2 = unique(LatOcc2[c('decimalLongitude', 'decimalLatitude', 'species'), ])###didn't work, unique rows
spp_to_model = unique(uniLoc$species)
dim(spp_to_model)
class(spp_to_model)
length(spp_to_model)###151
freqUnique = tabyl(uniLoc, species)##no unique occurrences points
barplot(freqUnique$n)
spp_to_model= freqUnique[freqUnique$n>40, ] #87 species after filtering by year
dim(spp_to_model) ####87
row.names(spp_to_model) ######## it only recognizes as vector
head(spp_to_model)
rownames(spp_to_model)=spp_to_model$species
row.names(spp_to_model) ####it is changed
dim(spp_to_model)
######### how to retain the species in species to model
rownames(uniLoc)####
head(uniLoc)
uniLoc2 = uniLoc[uniLoc$species%in%spp_to_model$species, ]
head(uniLoc2)
dim(uniLoc2)
m=tabyl(uniLoc2, species)
dim(m)
sapply(m, frequency)
barplot(m$n)
class(m)
m = as.data.frame(m)
write.csv(uniLoc2, "uniLoc2.csv")
uniLoc2 = read.csv("uniLoc2.csv")
##############
set.seed(42)
Abs_ls <- lapply(unique(uniLoc2$species), FUN = function(x){
  NotSpecies_df <- uniLoc2[uniLoc2$species != x, ]
  SampledAbs <- sample(1:nrow(NotSpecies_df), size = 1e5, replace = FALSE)
  Report_df <- NotSpecies_df[SampledAbs,c("decimalLongitude", "decimalLatitude")]
})
names(Abs_ls) <- unique(uniLoc2$species)
head(Abs_ls,3)
#######presenting names(Abs_ls) as list didn't help
m=l2<-lapply(Abs_ls, function(x) 
  cbind(x, Species = names(Abs_ls)))
absss=cbind(names(Abs_ls$species), Abs_ls)
species = names(Abs_ls)
Abs_ls_2 = cbind(species,Abs_ls)
nn=Abs_ls$`Lathyrus linifolius`
dim(nn)
head(nn)
names(nn)
nrow(Abs_ls[[1]])
Abs_ls$species = names(Abs_ls)
cbind(Abs_ls$species, Abs_ls)
head(mm)
#########
#zz=unlist(Abs_ls) #it doesn't work
y=cbind(Abs_ls,presence)
head(y,3)
uniLoc2presence =y
head(uniLoc2presence)
##################
colnames(Abs_ls[[1]])
#####################
for (i in Abs_ls[[i]]){
 abs2=cbind(Abs_ls[[i]],presence) 
}
####################
###################### aabsence  ##########
Abs_ls_presence=mapply(cbind, Abs_ls, "presence"=0, SIMPLIFY=F)
dim(Abs_ls_presence[[1]])
head(Abs_ls_presence)
class(Abs_ls_presence)
head(Abs_ls_presence,3)
names(Abs_ls_presence)
################### presence ##############
head(uniLoc2)
dim(uniLoc2)
class(uniLoc2)
uniLoc2_presence=cbind(uniLoc2,"presence"=1)
head(uniLoc2_presence)
names(uniLoc2_presence)=unique(uniLoc2$species)
############
uniLoc2_presence = uniLoc2_presence %>%
  group_by(species)%>%
  list()
head(uniLoc2_presence[[1]])
head(uniLoc2_presence)
class(uniLoc2_presence)
################################

########### downloading bioclim #################
#bio = raster:: getData('worldclim', var='bio', res=2.5, lon=12, lat=38)
#bb=getData("SRTM", lon >= -180 & lon <=180)###didn't work
bio2.5 = getData('worldclim', var='bio', res=2.5)
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
###########vifcor removes the one with the highest with vif 0.9 threshold default######
vifcor(bio)
###############but better to extract it at each of the presence points###
ex = raster::extract(bioc,gspg)
head(ex)
v = vifcor(ex, th=0.7)
############ No need to do this ##############
biomod = exclude(bioc,v)
biomod
##########################
d = sdmData(species ~ ., gspg, predictors = bioc, 
            bg = list (n = 100000, method = "gRandom")) ##need to be substituted by the points we have generated
d
##################
######## use getmethodsnames to see the methods######
install.packages("parallel")
library(parallel)
getmethodNames()
m = sdm(species~., d, methods = c ("maxent","gbm","GAM"), 
        replications = c("sub", "boot"), test.p = 25, n = 3,
        parallelSetting = list(ncore = 4, method = "parallel"))
#########################################
###### prediction ########filename is to identify the output####
####### you can run the model on smaller scale and then project on larger one###
p = predict(m, bio, filename = "pr.img", overwrite = TRUE)
#################################################
ev = getEvaluation(m, stat = c("TSS","threshold"))
#########################
############saving #######
save(p,file="p.Rdata")
save(ev,file="ev.Rdata")
#####################################

