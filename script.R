## Load Required Packages ##
library(maptools)
library(spatstat) 
library(rgdal)

## Read Dataset ##

data.Ext <- read.csv('./dataset.csv')
data.Ext <- data.Ext[order(data.Ext$ID_SITE),]
data.Core <- subset(data.Ext, AreaNoMarriedExcl != 'n') #generate a data.frame without unmarried sons

data.CoreDiwan <- subset(data.Core, ID_SITE!='ALO_07/1' & ID_SITE!='SUG_07/1' & ID_SITE!='TIB_07/1') #remove sites without diwan
data.CoreDiwan <- subset(data.CoreDiwan, DiwanVSDomestic != 'n') # remove 'n' in the DiwanVSDomestic column

##################
## Sample Sizes ##
##################

## Number of villages:
length(unique(data.Ext$ID_SITE)) #n=11
length(unique(data.Core$ID_SITE)) #n=11
length(unique(data.CoreDiwan$ID_SITE)) #n=8
Nsites=11
NsitesDiwan=8

## Total number of features:
nrow(data.Ext) #n=191
table(data.Ext$category) #n_d=141, n_l=50
nrow(data.Core) #n=171
table(data.Core$category) #n_d=121, n_l=49
table(data.Core$simbol) #n_diwan=8
table(data.Ext$simbol) #n_diwan=8
nrow(data.CoreDiwan) #n=94, of which 8 diwans

# convert to sp object
coordinates(data.Ext) <- c("long","lat")
proj4string(data.Ext) <- CRS("+proj=longlat +datum=WGS84")
data.Ext <- spTransform(data.Ext,CRS("+proj=utm +zone=32 +datum=WGS84"))

coordinates(data.Core) <- c("long","lat")
proj4string(data.Core) <- CRS("+proj=longlat +datum=WGS84")
data.Core <- spTransform(data.Core,CRS("+proj=utm +zone=32 +datum=WGS84"))

coordinates(data.CoreDiwan) <- c("long","lat")
proj4string(data.CoreDiwan) <- CRS("+proj=longlat +datum=WGS84")
data.CoreDiwan <- spTransform(data.CoreDiwan,CRS("+proj=utm +zone=32 +datum=WGS84"))



## Create a list of sites
CoreList <- vector("list",length=Nsites)
ExtList <- vector("list",length=Nsites)
CoreDiwanList <- vector("list",length=NsitesDiwan)


for (x in 1:Nsites)
{
    CoreList[[x]]<-subset(data.Core, ID_SITE==unique(data.Core$ID_SITE)[x])
    ExtList[[x]]<-subset(data.Ext, ID_SITE==unique(data.Ext$ID_SITE)[x])
}

for (x in 1:NsitesDiwan)
{
    CoreDiwanList[[x]]<-subset(data.CoreDiwan, ID_SITE==unique(data.CoreDiwan$ID_SITE)[x])
}




names(CoreList)=unique(data.Core$ID_SITE)
names(ExtList)=unique(data.Ext$ID_SITE)



## Export to ESRI Shapefile
for (i in 1:Nsites)
{
    writeOGR(CoreList[[i]],dsn='./shapefiles',layer=paste('coresite', i, sep='_'), driver='ESRI Shapefile')
    writeOGR(ExtList[[i]],dsn='./shapefiles',layer=paste('extsite', i, sep='_'), driver='ESRI Shapefile')
}


######################################################
## Evaluate Window Size using Ripas-Ripley Estimate ##
######################################################

windowListCore <- vector("list",length=length(unique(data.Core$ID_SITE)))
windowListExt <- vector("list",length=length(unique(data.Ext$ID_SITE)))
names(windowListCore)=unique(data.Core$ID_SITE)
names(windowListExt)=unique(data.Ext$ID_SITE)


for (x in 1:Nsites)
{
    windowListCore[[x]]<-ripras(coordinates(CoreList[[x]])[,1],coordinates(CoreList[[x]])[,2])
    windowListExt[[x]]<-ripras(coordinates(ExtList[[x]])[,1],coordinates(ExtList[[x]])[,2])
}


SiteAreas=data.frame(sitename=unique(data.Core$ID_SITE),
    CoreArea=unlist(lapply(windowListCore,area)),
    ExtArea=unlist(lapply(windowListExt,area)))
write.csv(SiteAreas,row.names=FALSE,"./tables/siteareas.csv")

## Export windows to ESRI Shapefile
for (i in 1:Nsites)
{
    auxCore <- as(windowListCore[[i]], "SpatialPolygons")
    auxExt <- as(windowListExt[[i]], "SpatialPolygons")
    auxCore <- SpatialPolygonsDataFrame(auxCore, data.frame(ID=i))
    auxExt <- SpatialPolygonsDataFrame(auxExt, data.frame(ID=i))
    writeOGR(auxCore, dsn = "./shapefiles", layer=paste('window_Core', i, sep='_'), driver = "ESRI Shapefile")  
    writeOGR(auxExt, dsn = "./shapefiles", layer=paste('window_Ext', i, sep='_'), driver = "ESRI Shapefile")  
}


###################################################################
## Assess Segregation/Aggregation Domestic vs Livestock Features ##
###################################################################


## Remove All columns except 'category' ##
CoreList.ppp = lapply(CoreList,function(x){x@data=data.frame(category=x@data$category);return(x)})


for (x in 1:Nsites)
{
    CoreList.ppp[[x]]<-as(CoreList.ppp[[x]],"ppp") #convert to ppp class object
    CoreList.ppp[[x]]$window=windowListCore[[x]]  #assign window
    CoreList.ppp[[x]]=rjitter(CoreList.ppp[[x]],radius=0.5) #jitter locations for matching coordinates
}

## Create a union of windows and point data
mergedWindow=Reduce(union.owin, windowListCore)
AggregatedPointPattern=unique.ppp(Reduce(function(...) superimpose(...,W= mergedWindow),CoreList.ppp))

## Define number of simulations 
nsim=9999

### Permutation Approach 1 (Complete Random Label) ###

result1=envelope(AggregatedPointPattern,Lcross,i="D",j="L",r=seq(0,100,1),correction="iso",nsim=nsim,simulate=expression(rlabel(AggregatedPointPattern)))


### Permutation Approach 2 (Stratified Random Label)###

simulatedPoints<-vector("list",length=nsim)
for (s in 1:nsim)
{
    tmpList=lapply(CoreList.ppp,rlabel)	
    simulatedPoints[[s]]=unique.ppp(Reduce(function(...) superimpose(...,W= mergedWindow), tmpList))	
}

result2=envelope(AggregatedPointPattern,Lcross,i="D",j="L",r=seq(0,100,1),correction="iso",nsim=nsim, simulate=simulatedPoints)

### Plot Results ###

par(mfrow=c(1,2))
plot(result1$r,result1$hi,type="n",main="Standard Permutation",xlab="Distance (in meters)",ylab="L function")
segregation=result1$r[which(result1$obs<result1$lo)]
rect(xleft=min(segregation),xright=max(segregation),ybottom=-1000,ytop=1000,col=rgb(1,0,0,0.1),border=NA)
polygon(x=c(result1$r,rev(result1$r)),y=c(result1$hi,rev(result1$lo)),col="lightgrey",border=NA)
lines(result1$r,result1$mmean,lty=2,col=2,lwd=1.5)
lines(result1$r,result1$obs,lty=1,col=1,lwd=1.5)
legend("bottomright",legend=c("Observed Statistic","Expected Statistic","95% Confidence Envelope","Significant Segregation"),lty=c(1,2,NA,NA),col=c(1,2,"lightgrey",rgb(1,0,0,0.1)),lwd=c(1.5,1.5,NA,NA),pch=c(NA,NA,15,15),cex=0.8,pt.cex=1.5)

plot(result2$r,result2$hi,type="n",main="Stratified Permutation",xlab="Distance (in meters)",ylab="L function")
segregation=result2$r[which(result2$obs<result2$lo)]
rect(xleft=min(segregation),xright=max(segregation),ybottom=-1000,ytop=1000,col=rgb(1,0,0,0.1),border=NA)
polygon(x=c(result2$r,rev(result2$r)),y=c(result2$hi,rev(result2$lo)),col="lightgrey",border=NA)
lines(result2$r,result2$mmean,lty=2,col=2,lwd=1.5)
lines(result2$r,result2$obs,lty=1,col=1,lwd=1.5)

dev.print(device=pdf,"./figures/figure3.pdf")

###################################################################
## Assess Segregation/Aggregation Domestic vs Livestock Features ##
###################################################################


## Create Z-transformed distance matrices

distMatrices <- list()
diwan <- list()
allHuts <- list()

for (i in 1:NsitesDiwan)
{
    # Z-scores
    b <- as.matrix(dist(CoreDiwanList[[i]]@coords)) # observed interdistance matrix
    distMatrices[[i]] <- (b - mean(b)) / sd(b) # z-score & store
    a <- CoreDiwanList[[i]]@data
    allHuts[[i]] <- c(which(a$HvsDnoUnmarried == 'D'), which(a$HvsDnoUnmarried == 'H')) # Index All Huts
    diwan[[i]] <- which(a$HvsDnoUnmarried == 'D') # # Index diwans
}


## Retrieve All permutations
perm <- expand.grid(allHuts) ## All permutations with one index value of allHuts randomly assigned as diwan

## compute all possible grand-mean distances from randomly assigned diwan
grandMean <- c()
for (i in 1:nrow(perm))
{
    MeansVector <- c()
    for (j in 1:ncol(perm))
    {
        MeansVector[j] <- mean(distMatrices[[j]][perm[i,j],]) # take the row in the distance matrix and compute the mean
    }
    grandMean[i] <- mean(MeansVector) # mean of means
}

## Compute diwan's grand mean distance
MeanDiwans <- c()
for (i in 1:length(diwan))
{
    MeanDiwans[i] <- mean(distMatrices[[i]][diwan[[i]],])
}

## Compute p-value
pval=(sum(grandMean>=mean(MeanDiwans)))/nrow(perm)
