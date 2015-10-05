library(adegenet)
library(maps)
library(mapdata)
library(mapplots)

### This script runs the DAPC analysis used in the Crucian status MS and plots
### the cluster memberships on a map.


#----------------------------------------------------------------------------------


body <- read.fstat("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Data/Purecruchecked.DAT") ## Load file - converts to GENIND object.

UK <- body[c(1:73, 176:196, 311:473)] ## Subset just UK, Belgian and German samples

UK$pop.names ## check


#----------------------------------------------------------------------------------
  ############################## Running DAPC #################################

## As mentioned in text, although the BIC scores suggest 10 clusters, we used 4 to simplify structure

grp <- find.clusters (UK, max.n.clust=100) ## All PCs, 4 clusters

table.value(table(pop(UK), grp$grp), col.lab=paste("inf", 1:12), row.lab=paste("ori", UK$pop.names)) ## look at cluster assignments

dapc2 <- dapc(UK, n.da=100, n.pca=50) ## Run quick DAPC test
temp <- optim.a.score(dapc2) ## check optimal number of PCs to retain

dapc1 <- dapc(UK, grp$grp) ## Run final DAPC

scatter(dapc1, col = mycol) ## Plot clusters on LD scatter

write.csv(dapc1$posterior, "C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Outputs/UK 4 clusters/UKpopclustmembs4.csv") ## write the cluster memberships to a file to be used for plotting on maps

body <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Outputs/UK 4 clusters/UKpopclustmembs4.csv", header = T) ## read in file if starting from here



UKmeanclustmembs <- sapply(split(body[2:7],body$population),colMeans) ## make population means for each cluster


write.csv(UKmeanclustmembs, "C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Outputs/UK 4 clusters/UKmeanclustmembs.csv") ## write to file

pies <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/R stats/Mapplots/UKmeanclustmembs.csv") ## read in cluster means if starting from here

cords <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Mantelcoordinates.csv", header = T) ##load cordinates for population points

keypie <- read.csv("C:/Users/Dan/Dropbox/PhD/Dan's PhD (Shared)/Data/Microsatellites/DAPC/Outputs/All data 4 clusters/keypie.csv", header = T) ## made a file for the key


#--------------------------------------------------------------------------------------------------------------------------------------

##############################    Plotting cluster memberships on the map ################################################



map("worldHires", xlim=c(-10, 10), ylim=c(48,57), col="gray80", fill=TRUE) ## UK zoomed map

points(cords$lat, cords$lon, pch=19, col="black", cex=0.7)  ## Plots my sample locations on the map. But dont know how to lable, yet!


## UK Pies ## Radiuses are AR/2


add.pie(pies$HOLT,x=-3.5,y=52,labels="",radius=0.847,edges=200,clockwise=T, col = mycol)

add.pie(pies$CAKE,x=?-3.5,y=53.5,labels="",radius=0.628833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$CCS,x=-2,y=49.7?,labels="",radius=0.664541667,edges=200,clockwise=T, col = mycol)

add.pie(pies$FFF,x=-3,y=50.7,labels="",radius=0.5,edges=200,clockwise=T, col = mycol)

add.pie(pies$BF,x=1.65,y=50.5,labels="",radius=0.715833334,edges=200,clockwise=T, col = mycol)

add.pie(pies$RM,x=3.2,y=51.8,labels="",radius=0.776166667,edges=200,clockwise=T, col = mycol)

add.pie(pies$PRIM,x=-1.5,y=54.8,labels="",radius=0.64975,edges=200,clockwise=T, col = mycol)

add.pie(pies$OTOM,x=0.8,y=55.3?,labels="",radius=0.635333334,edges=200,clockwise=T, col = mycol)

add.pie(pies$GFP,x=3.2,y=55, labels="",radius=0.737625,edges=200,clockwise=T, col = mycol)

add.pie(pies$RAIL,x=5.3,y=54.3, labels="",radius=0.77325,edges=200,clockwise=T, col = mycol)

add.pie(pies$MOAT,x=4.9,y=52.8?,labels="",radius=0.719875,edges=200,clockwise=T, col = mycol)



## Belgian Pies ##

add.pie(pies$BOK,x=7.5,y=50,labels="",radius=0.711208334,edges=200,clockwise=T, col = mycol)

add.pie(pies$MVW,x=2.8,y=49.2,labels="",radius=0.739416667,edges=200,clockwise=T, col = mycol)

add.pie(pies$MVWZ,x=5.2,y=49.6,labels="",radius=0.737166667,edges=200,clockwise=T, col = mycol)

## German Pie ##

add.pie(pies	$	FFG	,	x	=	7.558594	,	y	=	51.890053	,	labels	=	""	,	radius	=	1.186083334	,	edges	=	200,	clockwise	=	T, col = mycol)																																			

## Text labels ##
textxy(	-3.5	,	52	,	"GBR10"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	?-3.5	,	53.5	,	"GBR4"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	-2	,	49.7?	,	"GBR1"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	-3	,	50.7	,	"GBR2"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	1.65	,	50.5	,	"GBR8"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	3.2	,	51.8	,	"GBR6"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	-1.5	,	54.8	,	"GBR5"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	0.8	,	55.3?	,	"GBR9"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	3.2	,	55	,	"GBR3"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	5.3	,	54.3	,	"GBR11"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	4.9	,	52.8?	,	"GBR7"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	7.5	,	50	,	"BEL1"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	2.8	,	49.2	,	"BEL2"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	5.2	,	49.6	,	"BEL3"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)
textxy(	7.558594	,	51.890053	,	"GER2"	,	cex	=	0.8	,	col	=	"white",	offset	=	0)


## Add key ##


add.pie(keypie$K,x=-7.5,y=55.5,labels="",radius=1.2, edges=200,clockwise=T, col = "gray80")
add.pie(keypie$K,x=-7.5,y=55.5,labels="",radius=0.8, edges=200,clockwise=T, col = "gray80")
add.pie(keypie$K,x=-7.5,y=55.5,labels="",radius=0.4, edges=200,clockwise=T, col = "gray80")


text(-7.5,54, "Cluster 1A", col = "blue")
text(-7.5,53.7, "Cluster 2A", col = "red3")
text(-6,54, "Cluster 3A", col = "darkgreen")
text(-6,53.7, "Cluster 4A", col = "black")


text(-7.5, 56.5, "AR = 1.2", cex = 0.8)
text(-7.5, 56, "= 0.8", cex = 0.8)
text(-7.5, 55.5, "= 0.4", cex = 0.8)
