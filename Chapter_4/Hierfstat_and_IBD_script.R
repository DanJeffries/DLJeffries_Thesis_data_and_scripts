library(parallel)
library(MASS)
library(adegenet)
?read.fstat

body <- read.fstat("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Data/PurecruplusManc_numbered.DAT") ## Load file - converts to GENIND object.

pops <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Data/Fstat_pop_code_labels.txt", header = F)
pop(body) <- pops$V2

library(hierfstat) ## load this in after reading file so as to use the adegenet read.fstat function, not the hierfstat version

UK <- body[c(1:73, 176:196, 311:473)] ## Subset just UK, Belgian and German samples
## 257 individuals
table(pop(UK))
as.data.frame(pop(UK))

### Hierachal fstats in Hierfstat ###

## first need to make the levels column

UKheirFormat <- genind2hierfstat(UK) ## make data frame format
row.names(UKheirFormat) <- seq(257)

poplevel <- UKheirFormat$pop

pool_level <- c(rep(3, 13), rep(1, 27), rep(5, 33), rep(4,21), rep(2,68), rep(1, 27), rep(2,20), rep(1,14), rep(2, 20), rep(1,14))

country_level <- c(rep(1, 40), rep(2,33), rep(3, 21), rep(1,163))

Geolog_level <- c(rep(1,40), rep(2, 54), rep(1, 163))

levs <- data.frame(Geolog_level, country_level, pool_level, poplevel)

data.frame(pop(UK), pool_level) ## check

UKdata <- UKheirFormat[, c(2:14)] ## put the genetic data into a separate data frame

varCOMP <- varcomp.glob(levs,UKdata) ## print varCOMP to get results

### Below are various significance tests.

Geolog_test_temp <- test.between(UKdata, rand.unit=country_level, test=Geolog_level, nperm=1000)
## Note the above test permutes indidividuals between the geological classess but within countries - not that useful I don't think

## So I am going to permute across geological classes but within populations
Geolog_test_final <- test.between(UKdata, rand.unit=poplevel, test=Geolog_level, nperm=1000)

country_test_final <- test.between(UKdata, rand.unit=poplevel, test=country_level, nperm=1000)
country_test_final

pop_level_test <- test.within(UKdata, test=poplevel, within=country_level, nperm=1000)
pop_level_test

pool_level_test <- test.within(UKdata, test=pool_level, within=country_level, nperm=1000)
pool_level_test

## test pop variation between 

pop_in_pool_level_test <- test.within(UKdata, test=poplevel, within=pool_level, nperm=1000)

Geolog_test_temp
Geolog_test_final


### -----------------------------------------------------------------------------------------------------------

### IBD ########

coordinates <- read.csv("~/../Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/R stats/Adegenet IBD/Mantelcoordinates.csv",header = TRUE, row.names = 1) ##table containing xy coordinates
data.frame(row.names(coordinates))
NWEUcoords <- coordinates[c(1:6,11,18:25),]
NWEUcoords ## 15 pops
row.names(NWEUcoords)[9] <- "GFP2" ## ammend one of the pop names

UK$other <- NWEUcoords


UKgenp <- genind2genpop(UK)

Dgen <- dist.genpop(UKgenp,method=2) ##Generates the genetic distance semi-matrix

Dgeo <- dist(UKgenp$other) ## Generates the Euclidian distance semi matrix. Note- vignette says to use "dist(body$other$xy),I dont know why but this wont work for me!
Dgen

ibd <- mantel.randtest(Dgen,Dgeo)
ibd  ## gives the output for the test  NOTE that the "observation" is the R value. Sqaure to get R^2 (obviously)
R2 = ibd$obs^2
R2
## Density plot ---------------------------------------------------------------------

dens <- kde2d(Dgeo,Dgen, n=300, lims=c(-3, 50,-.1,1))  ##Uses a 2-dimensional kernel density estimation (kde2d) to create an object tht displays desinsty of dots

myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))  ## I think sets the colour palatte

plot(Dgeo, Dgen, pch=20,cex=.5) ## Plots the scatter with smaller dots

image(dens, col=transp(myPal(300),.7), add=TRUE)

abline(lm(Dgen~Dgeo))
text(4,0.1, "R^2 = 0.248***")



## --- PW Fsts ---------------------------------------------------------------------------

body <- read.fstat("~/../Dropbox/PhD/Dans_PhD_Shared/Data/Microsatellites/DAPC/Data/Purecruchecked.DAT")

UK <- body[c(1:73, 176:196, 311:473)] ## Subset just UK, Belgian and German samples
## note can do this with seppop and repool better



UKnoGF29 <- UK[,loc=c("L01", "L02", "L04", "L05", "L06", "L07", "L08", "L09", "L10", "L11", "L12", "L13")]

pop(UK)


library(adegenet)
PWfsts <- pairwise.fst(UKnoGF29, res.type = "matrix")

write.table(PWfsts, "~/../Dropbox/PhD/Dans_PhD_Shared/Papers/Crucian status paper/MS figs and Tabs/PW_fsts.tsv", sep = '\t')
