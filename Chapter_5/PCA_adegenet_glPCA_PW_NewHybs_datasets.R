library(adegenet)
library(scatterplot3d)


setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples_M4_m4_N4/NewHybrids_PW_datafiles/")


## Full datasets (strict Ustacks and Cstacks params) ----------------------------------------------------------------------------------------------

# Populations params: r07, p27/31 , m8

All_data <- read.PLINK('../populations_r07_p27_m8/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
All_pops <- read.delim("../populations_r07_p27_m8/pop.names", header = F)
pop(All_data) <- All_pops$V1
glPlot(All_data, posi="topleft", yaxt = 'n') ## Not too much missing data!
length(All_data$loc.names) ## 16048 loci


All_data_PCA <- glPca(All_data, nf = 5)
scatter.glPca(All_data_PCA , xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## PC1 = Cru -> Gold/Gib, PC2 = mtDNA lin1 -> Lin2
scatter.glPca(All_data_PCA , xax = 3, yax = 4, clabel = 0.7, posi = "bottomright") ## PC3 = Cru/Gold/Gib -> Common, ## PC4 = Gold -> gib ## PC4 = within Gib variation
scatter.glPca(All_data_PCA , xax = 4, yax = 5, clabel = 0.7, posi = "bottomright") ## PC5 = within Gib variation


# Populations params: r07, p26/30 (noSD), m8

All_noSD_p26 <- read.PLINK('../populations_r07_p26_m8_noSD/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
All_noSD_p26_names <- read.delim("../populations_r07_p26_m8_noSD/pop.names", header = F)
pop(All_noSD_p26) <- All_noSD_p26_names$V1
glPlot(All_noSD_p26, posi="topleft", yaxt = 'n') ## Less missing data than with SD in there, but still a lot falling out for common carp
length(All_noSD_p26$loc.names) ## 12117 Loci


All_noSD_p26_PCA <- glPca(All_noSD_p26, nf = 5)
scatter.glPca(All_noSD_p26_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topright") ## PC1 = Cru -> Gold/Gib, PC2 = mtDNA lin1 -> Lin2
scatter.glPca(All_noSD_p26_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## PC3 = Cru/Gold/Gib -> Common, ## PC4 = within Gib variation
scatter.glPca(All_noSD_p26_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "topright") ## PC5 = Gold -> gib


## Populations params: r0.7, p31/31, m8
All_r07_p31 <- read.PLINK("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/All_samples_M8_N5_m8/populations_r0.7_pall_m8/plink.raw", , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(All_r07_p31, posi="topleft", yaxt = 'n')

length(All_r07_p31$loc.names) ## 299 loci

All_r07_p31_PCA <- glPca(All_r07_p31, nf = 5) ## <- this is a weird one - probably don't use
scatter.glPca(All_r07_p31_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topright") ## PC1 - Pulls out the Cru-gold Common carp variation, PC2 - pulls out cru mtDNA lin1 -> lin2 variation
scatter.glPca(All_r07_p31_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "topleft") ## PC3 - Goldfish -> Gib -> SWED crucian var (a bit weird)
scatter.glPca(All_r07_p31_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "topleft") ## PC4 - Gold -> cru, PC5 pulls out north -> south in mtDNA lin1 crucian


## Populations params: r07, p30/30 (noSD), m8

All_data_p30_noSD <- read.PLINK('../populations_r07_p30_all_noSD/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

All_data_p30_noSD_names <- read.delim("../populations_r07_p30_all_noSD/pop.names", header = F)
pop(All_data_p30_noSD) <- All_data_p30_noSD_names$V1
glPlot(All_data_p30_noSD, posi="topleft", yaxt = 'n')
length(All_data_p30_noSD$loc.names) ## 585 loci

All_data_p30_noSD_PCA <- glPca(All_data_p30_noSD, nf = 5) ## <- this is a weird one - probably don't use
scatter.glPca(All_data_p30_noSD_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topright") ## PC1 - Pulls out the Cru-gold Common carp variation, PC2 - pulls out cru mtDNA lin1 -> lin2 variation
scatter.glPca(All_data_p30_noSD_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "topright") ## PC3 - Goldfish -> Gib -> Don lineage var (a bit weird)
scatter.glPca(All_data_p30_noSD_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "topleft") ## PC4 - similar to PC3, PC5 pulls out Gib -> crucian



## Populations r0.7, p31/31 where common carp are separated into their own population in the pop_codes map. 
## This should limit SNPs to only those present in both common carp and all other spp. 

All_noSD_common_sep <- read.PLINK("../populations_r07_p31_noSD_common_labelled_m8/plink.raw", , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
All_noSD_common_sep_names <- read.delim("../populations_r07_p31_noSD_common_labelled_m8/pop.names", header = F)
pop(All_noSD_common_sep) <- All_noSD_common_sep_names$V1
length(All_noSD_common_sep$loc.names) ## 316

All_noSD_common_sep_PCA <- glPca(All_noSD_common_sep, nf = 5) # <-  again this is a weird one - probably don't use
scatter.glPca(All_noSD_common_sep_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topright")
scatter.glPca(All_noSD_common_sep_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomleft")
scatter.glPca(All_noSD_common_sep_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft")


## Pooled invasive species - here all goldfish are labelled as such regardless of population, the same for gib and common.
# But any amigous samples are left in their origonal populations

# Populations params: r_0.7, p24/27 (using the relaxed Ustacks and Cstacks params)

pooled_invasives <- read.PLINK("/home/dan/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter 4 Hybridisation and introgression/data/RAD/plink.raw", , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
pooled_invasives$n.loc # = 985
glPlot(pooled_invasives, posi="topleft", yaxt = 'n')

pooled_invasives_PCA <- glPca(pooled_invasives, nf = 5) # <- the variation is much better spread across the PCs here (looking at eigenvalues)
scatter.glPca(pooled_invasives_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topleft") ## PC1 - Common -> Cru, PC2 = Gib -> Gold -> Cru
scatter.glPca(pooled_invasives_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomright") ## PC3 = Cru mtDNA lin1 - lin2
scatter.glPca(pooled_invasives_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "topleft") ## PC4 = gib -> gold
scatter.glPca(pooled_invasives_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "topleft") ## same as PC4


################################################################################

# I think overall - using all-datasets is not a good way to go about delimiting the species. 
# The Pw PCAs below allow me to be more stringent with what population criteria to use for SNPs. (note that I could have done this in PyRAD as well - which is more flexible for the presence absence of data in groups of populations.)

# Unless the PCA works better with the ref aligned data

##Reference aligned data  ----------------------------------------------------------------------------------------------

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_all/populations_r07_p_30_m8/")

# Populations params: r07, p27/31 , m8

All_data <- read.PLINK('./plink.raw', chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
All_pops <- read.delim("pops.txt", header = F)
All_samples <- read.delim("sample_names.txt", header = F)

pop(All_data) <- All_pops$V1
glPlot(All_data, posi="topleft", yaxt = 'n') ## Not too much missing data!
length(All_data$loc.names) #1296 loci


All_data_PCA <- glPca(All_data, nf = 5)


## Plotting PCA ----------------------------------------------

assignments <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/RAD_Figs_FINAL/All_spp_assignments", header = F)
table(assignments$V2)

length(assignments$V2)
length(All_data_PCA$scores[,1])

## Colours --------------------

C.car = "deepskyblue3"
C.a.aur = "firebrick4"
C.a.gib = "darkorange"
C.caxC.a.au = "darkorchid4"
C.caxC.a.gi = "deeppink"
C.caxC.carp = "grey54"
C.carp = "black"
CaxGi_BKX = "green"
c.a.aur_x_C.carp = "yellow"

spp.cols <- c( C.a.aur, c.a.aur_x_C.carp,  C.a.gib, C.car, C.caxC.a.au, C.caxC.a.gi, CaxGi_BKX, C.caxC.carp, C.carp)

## Get variation explained
total_eigs <- sum(All_data_PCA$eig)
PC1 <- paste("(",round((All_data_PCA$eig[1]/total_eigs)*100, digits = 2), "%)", sep = "")
PC2 <- paste("(",round((All_data_PCA$eig[2]/total_eigs)*100, digits = 2), "%)", sep = "")
PC3 <- paste("(",round((All_data_PCA$eig[3]/total_eigs)*100, digits = 2), "%)", sep = "")
PC4 <- paste("(",round((All_data_PCA$eig[4]/total_eigs)*100, digits = 2), "%)", sep = "")
PC5 <- paste("(",round((All_data_PCA$eig[5]/total_eigs)*100, digits = 2), "%)", sep = "")


## 3D Scatterplot

assignments <- read.delim("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/RAD_Figs_FINAL/All_spp_assignments", header = F)

scat3d <- scatterplot3d(All_data_PCA$scores[,3], All_data_PCA$scores[,5],All_data_PCA$scores[,1], 
                        color = "black", 
                        bg = spp.cols[assignments$V2], 
                        pch = 21, 
                        cex.symbols = 2, 
                        xlab = paste("PC3", PC3), 
                        ylab = paste("PC5", PC5),
                        zlab = paste("PC1", PC1),
                        angle = 60, 
                        type = "h",
                        scale.y = 1.5, 
                        lwd= 1)



scat3d <- scatterplot3d(All_data_PCA$scores[,3], All_data_PCA$scores[,5],All_data_PCA$scores[,1], 
                        color = "black", 
                        bg = spp.cols[assignments$V2], 
                        pch = 21, 
                        cex.symbols = 2, 
                        xlab = paste("PC3", PC3), 
                        ylab = paste("PC5", PC5),
                        zlab = paste("PC1", PC1),
                        angle = 80, 
                        type = "h",
                        scale.y = 2, 
                        lwd= 1)




#scat3d$points3d(All_data_PCA$scores[,1], All_data_PCA$scores[,5],All_data_PCA$scores[,2], col = spp.cols[assignments$V2], type = "h", cex.symbol = 0.1)

par(mar = c(4,4,4,4))
## PC2 vs 4
plot(All_data_PCA$scores[,2], All_data_PCA$scores[,4], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC2", PC2), ylab = paste("PC4", PC4))
text(All_data_PCA$scores[,2], All_data_PCA$scores[,4], All_samples$V1, cex = 1, col = "black")


# Key --------------------------------------------------------
points(17,6, bg = C.car_col, pch = 21, cex = 2)
points(17,5.8, bg = C.a.aur, pch = 21, cex = 2)
points(17,5.6, bg = C.a.gib, pch = 21, cex = 2)
points(17,5.4, bg = C.carp, pch = 21, cex = 2)
points(17,5.2, bg = C.caxC.a.au, pch = 21, cex = 2)
points(17,5.0, bg = C.caxC.a.gi, pch = 21, cex = 2)
points(17,4.8, bg = C.caxC.carp, pch = 21, cex = 2)
points(17,4.6, bg = CaxGi_BKX, pch = 21, cex = 2)
points(17,4.4, bg = c.a.aur_x_C.carp, pch = 21, cex = 2)



text(17,6, expression(italic("C. carassius")), pos = 4, cex = 0.9)
text(17,5.8, expression(italic("C. a. auratus")), pos = 4, cex = 0.9)
text(17,5.6, expression(italic("C. a. gibelio")),pos = 4, cex = 0.9)
text(17,5.4, expression(italic("C. carpio")),pos = 4, cex = 0.9)
text(17,5.2, expression(italic("C. carassius x C. a. auratus")),pos = 4, cex = 0.9)
text(17,5.0, expression(italic("C. carassius x C. a. gibelio")),pos = 4, cex = 0.9)
text(17,4.8, expression(italic("C. carassius x C. carpio")),pos = 4, cex = 0.9)
text(17,4.6, expression(italic("C. carassius x (C.carassius x C. carpio)")),pos = 4, cex = 0.9)
text(17,4.4, expression(italic("C. a. auratus x C. carpio")),pos = 4, cex = 0.9)


## PCs 3 and 4

plot(All_data_PCA$scores[,3], All_data_PCA$scores[,4], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC3", PC3), ylab = paste("PC4", PC4))
#text(All_data_PCA$scores[,3], All_data_PCA$scores[,4], All_samples$V1, cex = 0.7, col = "black")


# Key --------------------------------------------------------
points(2.2,-2, bg = C.car_col, pch = 21, cex = 2)
points(2.2,-2.1, bg = C.a.aur, pch = 21, cex = 2)
points(2.2,-2.2, bg = C.a.gib, pch = 21, cex = 2)
points(2.2,-2.3, bg = C.carp, pch = 21, cex = 2)
points(2.2,-2.4, bg = C.caxC.a.au, pch = 21, cex = 2)
points(2.2,-2.5, bg = C.caxC.a.gi, pch = 21, cex = 2)
points(2.2,-2.6, bg = C.caxC.carp, pch = 21, cex = 2)
points(2.2,-2.7, bg = CaxGi_BKX, pch = 21, cex = 2)
points(2.2,-2.8, bg = c.a.aur_x_C.carp, pch = 21, cex = 2)



text(2.2,-2, expression(italic("C. carassius")), pos = 4, cex = 0.9)
text(2.2,-2.1, expression(italic("C. a. auratus")), pos = 4, cex = 0.9)
text(2.2,-2.2, expression(italic("C. a. gibelio")),pos = 4, cex = 0.9)
text(2.2,-2.3, expression(italic("C. carpio")),pos = 4, cex = 0.9)
text(2.2,-2.4, expression(italic("C. carassius x C. a. auratus")),pos = 4, cex = 0.9)
text(2.2,-2.5, expression(italic("C. carassius x C. a. gibelio")),pos = 4, cex = 0.9)
text(2.2,-2.6, expression(italic("C. carassius x C. carpio")),pos = 4, cex = 0.9)
text(2.2,-2.7, expression(italic("C. carassius x (C.carassius x C. carpio)")),pos = 4, cex = 0.9)
text(2.2,-2.8, expression(italic("C. a. auratus x C. carpio")),pos = 4, cex = 0.9)




## PC1 vs PC5 most informative for species delimitation

plot(All_data_PCA$scores[,1], All_data_PCA$scores[,5], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC5", PC5))
text(All_data_PCA$scores[,1], All_data_PCA$scores[,5], All_samples$V1, cex = 0.7, col = "black")



# Key --------------------------------------------------------
points(3,2, bg = C.car_col, pch = 21, cex = 2)
points(3,1.87, bg = C.a.aur, pch = 21, cex = 2)
points(3,1.74, bg = C.a.gib, pch = 21, cex = 2)
points(3,1.61, bg = C.carp, pch = 21, cex = 2)
points(3,1.48, bg = C.caxC.a.au, pch = 21, cex = 2)
points(3,1.35, bg = C.caxC.a.gi, pch = 21, cex = 2)
points(3,1.22, bg = C.caxC.carp, pch = 21, cex = 2)
points(3,1.09, bg = CaxGi_BKX, pch = 21, cex = 2)
points(3,0.96, bg = c.a.aur_x_C.carp, pch = 21, cex = 2)



text(3,2, expression(italic("C. carassius")), pos = 4, cex = 0.9)
text(3,1.87, expression(italic("C. a. auratus")), pos = 4, cex = 0.9)
text(3,1.74, expression(italic("C. a. gibelio")),pos = 4, cex = 0.9)
text(3,1.61, expression(italic("C. carpio")),pos = 4, cex = 0.9)
text(3,1.48, expression(italic("C. carassius x C. a. auratus")),pos = 4, cex = 0.9)
text(3,1.35, expression(italic("C. carassius x C. a. gibelio")),pos = 4, cex = 0.9)
text(3,1.22, expression(italic("C. carassius x C. carpio")),pos = 4, cex = 0.9)
text(3,1.09, expression(italic("C. carassius x (C.carassius x C. carpio)")),pos = 4, cex = 0.9)
text(3,0.96, expression(italic("C. a. auratus x C. carpio")),pos = 4, cex = 0.9)






row.names(All_data_PCA$scores) <- All_samples$V1

head(pca_df)
pca_df <- as.data.frame(All_data_PCA$scores)
pca_df$pop <- All_pops$V1

pca <- qplot(PC1, PC2, data = pca_df, colour = pop, cex = 2)
?qplot
pca + annotate("text", x = pca_df$PC1, y = pca_df$PC2, label = row.names(pca_df), cex = 4)


# - GETTING LOADINGS ---------------------------------------------------------------------------------------------------------------

#PC1 - Carassius - cyprinus
loadingplot(All_data_PCA$loadings[,1], threshold = quantile(abs(All_data_PCA$loadings[,1]),0.95))
PC1_95th<- quantile(abs(All_data_PCA$loadings[,1]),0.95) ## this is the value for the 75th percentile
PC1_95th
loadings <- data.frame(Stacks_ID= All_data$loc.names, loading = abs(All_data_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_PC1 <- subset(loadings, loading >= PC1_95th) ## subset for only those above the 75th %ile
length(Diagnostic_PC1$loading) # = 66

#PC2 - between Carassius carassius lineages 1 and 2

loadingplot(All_data_PCA$loadings[,2], threshold = quantile(abs(All_data_PCA$loadings[,2]),0.95))
PC2_95th<- quantile(abs(All_data_PCA$loadings[,2]),0.95) ## this is the value for the 75th percentile
PC2_95th
loadings <- data.frame(Stacks_ID= All_data$loc.names, loading = abs(All_data_PCA$loadings[,2])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_PC2 <- subset(loadings, loading >= PC2_95th) ## subset for only those above the 75th %ile
length(Diagnostic_PC2$loading) # = 65


#PC3 - Carassius carassius - Carassius auratus complex

loadingplot(All_data_PCA$loadings[,3], threshold = quantile(abs(All_data_PCA$loadings[,3]),0.95))
PC3_95th<- quantile(abs(All_data_PCA$loadings[,3]),0.95) ## this is the value for the 75th percentile
PC3_95th
loadings <- data.frame(Stacks_ID= All_data$loc.names, loading = abs(All_data_PCA$loadings[,3])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_PC3 <- subset(loadings, loading >= PC3_95th) ## subset for only those above the 75th %ile
length(Diagnostic_PC3$loading) # = 65

#PC4 - Carassius carassius - Carassius auratus complex

loadingplot(All_data_PCA$loadings[,4], threshold = quantile(abs(All_data_PCA$loadings[,4]),0.95))
PC4_95th<- quantile(abs(All_data_PCA$loadings[,4]),0.95) ## this is the value for the 75th percentile
PC4_95th
loadings <- data.frame(Stacks_ID= All_data$loc.names, loading = abs(All_data_PCA$loadings[,4])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_PC4 <- subset(loadings, loading >= PC4_95th) ## subset for only those above the 75th %ile
length(Diagnostic_PC4$loading) # = 65



## --------Pairwise datasets and PCAs -----------------------------------------------------------------------------------------------

# PW DATASETS 

## Crucian common
Cru_comm <- read.PLINK('Cru_comm_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_comm_names <- read.delim("Cru_comm_pop.names", header = F)
Cru_comm$ind.names
pop(Cru_comm) <- Cru_comm_names$V1
glPlot(Cru_comm, posi="topleft", yaxt = 'n')
Cru_comm$n.loc # = 15054 - a lot of data and a lot missing!

## Common carp with common labelled as separate pop. All pops
Cru_comm_sep <- read.PLINK("../NewHybrids_PW_datafiles//Cru_comm_populations_comm_labelled/plink.raw", , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(Cru_comm_sep, posi="topleft", yaxt = 'n') ## Data looks very good! A couple of fish in there that may be c.auratus hybs though
pop(Cru_comm_sep)
Cru_comm_sep$n.loc # 314  - thats good enough I think

Cru_comm_sep_PCA <- glPca(Cru_comm_sep, nf = 5) ## <- looks good
scatter.glPca(Cru_comm_sep_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## PC1 = Cru -> Comm - very good. PC2 mtDNA lin1 -> lin2
scatter.glPca(Cru_comm_sep_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomright") ## PC3 = Within crucian mtDNA lin1 only
scatter.glPca(Cru_comm_sep_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## PC4 not very informative



## Crucian goldfish only
Cru_gold <- read.PLINK('Cru_gold_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_gold_names <- read.delim("Cru_gold_pop.names", header = F)
pop(Cru_gold) <- Cru_gold_names$V1
glPlot(Cru_gold, posi="topleft", yaxt = 'n') ## Data looks very good!
Cru_gold$n.loc # 1216


Cru_gold_PCA <- glPca(Cru_gold, nf = 5) ## <- looks good
scatter.glPca(Cru_gold_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## PC1 = Cru -> Goldfish - very good. PC2
scatter.glPca(Cru_gold_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomleft") ## PC3 = Within crucian mtDNA lin1 only
scatter.glPca(Cru_gold_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## PC4 EP01 - common carp hybrid

## Crucian gibel only 
Cru_gib <- read.PLINK('Cru_gib_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_gib_names <- read.delim("Cru_gib_pop.names", header = F)
pop(Cru_gib) <- Cru_gib_names$V1
glPlot(Cru_gib, posi="topleft", yaxt = 'n') ## Data looks very good!
Cru_gib$n.loc # 872
Cru_gib$ind.names

Cru_gib_PCA <- glPca(Cru_gib, nf = 5)

scatter.glPca(Cru_gib_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## PC1 = Cru - gibel, PC2 Cru mtDNA lin1 - lin2
scatter.glPca(Cru_gib_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomleft") ## PC3 = Within Cru mtDNA lin 1
scatter.glPca(Cru_gib_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomright") ## EP01 - common carp hybrid.
scatter.glPca(Cru_gib_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "bottomright") ## I think the var in the gibel lineages.


### Getting PC loadings for pairwise spp diagnostic SNPs ------------------------------------------------------------------------------------

## Cru-Comm

loadingplot(Cru_comm_sep_PCA$loadings[,1], threshold = quantile(abs(Cru_comm_sep_PCA$loadings[,1]),0.75))
Cru_Comm_Q_75 <- quantile(abs(Cru_comm_sep_PCA$loadings[,1]),0.75) ## this is the value for the 75th percentile

loadings <- data.frame(Stacks_ID= Cru_comm_sep$loc.names, loading = abs(Cru_comm_sep_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_Cru_comm <- subset(loadings, loading >= Cru_Comm_Q_75) ## subset for only those above the 75th %ile
length(Diagnostic_Cru_comm$loading) # = 84 loci 
Diagnostic_Cru_comm

write.table(Diagnostic_Cru_comm$Stacks_ID, "./Cru_comm_populations_comm_labelled/Diagnostic_snps.txt") ## writing to file 

## Cru-gold

loadingplot(Cru_gold_PCA$loadings[,1], threshold = quantile(abs(Cru_gold_PCA$loadings[,1]),0.95))
Cru_gold_Q_75 <- quantile(abs(Cru_gold_PCA$loadings[,1]),0.75)
Cru_gold_Q_75

loadings <- data.frame(Stacks_ID= Cru_gold$loc.names, loading = abs(Cru_gold_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_Cru_gold <- subset(loadings, loading >= Cru_gold_Q_75) ## subset for only those above the 75th %ile
length(Diagnostic_Cru_gold$loading) #  = 304
Diagnostic_Cru_gold

write.table(Diagnostic_Cru_gold$Stacks_ID, "./Cru_gold_populations/Diagnostic_snps.txt") ## writing to file 

## Cru - gib

loadingplot(Cru_gib_PCA$loadings[,1], threshold = quantile(abs(Cru_gib_PCA$loadings[,1]),0.75))
Cru_gib_Q_75 <- quantile(abs(Cru_gib_PCA$loadings[,1]),0.75)

loadings <- data.frame(Stacks_ID= Cru_gib$loc.names, loading = abs(Cru_gib_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_Cru_gib <- subset(loadings, loading >= Cru_gib_Q_75) ## subset for only those above the 75th %ile
length(Diagnostic_Cru_gib$loading) # = 218
Diagnostic_Cru_gib

write.table(Diagnostic_Cru_gib$Stacks_ID, "./Cru_gib_populations/Diagnostic_snps.txt") ## writing to file 


## These Diagnostic loci have been written to a file in the folder for that comparison (see above) and 
## 


## High FST Diagnostic SNPs



Cru_comm_DIAG <- read.PLINK("../NewHybrids_PW_datafiles/Cru_comm_populations_comm_labelled/plink.raw", chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
glPlot(Cru_comm_DIAG, posi="topleft", yaxt = 'n') ## Data looks very good! A couple of fish in there that may be c.auratus hybs though
pop(Cru_comm_DIAG)
Cru_comm_DIAG$n.loc # 56

##Check on PCA
Cru_comm_DIAG_PCA <- glPca(Cru_comm_DIAG, nf = 5) ## <- looks good
scatter.glPca(Cru_comm_DIAG_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## 
scatter.glPca(Cru_comm_DIAG_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomright") ## 
scatter.glPca(Cru_comm_DIAG_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## 


## Crucian goldfish only
Cru_gold_DIAG <- read.PLINK('../NewHybrids_PW_datafiles/Cru_gold_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

pop(Cru_gold_DIAG) <- Cru_gold_names$V1
glPlot(Cru_gold_DIAG, posi="topleft", yaxt = 'n') ## Data looks very good!
Cru_gold_DIAG$n.loc # 124


Cru_gold_DIAG_PCA <- glPca(Cru_gold_DIAG, nf = 5) ## <- looks good
scatter.glPca(Cru_gold_DIAG_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## 
scatter.glPca(Cru_gold_DIAG_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomleft") ## 
scatter.glPca(Cru_gold_DIAG_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") 

## Crucian gibel only 
Cru_gib_DIAG <- read.PLINK('../NewHybrids_PW_datafiles/Cru_gib_populations/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)

glPlot(Cru_gib_DIAG, posi="topleft", yaxt = 'n') ## Data looks very good!
Cru_gib_DIAG$n.loc # 60
Cru_gib_DIAG$ind.names

Cru_gib_DIAG_PCA <- glPca(Cru_gib_DIAG, nf = 5)

scatter.glPca(Cru_gib_DIAG_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## PC1 = Cru - gibel, PC2 Cru mtDNA lin1 - lin2
scatter.glPca(Cru_gib_DIAG_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomleft") ## PC3 = Within Cru mtDNA lin 1
scatter.glPca(Cru_gib_DIAG_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomright") ## EP01 - common carp hybrid.
scatter.glPca(Cru_gib_DIAG_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "bottomright") ## I think the var in the gibel lineages.



## Reference aligned ----------------------------------------------------------------------------------------------

setwd("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_all/NewHybrids_populations/")

cruComm_test <- read.PLINK("/media/dan/34D5D1CE642D7E36/2013076_Hanfling_Bernd/Stacks/Stacks_analyses_TRIMMED/With_reference/Data_links/pstacks_all/NewHybrids_populations/Cru_comm/Cru_comm_pop_codes.txt_dir/plink.raw")
cruComm_test$n.loc
glPlot(cruComm_test, posi="topleft", yaxt = 'n') ## Data looks very good!

## Cru_comm ---------------------

Cru_comm <- read.PLINK('./Cru_comm/Cru_comm_pop_codes.txt_dir/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_comm$n.loc # = 3853

glPlot(Cru_comm, posi="topleft", yaxt = 'n') ## Data looks very good!

Cru_comm$ind.names

Cru_comm_PCA <- glPca(Cru_comm, nf = 5)

scatter.glPca(Cru_comm_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "bottomright") ## PC1 = Cru - gibel, PC2 Cru mtDNA lin1 - lin2
scatter.glPca(Cru_comm_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomright") ## PC3 = Within Cru mtDNA lin 1
scatter.glPca(Cru_comm_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## EP01 - common carp hybrid.
scatter.glPca(Cru_gib_DIAG_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "bottomright") ## I think the var in the gibel lineages.


## Cru_gold ---------------------

Cru_gold <- read.PLINK('./Cru_gold/Cru_gold_pop_codes_noEP01.txt_dir/ plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_gold$n.loc # = 2728

glPlot(Cru_gold, posi="topleft", yaxt = 'n') ## Data looks very good!

Cru_gold$ind.names

Cru_gold_PCA <- glPca(Cru_gold, nf = 5)

scatter.glPca(Cru_gold_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topright") ## PC1 = Cru - gibel, PC2 Cru mtDNA lin1 - lin2
scatter.glPca(Cru_gold_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomright") ## PC3 = Within Cru mtDNA lin 1
scatter.glPca(Cru_gold_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## EP01 - common carp hybrid.
scatter.glPca(Cru_gold_DIAG_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "bottomright") ## I think the var in the gibel lineages.


## Cru_gib ---------------------

Cru_gib <- read.PLINK('./Cru_gib/Cru_gib_pop_codes_noEP01_newnames.txt_dir/plink.raw', , chunkSize=1000, parallel = TRUE, n.cores=3, saveNbAlleles=T)
Cru_gib$n.loc # = 951

glPlot(Cru_gib, posi="topleft", yaxt = 'n') ## Data looks very good!

Cru_gib$ind.names

Cru_gib_PCA <- glPca(Cru_gib, nf = 5)

scatter.glPca(Cru_gib_PCA, xax = 1, yax = 2, clabel = 0.7, posi = "topright") ## PC1 = Cru - gibel, PC2 Cru mtDNA lin1 - lin2
scatter.glPca(Cru_gib_PCA, xax = 2, yax = 3, clabel = 0.7, posi = "bottomright") ## PC3 = Within Cru mtDNA lin 1
scatter.glPca(Cru_gib_PCA, xax = 3, yax = 4, clabel = 0.7, posi = "bottomleft") ## EP01 - common carp hybrid.
scatter.glPca(Cru_gib_DIAG_PCA, xax = 4, yax = 5, clabel = 0.7, posi = "bottomright") ## I think the var in the gibel lineages.


#### Finding diagnostic SNPs from PCA ---------------------------------------------------------------------------------------------

## Cru-Comm

loadingplot(Cru_comm_PCA$loadings[,1], threshold = quantile(abs(Cru_comm_PCA$loadings[,1]),0.75))
Cru_Comm_Q_75 <- quantile(abs(Cru_comm_PCA$loadings[,1]),0.75) ## this is the value for the 75th percentile

loadings <- data.frame(Stacks_ID=Cru_comm$loc.names, loading = abs(Cru_comm_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_Cru_comm <- subset(loadings, loading >= Cru_Comm_Q_75) ## subset for only those above the 75th %ile
length(Diagnostic_Cru_comm$loading) # = 970 loci 
Diagnostic_Cru_comm

write.table(Diagnostic_Cru_comm$Stacks_ID, "./Cru_comm//Cru_comm_pop_codes.txt_dir/Diagnostic_snps.txt") ## writing to file 

## Cru-gold

loadingplot(Cru_gold_PCA$loadings[,1], threshold = quantile(abs(Cru_gold_PCA$loadings[,1]),0.75))
Cru_gold_Q_75 <- quantile(abs(Cru_gold_PCA$loadings[,1]),0.75)
Cru_gold_Q_75

loadings <- data.frame(Stacks_ID= Cru_gold$loc.names, loading = abs(Cru_gold_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_Cru_gold <- subset(loadings, loading >= Cru_gold_Q_75) ## subset for only those above the 75th %ile
length(Diagnostic_Cru_gold$loading) #  = 624
Diagnostic_Cru_gold

write.table(Diagnostic_Cru_gold$Stacks_ID, "./Cru_gold//Cru_gold_pop_codes_noEP01.txt_dir/Diagnostic_snps.txt") ## writing to file 

## Cru - gib

loadingplot(Cru_gib_PCA$loadings[,1], threshold = quantile(abs(Cru_gib_PCA$loadings[,1]),0.75))
Cru_gib_Q_75 <- quantile(abs(Cru_gib_PCA$loadings[,1]),0.75)

loadings <- data.frame(Stacks_ID= Cru_gib$loc.names, loading = abs(Cru_gib_PCA$loadings[,1])) ## Put loadings into a data frame with the stacks ID - order is retained

Diagnostic_Cru_gib <- subset(loadings, loading >= Cru_gib_Q_75) ## subset for only those above the 75th %ile
length(Diagnostic_Cru_gib$loading) # = 238
Diagnostic_Cru_gib

write.table(Diagnostic_Cru_gib$Stacks_ID, "./Cru_gib/Cru_gib_pop_codes_noEP01_newnames.txt_dir/Diagnostic_snps.txt") ## writing to file 


## These Diagnostic loci have been written to a file in the folder for that comparison (see above) and 
## 





