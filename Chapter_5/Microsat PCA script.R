library(adegenet)

setwd("~/../Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats")
setwd("~/Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats")

## Get data

micro_data <- read.fstat("Car_refined_final.dat", missing = "mean") ## replace NAs with mean genotype as advised in adegenet manual for centered PCA

assignments <- read.delim("../../data/Microsat_figs_final/All_spp_assignments_newcodes.txt", header= F) # Assignments based on Newhybs analysis

table(assignments$V2)
indNames(micro_data) <- assignments$V1


## All data PCA ---------------------------------------------------------------------------------------------------
X <- scaleGen(micro_data, missing="mean") # scale data and deal with missing data


micro_pca <- dudi.pca(X, scale = F, center = T, nf = 3) # Perform centered pca (Retained 4 PCs)

## Get variation explained
total_eigs <- sum(micro_pca$eig)
PC1 <- paste("(",round((micro_pca$eig[1]/total_eigs)*100, digits = 2), "%)", sep = "")
PC2 <- paste("(",round((micro_pca$eig[2]/total_eigs)*100, digits = 2), "%)", sep = "")
PC3 <- paste("(",round((micro_pca$eig[3]/total_eigs)*100, digits = 2), "%)", sep = "")
PC4 <- paste("(",round((micro_pca$eig[4]/total_eigs)*100, digits = 2), "%)", sep = "")


## Spp. colours -----------
C.car_col = "deepskyblue3"
C.car_dan = "darkblue"
C.car_don = "lightblue"
C.a.aur = "firebrick4"
C.a.gib = "darkorange"
C.caxC.a.au = "darkorchid4"
C.caxC.a.au_BKX = "white"
C.caxC.a.au_F2 = "darkgreen"
C.caxC.a.gi = "deeppink"
C.caxC.carp = "grey54"
C.carp = "black"
CaGy_BKx = "red"
CaGY_F2 ="green"
spp.cols <- c(C.a.aur, C.a.gib, C.car_col, C.car_dan, C.car_don, C.caxC.a.au, C.caxC.a.au_F2, C.caxC.a.gi, CaGY_F2, C.caxC.carp, C.carp )

## Plotting PCA -----------

## PC1 vs 2
assignments <- read.delim("../../data/Microsat_figs_final/All_spp_assignments_newcodes.txt", header= F) # Assignments based on Newhybs analysis
plot(micro_pca$li[,1], micro_pca$li[,2], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC2", PC2)) 
#text(micro_pca$li[,1], micro_pca$li[,2], assignments$V1, cex = 0.1)


points(-12,22, bg = C.car_col, pch = 21, cex = 2)
points(-12,21, bg = C.a.aur, pch = 21, cex = 2)
points(-12,20, bg = C.a.gib, pch = 21, cex = 2)
points(-12,19, bg = C.carp, pch = 21, cex = 2)
points(-12,18, bg = C.caxC.a.au, pch = 21, cex = 2)
points(-12,17, bg = C.caxC.a.gi, pch = 21, cex = 2)
points(-12,16, bg = C.caxC.carp, pch = 21, cex = 2)
points(-12,15, bg = CaGY_F2, pch = 21, cex = 2)
points(-12,14, bg = C.caxC.a.au_F2, pch = 21, cex = 2)

text(-12,22, expression(italic("C. carassius")), pos = 4, cex = 0.5)
text(-12,21, expression(italic("C. a. auratus")), pos = 4, cex = 0.5)
text(-12,20, expression(italic("C. a. gibelio")),pos = 4, cex = 0.5)
text(-12,19, expression(italic("C. carpio")),pos = 4, cex = 0.5)
text(-12,18, expression(italic("C. carassius x C. a. auratus")),pos = 4, cex = 0.5)
text(-12,17, expression(italic("C. carassius x C. a. gibelio")),pos = 4, cex = 0.5)
text(-12,16, expression(italic("C. carassius x C. carpio")),pos = 4, cex = 0.5)
text(-12,15, expression(italic("C. carassius x C. a. gibelio F2")),pos = 4, cex = 0.5)
text(-12,14, expression(italic("C. carassius x C. a. auratus F2")),pos = 4, cex = 0.5)

## PC2 vs 3
plot(micro_pca$li[,2], micro_pca$li[,3], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC2", PC2)) 
text(micro_pca$li[,2], micro_pca$li[,3], assignments$V1, cex = 0.5)

## PC3 vs 4
plot(micro_pca$li[,3], micro_pca$li[,4], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC2", PC2)) 
text(micro_pca$li[,3], micro_pca$li[,4],assignments$V1, cex = 0.5)

plot(micro_pca$li[,4], micro_pca$li[,5], col = "black", bg = spp.cols[assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC2", PC2)) 



## The DON and DANUBE clusters get in the way of the spp delimitation here - so trying without them below

## No Danube ---------------------------------------------------------------------------------------------------

## Getting rid of Danubian and DON pops 

micro_sep <- seppop(micro_data) # First separate the genind object into individual pops

NEU_assignments <- assignments[c(1:244,254:622,640:1104, 1114:1146,1157:1201,1237:1327),] # Get the associated lines of the assignments file

pop_names <- levels(pop(micro_data))
NEU_pop_names <- pop_names[c(1:13,15:25,28:59,61:65,69:75)] ## get subset of pop names

# paste together and copied and pasted them into line noted below (a bit hacky but quick)
sub_names <- paste("micro_sep$", NEU_pop_names, ", ", sep = "")
sub_names_pasted <- paste(sub_names, collapse = "")
sub_names_pasted
# line copied and pasted into:
micro_NEU <- repool(micro_sep$BOK, micro_sep$MVW, micro_sep$MVWZ, micro_sep$OKB, micro_sep$c2016, micro_sep$CCS, micro_sep$FFF, micro_sep$TL, micro_sep$NLP, micro_sep$MY20, micro_sep$GFP, micro_sep$VIIKCA, micro_sep$VIIKHY, micro_sep$CAKE, micro_sep$CALK, micro_sep$FFG, micro_sep$GD, micro_sep$PRIM, micro_sep$ROPE, micro_sep$SD, micro_sep$LUN, micro_sep$GR, micro_sep$POLEN, micro_sep$RM, micro_sep$SWED, micro_sep$MOAT, micro_sep$OST, micro_sep$BF, micro_sep$OTOM, micro_sep$UMCA, micro_sep$LMCA, micro_sep$JAR, micro_sep$H, micro_sep$M, micro_sep$FM, micro_sep$FCA, micro_sep$FCC, micro_sep$ECA, micro_sep$ECC, micro_sep$EC, micro_sep$EP, micro_sep$AL, micro_sep$EK, micro_sep$EST, micro_sep$SA, micro_sep$SK, micro_sep$TU, micro_sep$STEC, micro_sep$STYV, micro_sep$KAP, micro_sep$UKR, micro_sep$EST2, micro_sep$GFP2, micro_sep$HOLT, micro_sep$OU, micro_sep$RAIL, micro_sep$GOT, micro_sep$GOTH, micro_sep$KRA, micro_sep$LAS, micro_sep$COP, micro_sep$OBY, micro_sep$PED, micro_sep$TROM, micro_sep$WEN, micro_sep$GAM, micro_sep$BOR)

NEU_assignments$V1


### PCA ###----------------------------

Y <- scaleGen(micro_NEU, missing="mean")

micro_NEU_pca <- dudi.pca(Y, scale = F, center = T, nf = 3) # retained 4 PCs



## Species Colours

C.car_col = "deepskyblue3"
C.a.aur = "firebrick4"
C.a.gib = "darkorange"
C.caxC.a.au = "darkorchid4"
C.caxC.a.au_BKX = "white"
C.caxC.a.au_F2 = "darkgreen"
C.caxC.a.gi = "deeppink"
C.caxC.carp = "grey54"
C.carp = "black"
CaGy_BKx = "red"
CaGY_F2 ="green"
spp.cols <- c(C.a.aur, C.a.gib, C.car_col, C.caxC.a.au, C.caxC.a.gi, CaGY_F2, C.caxC.carp, C.carp )
table(NEU_assignments$V2)

## Get variation explained
total_eigs <- sum(micro_NEU_pca$eig)
PC1 <- paste("(",round((micro_NEU_pca$eig[1]/total_eigs)*100, digits = 2), "%)", sep = "")
PC2 <- paste("(",round((micro_NEU_pca$eig[2]/total_eigs)*100, digits = 2), "%)", sep = "")
PC3 <- paste("(",round((micro_NEU_pca$eig[3]/total_eigs)*100, digits = 2), "%)", sep = "")
PC4 <- paste("(",round((micro_NEU_pca$eig[4]/total_eigs)*100, digits = 2), "%)", sep = "")


## PC1 vs 2
plot(micro_NEU_pca$li[,1], micro_NEU_pca$li[,2], col = "black", bg = spp.cols[NEU_assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC2", PC2)) 
#text(micro_NEU_pca$li[,1], micro_NEU_pca$li[,2], indNames(micro_NEU), cex = .2, col = "black")

# Key --------------------------------------------------------
points(-13,24, bg = C.car_col, pch = 21, cex = 2)
points(-13,23, bg = C.a.aur, pch = 21, cex = 2)
points(-13,22, bg = C.a.gib, pch = 21, cex = 2)
points(-13,21, bg = C.carp, pch = 21, cex = 2)
points(-13,20, bg = C.caxC.a.au, pch = 21, cex = 2)
points(-13,19, bg = C.caxC.a.gi, pch = 21, cex = 2)
points(-13,18, bg = C.caxC.carp, pch = 21, cex = 2)
points(-13,17, bg = CaGY_F2, pch = 21, cex = 2)


text(-13,24, expression(italic("C. carassius")), pos = 4, cex = 0.5)
text(-13,23, expression(italic("C. a. auratus")), pos = 4, cex = 0.5)
text(-13,22, expression(italic("C. a. gibelio")),pos = 4, cex = 0.5)
text(-13,21, expression(italic("C. carpio")),pos = 4, cex = 0.5)
text(-13,20, expression(italic("C. carassius x C. a. auratus")),pos = 4, cex = 0.5)
text(-13,19, expression(italic("C. carassius x C. a. gibelio")),pos = 4, cex = 0.5)
text(-13,18, expression(italic("C. carassius x C. carpio")),pos = 4, cex = 0.5)
text(-13,17, expression(italic("C. carassius x C. gibeliio F2")),pos = 4, cex = 0.5)




## PC2 vs 3
plot(micro_NEU_pca$li[,2], micro_NEU_pca$li[,3], col = "black", bg = spp.cols[NEU_assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC1), ylab = paste("PC2", PC2)) 
#text(micro_NEU_pca$li[,2], micro_NEU_pca$li[,3], indNames(micro_NEU), cex = 0.1)


## PC3 vs 4
plot(micro_NEU_pca$li[,3], micro_NEU_pca$li[,4], col = "black", bg = spp.cols[NEU_assignments$V2], pch = 21, cex = 2, xlab = paste("PC1", PC3), ylab = paste("PC2", PC4)) 
#text(micro_NEU_pca$li[,3], micro_NEU_pca$li[,4], row.names(micro_NEU_pca$li), cex = 0.7)

# Key --------------------------------------------------------
points(6,-10, bg = C.car_col, pch = 21, cex = 2)
points(6,-11, bg = C.a.aur, pch = 21, cex = 2)
points(6,-12, bg = C.a.gib, pch = 21, cex = 2)
points(6,-13, bg = C.carp, pch = 21, cex = 2)
points(6,-14, bg = C.caxC.a.au, pch = 21, cex = 2)
points(6,-15, bg = C.caxC.a.gi, pch = 21, cex = 2)
points(6,-16, bg = C.caxC.carp, pch = 21, cex = 2)
points(6,-17, bg = CaGY_F2, pch = 21, cex = 2)


text(6,-10, expression(italic("C. carassius")), pos = 4, cex = 0.5)
text(6,-11, expression(italic("C. a. auratus")), pos = 4, cex = 0.5)
text(6,-12, expression(italic("C. a. gibelio")),pos = 4, cex = 0.5)
text(6,-13, expression(italic("C. carpio")),pos = 4, cex = 0.5)
text(6,-14, expression(italic("C. carassius x C. a. auratus")),pos = 4, cex = 0.5)
text(6,-15, expression(italic("C. carassius x C. a. gibelio")),pos = 4, cex = 0.5)
text(6,-16, expression(italic("C. carassius x C. carpio")),pos = 4, cex = 0.5)
text(6,-17, expression(italic("C. carassius x C. gibelio F2")),pos = 4, cex = 0.5)


# Subsetting for the NewHybs analysis --------------------------------------------------------------------------------------------

## This section was done prior to newhybs analysis usedto colour the above PCAs
## Cut-offs chosen on the basis of clustering and a priori morphological ID information

#### Crucian vs auratus spp.complex cut-off was taken at the value of SD08 (index )
Cru_vc_aur_PC2_cut_off = (micro_NEU_pca$li[442,2]+0.01) ## + 0.01 to keep SD08 in
Cru_vc_aur <- micro_NEU_pca[which(micro_NEU_pca$li[2] < Cru_vc_aur_PC2_cut_off)]

Cru_vs_aur_names <- row.names(Cru_vc_aur$li) # 1002 fish

# plot just these
plot(Cru_vc_aur$li[,1], Cru_vc_aur$li[,2], col = pop(micro_NEU), pch = 16) 
text(Cru_vc_aur$li[,1], Cru_vc_aur$li[,2], row.names(Cru_vc_aur$li), cex = 0.4)


##### Crucian vs cyp cut-off 
Cru_vs_cyp_PC1_cut_off = (micro_NEU_pca$li[1070,1]+0.01) ## + 0.01 to keep GODO2 OUT!
Cru_vs_cyp_PC2_cut_off = (micro_NEU_pca$li[442,2]+0.01) ## + 0.01 to keep SD08 in

Cru_vs_cyp <- micro_NEU_pca[which(micro_NEU_pca$li$Axis1 > Cru_vs_cyp_PC1_cut_off | micro_NEU_pca$li$Axis2 > Cru_vs_cyp_PC2_cut_off)]
Cru_vs_cyp_PCAcoords_final <- Cru_vs_cyp$li[-617,] ## remove the weird gibel carp

Cru_vs_cyp_names <- row.names(Cru_vs_cyp_PCAcoords_final) # 1029 fish

## plot just these
plot(Cru_vs_cyp_PCAcoords_final$Axis1, Cru_vs_cyp_PCAcoords_final$Axis2, col = pop(micro_NEU), pch = 16)  # colours dont apply properly now
text(Cru_vs_cyp$li[,1], Cru_vs_cyp$li[,2], row.names(Cru_vs_cyp$li), cex = 0.4)

## ---------------------------------------------------------------------------------------------------------------------------------

## Write name to file for NewHybs analysis

write(Cru_vs_cyp_names, "Cru_Cyp_names.txt")
write(Cru_vs_aur_names, "Cru_Aur_names.txt")

#--------------------------------------------------------------------------------------------------------------------------------------

## Getting Species specific allele frequencies for NewHybs


micro_data <- read.fstat("CRU_AU_NEWHYBS_DCODED.dat", missing = "mean") ## this is the NewHybrids coded file in Fstat format
cru_au_names <- read.delim("../../Cru_Aur_names.txt", header = F)

pop(micro_data)

indNames(micro_data) <- cru_au_names$V1
as.data.frame(indNames(micro_data)) # Use these to slice the dataset

Crucian_only <- micro_data[c(245:254, 274:283, 543:553, 657:659,661:663,933:942,1057:1065),]
pop(Crucian_only) <- rep("Crucian", 56)
indNames(Crucian_only)


AU_spp_combined <- micro_data[c(823,828,819,827,824,818,837,839,842,843,844,783,798,597,595,598,596,244),]
indNames(AU_spp_combined) ## check - looks good!
pop(AU_spp_combined) <- rep("AU_spp", 18) ## Change pop names so they all belong to the "same pop" for calculating allele frequencies


#Goldfish_only <- micro_data[c(868,853,907,909,912,914,913,893,898,889,897,894,888),] ## goldfish ID'd in the PCA
#indNames(Goldfish_only) ## check - looks good!
#pop(Goldfish_only) <- rep("GOLD", 13) ## Change pop names so they all belong to the "same pop" for calculating allele frequencies#

#Gibel_only <- micro_data[c(664,665,666,668,246),] ## goldfish ID'd in the PCA
#indNames(Gibel_only) ## check - looks good!
#pop(Gibel_only) <- rep("GIBEL", 5) ## Change pop names so they all belong to the "same pop" for calculating allele frequencies



Crucian_GP <- genind2genpop(Crucian_only)  # convert subsets to genpop objects
Goldfish_GP <- genind2genpop(Goldfish_only)  
Gibel_GP <- genind2genpop(Gibel_only)
AU_spp_GP <- genind2genpop(AU_spp_combined)

## Get allele frequencies per pop
Crucian_GP$pop.names
Crucian_freqs <- makefreq(Crucian_GP, missing = NA)
Crucian_freqs$tab

Goldfish_GP$pop.names
Gold_freqs <- makefreq(Goldfish_GP, missing = NA)
Gold_freqs$tab

Gibel_GP$pop.names
Gib_freqs <- makefreq(Gibel_GP, truenames = T, missing = NA)
Gib_freqs$tab

AU_spp_GP$pop.names
AU_spp_freqs <- makefreq(AU_spp_GP, missing = NA)
AU_spp_freqs$tab

Cru_gib<- cbind(t(Crucian_freqs$tab), t(Gib_freqs$tab))
Cru_gold<- cbind(t(Crucian_freqs$tab), t(Gold_freqs$tab))
Cru_AU_spp<- cbind(t(Crucian_freqs$tab), t(AU_spp_freqs$tab))



# subset - PCA-ID'd goldfish only

freqs <- makefreq(micro_GP) ## 

head(freqs)

micro_data$ind.names <- names$V1








