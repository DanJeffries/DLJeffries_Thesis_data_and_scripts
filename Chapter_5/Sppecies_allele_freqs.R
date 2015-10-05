library(adegenet)

setwd("~/../Dropbox/PhD/Dans_PhD_Shared/Thesis/Chapter_4_Hybridisation_and_introgression/data/Microsats")

micro_data <- read.fstat("Car_refined_final.dat", missing = "mean") ## replace NAs with mean genotype as advised in adegenet manual for centered PCA
assignments <- read.delim("../../MS Figs/Microsat_figs_final/All_samples_assignments.txt", header= F) # Assignments based on Newhybs analysis

pop(micro_data) <- assignments$V2

micro_gp<- genind2genpop(micro_data)
freqs <- makefreq(micro_gp)

write.table(t(freqs$tab), "Allele_freqs.tsv", sep = "\t", quote = F)
