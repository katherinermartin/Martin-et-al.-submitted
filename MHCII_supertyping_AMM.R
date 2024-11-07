# Martin et al. 2024: Supertyping analysis of MHC class II alleles (combined mono/chr 1 and polymorphic/chr14) from sea turtles
# this script includes k-means clustering, 2-step cross validation, and evaluation of k 3 through 8 at optimized values for dapc

library(adegenet)
library(ggplot2)
library(reshape2)
library(scales)
library(unikn)

setwd("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/combined_new_19Nov23/supertyping/combined/")

# Data frame of amino acid residue values for 308 MHC alleles across 20 positively selected sites

MHCII_matrix <- read.csv("MHCII_combined_supertyping_input.csv") # load in data matrix; electrochemical values for each amino acid in the full peptide binding region, based on 
MHCII_matrix <- MHCII_matrix[,-1] # remove first column (allele identifiers) so that it can be read into find.clusters (below)

lapply(MHCII_matrix,class) # all columns are numeric, too.

# k = 3
# Set up group membership
grp_3 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 3) # evaluate at k = 3 clusters

set.seed(1)
xval_3 <- xvalDapc(MHCII_matrix,
                   grp = grp_3$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.


# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.

xval_3[2:6] # x

# 3_a: evaluate x as center of optimal distribution:

set.seed(1)
xval_3_a <- xvalDapc(MHCII_matrix, grp = grp_3$grp, n.pca = 1:x, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to x PCAs, since x is the center.

xval_3_a[2:6] # x is optimal number of PCAs at 3 clusters for correctly predicting subsamples with the lowest error, some fluctuation

# Result: 3 clusters, x PCs

dapc_3 <- dapc(MHCII_matrix, grp_3$grp, n.pca = x) # retain x PCs automatically; retain 2 DFs

# Visualize DAPC results

pal3 <- seecol(pal_unikn_pref, n = 3) # set up palette

k3_scatter <- scatter(dapc_3, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal3)
dev.copy(svg, file="/MHCII_k3_scatter_combined_AMM_14Jul24.svg")
dev.off()

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # x

temp95 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.95))) # all those that have less than 95% probability of membership in a single cluster.
length(temp95) # x

temp99 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.99))) # all those that have less than 99% probability of membership in a single cluster.
length(temp99) # x

# Visualize:
compo3 <- compoplot(dapc_3, col = pal3)
dev.copy(svg, file="compo3_MHCHII_combined_14Jul24_AMM.svg")
dev.off()

supertypes3_MHCII_poly <- as.data.frame(dapc_3[["grp"]]) # gives supertype membership of each allele.
supertypes3_MHCII_poly

allele_names <- read.csv("MHCII_combined_supertyping_input.csv")
allele_names <- allele_names[1]

MHCII_combined_alleles_by_supertype3 <- cbind(supertypes3_MHCII_poly, allele_names)
MHCII_combined_alleles_by_supertype3

write.csv(MHCII_combined_alleles_by_supertype3, "MHCII_combined_alleles_by_supertype3_14Jul24_AMM.csv")