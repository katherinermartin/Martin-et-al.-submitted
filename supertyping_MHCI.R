# Martin et al. Supertyping analysis of MHCI alleles from sea turtles
# this script includes k-means clustering, 2-step cross validation, and evaluation of k x through xx at optimized values for dapc

library(adegenet)
library(ggplot2)
library(reshape2)
library(scales)
library(unikn)

setwd("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping")
# Data frame of amino acid residue values for 162 MHC alleles across 14 positively selected sites

MHCI_matrix <- read.csv("MHCI_supertype_input.csv") # load in data matrix; electrochemical values for each amino acid in the full peptide binding region, based on 
MHCI_matrix <- MHCI_matrix[,-1] # remove first column (allele identifiers) so that it can be read into find.clusters (below)

check_NAs <- is.na(MHCI_matrix) # no missing values
lapply(MHCI_matrix,class) # all columns are numeric, too.


# the following primer from the Grunwald lab was extremely helpful: http://grunwaldlab.github.io/Population_Genetics_in_R/

### K MEANS CLUSTERING

# "K" is the number of clusters to group the data into and so choosing a k is very important. First, perform K-means clustering over a number of values of K and repeat 10 times for each value to explore each k.

maxK <- 25 # argument for max.n.clust, which is an integer indicating the maximum number of clusters to be tried; find.clusters will evaluate 1 through k clusters (here, 25)

myMat <- matrix(nrow=120, ncol=maxK) # create empty matrix to be filled by the for loop below
colnames(myMat) <- 1:ncol(myMat) # give column names to matrix


for(i in 1:nrow(myMat)){
  grp <- find.clusters(MHCI_matrix, n.pca = 50, # retains 50 PCAs
                       choose.n.clust = FALSE, # FALSE = user doesn't choose cluster number
                       max.n.clust = maxK) # maximum number of clusters to be tried (defined above as 25 so it does 1-25)
  myMat[i,] <- grp$Kstat # fill matrix with Kstats from groups
}

# Visualizing k means clustering

my_df <- melt(myMat) # turn matrix into a df
colnames(my_df)[1:3] <- c("Group", "K", "BIC") # column headings for df
my_df$K <- as.factor(my_df$K) # make the K values a factor
head(my_df) # lists the groups per each K and the BIC

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1 # shows how BIC changes with the number of groups (K) that's chosen
#ggsave("figures/Kmeans_max25.pdf")

## There is no clear elbow point of minimizes BIC values, but 1-8 seems to encompass the majority of the change, so evaluate k = 1 thru 8.


### CROSS-VALIDATION
# This is a two-step validation procedure, where the number of principle components to retain is evaluated at each k 3:11
# Prior to validation, run find.clusters to get the group membership (required for the validation steps)
# Then, run a broad validation step using xvalDAPC, where 1-50 principle components are retained for 30 iterations. This will yield a principle component number that minimizes the mean square error (and is therefore a good choice for that k value)
# Finally, run a narrower validation step using xvalDAPC, but this time instead of retaining 1-50 principle components, center the calculation on the result from above. For example, if PC = 10 was the value that minimized MSE in the first validation, then run the second validation on PCs 1 to 20 (where 10 is the center of that distribution), and run for 1000 iterations.
# The result will be the number of PCs to retain at a particular k, which can then be used to run DAPC().

# NOTE: xval wouldn't run at 1 so starting at 2
# k = 2 
# Set up group membership
grp_2 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 2) # evaluate at k = 2 clusters

set.seed(1)
xval_2 <- xvalDapc(MHCI_matrix,
                   grp = grp_2$grp, # group membership should come from the "grp" variable in object grp_2
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_2[2:6] # 15, although 5 and 10 also had a RMSE of 0.

# 2_a: evaluate 15 as center of optimal distribution:

set.seed(1)
xval_2_a <- xvalDapc(MHCI_matrix, grp = grp_2$grp, n.pca = 1:30, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 30 PCAs, since 15 is the center.

xval_2_a[2:6] # 28 is optimal number of PCAs at 2 clusters for correctly predicting subsamples with the lowest error. Note: even with set.seed(), the optimal number has fluctuated between 26, 27, and 28 PCAs.

# Result: 2 clusters, 28 PCs

dapc_2 <- dapc(MHCI_matrix, grp_2$grp, n.pca = 28) # retain 28 PCs automatically; retain 1 DFs


# Visualize DAPC results

pal2 <- seecol(pal_unikn_pref, n = 2) # set up palette

k2_scatter <- scatter(dapc_2, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      col = pal2)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k2_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # none

temp99 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # none

# Visualize:
compo2 <- compoplot(dapc_2, col = pal2)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo2_3Jul24.svg")
dev.off()


##################################################
# k = 3
# Set up group membership
grp_3 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 3) # evaluate at k = 3 clusters

set.seed(1)
xval_3 <- xvalDapc(MHCI_matrix,
                   grp = grp_3$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.

xval_3[2:6] # 25 Note: even with set.seed(), the optimal number has fluctuated between 15-25 PCAs.

# 3_a: evaluate 25 as center of optimal distribution:

set.seed(1)
xval_3_a <- xvalDapc(MHCI_matrix, grp = grp_3$grp, n.pca = 1:50, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 50 PCAs, since 25 is the center.

xval_3_a[2:6] # 23 is optimal number of PCAs at 3 clusters for correctly predicting subsamples with the lowest error. Note: even with set.seed(), the optimal number fluctuates slightly.

# Result: 3 clusters, 23 PCs

dapc_3 <- dapc(MHCI_matrix, grp_3$grp, n.pca = 23) # retain 23 PCs automatically; retain 2 DFs

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
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k3_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # none

temp95 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 2

temp99 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 6

# Visualize:
compo3 <- compoplot(dapc_3, col = pal3)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo3_3Jul24.svg")
dev.off()


##################################################
# k = 4
# Set up group membership
grp_4 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 4) # evaluate at k = 4 clusters

set.seed(1)
xval_4 <- xvalDapc(MHCI_matrix,
                   grp = grp_4$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_4[2:6] # 35; fluctuates between 25 and 35.

# 4_a: evaluate 35 as center of optimal distribution:

set.seed(1)
xval_4_a <- xvalDapc(MHCI_matrix, grp = grp_4$grp, n.pca = 1:75, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 75 PCAs, since 35 is the center.

xval_4_a[2:6] # 26 is optimal number of PCAs at 4 clusters for correctly predicting subsamples with the lowest error; some fluctuation/variation.

# Result: 4 clusters, 26 PCs

dapc_4 <- dapc(MHCI_matrix, grp_4$grp, n.pca = 26) # retain 26 PCs automatically; retain 3 DFs

# Visualize DAPC results

pal4 <- seecol(pal_unikn_pref, n = 4) # set up palette

k4_scatter <- scatter(dapc_4, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      col = pal4)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k4_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 4

temp95 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 4

temp99 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 5

# Visualize:
compo4 <- compoplot(dapc_4, col = pal4)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo4_3Jul24.svg")
dev.off()

##################################################
# k = 5
# Set up group membership
grp_5 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 5) # evaluate at k = 5 clusters

set.seed(1)
xval_5 <- xvalDapc(MHCI_matrix,
                   grp = grp_5$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_5[2:6] # 25; fluctuates between 25 and 30

# 5_a: evaluate 30 as center of optimal distribution:

set.seed(1)
xval_5_a <- xvalDapc(MHCI_matrix, grp = grp_5$grp, n.pca = 1:60, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 60 PCAs, since 30 is the center.

xval_5_a[2:6] # 22 is optimal number of PCAs at 5 clusters for correctly predicting subsamples with the lowest error; some variation/fluctuation

# Result: 5 clusters, 22 PCs

dapc_5 <- dapc(MHCI_matrix, grp_5$grp, n.pca = 22) # retain 22 PCs automatically; retain 4 DFs

# Visualize DAPC results

pal5 <- seecol(pal_unikn_pref, n = 5) # set up palette

k5_scatter <- scatter(dapc_5, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      posi.da = "bottomright",
                      posi.pca = "topleft",
                      col = pal5)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k5_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 4

temp95 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 4

temp99 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 6

# Visualize:
compo5 <- compoplot(dapc_5, col = pal5)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo5_3Jul24.svg")
dev.off()

##################################################
# k = 6
# Set up group membership
grp_6 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 6) # evaluate at k = 6 clusters

set.seed(1)
xval_6 <- xvalDapc(MHCI_matrix,
                   grp = grp_6$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_6[2:6] # 25, some fluctuation

# 6_a: evaluate 25 as center of optimal distribution:

set.seed(1)
xval_6_a <- xvalDapc(MHCI_matrix, grp = grp_6$grp, n.pca = 1:50, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 50 PCAs, since 25 is the center.

xval_6_a[2:6] # 25 is optimal number of PCAs at 6 clusters for correctly predicting subsamples with the lowest error.

# Result: 6 clusters, 25 PCs

dapc_6 <- dapc(MHCI_matrix, grp_6$grp, n.pca = 25) # retain 25 PCs automatically; retain 5 DFs

# Visualize DAPC results

pal6 <- seecol(pal_unikn_pref, n = 6) # set up palette

k6_scatter <- scatter(dapc_6, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      posi.da = "bottomright",
                      posi.pca = "bottomleft",
                      col = pal6)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k6_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 1

temp95 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 1

temp99 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 1

# Visualize:
compo6 <- compoplot(dapc_6, col = pal6)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo6_3Jul24.svg")
dev.off()

##################################################
# k = 7
# Set up group membership
grp_7 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 7) # evaluate at k = 7 clusters

set.seed(1)
xval_7 <- xvalDapc(MHCI_matrix,
                   grp = grp_7$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_7[2:6] # 30

# 7_a: evaluate 30 as center of optimal distribution:

set.seed(1)
xval_7_a <- xvalDapc(MHCI_matrix, grp = grp_7$grp, n.pca = 1:60, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 60 PCAs, since 30 is the center.

xval_7_a[2:6] # 25 is optimal number of PCAs at 7 clusters for correctly predicting subsamples with the lowest error, some fluctuation

# Result: 7 clusters, 25 PCs

dapc_7 <- dapc(MHCI_matrix, grp_7$grp, n.pca = 25) # retain 25 PCs automatically; retain 6 DFs

# Visualize DAPC results

pal7 <- seecol(pal_unikn_pref, n = 7) # set up palette

k7_scatter <- scatter(dapc_7, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      posi.da = "bottomright",
                      posi.pca = "topleft",
                      col = pal7)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k7_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 1

temp95 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 1

temp99 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 1

# Visualize:
compo7 <- compoplot(dapc_7, col = pal7)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo7_3Jul24.svg")
dev.off()

##################################################
# k = 8
# Set up group membership
grp_8 <- find.clusters(MHCI_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 8) # evaluate at k = 8 clusters

set.seed(1)
xval_8 <- xvalDapc(MHCI_matrix,
                   grp = grp_8$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_8[2:6] # 20, some fluctuation though

# 8_a: evaluate 20 as center of optimal distribution:

set.seed(1)
xval_8_a <- xvalDapc(MHCI_matrix, grp = grp_8$grp, n.pca = 1:40, n.rep = 1000, parallel = "multicore", ncpus = 4L) # evaluate 1 to 40 PCAs, since 20 is the center.

xval_8_a[2:6] #  is optimal number of PCAs at 8 clusters for correctly predicting subsamples with the lowest error.

# Result: 8 clusters, 15 PCs, some fluctuation

dapc_8 <- dapc(MHCI_matrix, grp_8$grp, n.pca = 15) # retain 15 PCs automatically; retain 7 DFs

# Visualize DAPC results

pal8 <- seecol(pal_unikn_pref, n = 8) # set up palette

k8_scatter <- scatter(dapc_8, # data
                      bg = "white", # white background
                      pch = 20, # size and shape of points
                      cstar = 0, # cstar= lines btwn points, 0 for null
                      solid = 0.6,
                      cex = 1.5, #
                      clab = 0.375,
                      leg = TRUE,
                      scree.da = TRUE,
                      scree.pca = TRUE,
                      posi.da = "bottomright",
                      posi.pca = "bottomleft",
                      col = pal8)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/k8_scatter_3Jul24.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
temp90 # 11

temp95 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
temp95 # 18

temp99 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
temp99 # 32

# Visualize:
compo8 <- compoplot(dapc_8, col = pal8)
dev.copy(svg, file="/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/supertyping/figures/compo8_3Jul24.svg")
dev.off()


### K = 3 through 8 has been evaluated and 3 appears to be the optimal cluster number.
# Get a dataframe of each allele and its supertype
supertypes3_MHCI <- as.data.frame(dapc_3[["grp"]]) # gives supertype membership of each allele.
supertypes3_MHCI

allele_names <- read.csv("MHCI_supertype_input.csv")
allele_names <- allele_names[1]

MHCI_alleles_by_supertype3 <- cbind(supertypes3_MHCI, allele_names)
MHCI_alleles_by_supertype3

#write.csv(MHCI_alleles_by_supertype3, "MHCI_alleles_by_supertype_assuming_k3_30June2023.csv")