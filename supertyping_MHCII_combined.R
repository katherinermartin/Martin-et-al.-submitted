# Martin et al.: Supertyping analysis of MHC class II alleles (combined mono/chr 1 and polymorphic/chr14) from sea turtles
# this script includes k-means clustering, 2-step cross validation, and evaluation of k 3 through 8 at optimized values for dapc

library(adegenet)
library(ggplot2)
library(reshape2)
library(scales)
library(unikn)

# Data frame of amino acid residue values for 308 MHC alleles across 20 positively selected sites

MHCII_matrix <- read.csv("MHCII_chr14_chr1_supertyping_input.csv") # load in data matrix; electrochemical values for each amino acid in the full peptide binding region, based on 
MHCII_matrix <- MHCII_matrix[,-1] # remove first column (allele identifiers) so that it can be read into find.clusters (below)

lapply(MHCII_matrix,class) # all columns are numeric, too.


# the following primer from the Grunwald lab was extremely helpful: http://grunwaldlab.github.io/Population_Genetics_in_R/

### K MEANS CLUSTERING

# "K" is the number of clusters to group the data into and so choosing a k is very important. First, perform K-means clustering over a number of values of K and repeat 10 times for each value to explore each k.

maxK <- 25 # argument for max.n.clust, which is an integer indicating the maximum number of clusters to be tried; find.clusters will evaluate 1 through k clusters (here, 6; cannot do more than this because you'll get the error "more cluster centers than distinct data points")

myMat <- matrix(nrow=25, ncol=maxK) # create empty matrix to be filled by the for loop below
colnames(myMat) <- 1:ncol(myMat) # give column names to matrix


for(i in 1:nrow(myMat)){
  grp <- find.clusters(MHCII_matrix, n.pca = 50, # retains 50 PCAs
                       choose.n.clust = FALSE, # FALSE = user doesn't choose cluster number
                       max.n.clust = maxK) # maximum number of clusters to be tried (defined above as 14 so it does 1-14)
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
ggsave("MHCII_combined_Kmeans_max25_22Sept25.svg")

## There is no clear elbow point of minimizes BIC values, but 1-14 seems to encompass the majority of the change, so evaluate k = 2 through 5 to start


### CROSS-VALIDATION
# This is a two-step validation procedure, where the number of principle components to retain is evaluated at each k 3:8
# Prior to validation, run find.clusters to get the group membership (required for the validation steps)
# Then, run a broad validation step using xvalDAPC, where 1-50 principle components are retained for 30 iterations. This will yield a principle component number that minimizes the mean square error (and is therefore a good choice for that k value)
# Finally, run a narrower validation step using xvalDAPC, but this time instead of retaining 1-50 principle components, center the calculation on the result from above. For example, if PC = 10 was the value that minimized MSE in the first validation, then run the second validation on PCs 1 to 20 (where 10 is the center of that distribution), and run for 1000 iterations.
# The result will be the number of PCs to retain at a particular k, which can then be used to run DAPC().

# k = 2 
# Set up group membership
grp_2 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 2) # evaluate at k = 2 clusters

set.seed(1)
xval_2 <- xvalDapc(MHCII_matrix,
                   grp = grp_2$grp, # group membership should come from the "grp" variable in object grp_2
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_2[2:6] # 40 has lowest RMSE

# 2_a: evaluate 40 as center of optimal distribution:

set.seed(1)
xval_2_a <- xvalDapc(MHCII_matrix,
                     grp = grp_2$grp,
                     n.pca = 35:45,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_2_a[2:6] # 42 is optimal number of PCAs at 2 clusters for correctly predicting subsamples with the lowest error

# Result: 2 clusters, 42 PCs

dapc_2 <- dapc(MHCII_matrix,
               grp_2$grp,
               n.pca = 42) # retain 42 PCs automatically; retain 1 DFs

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
dev.copy(svg, file="MHCII_combined_k2_scatter_combined_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 4

temp95 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
length(temp95) # 6

temp99 <- which(apply(dapc_2$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
length(temp99) # 8

# Visualize:
compo2 <- compoplot(dapc_2, col = pal2)
dev.copy(svg, file="MHCII_combined_compo2_combined_22Sept25.svg")
dev.off()


##################################################
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

xval_3[2:6] # 5

# 3_a: evaluate 5 as center of optimal distribution:

set.seed(1)
xval_3_a <- xvalDapc(MHCII_matrix,
                     grp = grp_3$grp,
                     n.pca = 1:10,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_3_a[2:6] # 5 is optimal number of PCAs at 3 clusters for correctly predicting subsamples with the lowest error, some fluctuation

# Result: 3 clusters, 5 PCs

dapc_3 <- dapc(MHCII_matrix,
               grp_3$grp,
               n.pca = 5) # retain 5 PCs automatically; retain 2 DFs

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
                      posi.pca = "topleft",
                      col = pal3)
dev.copy(svg, file="MHCII_combined_k3_scatter_combined_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 3

temp95 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.95))) # all those that have less than 95% probability of membership in a single cluster.
length(temp95) # 3

temp99 <- which(apply(dapc_3$posterior,1, function(e) all(e<0.99))) # all those that have less than 99% probability of membership in a single cluster.
length(temp99) # 3

# Visualize:
compo3 <- compoplot(dapc_3, col = pal3)
dev.copy(svg, file="MHCII_combined_compo3_22Sept25.svg")
dev.off()

##################################################
# k = 4
# Set up group membership
grp_4 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 4) # evaluate at k = 4 clusters

set.seed(1)
xval_4 <- xvalDapc(MHCII_matrix,
                   grp = grp_4$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_4[2:6] # 5

# 4_a: evaluate 5 as center of optimal distribution:

set.seed(1)
xval_4_a <- xvalDapc(MHCII_matrix,
                     grp = grp_4$grp,
                     n.pca = 1:10,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_4_a[2:6] # 4 is optimal number of PCAs at 4 clusters for correctly predicting subsamples with the lowest error.

# Result: 4 clusters,  PCs

dapc_4 <- dapc(MHCII_matrix,
               grp_4$grp,
               n.pca = 4) # retain 4 PCs automatically; retain 3 DFs

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
                      posi.pca = "bottomleft",
                      posi.da = "topleft",
                      posi.leg = "bottomright",
                      col = pal4)
dev.copy(svg, file="MHCII_combined_k4_scatter_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to both. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 14

temp95 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
length(temp95) # 19

temp99 <- which(apply(dapc_4$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
length(temp99) # 38

# Visualize:
compo4 <- compoplot(dapc_4, col = pal4)
dev.copy(svg, file="MHCII_combined_compo4_22Sept25.svg")
dev.off()


##################################################
# k = 5
# Set up group membership
grp_5 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 5) # evaluate at k = 5 clusters

set.seed(1)
xval_5 <- xvalDapc(MHCII_matrix,
                   grp = grp_5$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_5[2:6] # 10

# 5_a: evaluate 10 as center of optimal distribution:

set.seed(1)
xval_5_a <- xvalDapc(MHCII_matrix,
                     grp = grp_5$grp,
                     n.pca = 5:15,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_5_a[2:6] # 8 is optimal number of PCAs at 5 clusters for correctly predicting subsamples with the lowest error, some fluctuation

# Result: 5 clusters, 8 PCs

dapc_5 <- dapc(MHCII_matrix,
               grp_5$grp,
               n.pca = 8) # retain 8 PCs automatically; retain 4 DFs

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
                      posi.pca = "topleft",
                      posi.da = "topright",
                      posi.leg = "bottomleft",
                      col = pal5)
dev.copy(svg, file="MHCII_combined_k5_scatter_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to multiple clusters. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 13

temp95 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
length(temp95) # 14

temp99 <- which(apply(dapc_5$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
length(temp99) # more than 28

# Visualize:
compo5 <- compoplot(dapc_5, col = pal5)
dev.copy(svg, file="MHCII_combined_compo5_22Sept25.svg")
dev.off()

##################################################
# k = 6
# Set up group membership
grp_6 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 6) # evaluate at k = 6 clusters

set.seed(1)
xval_6 <- xvalDapc(MHCII_matrix,
                   grp = grp_6$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.

xval_6[2:6] # 10

# 6_a: evaluate x as center of optimal distribution:

set.seed(1)
xval_6_a <- xvalDapc(MHCII_matrix,
                     grp = grp_6$grp,
                     n.pca = 5:15,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_6_a[2:6] # 8 is optimal number of PCAs at 6 clusters for correctly predicting subsamples with the lowest error

# Result: 6 clusters, 8 PCs

dapc_6 <- dapc(MHCII_matrix,
               grp_6$grp,
               n.pca = 8) # retain 8 PCs automatically; retain 5 DFs

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
                      posi.pca = "bottomleft",
                      col = pal6)
dev.copy(svg, file="MHCII_combined_k6_scatter_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to multiple clusters. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 7

temp95 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
length(temp95) # 9

temp99 <- which(apply(dapc_6$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
length(temp99) # 18

# Visualize:
compo6 <- compoplot(dapc_6, col = pal6)
dev.copy(svg, file="MHCII_combined_compo6_22Sept25.svg")
dev.off()

##################################################
# k = 7
# Set up group membership
grp_7 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 7) # evaluate at k = 7 clusters

set.seed(1)
xval_7 <- xvalDapc(MHCII_matrix,
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
xval_7_a <- xvalDapc(MHCII_matrix,
                     grp = grp_7$grp,
                     n.pca = 25:35,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_7_a[2:6] # x is optimal number of PCAs at 7 clusters for correctly predicting subsamples with the lowest error

# Result: 7 clusters, 31 PCs

dapc_7 <- dapc(MHCII_matrix,
               grp_7$grp,
               n.pca = 31) # retain 31 PCs automatically; retain 6 DFs

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
                      posi.pca = "bottomleft",
                      col = pal7)
dev.copy(svg, file="MHCII_combined_k7_scatter_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to multiple clusters. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 4

temp95 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
length(temp95) # 5

temp99 <- which(apply(dapc_7$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
length(temp99) # 9

# Visualize:
compox <- compoplot(dapc_7, col = pal7)
dev.copy(svg, file="MHCII_combined_compo7_22Sept25.svg")
dev.off()

##################################################
# k = 8
# Set up group membership
grp_8 <- find.clusters(MHCII_matrix, #input
                       n.pca = 50, # retain 50 PCs
                       n.iter = 100,
                       n.clust = 8) # evaluate at k = 8 clusters

set.seed(1)
xval_8 <- xvalDapc(MHCII_matrix,
                   grp = grp_8$grp, # group membership should come from the "grp" variable in object grp_1
                   training.set = 0.9, #the proportion of data (individuals) to be used for the training set; defaults to 0.9 if all groups have >= 10 members
                   result = "groupMean",
                   center = TRUE,
                   scale = FALSE,
                   n.pca = NULL,
                   n.rep = 100,
                   xval.plot = TRUE) # leaving the n.pca.max default so that it evaluates up to 162 PCs-- the number of alleles.

# xval gives several metrics on which number of PCs is best but these metrics aren't always congruent with one another. We used the number of PCs associated with the lowest RMSE as the "optimum" number of PCAs in the DAPC analysis.


xval_8[2:6] # 10

# 8_a: evaluate 10 as center of optimal distribution:

set.seed(1)
xval_8_a <- xvalDapc(MHCII_matrix,
                     grp = grp_8$grp,
                     n.pca = 5:15,
                     n.rep = 1000, # number of replicates at each level of PC retention
                     ncpus = 1)

xval_8_a[2:6] # 8 is optimal number of PCAs at 8 clusters for correctly predicting subsamples with the lowest error, some fluctuation

# Result: 8 clusters, 8 PCs

dapc_8 <- dapc(MHCII_matrix,
               grp_8$grp,
               n.pca = 8) # retain 8 PCs automatically; retain 7 DFs

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
                      posi.leg = "topright",
                      posi.pca = "bottomleft",
                      col = pal8)
dev.copy(svg, file="MHCII_combined_k8_scatter_22Sept25.svg")
dev.off()

# Compoplots shows us the alleles that have less than a certain probability of membership in a single cluster; i.e., those alleles that might belong to multiple clusters. For supertyping, we want to reduce this down to having NO alleles that are like this, since we're trying to break them down into the least divisible unit.

# Alleles with less than 100% probability of membership:
temp90 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.90))) # all those that have less than 90% probability of membership in a single cluster.
length(temp90) # 6

temp95 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.95))) # all those that have less than x% probability of membership in a single cluster.
length(temp95) # 15

temp99 <- which(apply(dapc_8$posterior,1, function(e) all(e<0.99))) # all those that have less than x% probability of membership in a single cluster.
length(temp99) # 24

# Visualize:
compo8 <- compoplot(dapc_8, col = pal8)
dev.copy(svg, file="compo8_combined_22Sept25.svg")
dev.off()


### K = 2 through 8 has been evaluated and 3 appears to be the optimal cluster number.
# Get a dataframe of each allele and its supertype assuming k = 3
supertypes3_MHCII_poly <- as.data.frame(dapc_3[["grp"]]) # gives supertype membership of each allele.
supertypes3_MHCII_poly

allele_names <- read.csv("MHCII_chr14_chr1_supertyping_input.csv")
allele_names <- allele_names[1]

MHCII_combined_alleles_by_supertype3 <- cbind(supertypes3_MHCII_poly, allele_names)
MHCII_combined_alleles_by_supertype3

write.csv(MHCII_combined_alleles_by_supertype3, "k3_MHCII_combined_alleles_by_supertype_23Sept25.csv")

