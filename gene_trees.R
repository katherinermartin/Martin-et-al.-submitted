# Martin et al. Gene tree of MHCI and MHCII across four species of sea turtles

library(treeio) # for read.mrbayes function
library(ggtree)
packageVersion("ggtree") # version 3.99.2. Using developer's version but should be updated Oct 29 during Bioconductor release
library(ggplot2) # plotting
library(tidyverse) # for general tidying of data
library(ape) # for various phylo functions
library(ggnewscale) # for heatmap
library(svglite)
library(ggpubr) # for faceting and annotating graphs

# read in dataframe of presence/absence of alleles in species
df_MHCI <- read.csv("allele_counts_by_species_supertypes_MHCI.csv")

# create supertype dataframe for heatmap
supertype_MHCI <- df_MHCI %>% select(c("allele", "supertype_k3"))

# make the haplotype column the name of the rows, otherwise ggtree won't be able to plot corresponding supertypes
supertype_MHCI <- supertype_MHCI %>% remove_rownames %>% column_to_rownames(var = "allele")
head(supertype_MHCI)

# remove supertype column and "type" column
df_MHCI <- df_MHCI[-c(7,8)]

# remove supertype column
df_MHCI <- df_MHCI[-7]

# create a list of lists for what species the alleles are found in
species_alleles_list_MHCI <- df_MHCI %>% 
  pivot_longer(cols = -allele) %>% 
  filter(value == 1) %>% 
  {split(.$allele, .$name)}

## Bayesian gene tree

### MHCI

# read in tree
tree_MHCI <- read.mrbayes("MHCI_infile.nex.con.tre")

tree_MHCI # this has 164 tips.

# drop the outgroup tips:
tree_MHCI <- treeio::drop.tip(tree_MHCI, c("MN514000.1", "AF156658.1"))

# read in the species/allele info with groupOTU
tree_MHCI <- groupOTU(tree_MHCI, .node = species_alleles_list_MHCI)

# Initialize tree: tree aesthetics
p1_MHCI <- ggtree(tree_MHCI, # tree read in
             ladderize = TRUE,
             right = TRUE) + 
  geom_tree(aes(color = group),
            lwd = 1.5) + # line thickness
  scale_color_manual(values = c(
    "#ff981a", # Cc
    "#4ca64c", # Cm
    "#ad95cf", # Dc
    "#95b7cf", #Lk
    "#ff1458")) + # multiple species 
  geom_treescale(fontsize = 2,
                 linesize = 1, # scale bar
                 offset = 4)
p1_MHCI

# add internal nodes colored by posterior probability
p2_MHCI <- p1_MHCI +
  geom_nodepoint(color = "black", fill = "black", aes(subset = prob >= 0.995), size = 4, shape = 21) + #posterior prob between 0.995- 1
  geom_nodepoint(color = "black", fill = "#bdbdbd", aes(subset = prob < 0.995 & prob >= 0.945), size = 4, shape = 21) + #posterior prob between 0.945-0.994
  geom_nodepoint(color = "black", fill = "#f0f0f0", aes(subset = prob < 0.945 & prob >= 0.895), size = 4, shape = 21) # posterior prob between 0.895-0.944
p2_MHCI

p3_MHCI <- gheatmap(p2_MHCI, # previous tree to use
               supertype_MHCI, # variable to create a heatmap from
               color = NA, # color of border around bins
               legend_title = "supertype", # legend title
               offset = 0.2,
               width = 0.075, # width of bars
               font.size = 1, # font size for labels
               colnames_position = "top", # If there were column labels, where they'd be
               colnames = FALSE, colnames_offset_y = 1) + # column names
  scale_fill_manual(values = c("#cae7e4",
                               "#84bab1",
                               "#457a71"),
                    na.value = "white",
                    name = "supertype k3")
p3_MHCI

# label the tips that are shared between species.
# to find which node numbers correspond to which allele, run:
tree_MHCI@phylo[["tip.label"]] # returns the alleles in their node order, then match up to which alleles you want to highlight.

p4_MHCI <- p3_MHCI + geom_tiplab(aes(subset = node %in% c(6, 9, 22, 24, 25, 29, 44, 50, 68, 73, 91, 110, 123, 126, 141, 143, 144, 152)), # tip labels for TSP alleles
                       color="black", # black color
                       size = 4, # font size
                       hjust= -0.05, # how far away from tips
                       align = TRUE)
p4_MHCI

# MHCII


# read in dataframe of presence/absence of alleles in species
df_MHCII <- read.csv("allele_counts_by_species_supertypes_MHCII.csv")

# create supertype dataframe for heatmap
supertype_MHCII <- df_MHCII %>% select(c("allele", "supertype_k3"))

# make the haplotype column the name of the rows, otherwise ggtree won't be able to plot corresponding supertypes
supertype_MHCII <- supertype_MHCII %>% remove_rownames %>% column_to_rownames(var = "allele")
head(supertype_MHCII)

# remove supertype column and "type" column
df_MHCII <- df_MHCII[-c(7,8)]

# remove supertype column
df_MHCII <- df_MHCII[-7]

# create a list of lists for what species the alleles are found in
species_alleles_list_MHCII <- df_MHCII %>% 
  pivot_longer(cols = -allele) %>% 
  filter(value == 1) %>% 
  {split(.$allele, .$name)}


# read in tree
tree_MHCII <- read.mrbayes("MHCII_chr14_chr1_infile.nex.con.tre")

tree_MHCII # this has 310 tips.

# drop the outgroup tips:
tree_MHCII <- treeio::drop.tip(tree_MHCII, c("OQ473796.1", "KR535940.1_trimmed"))

# read in the species/allele info with groupOTU
tree_MHCII <- groupOTU(tree_MHCII, .node = species_alleles_list_MHCII)

# Initialize tree: tree aesthetics
p1_MHCII <- ggtree(tree_MHCII, # tree read in
                   ladderize = TRUE,
                   right = TRUE) + 
  geom_tree(aes(color = group),
            lwd = 1.5) + # line thickness
  scale_color_manual(values = c(
    "#ff981a", # Cc
    "#4ca64c", # Cm
    "#ad95cf", # Dc
    "#95b7cf", #Lk
    "#ff1458")) + # multiple species 
  geom_treescale(fontsize = 2,
                 linesize = 1, # scale bar
                 offset = 4)
p1_MHCII

# add internal nodes colored by posterior probability
p2_MHCII <- p1_MHCII +
  geom_nodepoint(color = "black", fill = "black", aes(subset = prob >= 0.995), size = 4, shape = 21) + #posterior prob between 0.995- 1
  geom_nodepoint(color = "black", fill = "#bdbdbd", aes(subset = prob < 0.995 & prob >= 0.945), size = 4, shape = 21) + #posterior prob between 0.945-0.994
  geom_nodepoint(color = "black", fill = "#f0f0f0", aes(subset = prob < 0.945 & prob >= 0.895), size = 4, shape = 21) # posterior prob between 0.895-0.944
p2_MHCII

p3_MHCII <- gheatmap(p2_MHCII, # previous tree to use
                     supertype_MHCII, # variable to create a heatmap from
                     color = NA, # color of border around bins
                     legend_title = "supertype", # legend title
                     offset = 0.2,
                     width = 0.075, # width of bars
                     font.size = 1, # font size for labels
                     colnames_position = "top", # If there were column labels, where they'd be
                     colnames = FALSE, colnames_offset_y = 1) + # column names
  scale_fill_manual(values = c("#fccac5",
                               "#fa9fb5",
                               "#a64ca6"),
                    na.value = "white",
                    name = "supertype k3")
p3_MHCII

# label the tips that are shared between species.
# to find which node numbers correspond to which allele, run:
tree_MHCII@phylo[["tip.label"]] # returns the alleles in their node order, then match up to which alleles you want to highlight.

p4_MHCII <- p3_MHCII + geom_tiplab(aes(subset = node %in% c(2, 4, 5, 14, 29, 31, 84, 85, 90, 103, 232, 288)), # tip labels for TSP alleles
                                   color="black", # black color
                                   size = 4, # font size
                                   hjust= -0.05, # how far away from tips
                                   align = TRUE)
p4_MHCII

# facet the two trees

bayesian_gene_trees <- ggpubr::ggarrange(p4_MHCI,
                                         p4_MHCII,
                                         ncol = 2, nrow = 1)
bayesian_gene_trees

ggsave("figure3_genetrees_thinner.svg", dpi = 300, width = 30, height = 15)

## Maximum likelihood trees

### MHCI

## Maximum likelihood tree from CIPRES/IQTree
ml_tree_MHCI <- read.tree("MHCI_postreview_IQtree_output.contree")

# drop the outgroup tips:
ml_tree_MHCI <- treeio::drop.tip(ml_tree_MHCI, c("MN514000.1", "AF156658.1"))

# read in the species/allele info with groupOTU
ml_tree_MHCI <- groupOTU(ml_tree_MHCI, .node = species_alleles_list_MHCI)

# Initialize tree: tree aesthetics
p1_ml_MHCI <- ggtree(ml_tree_MHCI, # tree read in
                     ladderize = TRUE,
                     right = TRUE) + 
  geom_tree(aes(color = group),
            lwd = 1.5) + # line thickness
  scale_color_manual(values = c(
    "#ff981a", # Cc
    "#4ca64c", # Cm
    "#ad95cf", # Dc
    "#95b7cf", #Lk
    "#ff1458")) + # multiple species 
  geom_treescale(fontsize = 2,
                 linesize = 1, # scale bar
                 offset = 4)
p1_ml_MHCI

p2_ml_MHCI <- p1_ml_MHCI +
  geom_nodepoint(color = "black",
                 fill = "#000000",
                 aes(subset = label >= 90),
                 size = 7,
                 shape = 21)
p2_ml_MHCI # bootstrap labels above 90 labeled with solid black dot

p3_ml_MHCI <- gheatmap(p2_ml_MHCI, # previous tree to use
                    supertype_MHCI, # variable to create a heatmap from
                    color = NA, # color of border around bins
                    legend_title = "supertype", # legend title
                    offset = 0.2,
                    width = 0.075, # width of bars
                    font.size = 1, # font size for labels
                    colnames_position = "top", # If there were column labels, where they'd be
                    colnames = FALSE, colnames_offset_y = 1) + # column names
  scale_fill_manual(values = c("#cae7e4",
                               "#84bab1",
                               "#457a71"),
                    na.value = "white",
                    name = "supertype k3")
p3_ml_MHCI

### MHCII

## Maximum likelihood tree from CIPRES/IQTree
ml_tree_MHCII <- read.tree("MHCII_combined_postreview_IQtree_output.treefile")

# drop the outgroup tips:
ml_tree_MHCII <- treeio::drop.tip(ml_tree_MHCII, c("OQ473796.1", "KR535940.1_trimmed"))

# read in the species/allele info with groupOTU
ml_tree_MHCII <- groupOTU(ml_tree_MHCII, .node = species_alleles_list_MHCII)

# Initialize tree: tree aesthetics
p1_ml_MHCII <- ggtree(ml_tree_MHCII, # tree read in
                      ladderize = TRUE,
                      right = TRUE) + 
  geom_tree(aes(color = group),
            lwd = 1.5) + # line thickness
  scale_color_manual(values = c(
    "#ff981a", # Cc
    "#4ca64c", # Cm
    "#ad95cf", # Dc
    "#95b7cf", #Lk
    "#ff1458")) + # multiple species 
  geom_treescale(fontsize = 2,
                 linesize = 1, # scale bar
                 offset = 4)
p1_ml_MHCII

p2_ml_MHCII <- p1_ml_MHCII +
  geom_nodepoint(color = "black",
                 fill = "#000000",
                 aes(subset = label >= 90),
                 size = 7,
                 shape = 21)
p2_ml_MHCII # bootstrap labels above 90 labeled with solid black dot

p3_ml_MHCII <- gheatmap(p2_ml_MHCII, # previous tree to use
                        supertype_MHCII, # variable to create a heatmap from
                        color = NA, # color of border around bins
                        legend_title = "supertype", # legend title
                        offset = 0.2,
                        width = 0.075, # width of bars
                        font.size = 1, # font size for labels
                        colnames_position = "top", # If there were column labels, where they'd be
                        colnames = FALSE, colnames_offset_y = 1) + # column names
  scale_fill_manual(values = c("#fccac5",
                               "#fa9fb5",
                               "#a64ca6"),
                    na.value = "white",
                    name = "supertype k3")
p3_ml_MHCII

# facet the two trees

ml_gene_trees <- ggpubr::ggarrange(p3_ml_MHCI,
                                   p3_ml_MHCII,
                                   ncol = 2, nrow = 1)
ml_gene_trees

ggsave("figureS5_ML_genetrees.svg", dpi = 300, width = 45, height = 15)


