# Martin et al. Gene tree of MHCI across four species of sea turtles

library(treeio) # for read.mrbayes function
library(ggtree) # for phylogeny rendering
library(viridis) # for viridis palette
library(tidyverse) # for general tidying of data
library(ape) # for various phylo functions
library(ggnewscale) # for heatmap
library(svglite)

setwd("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/data/")

# read in dataframe of presence/absence of alleles in species
df <- read.csv("allele_counts_by_species_supertypes_MHCI.csv")

# create supertype dataframe for heatmap
supertype <- df %>% select(c("allele", "supertype_k3"))

# make the haplotype column the name of the rows, otherwise ggtree won't be able to plot corresponding supertypes
supertype <- supertype %>% remove_rownames %>% column_to_rownames(var="allele")
head(supertype)

# remove supertype column and "type" column
df <- df[-c(7,8)]

# remove supertype column
df <- df[-7]

# create a list of lists for what species the alleles are found in
species_alleles_list <- df %>% 
  pivot_longer(cols = -allele) %>% 
  filter(value == 1) %>% 
  {split(.$allele, .$name)}


# read in trree
tree <- read.mrbayes("../MrBayes/CIPRES_invariantgamma/infile.nex.con.tre")

tree # this has 164 tips.

# drop the outgroup tips:
tree <- treeio::drop.tip(tree, c("KF032390_Gallus_gallus", "KF466478_Tympanuchus_cupido"))

# read in the species/allele info with groupOTU
tree <- groupOTU(tree, .node = species_alleles_list)

# Initialize tree: tree aesthetics
p1 <- ggtree(tree, # tree read in
             ladderize = TRUE,
             right = TRUE, 
             size = 0.25) + 
  geom_tree(aes(color = group)) +
  scale_color_manual(values = c(
    "#ff981a", # Cc
    "#4ca64c", # Cm
    "#ad95cf", # Dc
    "#95b7cf", #Lk
    "#ff1458")) + # multiple species 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)
p1

# add internal nodes colored by posterior probability
p2 <- p1 +
  geom_nodepoint(color = "black", fill = "black", aes(subset = prob >= 0.995), size = 1, shape = 21) + #posterior prob between 0.995- 1
  geom_nodepoint(color = "black", fill = "#bdbdbd", aes(subset = prob < 0.995 & prob >= 0.945), size = 1, shape = 21) + #posterior prob between 0.945-0.994
  geom_nodepoint(color = "black", fill = "#f0f0f0", aes(subset = prob < 0.945 & prob >= 0.895), size = 1, shape = 21) # posterior prob between 0.895-0.944
p2

p3 <- gheatmap(p2, # previous tree to use
               supertype, # variable to create a heatmap from
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
p3

# label the tips that are shared between species
p4 <- p3 + geom_tiplab(aes(subset = node %in% c(1,	2,	11,	15,	18,	34,	39,	41,	42,	44,	46,	49,	58,	63,	73,	87,	98,	119, 122)), # tip labels for TSP alleles
                       color="black", # black color
                       size = 5, # font size
                       hjust= -0.05, # how far away from tips
                       align = TRUE)

p4

ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/phylo/final/MHCI_phylogeny_11June24_dots.svg", dpi = 300, width = 15, height = 15)

## Maximum likelihood tree from CIPRES/IQTree
ml_tree <- read.tree("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/IQTree/iqtree mhci 5 August 24/MHCI_IQTree.treefile")

# drop the outgroup tips:
ml_tree <- treeio::drop.tip(ml_tree, c("KF466478.1", "KF032390.1"))

# Initialize tree: tree aesthetics
ml_p1 <- ggtree(ml_tree,
                ladderize = TRUE,
                right = TRUE, 
                size = 0.25) + 
  geom_treescale(fontsize=2,
                 linesize=1,
                 offset=4)
ml_p1


ml_p1 <- ml_p1 +
  geom_nodepoint(color = "black", fill = "#000000", aes(subset = label >= 90), size = 1, shape = 21)
ml_p1 # bootstrap labels above 90 labeled with solid black dot

ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/phylo/final/MHCI_phylogeny_IQTree_05Aug24.svg", dpi = 300, width = 15, height = 15)

