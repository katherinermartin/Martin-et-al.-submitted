# Martin et al. Phylogeny for MHCII across four species of sea turtles, both chromosome 1 and 14

library(treeio) # for read.mrbayes function
library(ggtree) # for phylogeny rendering
library(viridis) # for viridis palette
library(tidyverse) # for general tidying of data
library(ape) # for various phylo functions
library(ggnewscale) # for heatmap
library(svglite)

setwd("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/data")

# read in dataframe of presence/absence of alleles in species
df <- read.csv("allele_counts_by_species_supertypes_MHCII.csv")

# create supertype dataframe for heatmap
supertype <- df %>% select(c("allele", "supertype"))

# make the haplotype column the name of the rows, otherwise ggtree won't be able to plot corresponding supertypes
supertype <- supertype %>% remove_rownames %>% column_to_rownames(var="allele")
head(supertype)

# remove supertype column and "type" column
df <- df[-c(7,8)]

# create a list of lists for what species the alleles are found in

species_alleles_list <- df %>% 
  pivot_longer(cols = -allele) %>% 
  filter(value == 1) %>% 
  {split(.$allele, .$name)}

# read in tree
tree <- read.mrbayes("/Volumes/DissertationData/MHC_evo/NGBW-JOB-MRBAYES_XSEDE-A82033EA3F944EFFA9D94C7B8C5E5D87/MHCII_chr14_chr1_infile.nex.con.tre")

tree # this has 310 tips.

# drop the outgroup tips:
tree <- treeio::drop.tip(tree, c("HQ203710.1_Gallus_gallus", "KX020876.1_Spheniscus_humboldti"))

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
               offset = 0.5, # how close the heatmap is
               legend_title = "supertype", # legend title
               width=0.075, # width of bars
               font.size = 1, # font size for labels
               colnames_position = "top", # If there were column labels
                                colnames = FALSE, colnames_offset_y = 1) + # column names
  scale_fill_manual(values = c("#fccac5","#fa9fb5", "#a64ca6"),
                    na.value = "white", name = "supertype")


# label the tips that are shared between species; I got these tip labels by running tree@phylo[["tip.label"]] and then the tips are returned in order
p4 <- p3 + geom_tiplab(aes(subset = node %in% c(8, 9, 10, 11, 60, 62, 83, 84, 89, 98, 191, 299)), # tip labels for TSP alleles
                       color="black", # black color
                       size = 5, # font size
                       hjust= -0.05, # how far away label is from tips
                       align = TRUE)
p4

ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/figures/phylo/final/MHCII_combined_phylogeny_15Sept24_dots.svg", dpi = 300, width = 15, height = 14)



## Maximum likelihood tree from CIPRES/IQTree
ml_tree <- read.tree("/Volumes/DissertationData/MHC_evo/combined_9Sept24/IQTree/iqtree mhcii 15 sept 24 trial 3/iqtree_MHCIIcombined_15Sept24.contree")

# drop the outgroup tips:
ml_tree <- treeio::drop.tip(ml_tree, c("KX020876.1", "HQ203710.1"))

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

ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/figures/phylo/final/MHCII_combined_phylogeny_IQTree_15Sept24.svg", dpi = 300, width = 15, height = 15)


