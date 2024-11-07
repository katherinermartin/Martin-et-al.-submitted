# Martin et al: population genetic diversity metrics on MHC class I, MHC class II monomorphic, and MHC class II polymorphic alleles in C. caretta, C. mydas, D. coriacea, and L. kempii.

library(tidyverse)
library(ggplot2) # for plotting
library(ggpubr) # for faceting and annotating graphs
library(ggbeeswarm) # for the jitter/dot plots

allelic_summary_stats <- read.csv("combined_MHCI_MHCII_allelic_div_stats_v4.csv")
head(allelic_summary_stats)
allelic_summary_stats$species <- as.factor(allelic_summary_stats$species)

# set factor level for species; necessary to
species_factor_levels <- c("Dc", "Lk", "Cc", "Cm")

alleles_controlled_for_sample <- ggplot(data = allelic_summary_stats,
                        aes(x = type,
                            y = (alleles_sampled/turtles_sampled),
                            fill = factor(species, levels = species_factor_levels))) +
  geom_bar(position = "dodge",
           stat = "identity") +
  theme_classic() +
  xlab("") +
  ylab("total alleles sampled / total individuals sampled") +
  theme(legend.position = "none",
        plot.title = element_text(size=18),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        text = element_text(size=14)) +
  scale_fill_manual(values = c("#ad95cf", # Dc
                               "#95b7cf", # Lk
                               "#ff981a", # Cc
                               "#4ca64c")) + # Cm
  coord_cartesian(ylim = c(0, 1)) + # force y limit to 1
  scale_y_continuous(breaks=c(0, 0.25, 0.50, 0.75, 1)) +
  labs(title = "C")
alleles_controlled_for_sample


nucleotide_div <- ggplot(data = allelic_summary_stats,
                         aes(x = type,
                             y = nucleotide_diversity,
                             fill = factor(species, levels = species_factor_levels))) +
  geom_bar(position = "dodge",
           stat = "identity") +
  theme_classic() +
  xlab("") +
  ylab("nucleotide diversity") +
  scale_fill_manual(values = c("#ad95cf", # Dc
                               "#95b7cf", # Lk
                               "#ff981a", # Cc
                               "#4ca64c")) + # Cm
  labs(title = "E") +
  theme(legend.position = "none",
        plot.title = element_text(size=18),
        axis.ticks.x = element_blank(),
        text = element_text(size = 14)) +
  geom_errorbar(aes(ymin = nucleotide_diversity - nuc_div_variance,
                    ymax= nucleotide_diversity + nuc_div_variance),
                width=.2,
                position = position_dodge(.9))
nucleotide_div

alleles_per_spp <- ggplot(data =
                            allelic_summary_stats,
                          aes(x = type,
                              y = alleles_sampled,
                              fill = factor(species, levels = species_factor_levels))) +
  geom_bar(position = "dodge",
           stat = "identity") +
  theme_classic() +
  xlab("") +
  ylab("total alleles sampled") +
  scale_fill_manual(values = c("#ad95cf", # Dc
                               "#95b7cf", # Lk
                               "#ff981a", # Cc
                               "#4ca64c")) + # Cm
  theme(legend.position = "none",
        plot.title = element_text(size=18),
        axis.text.x = element_blank(), # no x axis labels
        axis.ticks.x=element_blank(),
        text = element_text(size=14)) +
  labs(title = "B") +
  scale_y_continuous(breaks=c(0, 25, 50, 100, 150, 200))
alleles_per_spp


# read in a different dataframe; this is raw data on alleles per individual per spp, geom_box will calculate quartiles, avg, max, min
allele_count <- read.csv("combined_MHC_alleles_per_indiv_Cc_Cm_Dc_Lk.csv")
allele_count$species <- as.factor(allele_count$species)


set.seed(1234) # for figure repeatability; method = "pseudorandom" is well, random.
allele_count_swarm <- ggplot(data = allele_count,
                             aes(x = type,
                                 y = allele_count,
                                 fill = factor(species, levels = species_factor_levels))) +
  ggbeeswarm::geom_quasirandom(method ="pseudorandom",
                               dodge.width = 1, # to dodge/move horiz position
                               shape = 21,
                               color = "azure",
                               alpha = 0.8,
                               size = 3) +
  theme_classic() +
  xlab("") +
  ylab("alleles per individual") +
  scale_fill_manual(values = c("#ad95cf", # Dc
                               "#95b7cf", # Lk
                               "#ff981a", # Cc
                               "#4ca64c")) + # Cm
  labs(title = "D") +
  theme(legend.position = "none",
        plot.title = element_text(size=18),
        axis.ticks.x=element_blank(),
        text = element_text(size=14)) +
  scale_y_continuous(breaks=c(0, 1, 3, 5, 7, 9, 11))
allele_count_swarm


# now facet plots
allele_stats_figure <- ggpubr::ggarrange(alleles_per_spp,
                                         alleles_controlled_for_sample,
                                         allele_count_swarm,
                                         nucleotide_div,
                                         ncol = 2, nrow = 2)
allele_stats_figure

ggsave("facet_allele_stats_10Sept24.svg", dpi = 300)
