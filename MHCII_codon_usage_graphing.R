## Martin et al: codon usage distribution analyses, MHCII chr 14

library(tidyverse) # general data wrangling
library(ggplot2) # plotting
library(ggpubr) # for ggarrange function
setwd("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/combined_new_19Nov23/codon_usage/codon_usage_distributions_data")

# read in distributions. Columns are overall identical amino acids, then how many identical codons are expected under coancestry (CAd, for coancestry distribution) and convergent evolution (CEd for convergent evolution distribtuion)
Lk_Cc_dist <- read.csv("Lk_Cc_codon_usage_distributions_MHCII.csv")
Cm_Lk_dist <- read.csv("Cm_Lk_codon_usage_distributions_MHCII.csv")
Cc_Cm_dist <- read.csv("Cc_Cm_codon_usage_distributions_MHCII.csv")
Dc_Lk_dist <- read.csv("Dc_Lk_codon_usage_distributions_MHCII.csv")
Dc_Cc_dist <- read.csv("Dc_Cc_codon_usage_distributions_MHCII.csv")
Dc_Cm_dist <- read.csv("Dc_Cm_codon_usage_distributions_MHCII.csv")


# for each of these files, divide coancestry_distribution by identical amino acids and convergent evolution distribution by identical amino acids

Lk_Cc_dist$CAd_prop <- (Lk_Cc_dist$coancestry_distribution / Lk_Cc_dist$identical_amino_acids)
Lk_Cc_dist$CEd_prop <- (Lk_Cc_dist$convergent_evolution_distribution / Lk_Cc_dist$identical_amino_acids)

Cm_Lk_dist$CAd_prop <- (Cm_Lk_dist$coancestry_distribution / Cm_Lk_dist$identical_amino_acids)
Cm_Lk_dist$CEd_prop <- (Cm_Lk_dist$convergent_evolution_distribution / Cm_Lk_dist$identical_amino_acids)

Cc_Cm_dist$CAd_prop <- (Cc_Cm_dist$coancestry_distribution / Cc_Cm_dist$identical_amino_acids)
Cc_Cm_dist$CEd_prop <- (Cc_Cm_dist$convergent_evolution_distribution / Cc_Cm_dist$identical_amino_acids)

Dc_Lk_dist$CAd_prop <- (Dc_Lk_dist$coancestry_distribution / Dc_Lk_dist$identical_amino_acids)
Dc_Lk_dist$CEd_prop <- (Dc_Lk_dist$convergent_evolution_distribution / Dc_Lk_dist$identical_amino_acids)

Dc_Cc_dist$CAd_prop <- (Dc_Cc_dist$coancestry_distribution / Dc_Cc_dist$identical_amino_acids)
Dc_Cc_dist$CEd_prop <- (Dc_Cc_dist$convergent_evolution_distribution / Dc_Cc_dist$identical_amino_acids)

Dc_Cm_dist$CAd_prop <- (Dc_Cm_dist$coancestry_distribution / Dc_Cm_dist$identical_amino_acids)
Dc_Cm_dist$CEd_prop <- (Dc_Cm_dist$convergent_evolution_distribution / Dc_Cm_dist$identical_amino_acids)

# read in data on number of identical amino acids and number of identical codons coding for those amino acids, in each species pair

identical_aa_identical_codons <- read.csv("ratio_identical_codons_identical_amino_acids_per_spp_pair_MHCII.csv")
identical_aa_identical_codons

# create a column that is the proportion of identical codons out of total identical amino acids

identical_aa_identical_codons$prop_identical_codons <- (identical_aa_identical_codons$identical_codons) / (identical_aa_identical_codons$identical_amino_acids)

# create histograms
## need to make the dta from the distribution dataframes long in order to graph both distributions in the same plot as a bimodal histogram

# Lk_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Lk_Cc_CAd_prop <- Lk_Cc_dist %>% select(CAd_prop)
Lk_Cc_CEd_prop <- Lk_Cc_dist %>% select(CEd_prop)

# pivot to longer

Lk_Cc_CAd_prop <- Lk_Cc_CAd_prop %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Lk_Cc_CEd_prop <- Lk_Cc_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Lk_Cc_dist_simple <- rbind(Lk_Cc_CAd_prop, Lk_Cc_CEd_prop)


# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Lk_Cc_histo <- ggplot(Lk_Cc_dist_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons[1,4], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Lk and Cc
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Lk_Cc_histo


# Cm_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cm_Lk_CAd_prop <- Cm_Lk_dist %>% select(CAd_prop)
Cm_Lk_CEd_prop <- Cm_Lk_dist %>% select(CEd_prop)

# pivot to longer

Cm_Lk_CAd_prop <- Cm_Lk_CAd_prop %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cm_Lk_CEd_prop <- Cm_Lk_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cm_Lk_dist_simple <- rbind(Cm_Lk_CAd_prop, Cm_Lk_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cm_Lk_histo <- ggplot(Cm_Lk_dist_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Chelonia mydas - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11,face="italic"),
        axis.text.x = element_text(size = 10)) +
  
  geom_vline(xintercept = identical_aa_identical_codons[1,4], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Cm_Lk_histo

# Cc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Cm_CAd_prop <- Cc_Cm_dist %>% select(CAd_prop)
Cc_Cm_CEd_prop <- Cc_Cm_dist %>% select(CEd_prop)

# pivot to longer

Cc_Cm_CAd_prop <- Cc_Cm_CAd_prop %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Cm_CEd_prop <- Cc_Cm_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Cm_dist_simple <- rbind(Cc_Cm_CAd_prop, Cc_Cm_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Cm_histo <- ggplot(Cc_Cm_dist_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Chelonia mydas") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons[3,4], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Cc
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Cc_Cm_histo

# Dc_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Lk_CAd_prop <- Dc_Lk_dist %>% select(CAd_prop)
Dc_Lk_CEd_prop <- Dc_Lk_dist %>% select(CEd_prop)

# pivot to longer

Dc_Lk_CAd_prop <- Dc_Lk_CAd_prop %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Lk_CEd_prop <- Dc_Lk_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Lk_dist_simple <- rbind(Dc_Lk_CAd_prop, Dc_Lk_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Lk_histo <- ggplot(Dc_Lk_dist_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Dermochelys coriacea - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x = element_text(size = 10)) +
  geom_vline(xintercept = identical_aa_identical_codons[4,4], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Lk_histo

# Dc_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Cc_CAd_prop <- Dc_Cc_dist %>% select(CAd_prop)
Dc_Cc_CEd_prop <- Dc_Cc_dist %>% select(CEd_prop)

# pivot to longer

Dc_Cc_CAd_prop <- Dc_Cc_CAd_prop %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Cc_CEd_prop <- Dc_Cc_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Cc_dist_simple <- rbind(Dc_Cc_CAd_prop, Dc_Cc_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Cc_histo <- ggplot(Dc_Cc_dist_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Dermochelys coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  
  geom_vline(xintercept = identical_aa_identical_codons[5,4], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Cc_histo


# Dc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Cm_CAd_prop <- Dc_Cm_dist %>% select(CAd_prop)
Dc_Cm_CEd_prop <- Dc_Cm_dist %>% select(CEd_prop)

# pivot to longer

Dc_Cm_CAd_prop <- Dc_Cm_CAd_prop %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Cm_CEd_prop <- Dc_Cm_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Cm_dist_simple <- rbind(Dc_Cm_CAd_prop, Dc_Cm_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Cm_histo <- ggplot(Dc_Cm_dist_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Chelonia mydas - Dermochelys coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons[6,4], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Cm_histo

# if I want labels, just add this to the below:
# + xlab("proportion of identical codons") + ylab("frequency") +
Lk_Cc_histo
# ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/codon_usage/Lk_Cc_codon_usage_histogram.svg")

Cm_Lk_histo
#ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/codon_usage/Cm_Lk_codon_usage_histogram.svg")

Cc_Cm_histo
# ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/codon_usage/Cc_Cm_codon_usage_histogram.svg")

Dc_Lk_histo
#ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/codon_usage/Dc_Lk_codon_usage_histogram.svg")

Dc_Cc_histo
# ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/codon_usage/Dc_Cc_codon_usage_histogram.svg")

Dc_Cm_histo
# ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/figures/codon_usage/Dc_Cm_codon_usage_histogram.svg")

# combine plots into one
codon_usage_figure <- ggpubr::ggarrange(Cc_Cm_histo,
                                        Dc_Cc_histo,
                                        Lk_Cc_histo,
                                        Dc_Cm_histo,
                                        Cm_Lk_histo,
                                        Dc_Lk_histo,
                                        ncol = 2,
                                        nrow = 3)
annotate_figure(codon_usage_figure,
                bottom = text_grob("proportion identical codons",
                                  size = 15),
                left = text_grob("frequency", rot = 90, size = 15))
#ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/figures/codon_usage/MHCII_facet_codon_usage_figure_23Jan24.svg")


## one-tailed z-test to determine if the observed proportion is greater than 1) CEd mean proportion and 2) CAd mean proportion
# x: "successes", in this case, number of identical codons in each spp comparison
# n: total trials,in this case, number of identical amino acids in each spp comparison
# x and n are coming from "codon_usage_summary_MHCI.xlsx"
# p: proportion to test against, in this case, the mean proportion of either CEd or CAd

# convergent evolution scenarios
## Cc and Cm
Cc_Cm_CEd_res <- prop.test(x = 85345,
                           n = 89635,
                           p = mean(Cc_Cm_dist$CEd_prop), 
                           correct = FALSE)
Cc_Cm_CEd_res # significantly higher than CE


## Lk and Cc
Lk_Cc_CEd_res <- prop.test(x = 8554,
                           n = 9271,
                           p = mean(Lk_Cc_dist$CEd_prop), 
                           correct = FALSE)
Lk_Cc_CEd_res # significantly higher than CE

## Cm and Lk
Cm_Lk_CEd_res <- prop.test(x = 37155,
                           n = 40527,
                           p = mean(Cm_Lk_dist$CEd_prop), 
                           correct = FALSE)
Cm_Lk_CEd_res # significantly higher than CE

## Cc and Dc
Dc_Cc_CEd_res <- prop.test(x = 7940,
                           n = 8980,
                           p = mean(Dc_Cc_dist$CEd_prop), 
                           correct = FALSE)
Dc_Cc_CEd_res # significantly higher than CE

## Cm and Dc
Dc_Cm_CEd_res <- prop.test(x = 34253,
                           n = 37760,
                           p = mean(Dc_Cm_dist$CEd_prop), 
                           correct = FALSE)
Dc_Cm_CEd_res # significantly higher than CE

## Lk and Dc
Dc_Lk_CEd_res <- prop.test(x = 3081,
                           n = 3469,
                           p = mean(Dc_Lk_dist$CEd_prop), 
                           correct = FALSE)
Dc_Lk_CEd_res # significantly higher than CE

# coancestry scenario
## Cc and Cm
Cc_Cm_CAd_res <- prop.test(x = 85345,
                           n = 89635,
                           p = mean(Cc_Cm_dist$CAd_prop), 
                           correct = FALSE)
Cc_Cm_CAd_res # significantly higher than CA


## Lk and Cc
Lk_Cc_CAd_res <- prop.test(x = 8554,
                           n = 9271,
                           p = mean(Lk_Cc_dist$CAd_prop), 
                           correct = FALSE)
Lk_Cc_CAd_res # significantly lower than CA

## Cm and Lk
Cm_Lk_CAd_res <- prop.test(x = 37155,
                           n = 40527,
                           p = mean(Cm_Lk_dist$CAd_prop), 
                           correct = FALSE)
Cm_Lk_CAd_res # significantly lower than CA

## Cc and Dc
Dc_Cc_CAd_res <- prop.test(x = 7940,
                           n = 8980,
                           p = mean(Dc_Cc_dist$CAd_prop), 
                           correct = FALSE)
Dc_Cc_CAd_res # significantly lower than CA

## Cm and Dc
Dc_Cm_CAd_res <- prop.test(x = 34253,
                           n = 37760,
                           p = mean(Dc_Cm_dist$CAd_prop), 
                           correct = FALSE)
Dc_Cm_CAd_res # significantly lower than CA

## Lk and Dc
Dc_Lk_CAd_res <- prop.test(x = 3081,
                           n = 3469,
                           p = mean(Dc_Lk_dist$CAd_prop), 
                           correct = FALSE)
Dc_Lk_CAd_res # significantly lower than CA


# MHCII
# Get medians of each distribution for co-ancestry score analysis
## co-ancestry:
Lk_Cc_medianCA_MHCII <- median(Lk_Cc_dist$CAd_prop)
Lk_Cc_medianCA_MHCII # 0.9453133

Cc_Cm_medianCA_MHCII <- median(Cc_Cm_dist$CAd_prop)
Cc_Cm_medianCA_MHCII # 0.9374923

Cm_Lk_medianCA_MHCII <- median(Cm_Lk_dist$CAd_prop)
Cm_Lk_medianCA_MHCII # 0.9268636

Dc_Cc_medianCA_MHCII <- median(Dc_Cc_dist$CAd_prop)
Dc_Cc_medianCA_MHCII # 0.963196

Dc_Cm_medianCA_MHCII <- median(Dc_Cm_dist$CAd_prop)
Dc_Cm_medianCA_MHCII # 0.9317532

Dc_Lk_medianCA_MHCII <- median(Dc_Lk_dist$CAd_prop)
Dc_Lk_medianCA_MHCII # 0.9489767

## convergent evolution
Lk_Cc_medianCE_MHCII <- median(Lk_Cc_dist$CEd_prop)
Lk_Cc_medianCE_MHCII # 0.7231151

Cc_Cm_medianCE_MHCII <- median(Cc_Cm_dist$CEd_prop)
Cc_Cm_medianCE_MHCII # 0.7369104

Cm_Lk_medianCE_MHCII <- median(Cm_Lk_dist$CEd_prop)
Cm_Lk_medianCE_MHCII # 0.7143263

Dc_Cc_medianCE_MHCII <- median(Dc_Cc_dist$CEd_prop)
Dc_Cc_medianCE_MHCII # 0.7331849

Dc_Cm_medianCE_MHCII <- median(Dc_Cm_dist$CEd_prop)
Dc_Cm_medianCE_MHCII # 0.7518008

Dc_Lk_medianCE_MHCII <- median(Dc_Lk_dist$CEd_prop)
Dc_Lk_medianCE_MHCII # 0.7160565

