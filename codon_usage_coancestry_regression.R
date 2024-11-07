## Martin et al: codon usage distribution analyses and co-ancestry scores, MHCI and MHCII polymorphic

library(tidyverse) # general data wrangling
library(ggplot2) # plotting
library(ggpubr) # for ggarrange function
library(ape) # for cophenetic.phylo function
library(treeio) # for read.mrbayes function


# read in distributions. Columns are overall identical amino acids, then how many identical codons are expected under coancestry (CAd, for coancestry distribution) and convergent evolution (CEd for convergent evolution distribtuion)

Lk_Cc_dist_MHCI <- read.csv("Lk_Cc_codon_usage_distributions_MHCI.csv")
Cm_Lk_dist_MHCI <- read.csv("Lk_Cm_codon_usage_distributions_MHCI.csv")
Cc_Cm_dist_MHCI <- read.csv("Cc_Cm_codon_usage_distributions_MHCI.csv")
Dc_Lk_dist_MHCI <- read.csv("Dc_Lk_codon_usage_distributions_MHCI.csv")
Dc_Cc_dist_MHCI <- read.csv("Dc_Cc_codon_usage_distributions_MHCI.csv")
Dc_Cm_dist_MHCI <- read.csv("Dc_Cm_codon_usage_distributions_MHCI.csv")


# for each of these files, divide coancestry distribution by identical amino acids and convergent evolution distribution by identical amino acids to get the proportion of each

Lk_Cc_dist_MHCI$CAd_prop <- (Lk_Cc_dist_MHCI$coancestry_distribution / Lk_Cc_dist_MHCI$identical_amino_acids)
Lk_Cc_dist_MHCI$CEd_prop <- (Lk_Cc_dist_MHCI$convergent_evolution_distribution / Lk_Cc_dist_MHCI$identical_amino_acids)

Cm_Lk_dist_MHCI$CAd_prop <- (Cm_Lk_dist_MHCI$coancestry_distribution / Cm_Lk_dist_MHCI$identical_amino_acids)
Cm_Lk_dist_MHCI$CEd_prop <- (Cm_Lk_dist_MHCI$convergent_evolution_distribution / Cm_Lk_dist_MHCI$identical_amino_acids)

Cc_Cm_dist_MHCI$CAd_prop <- (Cc_Cm_dist_MHCI$coancestry_distribution / Cc_Cm_dist_MHCI$identical_amino_acids)
Cc_Cm_dist_MHCI$CEd_prop <- (Cc_Cm_dist_MHCI$convergent_evolution_distribution / Cc_Cm_dist_MHCI$identical_amino_acids)

Dc_Lk_dist_MHCI$CAd_prop <- (Dc_Lk_dist_MHCI$coancestry_distribution / Dc_Lk_dist_MHCI$identical_amino_acids)
Dc_Lk_dist_MHCI$CEd_prop <- (Dc_Lk_dist_MHCI$convergent_evolution_distribution / Dc_Lk_dist_MHCI$identical_amino_acids)

Dc_Cc_dist_MHCI$CAd_prop <- (Dc_Cc_dist_MHCI$coancestry_distribution / Dc_Cc_dist_MHCI$identical_amino_acids)
Dc_Cc_dist_MHCI$CEd_prop <- (Dc_Cc_dist_MHCI$convergent_evolution_distribution / Dc_Cc_dist_MHCI$identical_amino_acids)

Dc_Cm_dist_MHCI$CAd_prop <- (Dc_Cm_dist_MHCI$coancestry_distribution / Dc_Cm_dist_MHCI$identical_amino_acids)
Dc_Cm_dist_MHCI$CEd_prop <- (Dc_Cm_dist_MHCI$convergent_evolution_distribution / Dc_Cm_dist_MHCI$identical_amino_acids)

# read in data on number of identical amino acids and number of identical codons coding for those amino acids, in each species pair, from MCMC simulations (Lenz et al scripts)

identical_aa_identical_codons_MHCI <- read.csv("ratio_identical_codons_identical_amino_acids_per_spp_pair_MHCI.csv")
identical_aa_identical_codons_MHCI

# create a column that is the proportion of identical codons out of total identical amino acids

identical_aa_identical_codons_MHCI$prop_identical_codons <- (identical_aa_identical_codons_MHCI$identical_codons) / (identical_aa_identical_codons_MHCI$identical_amino_acids)


# create histograms
## need to make the dta from the distribution dataframes long in order to graph both distributions in the same plot as a bimodal histogram

# Lk_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Lk_Cc_CAd_prop_MHCI <- Lk_Cc_dist_MHCI %>% select(CAd_prop)
Lk_Cc_CEd_prop <- Lk_Cc_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Lk_Cc_CAd_prop_MHCI <- Lk_Cc_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Lk_Cc_CEd_prop <- Lk_Cc_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Lk_Cc_dist_MHCI_simple <- rbind(Lk_Cc_CAd_prop_MHCI, Lk_Cc_CEd_prop)
lapply(Lk_Cc_dist_MHCI_simple,class)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Lk_Cc_histo_MHCI <- ggplot(Lk_Cc_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[1,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Lk and Cc
  scale_x_continuous(limits = c(0.2, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Lk_Cc_histo_MHCI


# Cm_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cm_Lk_CAd_prop_MHCI <- Cm_Lk_dist_MHCI %>% select(CAd_prop)
Cm_Lk_CEd_prop_MHCI <- Cm_Lk_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Cm_Lk_CAd_prop_MHCI <- Cm_Lk_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cm_Lk_CEd_prop_MHCI <- Cm_Lk_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cm_Lk_dist_MHCI_simple <- rbind(Cm_Lk_CAd_prop_MHCI, Cm_Lk_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cm_Lk_histo_MHCI <- ggplot(Cm_Lk_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Chelonia mydas - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11,face="italic"),
        axis.text.x = element_text(size = 10)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[1,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.2, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Cm_Lk_histo_MHCI

# Cc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Cm_CAd_prop_MHCI <- Cc_Cm_dist_MHCI %>% select(CAd_prop)
Cc_Cm_CEd_prop_MHCI <- Cc_Cm_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Cc_Cm_CAd_prop_MHCI <- Cc_Cm_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Cm_CEd_prop_MHCI <- Cc_Cm_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Cm_dist_MHCI_simple <- rbind(Cc_Cm_CAd_prop_MHCI, Cc_Cm_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Cm_histo_MHCI <- ggplot(Cc_Cm_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Chelonia mydas") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[3,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.2, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Cc_Cm_histo_MHCI

# Dc_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Lk_CAd_prop_MHCI <- Dc_Lk_dist_MHCI %>% select(CAd_prop)
Dc_Lk_CEd_prop_MHCI <- Dc_Lk_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Dc_Lk_CAd_prop_MHCI <- Dc_Lk_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Lk_CEd_prop_MHCI <- Dc_Lk_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Lk_dist_MHCI_simple <- rbind(Dc_Lk_CAd_prop_MHCI, Dc_Lk_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Lk_histo_MHCI <- ggplot(Dc_Lk_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Dermochelys coriacea - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x = element_text(size = 10)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[4,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.20, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Lk_histo_MHCI

# Dc_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Cc_CAd_prop_MHCI <- Dc_Cc_dist_MHCI %>% select(CAd_prop)
Dc_Cc_CEd_prop_MHCI <- Dc_Cc_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Dc_Cc_CAd_prop_MHCI <- Dc_Cc_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Cc_CEd_prop_MHCI <- Dc_Cc_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Cc_dist_MHCI_simple <- rbind(Dc_Cc_CAd_prop_MHCI, Dc_Cc_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Cc_histo_MHCI <- ggplot(Dc_Cc_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Dermochelys coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[5,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.2, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Cc_histo_MHCI


# Dc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Cm_CAd_prop_MHCI <- Dc_Cm_dist_MHCI %>% select(CAd_prop)
Dc_Cm_CEd_prop_MHCI <- Dc_Cm_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Dc_Cm_CAd_prop_MHCI <- Dc_Cm_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Cm_CEd_prop_MHCI <- Dc_Cm_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Cm_dist_MHCI_simple <- rbind(Dc_Cm_CAd_prop_MHCI, Dc_Cm_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Cm_histo_MHCI <- ggplot(Dc_Cm_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Chelonia mydas - Dermochelys coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[6,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.2, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Cm_histo_MHCI


## one-tailed z-test to determine if the observed proportion is greater than 1) CEd mean proportion and 2) CAd mean proportion
# x: "successes", in this case, number of identical codons in each spp comparison
# n: total trials,in this case, number of identical amino acids in each spp comparison
# x and n are coming from "codon_usage_summary_MHCI.xlsx"
# p: proportion to test against, in this case, the mean proportion of either CEd or CAd

# convergent evolution scenarios
## Cc and Cm
Cc_Cm_CEd_res_MHCI <- prop.test(x = 24321,
                           n = 27726,
                           p = mean(Cc_Cm_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Cc_Cm_CEd_res_MHCI # significantly higher than CE


## Lk and Cc
Lk_Cc_CEd_res_MHCI <- prop.test(x = 5303,
                           n = 5794,
                           p = mean(Lk_Cc_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Lk_Cc_CEd_res_MHCI # significantly higher than CE

## Cm and Lk
Cm_Lk_CEd_res_MHCI <- prop.test(x = 14058,
                           n = 15537,
                           p = mean(Cm_Lk_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Cm_Lk_CEd_res_MHCI # significantly higher than CE

## Cc and Dc
Dc_Cc_CEd_res_MHCI <- prop.test(x = 1962,
                           n = 2098,
                           p = mean(Dc_Cc_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Dc_Cc_CEd_res_MHCI # significantly higher than CE

## Cm and Dc
Dc_Cm_CEd_res_MHCI <- prop.test(x = 6341,
                           n = 6483,
                           p = mean(Dc_Cm_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Dc_Cm_CEd_res_MHCI # significantly higher than CE

## Lk and Dc
Dc_Lk_CEd_res_MHCI <- prop.test(x = 1165,
                           n = 1197,
                           p = mean(Dc_Lk_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Dc_Lk_CEd_res_MHCI # significantly higher than CE


# coancestry scenarios
## Cc and Cm
Cc_Cm_CAd_res_MHCI <- prop.test(x = 24321,
                           n = 27726,
                           p = mean(Cc_Cm_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Cc_Cm_CAd_res_MHCI # significantly higher than CA


## Lk and Cc
Lk_Cc_CAd_res_MHCI <- prop.test(x = 5303,
                           n = 5794,
                           p = mean(Lk_Cc_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Lk_Cc_CAd_res_MHCI # significantly higher than CA

## Cm and Lk
Cm_Lk_CAd_res_MHCI <- prop.test(x = 14058,
                           n = 15537,
                           p = mean(Cm_Lk_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Cm_Lk_CAd_res_MHCI # significantly higher than CA

## Cc and Dc
Dc_Cc_CAd_res_MHCI <- prop.test(x = 1962,
                           n = 2098,
                           p = mean(Dc_Cc_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Dc_Cc_CAd_res_MHCI # significantly higher than CA

## Cm and Dc
Dc_Cm_CAd_res_MHCI <- prop.test(x = 6341,
                           n = 6483,
                           p = mean(Dc_Cm_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Dc_Cm_CAd_res_MHCI # significantly higher than CA

## Lk and Dc
Dc_Lk_CAd_res_MHCI <- prop.test(x = 1165,
                           n = 1197,
                           p = mean(Dc_Lk_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Dc_Lk_CAd_res_MHCI # significantly higher than CA

# MHCI
# Get medians of each distribution for co-ancestry score analysis
## co-ancestry:
Lk_Cc_medianCA_MHCI <- median(Lk_Cc_dist_MHCI$CAd_prop)
Lk_Cc_medianCA_MHCI # 0.8120469

Cc_Cm_medianCA_MHCI <- median(Cc_Cm_dist_MHCI$CAd_prop)
Cc_Cm_medianCA_MHCI # 0.8540359

Cm_Lk_medianCA_MHCI <- median(Cm_Lk_dist_MHCI$CAd_prop)
Cm_Lk_medianCA_MHCI # 0.859368

Dc_Cc_medianCA_MHCI <- median(Dc_Cc_dist_MHCI$CAd_prop)
Dc_Cc_medianCA_MHCI # 0.8098189

Dc_Cm_medianCA_MHCI <- median(Dc_Cm_dist_MHCI$CAd_prop)
Dc_Cm_medianCA_MHCI # 0.8654944

Dc_Lk_medianCA_MHCI <- median(Dc_Lk_dist_MHCI$CAd_prop)
Dc_Lk_medianCA_MHCI # 0.8638262

## convergent evolution
Lk_Cc_medianCE_MHCI <- median(Lk_Cc_dist_MHCI$CEd_prop)
Lk_Cc_medianCE_MHCI # 0.4837763

Cc_Cm_medianCE_MHCI <- median(Cc_Cm_dist_MHCI$CEd_prop)
Cc_Cm_medianCE_MHCI # 0.454952

Cm_Lk_medianCE_MHCI <- median(Cm_Lk_dist_MHCI$CEd_prop)
Cm_Lk_medianCE_MHCI # 0.4698462

Dc_Cc_medianCE_MHCI <- median(Dc_Cc_dist_MHCI$CEd_prop)
Dc_Cc_medianCE_MHCI # 0.2797903

Dc_Cm_medianCE_MHCI <- median(Dc_Cm_dist_MHCI$CEd_prop)
Dc_Cm_medianCE_MHCI # 0.2671603

Dc_Lk_medianCE_MHCI <- median(Dc_Lk_dist_MHCI$CEd_prop)
Dc_Lk_medianCE_MHCI # 0.267335

################### MHCII ####################

## Martin et al: codon usage distribution analyses, MHCII polymorphic

# read in distributions. Columns are overall identical amino acids, then how many identical codons are expected under coancestry (CAd, for coancestry distribution) and convergent evolution (CEd for convergent evolution distribtuion)
Lk_Cc_dist_MHCII <- read.csv("Cc_Lk_codon_usage_distributions_MHCII.csv")
Lk_Cm_dist_MHCII <- read.csv("Cm_Lk_codon_usage_distributions_MHCII.csv")
Cc_Cm_dist_MHCII <- read.csv("Cc_Cm_codon_usage_distributions_MHCII.csv")
Dc_Lk_dist_MHCII <- read.csv("Dc_Lk_codon_usage_distributions_MHCII.csv")
Dc_Cc_dist_MHCII <- read.csv("Dc_Cc_codon_usage_distributions_MHCII.csv")
Dc_Cm_dist_MHCII <- read.csv("Dc_Cm_codon_usage_distributions_MHCII.csv")


# for each of these files, divide coancestry_distribution by identical amino acids and convergent evolution distribution by identical amino acids

Lk_Cc_dist_MHCII$CAd_prop <- (Lk_Cc_dist_MHCII$coancestry_distribution / Lk_Cc_dist_MHCII$identical_amino_acids)
Lk_Cc_dist_MHCII$CEd_prop <- (Lk_Cc_dist_MHCII$convergent_evolution_distribution / Lk_Cc_dist_MHCII$identical_amino_acids)

Lk_Cm_dist_MHCII$CAd_prop <- (Lk_Cm_dist_MHCII$coancestry_distribution / Lk_Cm_dist_MHCII$identical_amino_acids)
Lk_Cm_dist_MHCII$CEd_prop <- (Lk_Cm_dist_MHCII$convergent_evolution_distribution / Lk_Cm_dist_MHCII$identical_amino_acids)

Cc_Cm_dist_MHCII$CAd_prop <- (Cc_Cm_dist_MHCII$coancestry_distribution / Cc_Cm_dist_MHCII$identical_amino_acids)
Cc_Cm_dist_MHCII$CEd_prop <- (Cc_Cm_dist_MHCII$convergent_evolution_distribution / Cc_Cm_dist_MHCII$identical_amino_acids)

Dc_Lk_dist_MHCII$CAd_prop <- (Dc_Lk_dist_MHCII$coancestry_distribution / Dc_Lk_dist_MHCII$identical_amino_acids)
Dc_Lk_dist_MHCII$CEd_prop <- (Dc_Lk_dist_MHCII$convergent_evolution_distribution / Dc_Lk_dist_MHCII$identical_amino_acids)

Dc_Cc_dist_MHCII$CAd_prop <- (Dc_Cc_dist_MHCII$coancestry_distribution / Dc_Cc_dist_MHCII$identical_amino_acids)
Dc_Cc_dist_MHCII$CEd_prop <- (Dc_Cc_dist_MHCII$convergent_evolution_distribution / Dc_Cc_dist_MHCII$identical_amino_acids)

Dc_Cm_dist_MHCII$CAd_prop <- (Dc_Cm_dist_MHCII$coancestry_distribution / Dc_Cm_dist_MHCII$identical_amino_acids)
Dc_Cm_dist_MHCII$CEd_prop <- (Dc_Cm_dist_MHCII$convergent_evolution_distribution / Dc_Cm_dist_MHCII$identical_amino_acids)

# read in data on number of identical amino acids and number of identical codons coding for those amino acids, in each species pair

identical_aa_identical_codons_MHCII <- read.csv("ratio_identical_codons_identical_amino_acids_per_spp_pair_MHCII.csv")
identical_aa_identical_codons_MHCII

# create a column that is the proportion of identical codons out of total identical amino acids

identical_aa_identical_codons_MHCII$prop_identical_codons <- (identical_aa_identical_codons_MHCII$identical_codons) / (identical_aa_identical_codons_MHCII$identical_amino_acids)

# create histograms
## need to make the dta from the distribution dataframes long in order to graph both distributions in the same plot as a bimodal histogram

# Lk_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Lk_Cc_CAd_prop_MHCII <- Lk_Cc_dist_MHCII %>% select(CAd_prop)
Lk_Cc_CEd_prop_MHCII <- Lk_Cc_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Lk_Cc_CAd_prop_MHCII <- Lk_Cc_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Lk_Cc_CEd_prop_MHCII <- Lk_Cc_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Lk_Cc_dist_MHCII_simple <- rbind(Lk_Cc_CAd_prop_MHCII, Lk_Cc_CEd_prop_MHCII)


# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Lk_Cc_histo_MHCII <- ggplot(Lk_Cc_dist_MHCII_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[1,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Lk and Cc
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Lk_Cc_histo_MHCII


# Cm_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cm_Lk_CAd_prop_MHCII <- Lk_Cm_dist_MHCII %>% select(CAd_prop)
Cm_Lk_CEd_prop_MHCII <- Lk_Cm_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Cm_Lk_CAd_prop_MHCII <- Cm_Lk_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cm_Lk_CEd_prop_MHCII <- Cm_Lk_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Lk_Cc_dist_MHCII_simple <- rbind(Cm_Lk_CAd_prop_MHCII, Cm_Lk_CEd_prop_MHCII)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cm_Lk_histo_MHCII <- ggplot(Lk_Cc_dist_MHCII_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Chelonia mydas - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11,face="italic"),
        axis.text.x = element_text(size = 10)) +
  
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[1,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Cm_Lk_histo_MHCII

# Cc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Cm_CAd_prop_MHCII <- Cc_Cm_dist_MHCII %>% select(CAd_prop)
Cc_Cm_CEd_prop_MHCII <- Cc_Cm_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Cc_Cm_CAd_prop_MHCII <- Cc_Cm_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Cm_CEd_prop_MHCII <- Cc_Cm_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Cm_dist_MHCII_simple <- rbind(Cc_Cm_CAd_prop_MHCII, Cc_Cm_CEd_prop_MHCII)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Cm_histo_MHCII <- ggplot(Cc_Cm_dist_MHCII_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Chelonia mydas") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[3,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Cc
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Cc_Cm_histo_MHCII

# Dc_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Lk_CAd_prop_MHCII <- Dc_Lk_dist_MHCII %>% select(CAd_prop)
Dc_Lk_CEd_prop_MHCII <- Dc_Lk_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Dc_Lk_CAd_prop_MHCII <- Dc_Lk_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Lk_CEd_prop_MHCII <- Dc_Lk_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Lk_dist_MHCII_simple <- rbind(Dc_Lk_CAd_prop_MHCII, Dc_Lk_CEd_prop_MHCII)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Lk_histo_MHCII <- ggplot(Dc_Lk_dist_MHCII_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Dermochelys coriacea - Lepidochelys kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x = element_text(size = 10)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[4,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(10)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Lk_histo_MHCII

# Dc_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Cc_CAd_prop_MHCII <- Dc_Cc_dist_MHCII %>% select(CAd_prop)
Dc_Cc_CEd_prop <- Dc_Cc_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Dc_Cc_CAd_prop_MHCII <- Dc_Cc_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Cc_CEd_prop <- Dc_Cc_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Cc_dist_MHCII_simple <- rbind(Dc_Cc_CAd_prop_MHCII, Dc_Cc_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Cc_histo_MHCII <- ggplot(Dc_Cc_dist_MHCII_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Caretta caretta - Dermochelys coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[5,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Cc_histo_MHCII


# Dc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Dc_Cm_CAd_prop_MHCII <- Dc_Cm_dist_MHCII %>% select(CAd_prop)
Dc_Cm_CEd_prop_MHCII <- Dc_Cm_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Dc_Cm_CAd_prop_MHCII <- Dc_Cm_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Dc_Cm_CEd_prop_MHCII <- Dc_Cm_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Dc_Cm_dist_MHCII_simple <- rbind(Dc_Cm_CAd_prop_MHCII, Dc_Cm_CEd_prop_MHCII)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Dc_Cm_histo_MHCII <- ggplot(Dc_Cm_dist_MHCII_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Chelonia mydas - Dermochelys coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size=11, face="italic"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[6,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  scale_x_continuous(limits = c(0.65, 1), breaks = scales::pretty_breaks(15)) +
  scale_y_continuous(limits = c(0, 500))
Dc_Cm_histo_MHCII

## one-tailed z-test to determine if the observed proportion is greater than 1) CEd mean proportion and 2) CAd mean proportion
# x: "successes", in this case, number of identical codons in each spp comparison
# n: total trials,in this case, number of identical amino acids in each spp comparison
# x and n are coming from "codon_usage_summary_MHCI.xlsx"
# p: proportion to test against, in this case, the mean proportion of either CEd or CAd

# convergent evolution scenarios
## Cc and Cm
Cc_Cm_CEd_res_MHCII <- prop.test(x = 85345,
                           n = 89635,
                           p = mean(Cc_Cm_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Cc_Cm_CEd_res_MHCII # significantly higher than CE


## Lk and Cc
Lk_Cc_CEd_res_MHCII <- prop.test(x = 8554,
                           n = 9271,
                           p = mean(Lk_Cc_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Lk_Cc_CEd_res_MHCII # significantly higher than CE

## Cm and Lk
Cm_Lk_CEd_res_MHCII <- prop.test(x = 37155,
                           n = 40527,
                           p = mean(Lk_Cc_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Cm_Lk_CEd_res_MHCII # significantly higher than CE

## Cc and Dc
Dc_Cc_CEd_res_MHCII <- prop.test(x = 7940,
                           n = 8980,
                           p = mean(Dc_Cc_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Dc_Cc_CEd_res_MHCII # significantly higher than CE

## Cm and Dc
Dc_Cm_CEd_res_MHCII <- prop.test(x = 34253,
                           n = 37760,
                           p = mean(Dc_Cm_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Dc_Cm_CEd_res_MHCII # significantly higher than CE

## Lk and Dc
Dc_Lk_CEd_res_MHCII <- prop.test(x = 3081,
                           n = 3469,
                           p = mean(Dc_Lk_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Dc_Lk_CEd_res_MHCII # significantly higher than CE

# coancestry scenario
## Cc and Cm
Cc_Cm_CAd_res_MHCII <- prop.test(x = 85345,
                           n = 89635,
                           p = mean(Cc_Cm_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Cc_Cm_CAd_res_MHCII # significantly higher than CA


## Lk and Cc
Lk_Cc_CAd_res_MHCII <- prop.test(x = 8554,
                           n = 9271,
                           p = mean(Lk_Cc_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Lk_Cc_CAd_res_MHCII # significantly lower than CA

## Cm and Lk
Cm_Lk_CAd_res_MHCII <- prop.test(x = 37155,
                           n = 40527,
                           p = mean(Lk_Cc_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Cm_Lk_CAd_res_MHCII # significantly lower than CA

## Cc and Dc
Dc_Cc_CAd_res_MHCII <- prop.test(x = 7940,
                           n = 8980,
                           p = mean(Dc_Cc_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Dc_Cc_CAd_res_MHCII # significantly lower than CA

## Cm and Dc
Dc_Cm_CAd_res_MHCII <- prop.test(x = 34253,
                           n = 37760,
                           p = mean(Dc_Cm_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Dc_Cm_CAd_res_MHCII # significantly lower than CA

## Lk and Dc
Dc_Lk_CAd_res_MHCII <- prop.test(x = 3081,
                           n = 3469,
                           p = mean(Dc_Lk_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Dc_Lk_CAd_res_MHCII # significantly lower than CA


# MHCII
# Get medians of each distribution for co-ancestry score analysis
## co-ancestry:
Lk_Cc_medianCA_MHCII <- median(Lk_Cc_dist_MHCII$CAd_prop)
Lk_Cc_medianCA_MHCII # 0.9453133

Cc_Cm_medianCA_MHCII <- median(Cc_Cm_dist_MHCII$CAd_prop)
Cc_Cm_medianCA_MHCII # 0.9374923

Cm_Lk_medianCA_MHCII <- median(Lk_Cm_dist_MHCII$CAd_prop)
Cm_Lk_medianCA_MHCII # 0.9268636

Dc_Cc_medianCA_MHCII <- median(Dc_Cc_dist_MHCII$CAd_prop)
Dc_Cc_medianCA_MHCII # 0.963196

Dc_Cm_medianCA_MHCII <- median(Dc_Cm_dist_MHCII$CAd_prop)
Dc_Cm_medianCA_MHCII # 0.9317532

Dc_Lk_medianCA_MHCII <- median(Dc_Lk_dist_MHCII$CAd_prop)
Dc_Lk_medianCA_MHCII # 0.9489767

## convergent evolution
Lk_Cc_medianCE_MHCII <- median(Lk_Cc_dist_MHCII$CEd_prop)
Lk_Cc_medianCE_MHCII # 0.7231151

Cc_Cm_medianCE_MHCII <- median(Cc_Cm_dist_MHCII$CEd_prop)
Cc_Cm_medianCE_MHCII # 0.7369104

Cm_Lk_medianCE_MHCII <- median(Lk_Cm_dist_MHCII$CEd_prop)
Cm_Lk_medianCE_MHCII # 0.7143263

Dc_Cc_medianCE_MHCII <- median(Dc_Cc_dist_MHCII$CEd_prop)
Dc_Cc_medianCE_MHCII # 0.7331849

Dc_Cm_medianCE_MHCII <- median(Dc_Cm_dist_MHCII$CEd_prop)
Dc_Cm_medianCE_MHCII # 0.7518008

Dc_Lk_medianCE_MHCII <- median(Dc_Lk_dist_MHCII$CEd_prop)
Dc_Lk_medianCE_MHCII # 0.7160565

# CO ANCESTRY SCORE AND REGRESSION FOR EACH LOCUS
# Co-ancestry scores for each species pair, Martin et al. 2024
## after Kaesler et al. 2017
## In a pairwise fashion, visualize the contribution of co-ancestry for the observed proportion of identical codons; the Lenz et al. analysis provides the distributions, this co-ancestry score controls it by phylogenetic relatedness.
## Co-ancestry score is calculated as: ((proportion of identical codons between species) - (median of the convergent evolution distribution)) / ((median of the coancestry distribution) - (median of the convergent evolution distribution))
## then see correlation between each co-ancestry score and phylogenetic distance


# calculate phylogenetic distances

# read in BEAST tree from Phillips et al. 2022
beasttree <- read.mrbayes("run2_treeannotator_Yule_uniform.txt")
beasttree <- as.phylo(beasttree)
phylo_dist <- as.data.frame(dist.nodes((beasttree))) # compute pairwise distances between pairs of tips, internal and terminal; this is a massive matrix and we don't need distances between every tip, just between the internal nodes of the species, so trim matrix to just those.

beast_phylo_dist_trial2 <- cophenetic(beasttree)

# values
CcCm <- phylo_dist[1062, 714]
CcCm

CcLk <- phylo_dist[1062, 1225]
CcLk

CmLk <- phylo_dist[714, 1225]
CmLk

CcDc <- phylo_dist[1062, 1391]
CcDc

CmDc <- phylo_dist[714, 1391]
CmDc

DcLk <- phylo_dist[1391, 1225]
DcLk

df <- read.csv("../codon_usage_co_ancestry/coancestry_score_codon_usage_table.csv")
df$species_pair <- as.factor(df$species_pair)

# compute co-ancestry score
df$coancestry_score <- ((df$proportion_identical_codons) - (df$median_CE)) / ((df$median_CA) - (df$median_CE))

df$locus <- as.factor(df$locus)

# MHCI
MHCI_df <- df %>% filter(locus == "MHCI")
MHCI_correlation <- cor.test(MHCI_df$coancestry_score, MHCI_df$phylo_dist_beast_distnode, method = 'pearson')
MHCI_correlation # correlation: -0.4041404; p-value = 0.4268


df_MHCI <- df %>% filter(locus == "MHCI")

plot_coancestry_MHCI <- ggplot(df_MHCI,
                 aes(x = coancestry_score,
                     y = phylo_dist_beast_distnode)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") + #lm method, default 95% confidence interval
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(text = element_text(size = 14))
plot_coancestry_MHCI

# MHCII
MHCII_df <- df %>% filter(locus == "MHCII poly")
MHCII_correlation <- cor.test(MHCII_df$coancestry_score, MHCII_df$phylo_dist_beast_distnode, method = 'pearson')
MHCII_correlation # correlation: -0.4561995; p-value = 0.3632


df_MHCII <- df %>% filter(locus == "MHCII poly")

plot_coancestry_MHCII <- ggplot(df_MHCII,
                  aes(x = coancestry_score,
                      y = phylo_dist_beast_distnode)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") + #lm method, default 95% confidence interval
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(text = element_text(size = 14))
plot_coancestry_MHCII

# facet coancestry score plots into one; I'm annotating first and then faceting so that the labels A and C are equivlanet in size and arrangement to the labels B and D for the coodn usage score
plot_coancestry_MHCI <- ggpubr::ggarrange(plot_coancestry_MHCI,
                                              ncol = 1,
                                              nrow = 1,
                                              labels = "A")
plot_coancestry_MHCI

plot_coancestry_MHCII <- ggpubr::ggarrange(plot_coancestry_MHCII,
                                          ncol = 1,
                                          nrow = 1,
                                          labels = "C")
plot_coancestry_MHCII

coancestry_scores_figure <- ggpubr::ggarrange(plot_coancestry_MHCI,
                                              plot_coancestry_MHCII,
                                              ncol = 1,
                                              nrow = 2)

# annotate coancestry score plots
coancestry_scores_figure <- annotate_figure(coancestry_scores_figure,
                bottom = text_grob("co-ancestry score",
                                   size = 15),
                left = text_grob("phylogenetic distance",
                                 rot = 90,
                                 size = 15))
coancestry_scores_figure

# facet codon usage plots into one for each locus
codon_usage_figure_MHCI <- ggpubr::ggarrange(Cc_Cm_histo_MHCI,
                                             Dc_Cc_histo_MHCI,
                                             Lk_Cc_histo_MHCI,
                                             Dc_Cm_histo_MHCI,
                                             Cm_Lk_histo_MHCI,
                                             Dc_Lk_histo_MHCI,
                                             ncol = 2,
                                             nrow = 3,
                                             labels = "B")
codon_usage_figure_MHCI

codon_usage_figure_MHCII <- ggpubr::ggarrange(Cc_Cm_histo_MHCII,
                                              Dc_Cc_histo_MHCII,
                                              Lk_Cc_histo_MHCII,
                                              Dc_Cm_histo_MHCII,
                                              Cm_Lk_histo_MHCII,
                                              Dc_Lk_histo_MHCII,
                                              ncol = 2,
                                              nrow = 3,
                                              labels = "D")
codon_usage_figure_MHCII

# facet codon usage figure together
codon_usage_figure <- ggpubr::ggarrange(codon_usage_figure_MHCI,
                                        codon_usage_figure_MHCII,
                                        ncol = 1,
                                        nrow = 2)
# annotate the faceted codon usage figure
codon_usage_figure <- annotate_figure(codon_usage_figure,
                   bottom = text_grob("proportion identical codons",
                                      size = 15),
                   left = text_grob("frequency",
                                    rot = 90,
                                    size = 15))
codon_usage_figure

# finally, facet together co-ancestry plots and codon usage plots; call the figures so that the elapsed time limit isn't reached

coancestry_scores_figure
codon_usage_figure
coancestry_codon_usage_facet <- ggpubr::ggarrange(coancestry_scores_figure,
                                                  codon_usage_figure,
                                                  ncol = 2,
                                                  nrow = 2)
coancestry_codon_usage_facet
ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/combined_analyses/codon_usage_co_ancestry/figures/coancestry_codon_usage_facet_3Feb24.svg", width = 45, height = 50, units = "cm", dpi = 300)


### ignore
# try saving them all separately then faceting after in inkscape

codon_usage_figure_MHCI
ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/combined_analyses/codon_usage_co_ancestry/figures/codon_usage_figure_MHCI_3Feb24.svg", width = 20, height = 20, units = "cm", dpi = 300)

codon_usage_figure_MHCII
ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/combined_analyses/codon_usage_co_ancestry/figures/codon_usage_figure_MHCII_3Feb24.svg", width = 20, height = 20, units = "cm", dpi = 300)

coancestry_scores_figure
ggsave("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/combined_analyses/codon_usage_co_ancestry/figures/coancestry_scores_figure_3Feb24.svg", width = 20, height = 20, units = "cm", dpi = 300)
