## Martin et al: codon usage distribution analyses and co-ancestry scores, MHCI and MHCII polymorphic

library(tidyverse) # general data wrangling
library(ggplot2)
library(ggpubr) # for ggarrange function
library(ape) # for cophenetic.phylo function
library(treeio) # for read.mrbayes function


# read in distributions. Columns are overall identical amino acids, then how many identical codons are expected under coancestry (CAd, for coancestry distribution) and convergent evolution (CEd for convergent evolution distribtuion)


Cc_Lk_dist_MHCI <- read.csv("Cc_Lk_codon_usage_distributions_MHCI.csv")
Cm_Lk_dist_MHCI <- read.csv("Cm_Lk_codon_usage_distributions_MHCI.csv")
Cc_Cm_dist_MHCI <- read.csv("Cc_Cm_codon_usage_distributions_MHCI.csv")
Dc_Lk_dist_MHCI <- read.csv("Dc_Lk_codon_usage_distributions_MHCI.csv")
Cc_Dc_dist_MHCI <- read.csv("Cc_Dc_codon_usage_distributions_MHCI.csv")
Cm_Dc_dist_MHCI <- read.csv("Cm_Dc_codon_usage_distributions_MHCI.csv")

# for each of these files, divide coancestry distribution by identical amino acids and convergent evolution distribution by identical amino acids to get the proportion of each

Cc_Lk_dist_MHCI$CAd_prop <- (Cc_Lk_dist_MHCI$coancestry_distribution / Cc_Lk_dist_MHCI$identical_amino_acids)
Cc_Lk_dist_MHCI$CEd_prop <- (Cc_Lk_dist_MHCI$convergent_evolution_distribution / Cc_Lk_dist_MHCI$identical_amino_acids)

Cm_Lk_dist_MHCI$CAd_prop <- (Cm_Lk_dist_MHCI$coancestry_distribution / Cm_Lk_dist_MHCI$identical_amino_acids)
Cm_Lk_dist_MHCI$CEd_prop <- (Cm_Lk_dist_MHCI$convergent_evolution_distribution / Cm_Lk_dist_MHCI$identical_amino_acids)

Cc_Cm_dist_MHCI$CAd_prop <- (Cc_Cm_dist_MHCI$coancestry_distribution / Cc_Cm_dist_MHCI$identical_amino_acids)
Cc_Cm_dist_MHCI$CEd_prop <- (Cc_Cm_dist_MHCI$convergent_evolution_distribution / Cc_Cm_dist_MHCI$identical_amino_acids)

Dc_Lk_dist_MHCI$CAd_prop <- (Dc_Lk_dist_MHCI$coancestry_distribution / Dc_Lk_dist_MHCI$identical_amino_acids)
Dc_Lk_dist_MHCI$CEd_prop <- (Dc_Lk_dist_MHCI$convergent_evolution_distribution / Dc_Lk_dist_MHCI$identical_amino_acids)

Cc_Dc_dist_MHCI$CAd_prop <- (Cc_Dc_dist_MHCI$coancestry_distribution / Cc_Dc_dist_MHCI$identical_amino_acids)
Cc_Dc_dist_MHCI$CEd_prop <- (Cc_Dc_dist_MHCI$convergent_evolution_distribution / Cc_Dc_dist_MHCI$identical_amino_acids)

Cm_Dc_dist_MHCI$CAd_prop <- (Cm_Dc_dist_MHCI$coancestry_distribution / Cm_Dc_dist_MHCI$identical_amino_acids)
Cm_Dc_dist_MHCI$CEd_prop <- (Cm_Dc_dist_MHCI$convergent_evolution_distribution / Cm_Dc_dist_MHCI$identical_amino_acids)

# read in data on number of identical amino acids and number of identical codons coding for those amino acids, in each species pair, from MCMC simulations (Lenz et al scripts)

identical_aa_identical_codons_MHCI <- read.csv("ratio_identical_codons_identical_amino_acids_per_spp_pair_MHCI.csv")
identical_aa_identical_codons_MHCI

# create a column that is the proportion of identical codons out of total identical amino acids

identical_aa_identical_codons_MHCI$prop_identical_codons <- (identical_aa_identical_codons_MHCI$identical_codons) / (identical_aa_identical_codons_MHCI$identical_amino_acids)


# create histograms
## need to make the dta from the distribution dataframes long in order to graph both distributions in the same plot as a bimodal histogram

# Lk_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Lk_CAd_prop_MHCI <- Cc_Lk_dist_MHCI %>% select(CAd_prop)
Cc_Lk_CEd_prop <- Cc_Lk_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Cc_Lk_CAd_prop_MHCI <- Cc_Lk_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Lk_CEd_prop <- Cc_Lk_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Lk_dist_MHCI_simple <- rbind(Cc_Lk_CAd_prop_MHCI, Cc_Lk_CEd_prop)
lapply(Cc_Lk_dist_MHCI_simple,class)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Lk_histo_MHCI <- ggplot(Cc_Lk_dist_MHCI_simple,
                           aes(x = count,
                               fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ca. caretta - L. kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 50),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[1,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Lk and Cc
  coord_cartesian(xlim = c(0.4,1.0),
                  ylim = c(0,300))
Cc_Lk_histo_MHCI

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
Cm_Lk_histo_MHCI <- ggplot(Cm_Lk_dist_MHCI_simple,
                           aes(x = count,
                               fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ch. mydas - L. kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_text(size = 50)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[1,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.4,1.0),
                  ylim = c(0,300))

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
Cc_Cm_histo_MHCI <- ggplot(Cc_Cm_dist_MHCI_simple,
                           aes(x = count,
                               fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ca. caretta - Ch. mydas") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50, face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[3,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.4,1.0),
                  ylim = c(0,300))

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
  labs(title = "D.coriacea - L. kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_text(size = 50)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[4,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.4,1.0),
                  ylim = c(0,300))

Dc_Lk_histo_MHCI

# Dc_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Dc_CAd_prop_MHCI <- Cc_Dc_dist_MHCI %>% select(CAd_prop)
Cc_Dc_CEd_prop_MHCI <- Cc_Dc_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Cc_Dc_CAd_prop_MHCI <- Cc_Dc_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Dc_CEd_prop_MHCI <- Cc_Dc_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Dc_dist_MHCI_simple <- rbind(Cc_Dc_CAd_prop_MHCI, Cc_Dc_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Dc_histo_MHCI <- ggplot(Cc_Dc_dist_MHCI_simple, aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ca. caretta - D. coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[5,5],
             linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons
  coord_cartesian(xlim = c(0.4,1.0),
                  ylim = c(0,300))

Cc_Dc_histo_MHCI


# Dc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cm_Dc_CAd_prop_MHCI <- Cm_Dc_dist_MHCI %>% select(CAd_prop)
Cm_Dc_CEd_prop_MHCI <- Cm_Dc_dist_MHCI %>% select(CEd_prop)

# pivot to longer

Cm_Dc_CAd_prop_MHCI <- Cm_Dc_CAd_prop_MHCI %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cm_Dc_CEd_prop_MHCI <- Cm_Dc_CEd_prop_MHCI %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cm_Dc_dist_MHCI_simple <- rbind(Cm_Dc_CAd_prop_MHCI, Cm_Dc_CEd_prop_MHCI)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cm_Dc_histo_MHCI <- ggplot(Cm_Dc_dist_MHCI_simple,
                           aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ch. mydas - D. coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50, face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCI[6,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.4,1.0),
                  ylim = c(0,300))

Cm_Dc_histo_MHCI


## one-tailed z-test to determine if the observed proportion is greater than 1) CEd mean proportion and 2) CAd mean proportion
# x: "successes", in this case, number of identical codons in each spp comparison
# n: total trials, in this case, number of identical amino acids in each spp comparison
# x and n are coming from "codon_usage_summary_MHCI.xlsx"
# p: proportion to test against, in this case, the mean proportion of either CEd or CAd

# convergent evolution scenarios
## Cc and Cm
Cc_Cm_CEd_res_MHCI <- prop.test(x = 19798,
                           n = 23097,
                           p = mean(Cc_Cm_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Cc_Cm_CEd_res_MHCI # significantly higher than CE
Cc_Cm_CEd_pvalue_MHCI <- Cc_Cm_CEd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# 
# data:  19798 out of 23097, null probability mean(Cc_Cm_dist_MHCI$CEd_prop)
# X-squared = 7969.1, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.5660462
# 95 percent confidence interval:
#   0.8525957 0.8616207
# sample estimates:
#   p 
# 0.8571676 

## Cc and Dc
Cc_Dc_CEd_res_MHCI <- prop.test(x = 1616,
                                n = 1752,
                                p = mean(Cc_Dc_dist_MHCI$CEd_prop), 
                                correct = FALSE)
Cc_Dc_CEd_res_MHCI # significantly higher than CE

Cc_Dc_CEd_pvalue_MHCI <- Cc_Dc_CEd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  1616 out of 1752, null probability mean(Cc_Dc_dist_MHCI$CEd_prop)
# X-squared = 1461.5, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.4667123
# 95 percent confidence interval:
#   0.9089004 0.9340003
# sample estimates:
#   p 
# 0.9223744 

## Cc and Lk
Cc_Lk_CEd_res_MHCI <- prop.test(x = 4288,
                           n = 4771,
                           p = mean(Cc_Lk_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Cc_Lk_CEd_res_MHCI # significantly higher than CE

Cc_Lk_CEd_pvalue_MHCI <- Cc_Lk_CEd_res_MHCI$p.value


# 1-sample proportions test without continuity correction
# data:  4288 out of 4771, null probability mean(Cc_Lk_dist_MHCI$CEd_prop)
# X-squared = 1790.2, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.5984842
# 95 percent confidence interval:
#   0.8898807 0.9070044
# sample estimates:
#   p 
# 0.8987634

## Cm and Dc
Cm_Dc_CEd_res_MHCI <- prop.test(x = 5080,
                                n = 5192,
                                p = mean(Cm_Dc_dist_MHCI$CEd_prop), 
                                correct = FALSE)
Cm_Dc_CEd_res_MHCI # significantly higher than CE

Cm_Dc_CEd_pvalue_MHCI <- Cm_Dc_CEd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  5080 out of 5192, null probability mean(Cm_Dc_dist_MHCI$CEd_prop)
# X-squared = 6468.4, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.4264193
# 95 percent confidence interval:
#   0.9741086 0.9820407
# sample estimates:
#   p 
# 0.9784284 

## Cm and Lk
Cm_Lk_CEd_res_MHCI <- prop.test(x = 11159,
                           n = 12583,
                           p = mean(Cm_Lk_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Cm_Lk_CEd_res_MHCI # significantly higher than CE

Cm_Lk_CEd_pvalue_MHCI <- Cm_Lk_CEd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  11159 out of 12583, null probability mean(Cm_Lk_dist_MHCI$CEd_prop)
# X-squared = 4587, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.5898626
# 95 percent confidence interval:
#   0.8811777 0.8922491
# sample estimates:
#   p 
# 0.8868314 


## Dc and Lk
Dc_Lk_CEd_res_MHCI <- prop.test(x = 936,
                           n = 968,
                           p = mean(Dc_Lk_dist_MHCI$CEd_prop), 
                           correct = FALSE)
Dc_Lk_CEd_res_MHCI # significantly higher than CE

Dc_Lk_CEd_pvalue_MHCI <- Dc_Lk_CEd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  936 out of 968, null probability mean(Dc_Lk_dist_MHCI$CEd_prop)
# X-squared = 1005.8, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.4589824
# 95 percent confidence interval:
#   0.9537054 0.9764875
# sample estimates:
#   p 
# 0.9669421 

# coancestry scenarios
## Cc and Cm
Cc_Cm_CAd_res_MHCI <- prop.test(x = 19798,
                           n = 23097,
                           p = mean(Cc_Cm_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Cc_Cm_CAd_res_MHCI # significantly higher than CA

Cc_Cm_CAd_pvalue_MHCI <- Cc_Cm_CAd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  19798 out of 23097, null probability mean(Cc_Cm_dist_MHCI$CAd_prop)
# X-squared = 172.86, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.8847886
# 95 percent confidence interval:
#   0.8525957 0.8616207
# sample estimates:
#   p 
# 0.8571676 

## Cc and Dc
Cc_Dc_CAd_res_MHCI <- prop.test(x = 1616,
                                n = 1752,
                                p = mean(Cc_Dc_dist_MHCI$CAd_prop), 
                                correct = FALSE)
Cc_Dc_CAd_res_MHCI # significantly higher than CA

Cc_Dc_CAd_pvalue_MHCI <- Cc_Dc_CAd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  1616 out of 1752, null probability mean(Cc_Dc_dist_MHCI$CAd_prop)
# X-squared = 104.76, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.830665
# 95 percent confidence interval:
#   0.9089004 0.9340003
# sample estimates:
#   p 
# 0.9223744 

## Cc and Lk
Cc_Lk_CAd_res_MHCI <- prop.test(x = 4288,
                           n = 4771,
                           p = mean(Cc_Lk_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Cc_Lk_CAd_res_MHCI # significantly higher than CA
Cc_Lk_CAd_pvalue_MHCI <- Cc_Lk_CAd_res_MHCI$p.value


# 1-sample proportions test without continuity correction
# data:  4288 out of 4771, null probability mean(Cc_Lk_dist_MHCI$CAd_prop)
# X-squared = 115.63, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.8419761
# 95 percent confidence interval:
#   0.8898807 0.9070044
# sample estimates:
#   p 
# 0.8987634 

## Cm and Dc
Cm_Dc_CAd_res_MHCI <- prop.test(x = 5080,
                                n = 5192,
                                p = mean(Cm_Dc_dist_MHCI$CAd_prop), 
                                correct = FALSE)
Cm_Dc_CAd_res_MHCI # significantly higher than CA

Cm_Dc_CAd_pvalue_MHCI <- Cm_Dc_CAd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  5080 out of 5192, null probability mean(Cm_Dc_dist_MHCI$CAd_prop)
# X-squared = 373.77, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.8968053
# 95 percent confidence interval:
#   0.9741086 0.9820407
# sample estimates:
#   p 
# 0.9784284 

## Cm and Lk
Cm_Lk_CAd_res_MHCI <- prop.test(x = 11159,
                           n = 12583,
                           p = mean(Cm_Lk_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Cm_Lk_CAd_res_MHCI # significantly higher than CA

Cm_Lk_CAd_pvalue_MHCI <- Cm_Lk_CAd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  11159 out of 12583, null probability mean(Cm_Lk_dist_MHCI$CAd_prop)
# X-squared = 5.1022, df = 1, p-value = 0.0239
# alternative hypothesis: true p is not equal to 0.8930545
# 95 percent confidence interval:
#   0.8811777 0.8922491
# sample estimates:
#   p 
# 0.8868314

## Dc and Lk
Dc_Lk_CAd_res_MHCI <- prop.test(x = 936,
                           n = 968,
                           p = mean(Dc_Lk_dist_MHCI$CAd_prop), 
                           correct = FALSE)
Dc_Lk_CAd_res_MHCI # significantly higher than CA

Dc_Lk_CAd_pvalue_MHCI <- Dc_Lk_CAd_res_MHCI$p.value

# 1-sample proportions test without continuity correction
# data:  936 out of 968, null probability mean(Dc_Lk_dist_MHCI$CAd_prop)
# X-squared = 61.453, df = 1, p-value = 4.534e-15
# alternative hypothesis: true p is not equal to 0.88725
# 95 percent confidence interval:
#   0.9537054 0.9764875
# sample estimates:
#   p 
# 0.9669421


# FDR correct the above pvalues

MHCI_pvalues <- c(Cc_Cm_CEd_pvalue_MHCI,
                       Cc_Dc_CEd_pvalue_MHCI,
                       Cc_Lk_CEd_pvalue_MHCI,
                       Cm_Dc_CEd_pvalue_MHCI,
                       Cm_Lk_CEd_pvalue_MHCI,
                       Dc_Lk_CEd_pvalue_MHCI,
                       Cc_Cm_CAd_pvalue_MHCI,
                       Cc_Dc_CAd_pvalue_MHCI,
                       Cc_Lk_CAd_pvalue_MHCI,
                       Cm_Dc_CAd_pvalue_MHCI,
                       Cm_Lk_CAd_pvalue_MHCI,
                       Dc_Lk_CAd_pvalue_MHCI
                       )
MHCI_pvalues

MHCI_fdrs <- p.adjust(MHCI_pvalues, method="fdr")

MHCI_fdrs


# Make a dataframe

species_pairs_distribution <- c("Cc_Cm_CEd",
                                "Cc_Dc_CEd",
                                "Cc_Lk_CEd",
                                "Cm_Dc_CEd",
                                "Cm_Lk_CEd",
                                "Dc_Lk_CEd",
                                "Cc_Cm_CAd",
                                "Cc_Dc_CAd",
                                "Cc_Lk_CAd",
                                "Cm_Dc_CAd",
                                "Cm_Lk_CAd",
                                "Dc_Lk_CAd")


MHCI_FDR_corrected_pvalue <- data.frame(species_pairs_distribution, MHCI_pvalues, MHCI_fdrs)
MHCI_FDR_corrected_pvalue
write.csv(MHCI_FDR_corrected_pvalue, "/Users/katiemartin/Documents/UCF/Research/MHC_species_evo/post_review_analysis/figures/supptableX_MHCI_CondonUsageAnalysis_pvalues.csv")

# MHCI
# Calculate medians of each distribution for co-ancestry score analysis; write to a csv
## co-ancestry:

Cc_Cm_medianCA_MHCI <- median(Cc_Cm_dist_MHCI$CAd_prop)
Cc_Cm_medianCA_MHCI # 0.8848768

Cc_Dc_medianCA_MHCI <- median(Cc_Dc_dist_MHCI$CAd_prop)
Cc_Dc_medianCA_MHCI # 0.8310502

Cc_Lk_medianCA_MHCI <- median(Cc_Lk_dist_MHCI$CAd_prop)
Cc_Lk_medianCA_MHCI # 0.8421715

Cm_Dc_medianCA_MHCI <- median(Cm_Dc_dist_MHCI$CAd_prop)
Cm_Dc_medianCA_MHCI # 0.8967643

Cm_Lk_medianCA_MHCI <- median(Cm_Lk_dist_MHCI$CAd_prop)
Cm_Lk_medianCA_MHCI # 0.8930303

Dc_Lk_medianCA_MHCI <- median(Dc_Lk_dist_MHCI$CAd_prop)
Dc_Lk_medianCA_MHCI # 0.8873967

## convergent evolution
Cc_Cm_medianCE_MHCI <- median(Cc_Cm_dist_MHCI$CEd_prop)
Cc_Cm_medianCE_MHCI # 0.5660908

Cc_Dc_medianCE_MHCI <- median(Cc_Dc_dist_MHCI$CEd_prop)
Cc_Dc_medianCE_MHCI # 0.466895

Cc_Lk_medianCE_MHCI <- median(Cc_Lk_dist_MHCI$CEd_prop)
Cc_Lk_medianCE_MHCI # 0.598407

Cm_Dc_medianCE_MHCI <- median(Cm_Dc_dist_MHCI$CEd_prop)
Cm_Dc_medianCE_MHCI # 0.4264253

Cm_Lk_medianCE_MHCI <- median(Cm_Lk_dist_MHCI$CEd_prop)
Cm_Lk_medianCE_MHCI # 0.5899229

Dc_Lk_medianCE_MHCI <- median(Dc_Lk_dist_MHCI$CEd_prop)
Dc_Lk_medianCE_MHCI # 0.4597107


################### MHCII ####################

## Martin et al: codon usage distribution analyses, MHCII chr 14

# read in distributions. Columns are overall identical amino acids, then how many identical codons are expected under coancestry (CAd, for coancestry distribution) and convergent evolution (CEd for convergent evolution distribution)


Cc_Lk_dist_MHCII <- read.csv("Cc_Lk_codon_usage_distributions_MHCII.csv")
Cm_Lk_dist_MHCII <- read.csv("Cm_Lk_codon_usage_distributions_MHCII.csv")
Cc_Cm_dist_MHCII <- read.csv("Cc_Cm_codon_usage_distributions_MHCII.csv")
Dc_Lk_dist_MHCII <- read.csv("Dc_Lk_codon_usage_distributions_MHCII.csv")
Cc_Dc_dist_MHCII <- read.csv("Cc_Dc_codon_usage_distributions_MHCII.csv")
Cm_Dc_dist_MHCII <- read.csv("Cm_Dc_codon_usage_distributions_MHCII.csv")


# for each of these files, divide coancestry_distribution by identical amino acids and convergent evolution distribution by identical amino acids

Cc_Lk_dist_MHCII$CAd_prop <- (Cc_Lk_dist_MHCII$coancestry_distribution / Cc_Lk_dist_MHCII$identical_amino_acids)
Cc_Lk_dist_MHCII$CEd_prop <- (Cc_Lk_dist_MHCII$convergent_evolution_distribution / Cc_Lk_dist_MHCII$identical_amino_acids)

Cm_Lk_dist_MHCII$CAd_prop <- (Cm_Lk_dist_MHCII$coancestry_distribution / Cm_Lk_dist_MHCII$identical_amino_acids)
Cm_Lk_dist_MHCII$CEd_prop <- (Cm_Lk_dist_MHCII$convergent_evolution_distribution / Cm_Lk_dist_MHCII$identical_amino_acids)

Cc_Cm_dist_MHCII$CAd_prop <- (Cc_Cm_dist_MHCII$coancestry_distribution / Cc_Cm_dist_MHCII$identical_amino_acids)
Cc_Cm_dist_MHCII$CEd_prop <- (Cc_Cm_dist_MHCII$convergent_evolution_distribution / Cc_Cm_dist_MHCII$identical_amino_acids)

Dc_Lk_dist_MHCII$CAd_prop <- (Dc_Lk_dist_MHCII$coancestry_distribution / Dc_Lk_dist_MHCII$identical_amino_acids)
Dc_Lk_dist_MHCII$CEd_prop <- (Dc_Lk_dist_MHCII$convergent_evolution_distribution / Dc_Lk_dist_MHCII$identical_amino_acids)

Cc_Dc_dist_MHCII$CAd_prop <- (Cc_Dc_dist_MHCII$coancestry_distribution / Cc_Dc_dist_MHCII$identical_amino_acids)
Cc_Dc_dist_MHCII$CEd_prop <- (Cc_Dc_dist_MHCII$convergent_evolution_distribution / Cc_Dc_dist_MHCII$identical_amino_acids)

Cm_Dc_dist_MHCII$CAd_prop <- (Cm_Dc_dist_MHCII$coancestry_distribution / Cm_Dc_dist_MHCII$identical_amino_acids)
Cm_Dc_dist_MHCII$CEd_prop <- (Cm_Dc_dist_MHCII$convergent_evolution_distribution / Cm_Dc_dist_MHCII$identical_amino_acids)

# read in data on number of identical amino acids and number of identical codons coding for those amino acids, in each species pair

identical_aa_identical_codons_MHCII <- read.csv("ratio_identical_codons_identical_amino_acids_per_spp_pair_MHCII.csv")
identical_aa_identical_codons_MHCII

# create a column that is the proportion of identical codons out of total identical amino acids

identical_aa_identical_codons_MHCII$prop_identical_codons <- (identical_aa_identical_codons_MHCII$identical_codons) / (identical_aa_identical_codons_MHCII$identical_amino_acids)

# create histograms
## need to make the dta from the distribution dataframes long in order to graph both distributions in the same plot as a bimodal histogram

# Lk_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Lk_CAd_prop_MHCII <- Cc_Lk_dist_MHCII %>% select(CAd_prop)
Cc_Lk_CEd_prop_MHCII <- Cc_Lk_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Cc_Lk_CAd_prop_MHCII <- Cc_Lk_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Lk_CEd_prop_MHCII <- Cc_Lk_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Lk_dist_MHCII_simple <- rbind(Cc_Lk_CAd_prop_MHCII, Cc_Lk_CEd_prop_MHCII)


# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Lk_histo_MHCII <- ggplot(Cc_Lk_dist_MHCII_simple,
                            aes(x = count,
                                fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ca. caretta - L. kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[1,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Lk and Cc
  coord_cartesian(xlim = c(0.7,1.0),
                  ylim = c(0, 500)) +
  scale_y_continuous(breaks = c(0, 100, 300, 500))
Cc_Lk_histo_MHCII


# Cm_Lk comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cm_Lk_CAd_prop_MHCII <- Cm_Lk_dist_MHCII %>% select(CAd_prop)
Cm_Lk_CEd_prop_MHCII <- Cm_Lk_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Cm_Lk_CAd_prop_MHCII <- Cm_Lk_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cm_Lk_CEd_prop_MHCII <- Cm_Lk_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Lk_dist_MHCII_simple <- rbind(Cm_Lk_CAd_prop_MHCII, Cm_Lk_CEd_prop_MHCII)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cm_Lk_histo_MHCII <- ggplot(Cc_Lk_dist_MHCII_simple,
                            aes(x = count, fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ch. mydas - L. kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_text(size = 50)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[1,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.7,1.0),
                  ylim = c(0, 500)) +
  scale_y_continuous(breaks = c(0, 100, 300, 500))
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
Cc_Cm_histo_MHCII <- ggplot(Cc_Cm_dist_MHCII_simple,
                            aes(x = count,
                                fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ca. caretta - Ch. mydas") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[3,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Cc
  coord_cartesian(xlim = c(0.7,1.0),
                  ylim = c(0, 500)) +
  scale_y_continuous(breaks = c(0, 100, 300, 500))
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
Dc_Lk_histo_MHCII <- ggplot(Dc_Lk_dist_MHCII_simple,
                            aes(x = count,
                                fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "D. coriacea - L. kempii") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50,
                                  face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_text(size = 50)) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[4,5],
             linetype = 8,
             linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.7,1.0),
                  ylim = c(0, 500)) +
  scale_y_continuous(breaks = c(0, 100, 300, 500))
Dc_Lk_histo_MHCII

# Dc_Cc comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cc_Dc_CAd_prop_MHCII <- Cc_Dc_dist_MHCII %>% select(CAd_prop)
Cc_Dc_CEd_prop <- Cc_Dc_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Cc_Dc_CAd_prop_MHCII <- Cc_Dc_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cc_Dc_CEd_prop <- Cc_Dc_CEd_prop %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cc_Dc_dist_MHCII_simple <- rbind(Cc_Dc_CAd_prop_MHCII, Cc_Dc_CEd_prop)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cc_Dc_histo_MHCII <- ggplot(Cc_Dc_dist_MHCII_simple,
                            aes(x = count,
                                fill = proportion)) +
  geom_histogram(position = "identity", bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ca. caretta - D. coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50, face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[5,5],
             linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.7,1.0),
                  ylim = c(0, 500)) +
  scale_y_continuous(breaks = c(0, 100, 300, 500))
Cc_Dc_histo_MHCII


# Dc_Cm comparison
## reduce dataframes to two separate dataframes: just the CAd_prop and CEd_prop

Cm_Dc_CAd_prop_MHCII <- Cm_Dc_dist_MHCII %>% select(CAd_prop)
Cm_Dc_CEd_prop_MHCII <- Cm_Dc_dist_MHCII %>% select(CEd_prop)

# pivot to longer

Cm_Dc_CAd_prop_MHCII <- Cm_Dc_CAd_prop_MHCII %>%
  pivot_longer(CAd_prop, names_to = "proportion", values_to = "count")

Cm_Dc_CEd_prop_MHCII <- Cm_Dc_CEd_prop_MHCII %>%
  pivot_longer(CEd_prop, names_to = "proportion", values_to = "count")

# Then merge
Cm_Dc_dist_MHCII_simple <- rbind(Cm_Dc_CAd_prop_MHCII, Cm_Dc_CEd_prop_MHCII)

# Map proportion type to fill, make the bars NOT stacked, and make them semitransparent
Cm_Dc_histo_MHCII <- ggplot(Cm_Dc_dist_MHCII_simple,
                            aes(x = count,
                                fill = proportion)) +
  geom_histogram(position = "identity",
                 bins = 250) + 
  theme_classic() +
  xlab("") +
  ylab("") +
  labs(title = "Ch. mydas - D. coriacea") +
  theme(legend.position = "none",
        plot.title = element_text(size = 50, face="italic"),
        axis.text.y = element_text(size = 50),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_vline(xintercept = identical_aa_identical_codons_MHCII[6,5], linetype = 8, linewidth = 1) + # intercept = calculated proportion of identical codons between Cm and Lk
  coord_cartesian(xlim = c(0.7,1.0),
                  ylim = c(0, 500)) +
  scale_y_continuous(breaks = c(0, 100, 300, 500))
Cm_Dc_histo_MHCII

## one-tailed z-test to determine if the observed proportion is greater than 1) CEd mean proportion and 2) CAd mean proportion
# x: "successes", in this case, number of identical codons in each spp comparison
# n: total trials,in this case, number of identical amino acids in each spp comparison
# x and n are coming from "codon_usage_summary_MHCI.xlsx"
# p: proportion to test against, in this case, the mean proportion of either CEd or CAd

# convergent evolution scenarios
## 1. Cc and Cm
Cc_Cm_CEd_res_MHCII <- prop.test(x = 101797,
                           n = 101797,
                           p = mean(Cc_Cm_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Cc_Cm_CEd_res_MHCII # significantly higher than CE


Cc_Cm_CEd_pvalue_MHCII <- Cc_Cm_CEd_res_MHCII$p.value

#1-sample proportions test without continuity correction
#data:  101797 out of 101797, null probability mean(Cc_Cm_dist_MHCII$CEd_prop)
#X-squared = 33295, df = 1, p-value < 2.2e-16
#alternative hypothesis: true p is not equal to 0.7535371
#95 percent confidence interval:
#  0.9999623 1.0000000
#sample estimates:
#  p 
# 1

## 2. Cc and Dc
Cc_Dc_CEd_res_MHCII <- prop.test(x = 7940,
                                 n = 8980,
                                 p = mean(Cc_Dc_dist_MHCII$CEd_prop), 
                                 correct = FALSE)
Cc_Dc_CEd_res_MHCII # significantly higher than CE

Cc_Dc_CEd_pvalue_MHCII <- Cc_Dc_CEd_res_MHCII$p.value

#1-sample proportions test without continuity correction

# data:  7940 out of 8980, null probability mean(Cc_Dc_dist_MHCII$CEd_prop)
# X-squared = 809.48, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.7550712
# 95 percent confidence interval:
#   0.8774037 0.8906419
# sample estimates:
#   p 
# 0.8841871 

## 3. Cc and Lk
Cc_Lk_CEd_res_MHCII <- prop.test(x = 10576,
                           n = 11305,
                           p = mean(Cc_Lk_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Cc_Lk_CEd_res_MHCII # significantly higher than CE

Cc_Lk_CEd_pvalue_MHCII <- Cc_Lk_CEd_res_MHCII$p.value

# 1-sample proportions test without continuity correction
# data:  10576 out of 11305, null probability mean(Cc_Lk_dist_MHCII$CEd_prop)
# X-squared = 2090.8, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.7490633
# 95 percent confidence interval:
#   0.9308381 0.9398966
# sample estimates:
#   p 
# 0.9355153 

## 4. Cm and Dc
Cm_Dc_CEd_res_MHCII <- prop.test(x = 34253,
                                 n = 37760,
                                 p = mean(Cm_Dc_dist_MHCII$CEd_prop), 
                                 correct = FALSE)
Cm_Dc_CEd_res_MHCII # significantly higher than CE

Cm_Dc_CEd_pvalue_MHCII <- Cm_Dc_CEd_res_MHCII$p.value

# 1-sample proportions test without continuity correction
# 
# data:  34253 out of 37760, null probability mean(Cm_Dc_dist_MHCII$CEd_prop)
# X-squared = 3919.8, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.7719368
# 95 percent confidence interval:
#   0.9041547 0.9100103
# sample estimates:
#   p 
# 0.9071239 


## 5. Cm and Lk
Cm_Lk_CEd_res_MHCII <- prop.test(x = 37155,
                           n = 40527,
                           p = mean(Cc_Lk_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Cm_Lk_CEd_res_MHCII # significantly higher than CE

Cm_Lk_CEd_pvalue_MHCII <- Cm_Lk_CEd_res_MHCII$p.value

# 1-sample proportions test without continuity correction
# 
# data:  37155 out of 40527, null probability mean(Cc_Lk_dist_MHCII$CEd_prop)
# X-squared = 6065.9, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.7490633
# 95 percent confidence interval:
#   0.9140676 0.9194458
# sample estimates:
#   p 
# 0.9167962 

## 6. Dc and Lk
Dc_Lk_CEd_res_MHCII <- prop.test(x = 3081,
                           n = 3469,
                           p = mean(Dc_Lk_dist_MHCII$CEd_prop), 
                           correct = FALSE)
Dc_Lk_CEd_res_MHCII # significantly higher than CE

Dc_Lk_CEd_pvalue_MHCII <- Dc_Lk_CEd_res_MHCII$p.value

# 1-sample proportions test without continuity correction
# 
# data:  3081 out of 3469, null probability mean(Dc_Lk_dist_MHCII$CEd_prop)
# X-squared = 354.42, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.7496885
# 95 percent confidence interval:
#   0.8772316 0.8982141
# sample estimates:
#   p 
# 0.8881522 

# coancestry scenario
## 1. Cc and Cm
Cc_Cm_CAd_res_MHCII <- prop.test(x = 97605,
                           n = 101797,
                           p = mean(Cc_Cm_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Cc_Cm_CAd_res_MHCII # significantly higher than CA

Cc_Cm_CAd_pvalue_MHCII <- Cc_Cm_CAd_res_MHCII$p.value

# 1-sample proportions test without continuity correction
# 
# data:  97605 out of 101797, null probability mean(Cc_Cm_dist_MHCII$CAd_prop)
# X-squared = 603.77, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.9406188
# 95 percent confidence interval:
#   0.9575819 0.9600234
# sample estimates:
#   p 
# 0.95882 

## 2. Cc and Dc
Cc_Dc_CAd_res_MHCII <- prop.test(x = 7940,
                                 n = 8980,
                                 p = mean(Cc_Dc_dist_MHCII$CAd_prop), 
                                 correct = FALSE)
Cc_Dc_CAd_res_MHCII # significantly lower than CA

Cc_Dc_CAd_pvalue_MHCII <- Cc_Dc_CAd_res_MHCII$p.value

# 1-sample proportions test without continuity correction
# 
# data:  7940 out of 8980, null probability mean(Cc_Dc_dist_MHCII$CAd_prop)
# X-squared = 1770.7, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.965374
# 95 percent confidence interval:
#   0.8774037 0.8906419
# sample estimates:
#   p 
# 0.8841871 

## 3. Cc and Lk
Cc_Lk_CAd_res_MHCII <- prop.test(x = 8554,
                           n = 9271,
                           p = mean(Cc_Lk_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Cc_Lk_CAd_res_MHCII # significantly lower than CA

Cc_Lk_CAd_pvalue_MHCII <- Cc_Lk_CAd_res_MHCII$p.value


# 1-sample proportions test without continuity correction
# 
# data:  8554 out of 9271, null probability mean(Cc_Lk_dist_MHCII$CAd_prop)
# X-squared = 154.27, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.9506122
# 95 percent confidence interval:
#   0.9170478 0.9279262
# sample estimates:
#   p 
# 0.9226621

## 4. Cm and Dc
Cm_Dc_CAd_res_MHCII <- prop.test(x = 34253,
                                 n = 37760,
                                 p = mean(Cm_Dc_dist_MHCII$CAd_prop), 
                                 correct = FALSE)
Cm_Dc_CAd_res_MHCII # significantly lower than CA

Cm_Dc_CAd_pvalue_MHCII <- Cm_Dc_CAd_res_MHCII$p.value


# 1-sample proportions test without continuity correction
# 
# data:  34253 out of 37760, null probability mean(Cm_Dc_dist_MHCII$CAd_prop)
# X-squared = 467.45, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.9346264
# 95 percent confidence interval:
#   0.9041547 0.9100103
# sample estimates:
#   p 
# 0.9071239 

## 5. Cm and Lk
Cm_Lk_CAd_res_MHCII <- prop.test(x = 37155,
                           n = 40527,
                           p = mean(Cc_Lk_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Cm_Lk_CAd_res_MHCII # significantly lower than CA

Cm_Lk_CAd_pvalue_MHCII <- Cm_Lk_CAd_res_MHCII$p.value


# 1-sample proportions test without continuity correction

# data:  37155 out of 40527, null probability mean(Cc_Lk_dist_MHCII$CAd_prop)
# X-squared = 987.11, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.9506122
# 95 percent confidence interval:
#   0.9140676 0.9194458
# sample estimates:
#   p 
# 0.9167962

## 6. Dc and Lk
Dc_Lk_CAd_res_MHCII <- prop.test(x = 3081,
                           n = 3469,
                           p = mean(Dc_Lk_dist_MHCII$CAd_prop), 
                           correct = FALSE)
Dc_Lk_CAd_res_MHCII # significantly lower than CA

Dc_Lk_CAd_pvalue_MHCII <- Dc_Lk_CAd_res_MHCII$p.value


# 1-sample proportions test without continuity correction
# 
# data:  3081 out of 3469, null probability mean(Dc_Lk_dist_MHCII$CAd_prop)
# X-squared = 442.7, df = 1, p-value < 2.2e-16
# alternative hypothesis: true p is not equal to 0.9589935
# 95 percent confidence interval:
#   0.8772316 0.8982141
# sample estimates:
#   p 
# 0.8881522 

# FDR correct the above pvalues

MHCII_pvalues <- c(Cc_Cm_CEd_pvalue_MHCII,
                   Cc_Dc_CEd_pvalue_MHCII,
                   Cc_Lk_CEd_pvalue_MHCII,
                   Cm_Dc_CEd_pvalue_MHCII,
                   Cm_Lk_CEd_pvalue_MHCII,
                   Dc_Lk_CEd_pvalue_MHCII,
                   Cc_Cm_CAd_pvalue_MHCII,
                   Cc_Dc_CAd_pvalue_MHCII,
                   Cc_Lk_CAd_pvalue_MHCII,
                   Cm_Dc_CAd_pvalue_MHCII,
                   Cm_Lk_CAd_pvalue_MHCII,
                   Dc_Lk_CAd_pvalue_MHCII
)
MHCII_fdrs <- p.adjust(MHCII_pvalues, method="fdr")

MHCII_fdrs

# Make a dataframe

MHCII_FDR_corrected_pvalue <- data.frame(species_pairs_distribution, MHCII_pvalues, MHCII_fdrs)

MHCII_FDR_corrected_pvalue
write.csv(MHCII_FDR_corrected_pvalue, "supptableX_MHCII_CondonUsageAnalysis_pvalues.csv")


# MHCII
# Calculate medians of each distribution for co-ancestry score analysis
## co-ancestry:
Cc_Cm_medianCA_MHCII <- median(Cc_Cm_dist_MHCII$CAd_prop)
Cc_Cm_medianCA_MHCII # 0.940622

Cc_Dc_medianCA_MHCII <- median(Cc_Dc_dist_MHCII$CAd_prop)
Cc_Dc_medianCA_MHCII # 0.9654917

Cc_Lk_medianCA_MHCII <- median(Cc_Lk_dist_MHCII$CAd_prop)
Cc_Lk_medianCA_MHCII # 0.9507298

Cm_Dc_medianCA_MHCII <- median(Cm_Dc_dist_MHCII$CAd_prop)
Cm_Dc_medianCA_MHCII # 0.9346478

Cm_Lk_medianCA_MHCII <- median(Cm_Lk_dist_MHCII$CAd_prop)
Cm_Lk_medianCA_MHCII # 0.932341

Dc_Lk_medianCA_MHCII <- median(Dc_Lk_dist_MHCII$CAd_prop)
Dc_Lk_medianCA_MHCII # 0.9590821

## convergent evolution
Cc_Cm_medianCE_MHCII <- median(Cc_Cm_dist_MHCII$CEd_prop)
Cc_Cm_medianCE_MHCII # 0.7534996

Cc_Dc_medianCE_MHCII <- median(Cc_Dc_dist_MHCII$CEd_prop)
Cc_Dc_medianCE_MHCII # 0.7550093

Cc_Lk_medianCE_MHCII <- median(Cc_Lk_dist_MHCII$CEd_prop)
Cc_Lk_medianCE_MHCII # 0.7490491

Cm_Dc_medianCE_MHCII <- median(Cm_Dc_dist_MHCII$CEd_prop)
Cm_Dc_medianCE_MHCII # 0.7720082

Cm_Lk_medianCE_MHCII <- median(Cm_Lk_dist_MHCII$CEd_prop)
Cm_Lk_medianCE_MHCII # 0.7431041

Dc_Lk_medianCE_MHCII <- median(Dc_Lk_dist_MHCII$CEd_prop)
Dc_Lk_medianCE_MHCII # 0.7497294

# CO ANCESTRY SCORE AND REGRESSION FOR EACH LOCUS
# Co-ancestry scores for each species pair
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

CcDc <- phylo_dist[1062, 1391]
CcDc

CcLk <- phylo_dist[1062, 1225]
CcLk

CmDc <- phylo_dist[714, 1391]
CmDc

CmLk <- phylo_dist[714, 1225]
CmLk

DcLk <- phylo_dist[1391, 1225]
DcLk


# Create dataframe of species pairs, proportion of identical codons, median CE and CA, locus, and phylogenetic distance based on beast:

species_pair <- c("Caretta caretta-Chelonia mydas",
                  "Caretta caretta-Dermochelys coriacea",
                  "Caretta caretta-Lepidochelys kempii",
                  "Chelonia mydas-Dermochelys coriacea",
                  "Chelonia mydas-Lepidochelys kempii",
                  "Dermochelys coriacea-Lepidochelys kempii",
                  "Caretta caretta-Chelonia mydas",
                  "Caretta caretta-Dermochelys coriacea",
                  "Caretta caretta-Lepidochelys kempii",
                  "Chelonia mydas-Dermochelys coriacea",
                  "Chelonia mydas-Lepidochelys kempii",
                  "Dermochelys coriacea-Lepidochelys kempii")

# from output of perl scripts (row 20 in the codon_usage_summary xlsx files); MHCI first then MHCII, and in the above order
proportion_identical_codons <- c(0.857167598,
                                 0.922374429,
                                 0.898763362,
                                 0.978428351,
                                 0.886831439,
                                 0.966942149,
                                 0.958820005,
                                 0.903525046,
                                 0.935515259,
                                 0.918235696,
                                 0.925955726,
                                 0.911669192
) 

median_CE <- c(Cc_Cm_medianCE_MHCI,
               Cc_Dc_medianCE_MHCI,
               Cc_Lk_medianCE_MHCI,
               Cm_Dc_medianCE_MHCI,
               Cm_Lk_medianCE_MHCI,
               Dc_Lk_medianCE_MHCI,
               Cc_Cm_medianCE_MHCII,
               Cc_Dc_medianCE_MHCII,
               Cc_Lk_medianCE_MHCII,
               Cm_Dc_medianCE_MHCII,
               Cm_Lk_medianCE_MHCII,
               Dc_Lk_medianCE_MHCII)

median_CA <- c(Cc_Cm_medianCA_MHCI,
               Cc_Dc_medianCA_MHCI,
               Cc_Lk_medianCA_MHCI,
               Cm_Dc_medianCA_MHCI,
               Cm_Lk_medianCA_MHCI,
               Dc_Lk_medianCA_MHCI,
               Cc_Cm_medianCA_MHCII,
               Cc_Dc_medianCA_MHCII,
               Cc_Lk_medianCA_MHCII,
               Cm_Dc_medianCA_MHCII,
               Cm_Lk_medianCA_MHCII,
               Dc_Lk_medianCA_MHCII)

locus <- c("MHCI",
           "MHCI",
           "MHCI",
           "MHCI",
           "MHCI",
           "MHCI",
           "MHCII_chr14",
           "MHCII_chr14",
           "MHCII_chr14",
           "MHCII_chr14",
           "MHCII_chr14",
           "MHCII_chr14")

phylo_dist_beast_distnode <- c(CcCm,
                               CcDc,
                               CcLk,
                               CmDc,
                               CmLk,
                               DcLk)

df <- data.frame(species_pair,
                 proportion_identical_codons,
                 median_CE,
                 median_CA,
                 locus,
                 phylo_dist_beast_distnode
                 )

# compute co-ancestry score
df$coancestry_score <- ((df$proportion_identical_codons) - (df$median_CE)) / ((df$median_CA) - (df$median_CE))


write.csv(df, "coancestry_score_codon_usage_table_2025.csv")

df$locus <- as.factor(df$locus)

# MHCI
MHCI_df <- df %>% filter(locus == "MHCI")
MHCI_correlation <- cor.test(MHCI_df$coancestry_score, MHCI_df$phylo_dist_beast_distnode, method = 'pearson')
MHCI_correlation

# Pearson's product-moment correlation
# 
# data:  MHCI_df$coancestry_score and MHCI_df$phylo_dist_beast_distnode
# t = 0.017826, df = 4, p-value = 0.9866
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.8084966  0.8145819
# sample estimates:
#         cor 
# 0.008912643 


df_MHCI <- df %>% filter(locus == "MHCI")

plot_coancestry_MHCI <- ggplot(df_MHCI,
                 aes(x = coancestry_score,
                     y = phylo_dist_beast_distnode)) +
  geom_point(size = 10) +
  geom_smooth(method = "lm") + #lm method, default 95% confidence interval
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.y = element_text(size = 50),
    axis.text.x = element_text(size = 50)) +
  labs(title = "A") +
  theme(text = element_text(size = 75))
plot_coancestry_MHCI


# MHCII
MHCII_df <- df %>% filter(locus == "MHCII_chr14")
MHCII_correlation <- cor.test(MHCII_df$coancestry_score, MHCII_df$phylo_dist_beast_distnode, method = 'pearson')
MHCII_correlation

# Pearson's product-moment correlation
# 
# data:  MHCII_df$coancestry_score and MHCII_df$phylo_dist_beast_distnode
# t = -1.0274, df = 4, p-value = 0.3623
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.9253470  0.5636498
# sample estimates:
#        cor 
# -0.4569261 


df_MHCII <- df %>% filter(locus == "MHCII_chr14")

plot_coancestry_MHCII <- ggplot(df_MHCII,
                  aes(x = coancestry_score,
                      y = phylo_dist_beast_distnode)) +
  geom_point(size = 10) +
  geom_smooth(method = "lm") + #lm method, default 95% confidence interval
  theme_classic() +
  xlab("") +
  ylab("") +
  theme(
    axis.text.y = element_text(size = 50),
    axis.text.x = element_text(size = 50)) +
  labs(title = "C") +
  theme(text = element_text(size = 75))
plot_coancestry_MHCII

# facet coancestry score plots into one; I'm annotating first and then faceting so that the labels A and C are equivlanet in size and arrangement to the labels B and D for the coodn usage score

coancestry_scores_figure <- ggpubr::ggarrange(plot_coancestry_MHCI,
                                              plot_coancestry_MHCII,
                                              ncol = 1,
                                              nrow = 2)
coancestry_scores_figure
# annotate coancestry score plots
coancestry_scores_figure <- annotate_figure(coancestry_scores_figure,
                bottom = text_grob("co-ancestry score",
                                   size = 75),
                left = text_grob("phylogenetic distance",
                                 rot = 90,
                                 size = 75))
coancestry_scores_figure


# facet codon usage plots into one for each locus
codon_usage_figure_MHCI <- ggpubr::ggarrange(Cc_Cm_histo_MHCI,
                                             Cc_Dc_histo_MHCI,
                                             Cc_Lk_histo_MHCI,
                                             Cm_Dc_histo_MHCI,
                                             Cm_Lk_histo_MHCI,
                                             Dc_Lk_histo_MHCI,
                                             ncol = 2,
                                             nrow = 3,
                                             labels = "B",
                                             font.label = list(size = 75))
codon_usage_figure_MHCI

codon_usage_figure_MHCII <- ggpubr::ggarrange(Cc_Cm_histo_MHCII,
                                              Cc_Dc_histo_MHCII,
                                              Cc_Lk_histo_MHCII,
                                              Cm_Dc_histo_MHCII,
                                              Cm_Lk_histo_MHCII,
                                              Dc_Lk_histo_MHCII,
                                              ncol = 2,
                                              nrow = 3,
                                              labels = "D",
                                              font.label = list(size = 75))
codon_usage_figure_MHCII

# facet codon usage figure together
codon_usage_figure <- ggpubr::ggarrange(codon_usage_figure_MHCI,
                                        codon_usage_figure_MHCII,
                                        ncol = 1,
                                        nrow = 2)
# annotate the faceted codon usage figure
codon_usage_figure <- annotate_figure(codon_usage_figure,
                   bottom = text_grob("proportion identical codons",
                                      size = 75),
                   left = text_grob("frequency",
                                    rot = 90,
                                    size = 75))
codon_usage_figure

# finally, facet together co-ancestry plots and codon usage plots; call the figures so that the elapsed time limit isn't reached

coancestry_scores_figure
codon_usage_figure
coancestry_codon_usage_facet <- ggpubr::ggarrange(coancestry_scores_figure,
                                                  codon_usage_figure,
                                                  ncol = 2,
                                                  nrow = 1)
coancestry_codon_usage_facet
ggsave("figure3_coancestry_codon_usage_facet.svg",
       width = 100,
       height = 75,
       units = "cm",
       limitsize = FALSE,
       dpi = 300)
