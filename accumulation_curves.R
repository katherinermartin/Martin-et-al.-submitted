## Martin et al: accumulation curve analysis of MHC alleles recovered
### To determine if we've sampled enough to give us an idea of how much more diversity could be out there.

### Running these analyses by species

library(vegan) # also loads vegan, for accumulation curves
library(tidyverse)
library(ggplotify) # for as.ggplot
# to export as a non Cairo SVG surface so that the text is edit-able

# loading data that has counts for alleles by species. Not including Hawaiian Cms just yet because I don't have those data.

## MHCI

df_MHCI <- read.csv("MHCI_allele_counts_for_rarefaction.csv")

df_MHCI <- df_MHCI %>% remove_rownames %>% column_to_rownames(var="allele")

alleleAbund_MHCI <- rowSums(df_MHCI) # gives the number of individuals "found" in each allele...
alleleAbund_MHCI

raremin_MHCI <- min(rowSums(df_MHCI))  #rarefaction uses the smallest number of individuals found per allele to extrapolate the expected number if all other alleles were found in only that number of individuals
raremin_MHCI # view smallest number of individuals

sRare_MHCI <- rarefy(df_MHCI, raremin_MHCI) # now use function rarefy
sRare_MHCI #gives an "expected"rarefied" number of alleles (not obs) if only 2 individuals were present

# use as.ggplot to save and get around it being a non-Cairo
as.ggplot(~rarecurve(df_MHCI,
          col = c("#ff981a",
                  "#4ca64c",
                  "#ad95cf",
                  "#95b7cf"),
          lty = 1, # solid line
          lwd = 3.5, # line thickness
          cex.axis = 1.3, # axis font size
          cex.lab = 1.3, # axis font size
          cex = 1.3, # overall font size
          xlab = "total alleles sampled",
          ylab = "unique alleles"))
ggsave("MHCI_rarefaction_28Dec23.svg")

## MHCII monomorphic
df_MHCII_mono <- read.csv("monomorphic_MHCII_allele_counts_for_rarefaction.csv")

df_MHCII_mono <- df_MHCII_mono %>% remove_rownames %>% column_to_rownames(var="allele")

alleleAbund_MHCII_mono <- rowSums(df_MHCII_mono) # gives the number of individuals "found" in each allele...
alleleAbund_MHCII_mono

raremin_MHCII_mono <- min(rowSums(df_MHCII_mono))  #rarefaction uses the smallest number of individuals found per allele to extrapolate the expected number if all other alleles were found in only that number of individuals
raremin_MHCII_mono # view smallest number of individuals

sRare_MHCII_mono <- rarefy(df_MHCII_mono, raremin_MHCII_mono) # now use function rarefy
sRare_MHCII_mono #gives an "expected"rarefied" number of alleles (not obs) if only 2 individuals were present

as.ggplot(~rarecurve(df_MHCII_mono,
          col = c("#ff981a",
                  "#4ca64c",
                  "#ad95cf",
                  "#95b7cf"),
          lty = 1, # solid line
          lwd = 3.5, # line thickness
          cex.axis = 1.3, # axis font size
          cex.lab = 1.3, # axis font size
          cex = 1.3, # overall font size
          xlab = "total alleles sampled",
          ylab = "unique alleles"))

ggsave("MHCII_mono_rarefaction_28Dec23.svg")

## MHCII polymorphic
df_MHCII_poly <- read.csv("polymorphic_MHCII_allele_counts_for_rarefaction.csv")

df_MHCII_poly <- df_MHCII_poly %>% remove_rownames %>% column_to_rownames(var="allele")

alleleAbund_MHCII_poly <- rowSums(df_MHCII_poly) # gives the number of individuals "found" in each allele...
alleleAbund_MHCII_poly

raremin_MHCII_poly <- min(rowSums(df_MHCII_poly))  #rarefaction uses the smallest number of individuals found per allele to extrapolate the expected number if all other alleles were found in only that number of individuals
raremin_MHCII_poly # view smallest number of individuals

sRare_MHCII_poly <- rarefy(df_MHCII_poly, raremin_MHCII_poly) # now use function rarefy
sRare_MHCII_poly #gives an "expected"rarefied" number of alleles (not obs) if only 2 individuals were present

as.ggplot(~rarecurve(df_MHCII_poly,
          col = c("#ff981a",
                  "#4ca64c",
                  "#ad95cf",
                  "#95b7cf"),
          lty = 1, # solid line,
          lwd = 3.5, # line thickness
          cex.axis = 1.3, # axis font size
          cex.lab = 1.3, # axis font size
          cex = 1.3, # overall font size
          xlab = "total alleles sampled",
          ylab = "unique alleles"))
ggsave("/MHCII_poly_rarefaction_28Dec23.svg")

