# Martin et al.: population genetic diversity metrics on MHC class I alleles in C. caretta, C. mydas, D. coriacea, and L. kempii, use as input to scripts to create Figure 1

library(tidyverse)
library(ape) # for manipulating sequence lists
library(pegas) # for nucleotide diversity calculation
library(ggplot2)
library(ggpubr)


# read in sequence list
all_seqs <- read.dna("MHCI_CcCmDcLk_May2023_clustalomega_alm_nooutgroups.fasta", format = "fasta") # must be an alignment
dim(all_seqs) # 162

# convert sequence list from DNAbin object to matrix
seqMat <- as.character(all_seqs)
dim(seqMat)

# subset to alleles found in C. caretta (FL only)
CcFL <- seqMat[c("Caca01",
                 "Carca09",
                 "Chmy07",
                 "Carca32",
                 "Carca15",
                 "Carca19",
                 "Chmy10",
                 "Chmy12",
                 "Caca02",
                 "Chmy17",
                 "Carca14",
                 "Carca11",
                 "Carca16",
                 "Carca17",
                 "Carca27",
                 "Carca13",
                 "Chmy26",
                 "Chmy31",
                 "Carca56",
                 "Caca03",
                 "Carca28",
                 "Chmy41",
                 "Carca109",
                 "Caca04",
                 "Chmy55",
                 "Caca05",
                 "Carca25",
                 "Carca81",
                 "Chmy66",
                 "Caca06",
                 "Caca07",
                 "Caca08",
                 "Caca09"),]

dim(CcFL) # 33 sequences

# Convert back to DNAbin
CcFL<-as.DNAbin(CcFL)

# Number of segregating sites:
CcFL_seg <- length(seg.sites(CcFL))
CcFL_seg # 95 segregating sites

# set length of nucleotide sequence
CcFL_length <- 162

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
CcFL_percent_segregating <- CcFL_seg/CcFL_length
CcFL_percent_segregating # 0.5864198 or 58.6% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(CcFL, TRUE) # pi = 0.23341049; variance = 0.01340812


# subset to alleles found in C. mydas in Florida
Cm_FL <- seqMat[c("Caca01",
                  "Chmy01",
                  "Chmy02",
                  "Carca09",
                  "Chmy03",
                  "Chmy04",
                  "Chmy05",
                  "Chmy06",
                  "Chmy07",
                  "Carca15",
                  "Carca19",
                  "Chmy08",
                  "Chmy09",
                  "Chmy10",
                  "Chmy11",
                  "Chmy12",
                  "Chmy13",
                  "Chmy14",
                  "Chmy15",
                  "Chmy16",
                  "Caca02",
                  "Chmy17",
                  "Chmy18",
                  "Chmy19",
                  "Chmy20",
                  "Chmy21",
                  "Chmy22",
                  "Chmy23",
                  "Chmy24",
                  "Chmy25",
                  "Chmy26",
                  "Chmy27",
                  "Chmy28",
                  "Chmy29",
                  "Chmy30",
                  "Chmy31",
                  "Chmy32",
                  "Chmy33",
                  "Chmy34",
                  "Chmy35",
                  "Chmy36",
                  "Chmy37",
                  "Chmy38",
                  "Chmy39",
                  "Chmy40",
                  "Chmy41",
                  "Chmy42",
                  "Chmy43",
                  "Chmy44",
                  "Chmy45",
                  "Chmy46",
                  "Chmy47",
                  "Chmy48",
                  "Chmy49",
                  "Chmy50",
                  "Chmy51",
                  "Chmy52",
                  "Chmy53",
                  "Chmy54",
                  "Chmy55",
                  "Chmy56",
                  "Chmy57",
                  "Chmy58",
                  "Chmy59",
                  "Chmy60",
                  "Chmy61",
                  "Chmy62",
                  "Chmy63",
                  "Chmy64",
                  "Chmy65",
                  "Chmy66",
                  "Chmy67",
                  "Chmy68",
                  "Chmy69",
                  "Chmy70",
                  "Chmy71",
                  "Chmy72",
                  "Chmy73",
                  "Chmy74",
                  "Chmy75",
                  "Chmy76",
                  "Chmy77",
                  "Chmy78",
                  "Chmy79",
                  "Chmy80",
                  "Chmy81",
                  "Chmy82",
                  "Chmy83",
                  "Chmy84",
                  "Chmy85",
                  "Chmy86",
                  "Chmy87",
                  "Chmy88",
                  "Chmy89",
                  "Chmy90",
                  "Chmy91",
                  "Chmy92",
                  "Chmy93"),]

dim(Cm_FL) # 98 sequences

# Convert back to DNAbin
Cm_FL <-as.DNAbin(Cm_FL)

# Number of segregating sites:
Cm_seg <- length(seg.sites(Cm_FL))
Cm_seg # 106 segregating sites

# set length of nucleotide sequence
Cm_length <- 162

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cm_percent_segregating <- Cm_seg/Cm_length
Cm_percent_segregating # 0.654321 or 65.4% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Cm_FL, TRUE) # pi = 0.201021837; variance = 0.009590092


# subset to alleles found in D. coriacea
Dc <- seqMat[c("Deco01",
               "Deco02",
               "Deco03",
               "Deco04",
               "Deco05",
               "Deco06",
               "Deco07",
               "Deco08"),]

dim(Dc) # 8 sequences

# Convert back to DNAbin
Dc <- as.DNAbin(Dc)

# Number of segregating sites:
Dc_seg <- length(seg.sites(Dc))
Dc_seg # 6 segregating sites

# set length of nucleotide sequence
Dc_length <- 162

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Dc_percent_segregating <- Dc_seg/Dc_length
Dc_percent_segregating # 0.03703704 or 3.7% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Dc, TRUE) # pi = 0.0178571429; variance = 0.000142145


# subset to alleles found in L. kempii
Lk <- seqMat[c("Chmy02",
               "Chmy09",
               "Chmy14",
               "Chmy17",
               "Chmy87",
               "Chmy90",
               "Leke01",
               "Leke02",
               "Leke03",
               "Leke04",
               "Leke05",
               "Leke06",
               "Leke07",
               "Leke08",
               "Leke09",
               "Leke10",
               "Leke11",
               "Leke12",
               "Leke13",
               "Leke14",
               "Leke15",
               "Leke16",
               "Leke17",
               "Leke18"),]

dim(Lk) # 24 sequences

# Convert back to DNAbin
Lk <- as.DNAbin(Lk)

# Number of segregating sites:
Lk_seg <- length(seg.sites(Lk))
Lk_seg # 92 segregating sites

# set length of nucleotide sequence
Lk_length <- 162

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Lk_percent_segregating <- Lk_seg/Lk_length
Lk_percent_segregating # 0.5679012 or 56.8% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Lk, TRUE) # pi = 0.2245706; variance = 0.0127448

## now graph these data
allelic_summary_stats <- read.csv("MHCI_summary_stats_allelic_div.csv")

indiv_per_spp <- ggplot(data = allelic_summary_stats, aes(x = species,
                                         y = turtles_sampled,
                                         fill = species)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("species") +
  ylab("individuals") +
  scale_fill_manual(values = c("#ff981a",
                               "#4ca64c",
                               "#ad95cf",
                               "#95b7cf")) +
  labs(title = "individuals sampled per species") +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  scale_y_continuous(limits = c(0, 300))
indiv_per_spp

alleles_per_spp <- ggplot(data = allelic_summary_stats, aes(x = species,
                                                            y = alleles_sampled,
                                                            fill = species)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("species") +
  ylab("alleles") +
  scale_fill_manual(values = c("#ff981a",
                               "#4ca64c",
                               "#ad95cf",
                               "#95b7cf")) +
  labs(title = "alleles sampled per species") +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  scale_y_continuous(limits = c(0, 125))
alleles_per_spp

nucleotide_div <- ggplot(data = allelic_summary_stats, aes(x = species,
                                                            y = nucleotide_diversity,
                                                            fill = species)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  xlab("species") +
  ylab("nucleotide diversity") +
  scale_fill_manual(values = c("#ff981a",
                               "#4ca64c",
                               "#ad95cf",
                               "#95b7cf")) +
  labs(title = "nucleotide diversity per species") +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  geom_errorbar(aes(ymin=nucleotide_diversity-nuc_div_variance, ymax=nucleotide_diversity+nuc_div_variance), width=.2,
                position=position_dodge(.9))
nucleotide_div


# read in a different dataframe; this is raw data on alleles per individual per spp, geom_box will calculate quartiles, avg, max, min
allele_count <- read.csv("MHCI_alleles_per_indiv_Cc_Cm_Dc_Lk.csv")

alleles <- ggplot(data = allele_count, aes(x = species,
                                                           y = allele_count,
                                                           fill = species)) +
  geom_boxplot() +
  theme_classic() +
  xlab("species") +
  ylab("alleles") +
  scale_fill_manual(values = c("#ff981a",
                               "#4ca64c",
                               "#ad95cf",
                               "#95b7cf")) +
  labs(title = "average alleles per individual species") +
  theme(legend.position = "none", plot.title = element_text(size=11)) +
  coord_cartesian(ylim = c(0, 8))
alleles




indiv_per_spp
ggsave("indiv_per_spp.svg")

alleles_per_spp
ggsave("allelic_summary_stats/alleles_per_spp.svg")

nucleotide_div
ggsave("allelic_summary_stats/nucleotide_div.svg")

alleles
ggsave("allelic_summary_stats/alleles.svg")

# now facet plots
allele_stats_figure <- ggpubr::ggarrange(indiv_per_spp,
                                         nucleotide_div,
                                         alleles_per_spp,
                                         alleles,
                                        ncol = 2, nrow = 2)
allele_stats_figure
