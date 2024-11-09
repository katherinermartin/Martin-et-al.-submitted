# Martin et al: population genetic diversity metrics on MHCII alleles in C. caretta, C. mydas, D. coriacea, and L. kempii, use as input to scripts to create Figure 1

library(tidyverse)
library(ape) # for manipulating sequence lists
library(pegas) # for nucleotide diversity calculation
library(ggplot2)
library(ggpubr)

setwd("/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_monomorphic_polymorphic_combined/combined_9Sept24/allelic_diversity/")

# read in sequence list
all_seqs <- read.dna("MHCII_all_alleles.fasta", format = "fasta") # must be an alignment

dim(all_seqs)
# convert sequence list from DNAbin object to matrix
seqMat <- as.character(all_seqs)
dim(seqMat) # 308 sequences, 272 bp

# just those sequences that are chr 1
mono_seqs <- read.dna("chr1_alm_for_nucdiv.fasta", format = "fasta") # must be an alignment

dim(mono_seqs)
# convert sequence list from DNAbin object to matrix
mono_seqMat <- as.character(mono_seqs)
dim(mono_seqMat) # 14 sequences, 269 bp

# just those sequences that are chr 14
poly_seqs <- read.dna("../../../MHCII_chr 1_chr 14_combined/allelic_diversity/chr14_alm_for_nucdiv.fasta", format = "fasta") # must be an alignment

dim(poly_seqs)
# convert sequence list from DNAbin object to matrix
poly_seqMat <- as.character(poly_seqs)
dim(poly_seqMat) # 294 sequences, 167 bp

# subset to alleles found in C. caretta
Cc <- seqMat[c("Caca001_MHCIIB", "Caca002_MHCIIB", "Caca003_MHCIIB", "Caca004_MHCIIB", "Caca005_MHCIIB", "Caca006_MHCIIB", "Caca007_MHCIIB", "Caca008_MHCIIB", "Caca009_MHCIIB", "Caca010_MHCIIB", "Caca011_MHCIIB", "Caca012_MHCIIB", "Caca013_MHCIIB", "Caca014_MHCIIB", "Caca015_MHCIIB", "Caca016_MHCIIB", "Caca017_MHCIIB", "Caca018_MHCIIB", "Caca019_MHCIIB", "Caca020_MHCIIB", "Caca021_MHCIIB", "Caca022_MHCIIB", "Caca023_MHCIIB", "Caca024_MHCIIB", "Caca025_MHCIIB", "Caca026_MHCIIB", "Caca027_MHCIIB", "Caca028_MHCIIB", "Caca029_MHCIIB", "Caca030_MHCIIB", "Caca031_MHCIIB", "Caca032_MHCIIB", "Caca033_MHCIIB", "Caca034_MHCIIB", "Caca035_MHCIIB", "Caca036_MHCIIB", "Caca037_MHCIIB", "Caca038_MHCIIB", "Caca039_MHCIIB", "Caca040_MHCIIB", "Caca041_MHCIIB", "Caca042_MHCIIB", "Caca043_MHCIIB", "Caca044_MHCIIB", "Caca045_MHCIIB", "Caca046_MHCIIB"),]

dim(Cc) #46 sequences, 272 bp

# Convert back to DNAbin
Cc <- as.DNAbin(Cc)

# Number of segregating sites:
Cc_seg <- length(seg.sites(Cc))
Cc_seg # 103 segregating sites

# set length of nucleotide sequence
Cc_length <- 272

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cc_percent_segregating <- Cc_seg/Cc_length
Cc_percent_segregating # 0.3786765 or 37.86% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Cc, TRUE) # pi = 0.143979027; variance = 0.005002361

# subset to chr 1 alleles found in C. caretta
Cc_mono <- mono_seqMat[c("Caca001_MHCIIB", "Caca002_MHCIIB"),]

dim(Cc_mono) #2 sequences

# Convert back to DNAbin
Cc_mono <- as.DNAbin(Cc_mono)

# Number of segregating sites:
Cc_mono_seg <- length(seg.sites(Cc_mono))
Cc_mono_seg # 1 segregating sites

# set length of nucleotide sequence
Cc_mono_length <- 269

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cc_mono_percent_segregating <- Cc_mono_seg/Cc_mono_length
Cc_mono_percent_segregating # 0.003717472 or 0.371% segregating sites

# Calculate nucleotide diversity or pi
pegas::nuc.div(Cc_mono, TRUE) # pi = 0.0037174721; variance = 0.0000276392

# subset to chr 14 alleles found in C. caretta
Cc_poly <- poly_seqMat[c("Caca003_MHCIIB", "Caca004_MHCIIB", "Caca005_MHCIIB", "Caca006_MHCIIB", "Caca007_MHCIIB", "Caca008_MHCIIB", "Caca009_MHCIIB", "Caca010_MHCIIB", "Caca011_MHCIIB", "Caca012_MHCIIB", "Caca013_MHCIIB", "Caca014_MHCIIB", "Caca015_MHCIIB", "Caca016_MHCIIB", "Caca017_MHCIIB", "Caca018_MHCIIB", "Caca019_MHCIIB", "Caca020_MHCIIB", "Caca021_MHCIIB", "Caca022_MHCIIB", "Caca023_MHCIIB", "Caca024_MHCIIB", "Caca025_MHCIIB", "Caca026_MHCIIB", "Caca027_MHCIIB", "Caca028_MHCIIB", "Caca029_MHCIIB", "Caca030_MHCIIB", "Caca031_MHCIIB", "Caca032_MHCIIB", "Caca033_MHCIIB", "Caca034_MHCIIB", "Caca035_MHCIIB", "Caca036_MHCIIB", "Caca037_MHCIIB", "Caca038_MHCIIB", "Caca039_MHCIIB", "Caca040_MHCIIB", "Caca041_MHCIIB", "Caca042_MHCIIB", "Caca043_MHCIIB", "Caca044_MHCIIB", "Caca045_MHCIIB", "Caca046_MHCIIB"),]

dim(Cc_poly) #44 sequences, 167 bp

# Convert back to DNAbin
Cc_poly <- as.DNAbin(Cc_poly)

# Number of segregating sites:
Cc_poly_seg <- length(seg.sites(Cc_poly))
Cc_poly_seg # 71 segregating sites

# set length of nucleotide sequence
Cc_poly_length <- 167

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cc_poly_percent_segregating <- Cc_poly_seg/Cc_poly_length
Cc_poly_percent_segregating # 0.4251497 or 42.5% segregating sites

# Calculate nucleotide diversity or pi
pegas::nuc.div(Cc_poly, TRUE) # pi = 0.125533289; variance = 0.003932565

# subset to alleles found in C. mydas
Cm <- seqMat[c("Chmy001_MHCIIB", "Chmy002_MHCIIB", "Chmy003_MHCIIB", "Chmy004_MHCIIB", "Chmy005_MHCIIB", "Chmy006_MHCIIB", "Chmy007_MHCIIB", "Chmy008_MHCIIB", "Chmy009_MHCIIB", "Chmy010_MHCIIB", "Chmy011_MHCIIB", "Chmy012_MHCIIB", "Chmy013_MHCIIB", "Chmy014_MHCIIB", "Chmy015_MHCIIB", "Chmy016_MHCIIB", "Chmy017_MHCIIB", "Chmy018_MHCIIB", "Chmy019_MHCIIB", "Chmy020_MHCIIB", "Chmy021_MHCIIB", "Chmy022_MHCIIB", "Chmy023_MHCIIB", "Chmy024_MHCIIB", "Chmy025_MHCIIB", "Chmy026_MHCIIB", "Chmy027_MHCIIB", "Chmy028_MHCIIB", "Chmy029_MHCIIB", "Chmy030_MHCIIB", "Chmy031_MHCIIB", "Chmy032_MHCIIB", "Chmy033_MHCIIB", "Chmy034_MHCIIB", "Chmy035_MHCIIB", "Chmy036_MHCIIB", "Chmy037_MHCIIB", "Chmy038_MHCIIB", "Chmy039_MHCIIB", "Chmy040_MHCIIB", "Chmy041_MHCIIB", "Chmy042_MHCIIB", "Chmy043_MHCIIB", "Chmy044_MHCIIB", "Chmy045_MHCIIB", "Chmy046_MHCIIB", "Chmy047_MHCIIB", "Chmy048_MHCIIB", "Chmy049_MHCIIB", "Chmy050_MHCIIB", "Chmy051_MHCIIB", "Chmy052_MHCIIB", "Chmy053_MHCIIB", "Chmy054_MHCIIB", "Chmy055_MHCIIB", "Chmy056_MHCIIB", "Chmy057_MHCIIB", "Chmy058_MHCIIB", "Chmy059_MHCIIB", "Chmy060_MHCIIB", "Chmy061_MHCIIB", "Chmy062_MHCIIB", "Chmy063_MHCIIB", "Chmy064_MHCIIB", "Chmy065_MHCIIB", "Chmy066_MHCIIB", "Chmy067_MHCIIB", "Chmy068_MHCIIB", "Chmy069_MHCIIB", "Chmy070_MHCIIB", "Chmy071_MHCIIB", "Chmy072_MHCIIB", "Chmy073_MHCIIB", "Chmy074_MHCIIB", "Chmy075_MHCIIB", "Chmy076_MHCIIB", "Chmy077_MHCIIB", "Chmy078_MHCIIB", "Chmy079_MHCIIB", "Chmy080_MHCIIB", "Chmy081_MHCIIB", "Chmy082_MHCIIB", "Chmy083_MHCIIB", "Chmy084_MHCIIB", "Chmy085_MHCIIB", "Chmy086_MHCIIB", "Chmy087_MHCIIB", "Chmy088_MHCIIB", "Chmy089_MHCIIB", "Chmy090_MHCIIB", "Chmy091_MHCIIB", "Chmy092_MHCIIB", "Chmy093_MHCIIB", "Chmy094_MHCIIB", "Chmy095_MHCIIB", "Chmy096_MHCIIB", "Chmy097_MHCIIB", "Chmy098_MHCIIB", "Chmy099_MHCIIB", "Chmy100_MHCIIB", "Chmy101_MHCIIB", "Chmy102_MHCIIB", "Chmy103_MHCIIB", "Chmy104_MHCIIB", "Chmy105_MHCIIB", "Chmy106_MHCIIB", "Chmy107_MHCIIB", "Chmy108_MHCIIB", "Chmy109_MHCIIB", "Chmy110_MHCIIB", "Chmy111_MHCIIB", "Chmy112_MHCIIB", "Chmy113_MHCIIB", "Chmy114_MHCIIB", "Chmy115_MHCIIB", "Chmy116_MHCIIB", "Chmy117_MHCIIB", "Chmy118_MHCIIB", "Chmy119_MHCIIB", "Chmy120_MHCIIB", "Chmy121_MHCIIB", "Chmy122_MHCIIB", "Chmy123_MHCIIB", "Chmy124_MHCIIB", "Chmy125_MHCIIB", "Chmy126_MHCIIB", "Chmy127_MHCIIB", "Chmy128_MHCIIB", "Chmy129_MHCIIB", "Chmy130_MHCIIB", "Chmy131_MHCIIB", "Chmy132_MHCIIB", "Chmy133_MHCIIB", "Chmy134_MHCIIB", "Chmy135_MHCIIB", "Chmy136_MHCIIB", "Chmy137_MHCIIB", "Chmy138_MHCIIB", "Chmy139_MHCIIB", "Chmy140_MHCIIB", "Chmy141_MHCIIB", "Chmy142_MHCIIB", "Chmy143_MHCIIB", "Chmy144_MHCIIB", "Chmy145_MHCIIB", "Chmy146_MHCIIB", "Chmy147_MHCIIB", "Chmy148_MHCIIB", "Chmy149_MHCIIB", "Chmy150_MHCIIB", "Chmy151_MHCIIB", "Chmy152_MHCIIB", "Chmy153_MHCIIB", "Chmy154_MHCIIB", "Chmy155_MHCIIB", "Chmy156_MHCIIB", "Chmy157_MHCIIB", "Chmy158_MHCIIB", "Chmy159_MHCIIB", "Chmy160_MHCIIB", "Chmy161_MHCIIB", "Chmy162_MHCIIB", "Chmy163_MHCIIB", "Chmy164_MHCIIB", "Chmy165_MHCIIB", "Chmy166_MHCIIB", "Chmy167_MHCIIB", "Chmy168_MHCIIB", "Chmy169_MHCIIB", "Chmy170_MHCIIB", "Chmy171_MHCIIB", "Chmy172_MHCIIB", "Chmy173_MHCIIB", "Chmy174_MHCIIB", "Chmy175_MHCIIB", "Chmy176_MHCIIB", "Chmy177_MHCIIB", "Chmy178_MHCIIB", "Chmy179_MHCIIB", "Chmy180_MHCIIB", "Chmy181_MHCIIB", "Chmy182_MHCIIB", "Chmy183_MHCIIB", "Chmy184_MHCIIB", "Chmy185_MHCIIB", "Chmy186_MHCIIB", "Chmy187_MHCIIB", "Chmy188_MHCIIB", "Chmy189_MHCIIB", "Chmy190_MHCIIB", "Chmy191_MHCIIB", "Chmy192_MHCIIB", "Chmy193_MHCIIB", "Chmy194_MHCIIB", "Chmy195_MHCIIB", "Chmy196_MHCIIB", "Chmy197_MHCIIB", "Chmy198_MHCIIB", "Chmy199_MHCIIB", "Chmy200_MHCIIB", "Chmy201_MHCIIB", "Chmy202_MHCIIB", "Chmy203_MHCIIB", "Chmy204_MHCIIB", "Chmy205_MHCIIB", "Chmy206_MHCIIB", "Chmy207_MHCIIB", "Chmy208_MHCIIB", "Chmy209_MHCIIB", "Chmy210_MHCIIB", "Chmy211_MHCIIB", "Chmy212_MHCIIB", "Chmy213_MHCIIB", "Chmy214_MHCIIB", "Chmy215_MHCIIB", "Chmy216_MHCIIB", "Chmy217_MHCIIB", "Chmy218_MHCIIB", "Chmy219_MHCIIB"),]

dim(Cm) # 219 sequences, 292 bp

# Convert back to DNAbin
Cm <- as.DNAbin(Cm)

# Number of segregating sites:
Cm_seg <- length(seg.sites(Cm))
Cm_seg # 110 segregating sites

# set length of nucleotide sequence
Cm_length <- 292

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cm_percent_segregating <- Cm_seg/Cm_length
Cm_percent_segregating # 0.3767123 or 37.7% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Cm, TRUE) # pi = 0.174329807; variance = 0.007016745

# subset to chr 1 alleles found in C. mydas
Cm_mono <- mono_seqMat[c("Chmy001_MHCIIB", "Chmy002_MHCIIB", "Chmy003_MHCIIB", "Chmy004_MHCIIB", "Chmy005_MHCIIB", "Chmy006_MHCIIB", "Chmy007_MHCIIB", "Chmy008_MHCIIB", "Chmy009_MHCIIB", "Chmy010_MHCIIB"),]

dim(Cm_mono) # 10 sequences, 269 bp

# Convert back to DNAbin
Cm_mono <- as.DNAbin(Cm_mono)

# Number of segregating sites:
Cm_mono_seg <- length(seg.sites(Cm_mono))
Cm_mono_seg # 9 segregating sites

# set length of nucleotide sequence
Cm_mono_length <- 289

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cm_mono_percent_segregating <- Cm_mono_seg/Cm_mono_length
Cm_mono_percent_segregating # 0.03114187 or 3.11% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Cm_mono, TRUE) # pi = 7.682776e-03; variance = 2.810448e-05

# subset to chr 14 alleles found in C. mydas
Cm_poly <- poly_seqMat[c("Chmy011_MHCIIB", "Chmy012_MHCIIB", "Chmy013_MHCIIB", "Chmy014_MHCIIB", "Chmy015_MHCIIB", "Chmy016_MHCIIB", "Chmy017_MHCIIB", "Chmy018_MHCIIB", "Chmy019_MHCIIB", "Chmy020_MHCIIB", "Chmy021_MHCIIB", "Chmy022_MHCIIB", "Chmy023_MHCIIB", "Chmy024_MHCIIB", "Chmy025_MHCIIB", "Chmy026_MHCIIB", "Chmy027_MHCIIB", "Chmy028_MHCIIB", "Chmy029_MHCIIB", "Chmy030_MHCIIB", "Chmy031_MHCIIB", "Chmy032_MHCIIB", "Chmy033_MHCIIB", "Chmy034_MHCIIB", "Chmy035_MHCIIB", "Chmy036_MHCIIB", "Chmy037_MHCIIB", "Chmy038_MHCIIB", "Chmy039_MHCIIB", "Chmy040_MHCIIB", "Chmy041_MHCIIB", "Chmy042_MHCIIB", "Chmy043_MHCIIB", "Chmy044_MHCIIB", "Chmy045_MHCIIB", "Chmy046_MHCIIB", "Chmy047_MHCIIB", "Chmy048_MHCIIB", "Chmy049_MHCIIB", "Chmy050_MHCIIB", "Chmy051_MHCIIB", "Chmy052_MHCIIB", "Chmy053_MHCIIB", "Chmy054_MHCIIB", "Chmy055_MHCIIB", "Chmy056_MHCIIB", "Chmy057_MHCIIB", "Chmy058_MHCIIB", "Chmy059_MHCIIB", "Chmy060_MHCIIB", "Chmy061_MHCIIB", "Chmy062_MHCIIB", "Chmy063_MHCIIB", "Chmy064_MHCIIB", "Chmy065_MHCIIB", "Chmy066_MHCIIB", "Chmy067_MHCIIB", "Chmy068_MHCIIB", "Chmy069_MHCIIB", "Chmy070_MHCIIB", "Chmy071_MHCIIB", "Chmy072_MHCIIB", "Chmy073_MHCIIB", "Chmy074_MHCIIB", "Chmy075_MHCIIB", "Chmy076_MHCIIB", "Chmy077_MHCIIB", "Chmy078_MHCIIB", "Chmy079_MHCIIB", "Chmy080_MHCIIB", "Chmy081_MHCIIB", "Chmy082_MHCIIB", "Chmy083_MHCIIB", "Chmy084_MHCIIB", "Chmy085_MHCIIB", "Chmy086_MHCIIB", "Chmy087_MHCIIB", "Chmy088_MHCIIB", "Chmy089_MHCIIB", "Chmy090_MHCIIB", "Chmy091_MHCIIB", "Chmy092_MHCIIB", "Chmy093_MHCIIB", "Chmy094_MHCIIB", "Chmy095_MHCIIB", "Chmy096_MHCIIB", "Chmy097_MHCIIB", "Chmy098_MHCIIB", "Chmy099_MHCIIB", "Chmy100_MHCIIB", "Chmy101_MHCIIB", "Chmy102_MHCIIB", "Chmy103_MHCIIB", "Chmy104_MHCIIB", "Chmy105_MHCIIB", "Chmy106_MHCIIB", "Chmy107_MHCIIB", "Chmy108_MHCIIB", "Chmy109_MHCIIB", "Chmy110_MHCIIB", "Chmy111_MHCIIB", "Chmy112_MHCIIB", "Chmy113_MHCIIB", "Chmy114_MHCIIB", "Chmy115_MHCIIB", "Chmy116_MHCIIB", "Chmy117_MHCIIB", "Chmy118_MHCIIB", "Chmy119_MHCIIB", "Chmy120_MHCIIB", "Chmy121_MHCIIB", "Chmy122_MHCIIB", "Chmy123_MHCIIB", "Chmy124_MHCIIB", "Chmy125_MHCIIB", "Chmy126_MHCIIB", "Chmy127_MHCIIB", "Chmy128_MHCIIB", "Chmy129_MHCIIB", "Chmy130_MHCIIB", "Chmy131_MHCIIB", "Chmy132_MHCIIB", "Chmy133_MHCIIB", "Chmy134_MHCIIB", "Chmy135_MHCIIB", "Chmy136_MHCIIB", "Chmy137_MHCIIB", "Chmy138_MHCIIB", "Chmy139_MHCIIB", "Chmy140_MHCIIB", "Chmy141_MHCIIB", "Chmy142_MHCIIB", "Chmy143_MHCIIB", "Chmy144_MHCIIB", "Chmy145_MHCIIB", "Chmy146_MHCIIB", "Chmy147_MHCIIB", "Chmy148_MHCIIB", "Chmy149_MHCIIB", "Chmy150_MHCIIB", "Chmy151_MHCIIB", "Chmy152_MHCIIB", "Chmy153_MHCIIB", "Chmy154_MHCIIB", "Chmy155_MHCIIB", "Chmy156_MHCIIB", "Chmy157_MHCIIB", "Chmy158_MHCIIB", "Chmy159_MHCIIB", "Chmy160_MHCIIB", "Chmy161_MHCIIB", "Chmy162_MHCIIB", "Chmy163_MHCIIB", "Chmy164_MHCIIB", "Chmy165_MHCIIB", "Chmy166_MHCIIB", "Chmy167_MHCIIB", "Chmy168_MHCIIB", "Chmy169_MHCIIB", "Chmy170_MHCIIB", "Chmy171_MHCIIB", "Chmy172_MHCIIB", "Chmy173_MHCIIB", "Chmy174_MHCIIB", "Chmy175_MHCIIB", "Chmy176_MHCIIB", "Chmy177_MHCIIB", "Chmy178_MHCIIB", "Chmy179_MHCIIB", "Chmy180_MHCIIB", "Chmy181_MHCIIB", "Chmy182_MHCIIB", "Chmy183_MHCIIB", "Chmy184_MHCIIB", "Chmy185_MHCIIB", "Chmy186_MHCIIB", "Chmy187_MHCIIB", "Chmy188_MHCIIB", "Chmy189_MHCIIB", "Chmy190_MHCIIB", "Chmy191_MHCIIB", "Chmy192_MHCIIB", "Chmy193_MHCIIB", "Chmy194_MHCIIB", "Chmy195_MHCIIB", "Chmy196_MHCIIB", "Chmy197_MHCIIB", "Chmy198_MHCIIB", "Chmy199_MHCIIB", "Chmy200_MHCIIB", "Chmy201_MHCIIB", "Chmy202_MHCIIB", "Chmy203_MHCIIB", "Chmy204_MHCIIB", "Chmy205_MHCIIB", "Chmy206_MHCIIB", "Chmy207_MHCIIB", "Chmy208_MHCIIB", "Chmy209_MHCIIB", "Chmy210_MHCIIB", "Chmy211_MHCIIB", "Chmy212_MHCIIB", "Chmy213_MHCIIB", "Chmy214_MHCIIB", "Chmy215_MHCIIB", "Chmy216_MHCIIB", "Chmy217_MHCIIB", "Chmy218_MHCIIB", "Chmy219_MHCIIB"),]

dim(Cm_poly) # 209 sequences, 167 bp

# Convert back to DNAbin
Cm_poly <- as.DNAbin(Cm_poly)

# Number of segregating sites:
Cm_poly_seg <- length(seg.sites(Cm_poly))
Cm_poly_seg # 92 segregating sites

# set length of nucleotide sequence
Cm_poly_length <- 167

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Cm_poly_percent_segregating <- Cm_poly_seg/Cm_poly_length
Cm_poly_percent_segregating # 0.5508982 or 55.08% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Cm_poly, TRUE) # pi = 0.159703596; variance = 0.006044555


# subset to alleles found in D. coriacea
Dc <- seqMat[c("Deco001_MHCIIB", "Deco002_MHCIIB", "Deco003_MHCIIB", "Deco004_MHCIIB", "Deco005_MHCIIB", "Deco006_MHCIIB", "Deco007_MHCIIB", "Deco008_MHCIIB", "Deco009_MHCIIB", "Deco010_MHCIIB", "Deco011_MHCIIB", "Deco012_MHCIIB", "Deco013_MHCIIB", "Deco014_MHCIIB", "Deco015_MHCIIB", "Deco016_MHCIIB", "Deco017_MHCIIB", "Deco018_MHCIIB", "Deco019_MHCIIB", "Deco020_MHCIIB", "Deco021_MHCIIB"),]

dim(Dc) # 21 sequences, 292 bp

# Convert back to DNAbin
Dc <- as.DNAbin(Dc)

# Number of segregating sites:
Dc_seg <- length(seg.sites(Dc))
Dc_seg # 78 segregating sites

# set length of nucleotide sequence
Dc_length <- 292

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Dc_percent_segregating <- Dc_seg/Dc_length
Dc_percent_segregating # 0.2671233 or 26.7% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Dc, TRUE) # pi = 0.068902439; variance = 0.001254568

# cannot do analyses on chr 1 alleles found in D. coriacea because there is only one

# subset to chr 14 alleles found in D. coriacea
Dc_poly <- poly_seqMat[c("Deco002_MHCIIB", "Deco003_MHCIIB", "Deco004_MHCIIB", "Deco005_MHCIIB", "Deco006_MHCIIB", "Deco007_MHCIIB", "Deco008_MHCIIB", "Deco009_MHCIIB", "Deco010_MHCIIB", "Deco011_MHCIIB", "Deco012_MHCIIB", "Deco013_MHCIIB", "Deco014_MHCIIB", "Deco015_MHCIIB", "Deco016_MHCIIB", "Deco017_MHCIIB", "Deco018_MHCIIB", "Deco019_MHCIIB", "Deco020_MHCIIB", "Deco021_MHCIIB"),]

dim(Dc_poly) # 20 sequences, 167 bp

# Convert back to DNAbin
Dc_poly <- as.DNAbin(Dc_poly)

# Number of segregating sites:
Dc_poly_seg <- length(seg.sites(Dc_poly))
Dc_poly_seg # 14 segregating sites

# set length of nucleotide sequence
Dc_poly_length <- 167

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Dc_poly_percent_segregating <- Dc_poly_seg/Dc_poly_length
Dc_poly_percent_segregating # 0.08383234 or 8.38% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Dc_poly, TRUE) # pi = 0.0338796092; variance = 0.0003586787


# subset to alleles found in L. kempii
Lk <- seqMat[c("Leke001_MHCIIB","Leke002_MHCIIB","Leke003_MHCIIB","Leke004_MHCIIB","Leke005_MHCIIB","Leke006_MHCIIB","Leke007_MHCIIB","Leke008_MHCIIB","Leke009_MHCIIB","Leke010_MHCIIB","Leke011_MHCIIB","Leke012_MHCIIB","Leke013_MHCIIB","Leke014_MHCIIB","Leke015_MHCIIB","Leke016_MHCIIB","Leke017_MHCIIB","Leke018_MHCIIB","Leke019_MHCIIB","Leke020_MHCIIB","Leke021_MHCIIB","Leke022_MHCIIB"),]

dim(Lk) # 22 sequences, 292 bp

# Convert back to DNAbin
Lk <- as.DNAbin(Lk)

# Number of segregating sites:
Lk_seg <- length(seg.sites(Lk))
Lk_seg # 97 segregating sites

# set length of nucleotide sequence
Lk_length <- 292

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Lk_percent_segregating <- Lk_seg/Lk_length
Lk_percent_segregating # 0.3321918 or 33.2% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Lk, TRUE) # pi = 0.157903073; variance = 0.006301842

# cannot do analyses on chr 1 alleles found in L. kempii because there is only one

# subset to chr 14 alleles found in L. kempii
Lk_poly <- poly_seqMat[c("Leke002_MHCIIB","Leke003_MHCIIB","Leke004_MHCIIB","Leke005_MHCIIB","Leke006_MHCIIB","Leke007_MHCIIB","Leke008_MHCIIB","Leke009_MHCIIB","Leke010_MHCIIB","Leke011_MHCIIB","Leke012_MHCIIB","Leke013_MHCIIB","Leke014_MHCIIB","Leke015_MHCIIB","Leke016_MHCIIB","Leke017_MHCIIB","Leke018_MHCIIB","Leke019_MHCIIB","Leke020_MHCIIB","Leke021_MHCIIB","Leke022_MHCIIB"),]

dim(Lk_poly) # 21 sequences, 167 bp

# Convert back to DNAbin
Lk_poly <- as.DNAbin(Lk_poly)

# Number of segregating sites:
Lk_poly_seg <- length(seg.sites(Lk_poly))
Lk_poly_seg # 64 segregating sites

# set length of nucleotide sequence
Lk_poly_length <- 167

# Percentage of segregating sites:
# divide segregating sites by sequence length for percentage segregating sites
Lk_poly_percent_segregating <- Lk_poly_seg/Lk_poly_length
Lk_poly_percent_segregating # 0.3832335 or 38.3% segregating sites

# Calculate nucleotide diversity or pi

pegas::nuc.div(Lk_poly, TRUE) # pi = 0.138323353; variance = 0.005011115
