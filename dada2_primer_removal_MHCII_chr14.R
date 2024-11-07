## dada2 of MHC reads to remove forward and reverse primers for MHCII chromosome 14

library(dada2)

### PRE ASSEMBLY
path <- "/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCII_chr1/miseq_data/sea_turtle_MHCII_chr1_fastqz/"

# read in fastq files
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Assign the filenames for the trimmed fastq files. This will append the R1 and R2 to the ends of the trimmed files and will create a subdirectory called "primers_trimmed"

trimFs <- file.path(path, "primers_trimmed", paste0(sample.names, "_R1_001_trim.fastq"))
trimRs <- file.path(path, "primers_trimmed", paste0(sample.names, "_R2_001_trim.fastq"))
names(trimFs) <- sample.names
names(trimRs) <- sample.names

# forward primer: TGTCTCTACACCAACGGCACC: 21 characters
# reverse primer: CGTAGTTGTGCCGGCAGAAC: 20 characters

out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs,
                     trimLeft = c(21,20), compress = FALSE) # very important to have compress there as false, otherwise it'll compress them as fastq.gz files but with a ".fastq" appendage and fastqc won't be happy

