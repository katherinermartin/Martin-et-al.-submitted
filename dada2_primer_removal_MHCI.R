## Martin et al: dada2 of MHCI reads to remove forward and reverse primers

library("dada2")

### UCF samples:

### PRE ASSEMBLY
path <- "/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCImiseq_data/11May2023/May23_MHC_run1-388182164/FASTQ_Generation_2023-05-12_17_22_36Z-670700131/fastqc_files/"

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

# forward primer: GATGTATGGGTGTGATCTCCGGG: 23 characters
# reverse primer: TTCACTCGATGCAGGTCDNCTCCAGGT: 27 characters

out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs,
                     trimLeft = c(23,27), compress = FALSE) # very important to have compress there as false, otherwise it'll compress them as fastq.gz files but with a ".fastq" appendage and fastqc won't be happy


### Massachusetts Kemp's ridley reads:
# Not every read in the fastq file has a primer sequence, so I cannot use the above approach, which removes the first 23 and last 27 characters from each read. This would result in loss of data.
# Instead, I'll use the ITS approach, detailed here: https://benjjneb.github.io/dada2/ITS_workflow.html#identify-primers
# This allows me to remove primers based on their nucleotide content, rather than relying on them being in a fixed place.

library(dada2)
library(ShortRead)
library(Biostrings)
path <- "/Users/KatieMartin/Documents/UCF/Research/MHC_species_evo/analysis/MHCI/Massachusetts_Kemps/Kemps_MHC_fastqfiles"  ## CHANGE ME to the directory containing the fastq files.
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_SS.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_SS.fastq.gz", full.names = TRUE))

# Identify primers
FWD <- "GATGTATGGGTGTGATCTCCGGG"  ## CHANGE ME to your forward primer sequence
REV <- "TTCACTCGATGCAGGTCDNCTCCAGGT"  ## CHANGE ME...

# verify the presence and orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# The presence of ambiguous bases (Ns) in the sequencing reads makes accurate mapping of short primer sequences difficult. Next we are going to “pre-filter” the sequences just to remove those with Ns, but perform no other filtering.

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# count the number of times primers appear in the forward and reverse reads
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                          primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                       fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

# Remove primers
cutadapt <- "/Users/KatieMartin/anaconda2/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R


# We now create output filenames for the cutadapt-ed files, and define the parameters we are going to give the cutadapt command. The critical parameters are the primers, and they need to be in the right orientation, i.e. the FWD primer should have been matching the forward-reads in its forward orientation, and the REV primer should have been matching the reverse-reads in its forward orientation. Warning: A lot of output will be written to the screen by cutadapt!


path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# Sanity check that primers were removed:
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), FWD.ReverseReads = sapply(FWD.orients,
                                                                                                        primerHits, fn = fnRs.cut[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits,
                                                                                                                                                                   fn = fnFs.cut[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Success!

# The primer-free sequence files are now ready to be analyzed through the DADA2 pipeline. Similar to the earlier steps of reading in FASTQ files, we read in the names of the cutadapt-ed FASTQ files and applying some string manipulation to get the matched lists of forward and reverse fastq files.

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_R1_SS.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_RR_SS.fastq.gz", full.names = TRUE))

# samples in directory "~/cutadapt" are ready for PEAR

