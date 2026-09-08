
# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)
library(metagenomeSeq)
library(ADAPT)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ConsensusMetaDA)

############### Figure 2: Differential abundance UpSet and Venn Diagram ##################

## Load data Marine plastics
biome_file <- "./Data/sed_plastic_rosato_ASVs_table.biom"

sample_table_file <-  "./Data/sed_plastic_rosato_metadata.tsv"

## Load data to create Phyloseq Object
sediment <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

## Marine plastics Differential abundance
sediment_DA <- OTUs_multi_DA(sediment)

# An example vignette is provided (https://github.com/kmanoharan01/ConsensusMetaDA/blob/main/inst/doc/ConsensusMetaDA_Manual.html)
