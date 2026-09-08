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

################### Figure 3 #####################
wget ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_150bp.subset_2k.rare_5000.biom or 

Download via https://ftp.microbio.me/emp/release1/otu_tables/deblur/ 

file name: emp_deblur_150bp.subset_2k.rare_5000.biom

or
unzip ../Data/emp_deblur_150bp.subset_2k.rare_5000.biom

## emp
## Load data
biome_file <- "./Data/emp_deblur_150bp.subset_2k.rare_5000.biom"

sample_table_file <-  "./Data/samples_table_emp.txt"


## Load data to create phyloseq object
emp <- build_OTU_counts_NoTaxa(biom = biome_file, sample_table = sample_table_file)

## Plot using phyloseq object
emp <- OTU_plots(emp)


# An example vignette is provided (https://github.com/kmanoharan01/ConsensusMetaDA/blob/main/inst/doc/ConsensusMetaDA_Manual.html)
