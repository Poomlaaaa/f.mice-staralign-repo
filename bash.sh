# Clean environment
rm(list = ls(all.names = TRUE))

# Set working directory to where your script and count files are located
setwd("~/starmm")

# Source the script
source("combine_counts.R")

install.packages("BiocManager")  # Install BiocManager if not installed
BiocManager::install("rtracklayer")
BiocManager::install("GenomicFeatures")
install.packages("dplyr")  # Install dplyr for data manipulation

# Load necessary libraries
library(rtracklayer)
library(dplyr)

# Define the GTF file location
gtf_file <- "gencode.vM34.primary_assembly.annotation.gtf.gz"

# Import the GTF file
gtf <- import(gtf_file)

# Convert to a data frame
gtf_df <- as.data.frame(gtf)

# Extract gene lengths
gene_lengths <- gtf_df %>%
  filter(type == "gene") %>%
  mutate(Length = end - start + 1) %>%
  select(GeneID = gene_id, Length)

# Write gene lengths to a file
write.table(gene_lengths, "gene_lengths.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Load gene length data
gene_length_file <- "gene_lengths.txt"
gene_lengths <- read.table(gene_length_file, header=TRUE, stringsAsFactors=FALSE)

# Load combined count data
combined_counts <- read.table("combined_counts.txt", header=TRUE, stringsAsFactors=FALSE)

# Merge the gene lengths with the count data
combined_data <- merge(combined_counts, gene_lengths, by="GeneID")

# Function to calculate FPKM
calculate_fpkm <- function(counts, lengths, total_counts) {
  fpkm <- (counts / lengths) / (total_counts / 1e6)
  return(fpkm)
}

# Calculate total counts for each sample
total_counts <- colSums(combined_data[, -c(1, ncol(combined_data))])

# Calculate FPKM for each sample
fpkm_data <- combined_data %>%
  rowwise() %>%
  mutate(across(starts_with("SRR"), ~ calculate_fpkm(.x, Length, total_counts)))

# Look at the combined data file
str(combined_data)
View(combined_data)

# Function to calculate FPKM
calculate_fpkm <- function(counts, lengths, total_counts) {
  fpkm <- (counts / lengths) / (total_counts / 1e6)
  return(fpkm)
}

# Calculate FPKM for each sample explicitly
fpkm_data <- combined_data %>%
  mutate(
    SRR13893218_fpkm = calculate_fpkm(SRR13893218, Length, total_counts),
    SRR13893219_fpkm = calculate_fpkm(SRR13893219, Length, total_counts),
    SRR13893220_fpkm = calculate_fpkm(SRR13893220, Length, total_counts),
    SRR13893221_fpkm = calculate_fpkm(SRR13893221, Length, total_counts),
    SRR13893222_fpkm = calculate_fpkm(SRR13893222, Length, total_counts),
    SRR13893223_fpkm = calculate_fpkm(SRR13893223, Length, total_counts),
  ) %>%
  select(-Length)  # Remove the Length column after computation

# Write the FPKM data to a file
write.table(fpkm_data, "fpkm_counts.txt", sep="\t", quote=FALSE, row.names=FALSE)












 


