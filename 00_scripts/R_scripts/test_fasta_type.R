# Load the RDS file
GCF_011617245.1 <- readRDS("01_inputs/04_fastas/GCF_011617245.1.rds")

# Get the accession name (assuming it's the first/only name in the list)
accession <- "GCF_011617245.1"

# Extract genome and protein sequences
# (These paths may need adjustment based on your actual structure)
genome_fasta <- GCF_011617245.1$genome  # or however genome is stored
protein_fasta <- GCF_011617245.1$protein  # or however protein is stored

# Write genome fasta to .fna file
writeLines(genome_fasta, paste0(accession, ".fna"))

# Write protein fasta to .faa file  
writeLines(protein_fasta, paste0(accession, ".faa"))
