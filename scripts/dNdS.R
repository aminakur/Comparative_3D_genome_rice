setwd('/scratch/ak8725/orthologr/')

# compute dN/dS values
library(orthologr, lib.loc = '/scratch/ak8725/R/x86_64-pc-linux-gnu/4.3')

CDS_dNdS <- 
dNdS(query_file = 
'Oryza_rufipogon.OR_W1943.cds.all.fa',
     subject_file = 
'Oryza_meridionalis.Oryza_meridionalis_v1.3.cds.all.fa',
     aligner = "blast",
     aligner_path = "/share/apps/blast+/2.15.0/bin/",	
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 16 )

# store result in Excel readable csv file
#install.packages("readr")
library(readr)
readr::write_excel_csv(CDS_dNdS, "oruf_omer_dNdS.csv")

#analyze the results
#1. NA dS values are ignored

#function for calculating mode
Mode <- function(x) {
  # Remove NA values from the vector
  x <- x[!is.na(x)]
  
  # Check if vector is empty
  if(length(x) == 0) {
    return(NA)
  }
  
  # Calculate mode
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  mode_index <- which.max(tab)
  
  # Return mode
  return(ux[mode_index])
}

# Calculate mean, median, and mode
mean_dS <- mean(CDS_dNdS$dS, na.rm = TRUE)
median_dS <- median(CDS_dNdS$dS, na.rm = TRUE)
mode_dS <- Mode(CDS_dNdS$dS)

#2.NA dS values are treated as 0s
# Replace NA values with 0 in the 'dNdS' column
dS_with_zeros <- ifelse(is.na(CDS_dNdS$dS), 0, CDS_dNdS$dS)

# Calculate mean, median, and mode
mean_dS_zeros <- mean(dS_with_zeros)
median_dS_zeros <- median(dS_with_zeros)
mode_dS_zeros <- Mode(dS_with_zeros)

# Count non-NA values in the 'dS' column before replacing NA values with 0
non_na_count_original <- sum(!is.na(CDS_dNdS$dS))

#I get 8 values as output of 1 and 2
cat("Number of orthologous gene pairs:", nrow(CDS_dNdS), "\n")
cat("Number of non-NA dS values (including 0):", 
non_na_count_original, "\n")
cat("Mean dS:", mean_dS, "\n")
cat("Median dS:", median_dS, "\n")
cat("Mode dS:", mode_dS, "\n")
cat("Mean dS with NAs as 0s:", mean_dS_zeros, "\n")
cat("Median dS with NAs as 0s:", median_dS_zeros, "\n")
cat("Mode dS with NAs as 0s:", mode_dS_zeros, "\n")
