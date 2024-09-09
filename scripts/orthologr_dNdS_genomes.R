setwd("/Users/aminosaurier/Downloads/dNdS")

# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# Install package dependencies
BiocManager::install(c(
  "Biostrings",
  "GenomicRanges",
  "GenomicFeatures",
  "Rsamtools",
  "rtracklayer"
))

# install CRAN dependencies
install.packages(c("doParallel", "foreach", "ape", "Rdpack", "benchmarkme", "devtools"))

# install BLAST dependency metablastr from GitHub
devtools::install_github("drostlab/metablastr")

# install DIAMOND dependency rdiamond from GitHub
devtools::install_github("drostlab/rdiamond")

# install orthologr from GitHub
devtools::install_github("drostlab/orthologr")

# compute dN/dS values
library(orthologr)
CDS_dNdS <- 
dNdS(query_file = '/Users/aminosaurier/Downloads/cds/Oryza_sativa.IRGSP-1.0.cds.all.fa',
     subject_file = '/Users/aminosaurier/Downloads/cds/Oryza_sativa_azucena.AzucenaRS1.cds.all.fa',
     aligner = "diamond",
     ortho_detection = "RBH", 
     aa_aln_type     = "pairwise",
     aa_aln_tool     = "NW", 
     codon_aln_tool  = "pal2nal", 
     dnds_est.method = "Comeron", 
     comp_cores      = 1 )

# store result in Excel readable csv file
#install.packages("readr")
library(readr)
readr::write_excel_csv(CDS_dNdS, "NPB_Az_dNdS.csv")

#analyze the results
CDS_dNdS <- read.csv("/Users/aminosaurier/Downloads/dNdS/Az_oruf_dNdS.csv")

#1. NA dNdS values are ignored

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

# Calculate mean, median, and mode dS
mean_dS <- mean(CDS_dNdS$dS, na.rm = TRUE)
median_dS <- median(CDS_dNdS$dS, na.rm = TRUE)
mode_dS <- Mode(CDS_dNdS$dS)

#2.NA dNdS values are treated as 0s
# Replace NA values with 0 in the 'dNdS' column
dS_with_zeros <- ifelse(is.na(CDS_dNdS$dS), 0, CDS_dNdS$dS)

# Calculate mean, median, and mode
mean_dS_zeros <- mean(dS_with_zeros)
median_dS_zeros <- median(dS_with_zeros)
mode_dS_zeros <- Mode(dS_with_zeros)

# Count non-NA values in the 'dNdS' column before replacing NA values with 0
non_na_count_original <- sum(!is.na(CDS_dNdS$dS))

#I get 8 values as output of 1 and 2
cat("Number of orthologous gene pairs:", nrow(CDS_dNdS), "\n")
cat("Number of non-NA dS values (including 0):", non_na_count_original, "\n")
cat("Mean dS:", mean_dS, "\n")
cat("Median dS:", median_dS, "\n")
cat("Mode dS:", mode_dS, "\n")
cat("Mean dS with NAs as 0s:", mean_dS_zeros, "\n")
cat("Median dS with NAs as 0s:", median_dS_zeros, "\n")
cat("Mode dS with NAs as 0s:", mode_dS_zeros, "\n")



#3.Plot the results
install.packages("openxlsx")
library(openxlsx)

getSheetNames("/Users/aminosaurier/Downloads/dNdS/pairwise.xlsx")

pairwise <- read.xlsx("/Users/aminosaurier/Downloads/dNdS/pairwise.xlsx",  sheet = "Sheet1", rowNames = FALSE)

# Extract the first column (variable names)
variable_names <- pairwise[[1]]

# Remove the first column from the data frame
pairwise <- pairwise[, -1]

# Set the column names of the data frame using the variable names
rownames(pairwise) <- variable_names

plot(pairwise[, c(1:9)])

#NA0 means NAs were replaced by 0s (not deleted)
plot(SSIM ~ mean_dS, data = pairwise)
plot(SSIM ~ median_dS, data = pairwise)
plot(SSIM ~ NA0_mean_dS, data = pairwise)
plot(SSIM ~ NA0_median_dS, data = pairwise)

#1
# Calculate the correlation coefficient between SSIM and mean_dS
correlation_coefficient <- cor(pairwise$SSIM, pairwise$mean_dS)

# Perform a correlation test to get the p-value
correlation_test <- cor.test(pairwise$SSIM, pairwise$mean_dS)

# Extract the p-value from the correlation test
p_value <- correlation_test$p.value

# Plot SSIM against mean_dS
plot(SSIM ~ mean_dS, data = pairwise)

# Add text annotation with the correlation coefficient and p-value
text(x = mean(pairwise$mean_dS), y = 0.4, 
     labels = paste("r:", round(correlation_coefficient, 2),
                    "\np-value:", round(p_value, 4)),
     pos = 4, cex = 0.8, col = "blue")


#2
# Calculate the correlation coefficient between SSIM and mean_dS
correlation_coefficient2 <- cor(pairwise$SSIM, pairwise$median_dS)

# Perform a correlation test to get the p-value
correlation_test2 <- cor.test(pairwise$SSIM, pairwise$median_dS)

# Extract the p-value from the correlation test
p_value <- correlation_test2$p.value

# Plot SSIM against mean_dS
plot(SSIM ~ median_dS, data = pairwise)

# Add text annotation with the correlation coefficient and p-value
text(x = mean(pairwise$median_dS), y = 0.4, 
     labels = paste("r:", round(correlation_coefficient2, 2),
                    "\np-value:", round(p_value, 4)),
     pos = 4, cex = 0.8, col = "green")




#3
# Calculate the correlation coefficient between SSIM and mean_dS
correlation_coefficient <- cor(pairwise$SSIM, pairwise$NA0_mean_dS)

# Perform a correlation test to get the p-value
correlation_test <- cor.test(pairwise$SSIM, pairwise$NA0_mean_dS)

# Extract the p-value from the correlation test
p_value <- correlation_test$p.value

# Plot SSIM against mean_dS
plot(SSIM ~ NA0_mean_dS, data = pairwise)

# Add text annotation with the correlation coefficient and p-value
text(x = mean(pairwise$NA0_mean_dS), y = 0.4, 
     labels = paste("r:", round(correlation_coefficient, 2),
                    "\np-value:", round(p_value, 4)),
     pos = 4, cex = 0.8, col = "red")


#4
# Calculate the correlation coefficient between SSIM and mean_dS
correlation_coefficient2 <- cor(pairwise$SSIM, pairwise$NA0_median_dS)

# Perform a correlation test to get the p-value
correlation_test2 <- cor.test(pairwise$SSIM, pairwise$NA0_median_dS)

# Extract the p-value from the correlation test
p_value <- correlation_test2$p.value

# Plot SSIM against mean_dS
plot(SSIM ~ NA0_median_dS, data = pairwise)

# Add text annotation with the correlation coefficient and p-value
text(x = mean(pairwise$NA0_median_dS), y = 0.4, 
     labels = paste("r:", round(correlation_coefficient2, 2),
                    "\np-value:", round(p_value, 4)),
     pos = 4, cex = 0.8, col = "purple")