source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/Load_libraries.r")
source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/dataload.r")
source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/marker_score_converter.r")
source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/chi_sq_mendel.r")
source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/significant_loci_mendel_test.r")
source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/percent_none_mendel.r")
source("C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/scripts/functions/merged_matrix.r")
# Usage load libraries
load_required_libraries()
library(dplyr)
library(package = "AlphaSimR")
library(package = "ade4")
library(package = "ggplot2")
library(package = "rrBLUP")
library(package = "qqman")
library(openxlsx)
library(lme4)
# Usage load datafiles 
your_map_object <- "C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/data_matrixes/input_quality/2024_gustavo_GMASS/R-files"
your_file_name <- "F1.BR.01300_one_brazil_2024.xlsx"
your_file_loc_name<-"Parents_Brazil_2024_loc.xlsx"
df  <- read_excel_to_df(your_map_object, your_file_name)
genome_loc<- read_excel_to_df(your_map_object, your_file_loc_name)

# Applying Min-Max Scaling by LG
# Retry the operation with explicit namespace prefix if necessary
genome_scaled <- genome_loc %>%
  dplyr::group_by(LG) %>%
  dplyr::mutate(Scaled_Position = (Position - min(Position)) / (max(Position) - min(Position)))
genome_scaled<-data.frame(genome_scaled)
# Using the explicit dplyr prefix
genMap_1 <- genome_scaled %>%
  dplyr::rename(marker = marker_ID, chromosome = LG, position = Scaled_Position) %>%
  dplyr::select(marker, chromosome, position) 

# Convert to a data frame if needed 
genMap_1 <- as.data.frame(genMap_1)

# Usage convert DNA
converted_score <- convert_to_converted_score(df)
#


# # Assign marker names to the haplotypes
# # These can be in any order, because the software will order them based
# # on the genetic map automatically
# colnames(converted_score) = genMap_1$marker

#############allinging the genetic map with the marker data

# Creating a vector of markers from genMap_1 that should be kept in haplo
markers_to_keep <- genMap_1$marker

# Subsetting the haplo data frame to keep only the columns that are present in the markers_to_keep vector
haplo_filtered <- converted_score [, colnames(converted_score ) %in% markers_to_keep]

# Convert matrix to data frame
haplo_filtered <- as.data.frame(haplo_filtered)
# Adjust values to be within the 0-255 range
haplo_filtered  <- pmin(pmax(haplo_filtered , 0), 255)

# Replace NA with median per column using dplyr
# Replace NA with median per column using dplyr
haplo_filtered <- haplo_filtered %>%
  dplyr::mutate(dplyr::across(everything(), ~ifelse(is.na(.), median(., na.rm = TRUE), .)))

####converting the data to a genotype file for AlphsimR
#### this means the genotype represents allele dosage 0, 1, 2

# Replace all 1s with "2m"
allele_dose<-haplo_filtered
allele_dose[allele_dose == 1] <- 2

# Replace all 0s with "1"
allele_dose[allele_dose == 0] <- 1

# Finally, replace all -1s with "0"
allele_dose[allele_dose == -1] <- 0

# Print the updated data frame
 
  
# Ensure the list contains only numeric vectors and then bind them into a matrix
if (is.list(allele_dose) && all(sapply(allele_dose, is.numeric))) {
    allele_dose_matrix <- do.call(cbind, allele_dose)
} else {
    stop("haplo_filtered is not a list or contains non-numeric elements or elements of varying lengths")
}

# Remove columns where all values are NA
#allele_dose_matrix <- allele_dose_matrix[, colSums(!is.na(allele_dose_matrix)) > 0]
#####handeling missing data
# Assuming your first matrix is named matrix1
columns_with_na <- colnames(allele_dose_matrix)[colSums(is.na(allele_dose_matrix)) > 0]
# Remove rows where the 'marker' value is in columns_with_na
genMap_no_NA <- genMap_1[!genMap_1$marker %in% columns_with_na, ]

# Remove columns from allele_dose_matrix that have names listed in columns_with_na
allele_dose_matrix_no_NA <- allele_dose_matrix[, !colnames(allele_dose_matrix) %in% columns_with_na]


# Create the founder population
# we have found an issue in our code for the data import loading of alphsimR
#####when we load a haplotype we expect two rows per individual that are sumed to make the haplotype
#### so we need to modify the haplo_matrix to a geno matrix
founderPop = importInbredGeno(
  geno  =  allele_dose_matrix_no_NA,
  genMap = genMap_no_NA,
 
)

SP = SimParam$new(founderPop)

# Accessing marker data
####this is how to pull all the makers
parents_geno<-pullSegSiteGeno(founderPop )

###set new population
pop = newPop(founderPop, simParam=SP)
# # Perform a full diallel cross among the population 
pop2 = hybridCross(pop, pop, simParam=SP)
#pull the genotypes of the alleles
F1_geno<-pullSegSiteGeno(pop2)

parents_geno<-data.frame(t(parents_geno))
F1_geno<-data.frame(t(F1_geno))

# Write 'converted_score' to an Excel file with row names
write.xlsx(parents_geno, "C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/data_matrixes/input_quality/2024_gustavo_GMASS/output/parents_geno_sim.xlsx",rowNames = TRUE)
write.xlsx(F1_geno, "C:/Users/joseg/OneDrive - Enza Zaden Beheer B.V/File_Platform/data_matrixes/input_quality/2024_gustavo_GMASS/output/F1_geno_sim.xlsx",rowNames = TRUE)



# ######################################################################setting up the genome
# # Force QTLs and SNPs to overlap
# SP$restrSegSites(overlap = T)

# # Create additive trait with 5 QTLs per chromosome
# SP$addTraitAG(nQtlPerChr = 5,
              # mean = 0,
              # var  = 1)

# # Add SNP-chip with 50 SNPs per chromosome, if we have 12 chromossomes we have 600 snps
# SP$addSnpChip(nSnpPerChr = 50)

# # Check that QTLs and SNPs overlap
# sum(getQtlMap(trait = 1)$id %in% getSnpMap(1)$id)
# founder_hap = pullSegSiteHaplo(founderPop)
# # Extract the genotypes of individuals
# founder_geno<-pullSegSiteGeno(founderPop)


# # Assuming founderPop is already generated

# # ---- Create two subpopulations ----


# # Create five subpopulations

# # Create a new population from the founder haplotypes
# pop = newPop(founderPop, simParam=SP)

# # Perform a full diallel cross among the population 
# pop2 = hybridCross(pop, pop, simParam=SP)

# # Assuming geno is already pulled for the full population
# #geno = pullQtlGeno(pop2)
# ####this wil get only the snps, 
# geno=pullSnpGeno(pop2, simParam=SP)

# # Accessing marker data
# ####this is how to pull all the makers
# hybrid_geno<-pullSegSiteGeno(pop2)

# PCA  = dudi.pca(df = hybrid_geno, center = TRUE, scale = FALSE, scannf = FALSE, nf = 5)
# (VAF = 100 * PCA$eig[1:5] / sum(PCA$eig)) # Variance explained
# # Assuming 'PCA' is already defined correctly
# # and 'VAF' is the vector of variance accounted for by each principal component

# # Create a data frame with PCA scores for the dataset
# df.PCA = data.frame(
  # "PC1" = PCA$l1[, 1],  # Accessing the first Principal Component scores
  # "PC2" = PCA$l1[, 2]   # Accessing the second Principal Component scores
# )

# # Plotting the PCA results
# library(ggplot2)  # Load ggplot2 for plotting

# ggplot(df.PCA, aes(x = PC1, y = PC2)) +
  # geom_point() +
  # ggtitle("PCA of Dataset") +
  # xlab(paste("PC1: ", round(VAF[1], 2), "% Variance Explained", sep = "")) +
  # ylab(paste("PC2: ", round(VAF[2], 2), "% Variance Explained", sep = ""))


# # geno1<-geno
# # geno2<-geno