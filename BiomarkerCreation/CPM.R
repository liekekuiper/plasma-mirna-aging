######## Recalculating the CPM

###################################################################################################
####################
# Folder locations
####################
rm(list=ls())
setwd("/Volumes/Biological_Age/miRNA")
library("haven")
library("dplyr")
library("readxl")
library("edgeR")

rna12 = read_spss("Data/591miRNAs_WellExpressed_2000_RS1&2.sav")
rnacounts12 = read_excel("Data/BeforeCPM/JPC-17-027_miRNA_RawCounts_Donor.xlsx")
rna25pct = read_excel("Data/miRNA normalization data.xlsx", sheet = "n2000 25pct cut 687 miRNA")
rna25pct$miRNA <- gsub("-", "", rna25pct$miRNA)
rnatr=data.frame(read.csv("Data/rnaallageJulia.csv"))
rnatr$ergoid_subset = paste0("R.", rnatr$ergoid)
rnacounts12 = select(rnacounts12, c(miRNA, all_of(rnatr$ergoid_subset))) #Remove outliers
rnacounts04 = read_excel("Data/BeforeCPM/JPC-17-046_miRNA_RawCounts_Donor.xlsx")


ApplyCPM <- function(df, wellexpressed, cutoff = NA) {
  # Remove hyphens from the 'miRNA' column
  df$miRNA <- gsub("-", "", df$miRNA)
  
  if(wellexpressed){
    if(cutoff == 50){
  # Filter rows based on 'miRNA' names present in 'rna12'
  rnaCPM <- as.data.frame(df[df$miRNA %in% names(rna12), ])}
    else if (cutoff == 25){
  # Filter rows based on 'miRNA' names present in 'rna12'
  rnaCPM <- as.data.frame(df[df$miRNA %in% rna25pct$miRNA, ])}    
    }
  else{
    rnaCPM <- as.data.frame(df[!grepl("HK_|CTRL_", df$miRNA),])   
  }
  
 # Set the 'miRNA' column as row names and remove it from the data frame
  row.names(rnaCPM) <- rnaCPM$miRNA
  rnaCPM <- select(rnaCPM, -miRNA)
  
  ## Normalisation by the TMM method (Trimmed Mean of M-value)
  dge <- DGEList(rnaCPM)                        # DGEList object created from the count data
  calc.normfact.out <- calcNormFactors(dge, method = "TMM")  
  
  # Calculate CPM values and set values less than 1 to 0
  cpm.tmm.log2 <- cpm(calc.normfact.out, log=TRUE, normalized.lib.sizes = TRUE)
  cpm.tmm.log2[cpm.tmm.log2 < 1.0] <- 0 
  
  # Transpose the data frame and convert it back to a data frame
  transposed_df <- as.data.frame(t(cpm.tmm.log2))
  
  # Set appropriate column and row names
  colnames(transposed_df) <- rownames(rnaCPM)
  rownames(transposed_df) <- colnames(rnaCPM)
  
  # Make ergoid by removing 'R.' from row names
  transposed_df$ergoid <- gsub("R.", "", rownames(transposed_df))
  # Remove row names
  rownames(transposed_df) <- NULL
  
  
  # Re-order dataframe
  transposed_df <- transposed_df[,c("ergoid", names(transposed_df)[1:(ncol(transposed_df)-1)])]
  
  return(transposed_df)
}

cpm12 = ApplyCPM(rnacounts12, T, 50)
cpm04 = ApplyCPM(rnacounts04, T, 50)
cpm12sens = ApplyCPM(rnacounts12, F)
cpm12sens25 = ApplyCPM(rnacounts12, T, 25)

write.csv(cpm12, "Data/cpm_rs12.csv")
write.csv(cpm04, "Data/cpm_rs4.csv")
write.csv(cpm12sens, "Data/cpm_rs12_sens.csv")
write.csv(cpm12sens25, "Data/cpm_rs12_sens25.csv")


rna12 = subset(rna12, ergoid %in% cpm12$ergoid)
rna12 = rna12[order(rna12$ergoid),]
# Check if the 'ergoid' column matches between the two dataframes
if (nrow(rna12[rna12$ergoid != cpm12$ergoid, ]) > 1) {
  stop("The 'ergoid' column does not match between the two dataframes.")
}

# Initialize a vector to store the correlation values and their column names
correlation_data <- data.frame(
  MiRNA = names(cpm12)[-1], # Exclude 'ergoid'
  Correlation = numeric(length = 591)  # Assuming 592 columns, one for 'ergoid'
)

# Calculate correlations for each pair of columns (excluding 'ergoid')
for (i in 2:592) {
  col_cpm <- cpm12[, i]
  col_rna <- rna12[rna12$ergoid %in% cpm12$ergoid, i]
  correlation_data$Correlation[i - 1] <- cor(col_cpm, col_rna)
}

# Find the column name with the lowest and highest correlation
lowest_corr_column <- correlation_data$MiRNA[which.min(correlation_data$Correlation)] #miR-4461
highest_corr_column <- correlation_data$MiRNA[which.max(correlation_data$Correlation)] #miR-6851-3p

write.table(correlation_data, file = "Results/CorrelationCPM.txt", row.names = F)
