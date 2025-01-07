rm(list=ls())
setwd("/Volumes/Biological_Age/miRNA")
library("haven")
library("dplyr")
library("readxl")
library("DESeq2")
library("EnhancedVolcano")

rna12 = read_spss("Data/591miRNAs_WellExpressed_2000_RS1&2.sav")
rnacounts12 = read_excel("Data/BeforeCPM/JPC-17-027_miRNA_RawCounts_Donor.xlsx")
rna25pct = read_excel("Data/miRNA normalization data.xlsx", sheet = "n2000 25pct cut 687 miRNA")
rna25pct$miRNA <- gsub("-", "", rna25pct$miRNA)
rnatr=data.frame(read.csv("Data/rnaallageJulia.csv"))
rnatr$ergoid_subset = paste0("R.", rnatr$ergoid)
rnacounts12 = dplyr::select(rnacounts12, c(miRNA, all_of(rnatr$ergoid_subset))) #Remove outliers
rnacounts04 = read_excel("Data/BeforeCPM/JPC-17-046_miRNA_RawCounts_Donor.xlsx")

# Covariates and outcomes
trainadj= data.table::fread("Data/train_standardadjustments.csv", data.table = F)
testadj = data.table::fread("Data/test_standardadjustments.csv", data.table = F)
extraadj= data.table::fread("Data/ex_standardadjustments.csv", data.table=F)
trainadj = left_join(trainadj, rnatr[,c("ergoid", "PhenoAge")])
testadj = left_join(testadj, rnatr[,c("ergoid", "PhenoAge")])

trainadj$AGE = scale(trainadj$AGE)
testadj$AGE  = scale(testadj$AGE)
extraadj$AGE = scale(extraadj$AGE)
trainadj$PhenoAge = scale(trainadj$PhenoAge)
testadj$PhenoAge  = scale(testadj$PhenoAge)
trainadj$event = factor(trainadj$event, levels = c(0, 1))
testadj$event  = factor(testadj$event, levels = c(0, 1))

rownames(trainadj)= paste0("R.", trainadj$ergoid)
rownames(testadj) = paste0("R.", testadj$ergoid)
rownames(extraadj)= paste0("R.", extraadj$ergoid)

# Extended function to calculate differential expression and extract test characteristics
calculate_differential_expression <- function(df, meta_data, wellexpressed = F, cutoff = NA, outcome_var, validation = F, previousdf = NA, outcome_type) {
 
  # Create a mapping for the original miRNA names
  miRNA_map <- setNames(df$miRNA, gsub("-", "", df$miRNA))
  
  # Remove hyphens from the 'miRNA' column
  df$miRNA <- gsub("-", "", df$miRNA)
  
  # Determine the dataset type for naming
  if (!validation) {
    dataset_type <- "Training"
  } else {
    if (4 %in% unique(meta_data$rs_cohort)) {
      dataset_type <- "Validation"
    } else {
      dataset_type <- "Test"
    }
  }

  if (wellexpressed) {
    if (cutoff == 50) {
      # Filter rows based on 'miRNA' names present in 'rna12'
      data <- df[df$miRNA %in% names(rna12), ]
    } else if (cutoff == 25) {
      # Filter rows based on 'miRNA' names present in 'rna25pct'
      data <- df[df$miRNA %in% rna25pct$miRNA, ]
    }
  } else if (validation) {
    previousdf_fdrsig = subset(previousdf, padj < 0.05)
    previousdf_fdrsig$miRNAkeep = gsub("-", "", previousdf_fdrsig$miRNA)
    data <- df[df$miRNA %in% previousdf_fdrsig$miRNAkeep, ]
  } else {
    data <- df[!grepl("HK_|CTRL_", df$miRNA), ]   
  }
  
  # Check covariates based on outcome and rs_cohort
  outcome_covariates <- if (outcome_var == "AGE") {
    c("sex", "lympho", "mono", "RBC", "border", "PlateNR")
  } else {
    c("sex", "AGE", "lympho", "mono", "RBC", "border", "PlateNR")
  }
  
  # Add rs_cohort to covariates if it has more than one level
  if (length(unique(meta_data$rs_cohort)) > 1) {
    outcome_covariates <- c(outcome_covariates, "rs_cohort")
  }
  
  meta_data = meta_data[complete.cases(meta_data[, c(outcome_var, outcome_covariates)]), ]
  meta_data$sex <- factor(meta_data$sex)
  meta_data$border <- factor(meta_data$border)
  meta_data$PlateNR<- factor(meta_data$PlateNR)
  meta_data$rs_cohort<- factor(meta_data$rs_cohort)
  meta_data$lympho <- scale(meta_data$lympho)
  meta_data$mono   <- scale(meta_data$mono)
  meta_data$RBC    <- scale(meta_data$RBC)
  
  rownames(data) = rownamessave = data$miRNA
  data = dplyr::select(data, rownames(meta_data))
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = data, colData = meta_data, design = as.formula(paste("~", paste(c(outcome_covariates, outcome_var), collapse = "+"))))
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Get results
  if (outcome_var == "event") {
    res <- results(dds, contrast = c("event", 1, 0)) 
  #  print(head(res))
  } else {
    res <- results(dds, name = outcome_var) 
    }
  
  
  rownames(res) = miRNA_map[rownamessave]
  
  # Adjust p-values
  res$padj <- p.adjust(res$pvalue, method = "fdr")
  res = as.data.frame(res)
  res$miRNA = rownames(res)
  res = dplyr::select(res, c(miRNA, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj))
  
  if(dataset_type == "Test" & outcome_type == "Standard50"){
    datatestdown = subset(res, log2FoldChange < 0 & padj < 0.05)
    datatestup   = subset(res, log2FoldChange > 0 & padj < 0.05)
    sigtest = full_join(datatestdown, datatestup)
    # Order by padj and keep only the top 20 rows with the lowest padj values if >20 sig in test set
    if(nrow(sigtest) > 20){
      sigtest <- sigtest %>%
        arrange(padj) %>%
        slice_head(n = 20)
      subtitle_plot = 'Vulcano plot in the training set\nAnnotated miRNAs represent the top 20 miRNAs in the test set'
    } else{
      subtitle_plot = paste('Vulcano plot in the training set\nAnnotated miRNAs represent the', nrow(sigtest),'miRNAs significantly associated in the test set')
      }
    selectlabstest = sigtest$miRNA
    
    # Determine plot title based on outcome_var
    plot_title <- switch(outcome_var,
                         "AGE" = "Chronological Age",
                         "event" = "10-Year Mortality",
                         "FI" = "the Frailty Index",
                         "PhenoAge" = "PhenoAge",
                         outcome_var)  # Default to outcome_var if no match
   
    # Assign colors and names for the legend
    keyvals.colour <- ifelse(
      previousdf$padj >= 0.05, "black",  # Not significant
      ifelse(
        previousdf$padj < 0.05 & !previousdf$miRNA %in% c(datatestdown$miRNA, datatestup$miRNA), "lightgrey",  # Significant only in training
        ifelse(
          previousdf$padj < 0.05 & previousdf$miRNA %in% datatestdown$miRNA, "blue",  # Down-regulated
          "red"  # Up-regulated
        )
      )
    )
    
    # Replace NAs with "black"
    keyvals.colour[is.na(keyvals.colour)] <- "black"
    
    # Assign names for the legend
    names(keyvals.colour)[keyvals.colour == "black"] <- "Not significant"
    names(keyvals.colour)[keyvals.colour == "lightgrey"] <- "Only significant in training set"
    names(keyvals.colour)[keyvals.colour == "blue"] <- "Down-regulated"
    names(keyvals.colour)[keyvals.colour == "red"] <- "Up-regulated"
    print(head(previousdf))
    volcano_plot <- EnhancedVolcano(previousdf, 
                                    lab = rownames(previousdf), 
                                    selectLab = selectlabstest,
                                    labCol = 'black',
                                    labFace = 'italic',
                                    boxedLabels = TRUE,
                                    x = 'log2FoldChange',  # No fold change
                                    y = 'padj', 
                                    title = paste('Differential Expression Analysis of', plot_title), 
                                    subtitle = subtitle_plot, 
                                    ylab = bquote(~-Log[10] ~ italic(pFDR)), 
                                    pCutoff = 0.05,
                                    pointSize = 3.0,
                                    labSize = 6.0,
                                    xlim = c(min(previousdf$log2FoldChange, na.rm = TRUE) - 0.5, max(previousdf$log2FoldChange, na.rm = TRUE) +
                                               0.5),
                                    ylim = c(0, max(-log10(previousdf$padj), na.rm = TRUE) + 1),
                                    colCustom = keyvals.colour,
                                    max.overlaps = Inf,
                                    maxoverlapsConnectors = Inf,
                                    drawConnectors = TRUE,
                                    typeConnectors = 'closed',
                                    widthConnectors = 0.75)
    save(volcano_plot, file = paste0("./Results/Rebuttal/DE_",dataset_type, "_",  outcome_var,"_", outcome_type, ".RData"))
    
    ggsave(paste0("./Results/Rebuttal/DE_",dataset_type, "_",  outcome_var,"_", outcome_type, ".jpg"), plot = volcano_plot, units = "in", width = 17, height = 11)
    
    
  }
  
  # Construct the filename based on input parameters
  dataset_name <- paste0("Results/Rebuttal/DifferentialExpression/DE_", 
                         outcome_var, 
                         ifelse(!is.null(outcome_type), paste0("_", outcome_type), ""),
                         "_", dataset_type,
                         ".csv")
    # Write results to CSV
  write.csv(res, dataset_name, row.names = F)
  # Convert to data frame and return
  return(res)
}

#### Age ####
# Age ~ 591 miRNA ~ standard set 50% cutoff
trainage <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "AGE", wellexpressed = T, cutoff = 50, outcome_type = "Standard50")
testage  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "AGE", validation = T, previousdf = trainage, outcome_type = "Standard50")
validage <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "AGE", validation = T, previousdf = testage, outcome_type = "Standard50")

# Age ~ 687 miRNA ~ sensitivity 25% cutoff
trainage25 <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "AGE", wellexpressed = T, cutoff = 25, outcome_type = "Sensitivity25")
testage25  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "AGE", validation = T, previousdf = trainage25, outcome_type = "Sensitivity25")
validage25 <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "AGE", validation = T, previousdf = testage25, outcome_type = "Sensitivity25")

# Age ~ 2083 miRNA ~ sensitivity no cutoff
trainageall <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "AGE", wellexpressed = F, outcome_type = "NoCutoff")
testageall  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "AGE", validation = T, previousdf = trainageall, outcome_type = "NoCutoff")
validageall <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "AGE", validation = T, previousdf = testageall, outcome_type = "NoCutoff")

#### PhenoAge ####
# PhenoAge ~ 591 miRNA ~ standard set 50% cutoff
trainPA <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "PhenoAge", wellexpressed = T, cutoff = 50, outcome_type = "Standard50")
testPA  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "PhenoAge", validation = T,  previousdf = trainPA, outcome_type = "Standard50")
#validPA <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "PhenoAge", validation = T,  previousdf = testPA, outcome_type = "Standard50")

# PhenoAge ~ 687 miRNA ~ sensitivity 25% cutoff
trainPA25 <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "PhenoAge", wellexpressed = T, cutoff = 25, outcome_type = "Sensitivity25")
testPA25  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "PhenoAge", validation = T, previousdf = trainPA25, outcome_type = "Sensitivity25")
#validPA25 <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "PhenoAge", validation = T, previousdf = testPA25, outcome_type = "Sensitivity25")

# PhenoAge ~ 2083 miRNA ~ sensitivity no cutoff
trainPAall <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "PhenoAge", wellexpressed = F, outcome_type = "NoCutoff")
testPAall  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "PhenoAge", validation = T, previousdf = trainPAall, outcome_type = "NoCutoff")
#validPAall <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "PhenoAge", validation = T, previousdf = testPAall, outcome_type = "NoCutoff")

#### Frailty index ####
# Frailty index ~ 591 miRNA ~ standard set 50% cutoff
trainFI <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "FI", wellexpressed = T, cutoff = 50, outcome_type = "Standard50")
testFI  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "FI", validation = T,  previousdf = trainFI, outcome_type = "Standard50")
#validFI <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "FI", validation = T,  previousdf = testFI, outcome_type = "Standard50")

# Frailty index ~ 687 miRNA ~ sensitivity 25% cutoff
trainFI25 <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "FI", wellexpressed = T, cutoff = 25, outcome_type = "Sensitivity25")
testFI25  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "FI", validation = T, previousdf = trainFI25, outcome_type = "Sensitivity25")
#validFI25 <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "FI", validation = T, previousdf = testFI25, outcome_type = "Sensitivity25")

# Frailty index ~ 2083 miRNA ~ sensitivity no cutoff
trainFIall <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "FI", wellexpressed = F, outcome_type = "NoCutoff")
testFIall  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "FI", validation = T, previousdf = trainFIall, outcome_type = "NoCutoff")
#validFIall <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "FI", validation = T, previousdf = testFIall, outcome_type = "NoCutoff")

#### Mortality ####
# Mortality ~ 591 miRNA ~ standard set 50% cutoff
trainmort <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "event", wellexpressed = T, cutoff = 50, outcome_type = "Standard50")
testmort  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "event", validation = T,  previousdf = trainmort, outcome_type = "Standard50")
#validmort <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj , outcome_var = "event", validation = T,  previousdf = testmort, outcome_type = "Standard50")

# Mortality ~ 687 miRNA ~ sensitivity 25% cutoff
trainmort25 <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "event", wellexpressed = T, cutoff = 25, outcome_type = "Sensitivity25")
testmort25  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "event", validation = T, previousdf = trainmort25, outcome_type = "Sensitivity25")
#validmort25 <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "event", validation = T, previousdf = testmort25, outcome_type = "Sensitivity25")

# Mortality ~ 2083 miRNA ~ sensitivity no cutoff
trainmortall <- calculate_differential_expression(df = rnacounts12, meta_data = trainadj, outcome_var = "event", wellexpressed = F, outcome_type = "NoCutoff")
testmortall  <- calculate_differential_expression(df = rnacounts12, meta_data = testadj , outcome_var = "event", validation = T, previousdf = trainmortall, outcome_type = "NoCutoff")
#validmortall <- calculate_differential_expression(df = rnacounts04, meta_data = extraadj, outcome_var = "event", validation = T, previousdf = testmortall, outcome_type = "NoCutoff")

#### Figures ####
library("ComplexUpset")

csv_files <- list.files(path = "./Results/Rebuttal/DifferentialExpression", pattern = "\\.csv$", full.names = TRUE)


# Read each CSV file and assign to a variable in the environment
lapply(csv_files, function(file) {
  data <- read.csv(file, row.names = 1)
  assign(sub(".csv$", "", basename(file)), data, envir = .GlobalEnv)
  
})

# Function to filter miRNAs with padj < 0.05 for a given dataframe
make_sig <- function(df, outcome) {
  df$miRNA = rownames(df)
  df = subset(df, padj < 0.05)
  df = df %>% select(., miRNA)
  df[[outcome]] = 1
  return(df)
}

# Apply function to each dataset and store results
sigage <- make_sig(DE_AGE_Standard50_Test, "Age")
sigPA <- make_sig(DE_PhenoAge_Standard50_Test, "PA")
sigFI <- make_sig(DE_FI_Standard50_Test, "FI")
sigmort <- make_sig(DE_event_Standard50_Test, "mort")


DE_sig = list(sigage, sigPA, sigFI, sigmort)
allupsetRuse = Reduce(function(x,y)full_join(x,y), DE_sig)
allupsetRuse = allupsetRuse %>% select(., -miRNA)
allupsetRuse[is.na(allupsetRuse)] <- 0
set_vars <- c("Age", "PhenoAge", "Frailty Index", "10-Year Mortality")
names(allupsetRuse) = set_vars

#Differential Expression
DE_upset = ComplexUpset::upset(allupsetRuse, 
                                set_vars,
                                width_ratio=0.4,
                                name = 'Combination',
                                group_by='degree',  
                                matrix=(
                                  intersection_matrix(geom=geom_point(shape='circle filled', size=5)) + 
                                    scale_color_manual(values=c('Age'="#ffb000", 'PhenoAge'='#fe6100', 'Frailty Index'='#dc267f',"10-Year Mortality" = "#361ae5" ),
                                                       guide=guide_legend(override.aes=list(shape='circle')))), 
                                queries=list(
                                  upset_query(set='Age', fill="#ffb000"),
                                  upset_query(set='PhenoAge', fill= "#fe6100"),
                                  upset_query(set='Frailty Index', fill="#dc267f"),
                                  upset_query(set="10-Year Mortality", fill= "#361ae5")),
                                themes=upset_default_themes(text=element_text(size=24)),
                                base_annotations=list('Number of miRNAs associated\nwith combination'=intersection_size(text = list(size=10), counts=TRUE)),
                                set_sizes = upset_set_size() + ylab('Number of miRNAs\nassociated with outcome'))

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}
plotage <- loadRData("Results/Rebuttal/DE_Test_AGE_Standard50.RData")
plotPA <- loadRData("Results/Rebuttal/DE_Test_PhenoAge_Standard50.RData")
plotFI <- loadRData("Results/Rebuttal/DE_Test_FI_Standard50.RData")
plotmort <- loadRData("Results/Rebuttal/DE_Test_event_Standard50.RData")

library("ggpubr")
fig2 = ggarrange(ggarrange(plotage, plotPA, plotFI, plotmort, nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"), font.label = list(size = 24, color = "black", face = "bold", family = NULL)), DE_upset, labels = c("", "e"), ncol = 2, font.label = list(size = 24, color = "black", face = "bold", family = NULL), widths = c(1.3,1))

ggsave(fig2, filename = "./Results/VulcanotrainAllPlots_Rebuttal.png", width = 33, height = 15, units='in', bg = "white")
ggsave(DE_upset, filename = "./Results/VennGraphicAbstr_Rebuttal.png", width = 8, height = 14, units='in', bg = "white")

#Biomarkers
library("ComplexUpset")
library("ggplot2")
biomarkers = read_excel("./Results/Results_miRNA_29Nov2024_EN.xlsx", sheet = "All")

biomarkersvenn <- biomarkers %>%
  mutate(across(-miRNA, ~ as.numeric(.))) %>%
  mutate(across(-miRNA, ~ ifelse(. != 0, 1, .))) %>%
  subset(., miRNA != "(Intercept)")%>%
  select(., -miRNA) %>%
  subset(., rowSums(.) > 0)
set_vars_en <- c("mirAge", "mirPA", "mirFI", "mirMort")
names(biomarkersvenn) = set_vars_en



upset_en = ComplexUpset::upset(biomarkersvenn, 
                               set_vars_en,
                                width_ratio=0.4,
                                name='Combination of biomarkers',
                                group_by='degree',  
                                matrix=(
                                  intersection_matrix(geom=geom_point(shape='circle filled', size=5)) + 
                                    scale_color_manual(values=c('mirAge'="#ffb000", 'mirPA'='#fe6100', 'mirFI'='#dc267f',"mirMort" = "#361ae5" ),
                                                       guide=guide_legend(override.aes=list(shape='circle')))), 
                                queries=list(
                                  upset_query(set='mirAge', fill="#ffb000"),
                                  upset_query(set='mirPA', fill= "#fe6100"),
                                  upset_query(set='mirFI', fill="#dc267f"),
                                  upset_query(set='mirMort', fill= "#361ae5")),
                                themes=upset_default_themes(text=element_text(size=24)),
                                base_annotations=list('Number of miRNAs\nincluded in combination'=intersection_size(text = list(size=10), counts=TRUE)),
                                set_sizes = upset_set_size() + ylab('Number of miRNAs\nincluded in biomarker'))

ggsave(upset_en, filename = "Results/EN_Upsetplot.png", width = 13, height = 9, units = "in", dpi = 600, limitsize = F)
