rm(list=ls())
setwd("/Volumes/Biological_Age/miRNA/Results/Rebuttal/DifferentialExpression/")
library("haven")
library("dplyr")
library("readxl")
library("DESeq2")
library("VennDiagram")
library("enrichplot")
library("ggpubr")
library("clusterProfiler")
library("org.Hs.eg.db") 
library("seqinr")

# List all CSV files in the directory
csv_files <- list.files(pattern = "*.csv")

# Read each CSV file and assign to a variable in the environment
lapply(csv_files, function(file) {
  data <- read.csv(file, row.names = 1)
  assign(sub(".csv", "", file), data, envir = .GlobalEnv)
})

#Get file with full miRNA names
rnacounts12 = read_excel("../../../Data/BeforeCPM/JPC-17-027_miRNA_RawCounts_Donor.xlsx")
rnacounts12$miRNA2 = rnacounts12$miRNA
rnacounts12$miRNA =  gsub("-", "", rnacounts12$miRNA)


# Function to filter miRNAs with padj < 0.05 for a given dataframe
select_significant_miRNAs <- function(df) {
  df$miRNA2 = rownames(df)
  df = subset(df, padj < 0.05)
  df = left_join(df, rnacounts12[,c("miRNA", "miRNA2")])
}

all591mirnas = DE_AGE_Standard50_Training
all591mirnas$miRNA = rownames(DE_AGE_Standard50_Training)

# Apply function to each dataset and store results
significant_testage_miRNAs <- select_significant_miRNAs(DE_AGE_Standard50_Test)
significant_testPA_miRNAs <- select_significant_miRNAs(DE_PhenoAge_Standard50_Test)
significant_testFI_miRNAs <- select_significant_miRNAs(DE_FI_Standard50_Test)
significant_testmort_miRNAs <- select_significant_miRNAs(DE_event_Standard50_Test)

# Find the number of common miRNAs across all datasets
common_miRNAs <- Reduce(intersect, list(
  significant_testage_miRNAs$miRNA2,
  significant_testPA_miRNAs$miRNA2,
  significant_testFI_miRNAs$miRNA2,
  significant_testmort_miRNAs$miRNA2
))

# Print the count of common miRNAs
length(common_miRNAs)

# Select and rename only the required columns
age_subset <- significant_testage_miRNAs[significant_testage_miRNAs$miRNA2 %in% common_miRNAs, 
                                         c("miRNA2", "log2FoldChange", "padj")]
colnames(age_subset) <- c("miRNA2", "log2FoldChange_age", "padj_age")

PA_subset <- significant_testPA_miRNAs[significant_testPA_miRNAs$miRNA2 %in% common_miRNAs, 
                                       c("miRNA2", "log2FoldChange", "padj")]
colnames(PA_subset) <- c("miRNA2", "log2FoldChange_PA", "padj_PA")

FI_subset <- significant_testFI_miRNAs[significant_testFI_miRNAs$miRNA2 %in% common_miRNAs, 
                                       c("miRNA2", "log2FoldChange", "padj")]
colnames(FI_subset) <- c("miRNA2", "log2FoldChange_FI", "padj_FI")

mort_subset <- significant_testmort_miRNAs[significant_testmort_miRNAs$miRNA2 %in% common_miRNAs, 
                                           c("miRNA2", "log2FoldChange", "padj")]
colnames(mort_subset) <- c("miRNA2", "log2FoldChange_mort", "padj_mort")

# Merge all datasets by miRNA2
merged_data <- Reduce(function(x, y) merge(x, y, by = "miRNA2", all = TRUE), 
                      list(age_subset, PA_subset, FI_subset, mort_subset))

##Elastic net results

# Weights of biomarkers
addfile12 <- read_excel("../../Results_miRNA_29Nov2024_EN.xlsx", sheet = "All") %>%
  mutate(across(-miRNA, as.numeric)) %>%
  left_join(rnacounts12[, c("miRNA", "miRNA2")], by = "miRNA") %>%
  mutate(miRNA2 = case_when(is.na(miRNA2) ~ "Intercept", TRUE ~ miRNA2)) %>%
  dplyr::select(miRNA2, age, PhenoAge, FI, mortality) %>%
  rename("miRNA" = "miRNA2", 
         "mirAge" = "age", 
         "mirPA" = "PhenoAge", 
         "mirFI" = "FI",
         "mirMort" = "mortality")
openxlsx::write.xlsx(addfile12, "../../AdditionalFile12.xlsx")

# Venn diagram
# Load data and filter non-zero weights
load_data <- function(sheet_name) {
  read_excel("../../Results_miRNA_29Nov2024_EN.xlsx", sheet = sheet_name) %>%
    mutate(Weight = as.numeric(Weight)) %>%
    filter(Weight != 0) %>%
    left_join(rnacounts12[, c("miRNA", "miRNA2")])
}

enpredage <- load_data("age")
enpredPA <- load_data("PhenoAge")
enpredFI <- load_data("FI")
enpredmort <- load_data("mortality")

# Extract miRNA2 lists for Venn diagram
get_miRNA2 <- function(data) {
  data %>% filter(miRNA != "(Intercept)") %>% pull(miRNA2)
}

agevenn <- get_miRNA2(enpredage)
pavenn <- get_miRNA2(enpredPA)
fivenn <- get_miRNA2(enpredFI)
mortvenn <- get_miRNA2(enpredmort)

# Combine data and add direction labels
combined_data <- bind_rows(
  enpredage %>% mutate(Source = "mirAge"),
  enpredPA %>% mutate(Source = "mirPA"),
  enpredFI %>% mutate(Source = "mirFI"),
  enpredmort %>% mutate(Source = "mirMort")
) %>%
  mutate(Direction = ifelse(Weight > 0, "+", ifelse(Weight < 0, "-", NA)))

# Create contingency table
contingency_table <- combined_data %>%
  dplyr::select(miRNA2, Source, Weight) %>%
  pivot_wider(names_from = Source, values_from = Weight, values_fill = NA) %>%
  filter(!is.na(miRNA2)) %>%
  mutate(included_in = rowSums(!is.na(across(where(is.numeric)))))

# Identify miRNAs in >1 biomarker or with inconsistent direction
morethanone <- contingency_table %>% filter(included_in > 1)
not_consistent <- contingency_table %>%
  rowwise() %>%
  filter(all(c(-1, 1) %in% sign(across(where(is.numeric))))) %>%
  ungroup()

# Define colors and data for Venn diagram
colorblind <- c("#361ae5", "#dc267f", "#fe6100", "#ffb000")
venn_data <- list(
  "MiRNA in mirAge" = agevenn,
  "MiRNA in mirPA" = pavenn,
  "MiRNA in mirFI" = fivenn,
  "MiRNA in mirMort" = mortvenn
)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = venn_data,
  imagetype = 'png',
  height = 5, width = 8, units = "in",
  category.names = names(venn_data),
  col = colorblind, fill = colorblind, alpha = 0.4,
  cat.cex = 1, cat.fontfamily = "sans",
  main.fontfamily = "sans", fontfamily = "sans",
  main = "Venn Diagram of MiRNAs selected in aging biomarkers",
  main.cex = 1.5,
  filename = "../Venn_aging_biomarkers.png",
  resolution = 1200
)

miRNAs = unique(c(paste0("hsa-",significant_testage_miRNAs$miRNA2),
           paste0("hsa-",significant_testPA_miRNAs$miRNA2),
           paste0("hsa-",significant_testFI_miRNAs$miRNA2),
           paste0("hsa-",significant_testmort_miRNAs$miRNA2)))

#Genes
predicted = data.table::fread("../../../Data/BiologicalInterpretation/miRDB_v6.0_prediction_result.txt") #mirDB
predictedkeep = subset(predicted, V3 > 80 & V1 %in% paste0("hsa-", all591mirnas$miRNA))
names(predictedkeep) = c("miRNA", "Gene", "certainty")
predictedkeep = predictedkeep %>%
  mutate(Gene = case_when(Gene == "NM_001017973" ~ "NM_001365679",
                          Gene == "NM_001286549" ~ "NR_149155",
                          Gene == "NM_001324033" ~ "NR_172049",
                          Gene == "NM_001321382" ~ "NR_169872",
                          Gene == "NM_001317989" ~ "NR_161375",
                          Gene == "NM_001319291" ~ "NM_197962",
                          Gene == "NM_032374" ~ "NM_001370595",
                          Gene == "NM_001292023" ~ "NR_158984",
                          Gene == "NM_001242804" ~ "NR_164148",
                          Gene == "NM_001009612" ~ "NM_001394958",
                          Gene == "NM_020732" ~ "NM_001374820", 
                          Gene == "NM_001277348" ~ "NR_161235",
                          Gene == "NM_001243552" ~ "NR_164143",
                          Gene == "NM_207422" ~ NA_character_, #NM_207422.3: This RefSeq was removed because it is now thought that this gene does not encode a protein.
                          Gene == "NM_001242671" ~ "NR_161342",
                          Gene == "NM_001242502" ~ "NR_158991",
                          Gene == "NM_001302818" ~ "NR_161343",
                          Gene == "NM_001257177" ~ "NR_161333",
                          Gene == "NM_001302813" ~ "NR_161292",
                          Gene == "NM_001302815" ~ "NR_161293",
                          Gene == "NM_173644" ~ "NR_161291",
                          Gene == "NM_001305396" ~ "NR_161263",
                          Gene == "NM_001305397" ~ "NR_161262",
                          Gene == "NM_001305395" ~ "NR_161261",
                          Gene == "NM_173667" ~ "NR_161251",
                          Gene == "NM_001242575" ~ "NR_161186",
                          Gene == "NM_001001709" ~ "NR_161226",
                          Gene == "NM_182626" ~ NA_character_, #NR_161344.1: This RefSeq was removed because currently there is insufficient support for the transcript.
                          Gene ==  "NM_175918" ~ NA_character_,#NM_175918.3: This RefSeq was removed because the gene was discontinued.
                          Gene == "NM_000060" ~ "NM_001370658",
                          Gene == "NM_182626" ~ NA_character_, #NR_161344.1: This RefSeq was removed because currently there is insufficient support for the transcript.
                          Gene == "NM_001101337" ~ "NR_160966",
                          Gene == "NM_001256368" ~ "NR_161181",
                          Gene == "NM_001242780" ~ "NR_164141",
                          Gene == "NM_001322103" ~ "NR_161199",
                          Gene == "NM_006542" ~ NA_character_, #NM_006542.3: This RefSeq was removed because the gene was discontinued.
                          Gene == "NM_001304804" ~ "NR_158221",
                          Gene == "NM_033517" ~ "NM_001372044",
                          Gene == "NM_001291905" ~ "NR_160286",
                          Gene == "NM_001317033" ~ "NR_172094",
                          Gene == "NM_001291904" ~ "NR_160285",
                          Gene == "NM_001291379" ~ "NR_171664",
                          Gene == "NM_001256795" ~ "NR_164124",
                          Gene == "NM_182832" ~ "NR_148920",
                          Gene == "NM_001001891" ~ "NM_001370694",
                          Gene == "NM_001294306" ~ "NR_164115",
                          Gene == "NM_004709" ~ NA_character_, #NM_004709.2: This RefSeq was removed because it is 3' UTR sequence and is not a distinct protein coding transcript.
                          Gene == "NM_001014442" ~ "NR_172510",
                          Gene == "NM_001135816" ~ "NR_172511",
                          Gene == "NM_001242791" ~ "NR_164135",
                          Gene == "NM_001242901" ~ "NR_164163",
                          Gene == "NM_001281727" ~ "NR_168132",
                          Gene == "NM_001281728" ~ "NR_168133",
                          Gene == "NM_177953" ~ "NR_168131",
                          Gene == "NM_001281729" ~ "NR_168134",
                          TRUE ~ Gene)) %>% subset(., !is.na(Gene))

length(unique(predictedkeep$miRNA)) #571

predicted_age_genes = subset(predictedkeep, miRNA %in% paste0("hsa-",significant_testage_miRNAs$miRNA2))
predicted_PA_genes = subset(predictedkeep, miRNA %in% paste0("hsa-",significant_testPA_miRNAs$miRNA2))
predicted_FI_genes = subset(predictedkeep, miRNA %in% paste0("hsa-",significant_testFI_miRNAs$miRNA2))
predicted_mort_genes = subset(predictedkeep, miRNA %in% paste0("hsa-",significant_testmort_miRNAs$miRNA2))

length(unique(predicted_age_genes$miRNA)) #176
length(unique(predicted_PA_genes$miRNA))  #216
length(unique(predicted_FI_genes$miRNA))  #59
length(unique(predicted_mort_genes$miRNA))#15

# Function for OR with ID conversion from RefSeq to Entrez
perform_OR_with_conversion <- function(significant_mirna_lists, all_mirna_targets, organism = "hsa") {
  # significant_mirna_lists: Named list of data frames with columns "miRNA" and "Gene" (RefSeq IDs)
  # all_mirna_targets: Data frame with "miRNA" and "Gene" (RefSeq IDs) for all miRNA-gene mappings
  
  results <- list()
  
  # Helper function to convert RefSeq to Entrez
  convert_refseq_to_entrez <- function(refseq_ids) {
    mapped_ids <- bitr(refseq_ids, fromType = "REFSEQ", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    mapped_ids$ENTREZID
  }
  
  # Convert all background genes to Entrez IDs
  all_background_genes <- convert_refseq_to_entrez(all_mirna_targets$Gene)
  
  for (outcome in names(significant_mirna_lists)) {
    sig_mirna_genes <- significant_mirna_lists[[outcome]]
    
    # Convert significant genes to Entrez IDs, preserving duplicates
    significant_genes <- convert_refseq_to_entrez(sig_mirna_genes$Gene)
    
    # Perform GO enrichment
    go_enrichment <- enrichGO(
      gene = significant_genes,
      universe = all_background_genes,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )
    go_enrichment@result$Outcome = outcome
    a = dotplot(go_enrichment, showCategory=20, font.size = 15)
    enrichres2 = pairwise_termsim(go_enrichment, method = "JC", semData = NULL, showCategory = 200)
   # b = treeplot(enrichres2)
    b = enrichplot::cnetplot(
      enrichres2,
      #layout = 'kk',
      color.params = list(
        foldChange = NULL, 
        edge = TRUE, 
        category = "#882255", 
        gene = "#332288"
      ), 
      max.overlaps = 1000, 
      cex.params = list(
        category_label = 1, 
        gene_label = 0.5
      )
    )
    # Perform KEGG enrichment
    kegg_enrichment <- enrichKEGG(
      gene = significant_genes,
      universe = all_background_genes,
      organism = organism,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    )
    kegg_enrichment@result$Outcome = outcome
    if(sum(kegg_enrichment@result$p.adjust < 0.05) != 0){
    c = dotplot(kegg_enrichment, showCategory=20, font.size = 15)
    enrichreskegg = pairwise_termsim(kegg_enrichment, method = "JC", semData = NULL, showCategory = 200)
    enrichreskeggx <- setReadable(enrichreskegg, 'org.Hs.eg.db', 'ENTREZID')
    d <- enrichplot::cnetplot(
      enrichreskeggx,
     #layout = 'kk',
      color.params = list(
        foldChange = NULL, 
        edge = TRUE, 
        category = "#882255", 
        gene = "#332288"
      ), 
      max.overlaps = 1000, 
      cex.params = list(
        category_label = 1, 
        gene_label = 0.6
      )
    )
    allplots = ggarrange(a, b, c, d, ncol = 2, nrow = 2, align = "h", widths = c(1, 1.8), labels = c("a", "b", "c", "d"), font.label = list(size = 20, color = "black", face = "bold", family = NULL))
    ggsave(plot = allplots, filename = paste0("../pathwayplots", outcome, ".png"), width = 20, height = 16, units = "in", limitsize = F, bg = "white", dpi = 600)
      }
    else{
  
    allplots = ggarrange(a, b, ncol = 2, nrow = 1, align = "h", widths =  c(1, 1.8), labels = c("a", "b"), font.label = list(size = 20, color = "black", face = "bold", family = NULL))
    ggsave(plot = allplots, filename = paste0("../pathwayplots", outcome, ".png"), width = 20, height = 8, units = "in", limitsize = F, bg = "white", dpi = 600)
    
    }
    
    # Store results for both GO and KEGG
    results[[outcome]] <- list(GO = go_enrichment, KEGG = kegg_enrichment)
  }
  
  return(results)
}


significant_mirna_lists = list("age" = predicted_age_genes, "PA" = predicted_PA_genes, "FI" = predicted_FI_genes, "mort" =predicted_mort_genes)
# Run the enrichment function 
results <- perform_OR_with_conversion(significant_mirna_lists, predictedkeep)
print(results)

#Save results per outcome
age_go = results$age$GO@result
PA_go = results$PA$GO@result
FI_go = results$FI$GO@result
mort_go = results$mort$GO@result
results_go = list(age_go, PA_go, FI_go, mort_go)
results_go_df = Reduce(function(x,y) full_join(x,y), results_go)
results_go_df = dplyr::select(results_go_df, c(Outcome, ID, Description, GeneRatio, BgRatio, FoldEnrichment, p.adjust, Count))

# Make file of GO results
results_go_df = results_go_df %>% mutate(
  Outcome = case_when(Outcome == "age" ~ "Chronological age",
                      Outcome == "PA" ~ "PhenoAge",
                      Outcome == "FI" ~ "Frailty Index",
                      Outcome == "mort" ~ "Ten-Year Mortality"),
  p.adjust = case_when(p.adjust < 0.01 ~ sub("e","x10^",sprintf(p.adjust, fmt="%.2e")),
                   p.adjust>= 0.01 ~ sprintf(p.adjust, fmt = "%.2f"))
) %>% rename(., "GO Pathway Identifier" = "ID",
                "Proportion of input genes mapped to the pathway" = "GeneRatio",
                "Proportion of background genes mapped to the pathway" = "BgRatio",
                "Fold Enrichment" = "FoldEnrichment",
                "p-value after FDR-correction" = "p.adjust") 

openxlsx::write.xlsx(results_go_df, "../../AdditionalFile9.xlsx")

# Make file of KEGG results
age_kegg = results$age$KEGG@result
PA_kegg = results$PA$KEGG@result
FI_kegg = results$FI$KEGG@result
mort_kegg = results$mort$KEGG@result
results_kegg = list(age_kegg, PA_kegg, FI_kegg, mort_kegg)
results_kegg_df = Reduce(function(x,y) full_join(x,y), results_kegg)
results_kegg_df = dplyr::select(results_kegg_df, c(Outcome, ID, Description, subcategory, category, GeneRatio, BgRatio, FoldEnrichment, p.adjust, Count))
results_kegg_df = results_kegg_df %>% mutate(
  Outcome = case_when(Outcome == "age" ~ "Chronological age",
                      Outcome == "PA" ~ "PhenoAge",
                      Outcome == "FI" ~ "Frailty Index",
                      Outcome == "mort" ~ "Ten-Year Mortality"),
  p.adjust = case_when(p.adjust < 0.01 ~ sub("e","x10^",sprintf(p.adjust, fmt="%.2e")),
                       p.adjust>= 0.01 ~ sprintf(p.adjust, fmt = "%.2f"))
) %>% rename(., "KEGG Pathway Identifier" = "ID",
             "Proportion of input genes mapped to the pathway" = "GeneRatio",
             "Proportion of background genes mapped to the pathway" = "BgRatio",
             "Fold Enrichment" = "FoldEnrichment",
             "p-value after FDR-correction" = "p.adjust",
             "Subcategory" = "subcategory",
             "Category" = "category") 
#results_kegg_df$Outcome <- with(results_kegg_df, ifelse(duplicated(Outcome), NA, Outcome))
# Save the table to an Excel file
openxlsx::write.xlsx(results_kegg_df, "../../AdditionalFile10.xlsx")

#Save results per outcome ~ en results
age_en_go = results_en$en_age$GO@result
PA_en_go = results_en$en_PA$GO@result
FI_en_go = results_en$en_FI$GO@result
mort_en_go = results_en$en_mort$GO@result

age_en_kegg = results_en$en_age$KEGG@result
PA_en_kegg = results_en$en_PA$KEGG@result
FI_en_kegg = results_en$en_FI$KEGG@result
mort_en_kegg = results_en$en_mort$KEGG@result

#Tissue
expression = data.table::fread("../../../Data/BiologicalInterpretation/hsa_snc_expression.csv", data.table=F) #Tissue Atlas v2
mirexpression = subset(expression, type == "mirna" & acc %in% paste0("hsa-", all591mirnas$miRNA))
# Calculate highest expression

highest_expressed_organ <- mirexpression %>%
  group_by(acc) %>%
  slice_max(expression, with_ties = FALSE) %>%
  dplyr::select(acc, organ, expression)


# Overexpression
mirexpression_keep_age = subset(highest_expressed_organ, acc %in% paste0("hsa-",significant_testage_miRNAs$miRNA) & expression > 0)
mirexpression_keep_PA = subset(highest_expressed_organ, acc %in% paste0("hsa-",significant_testPA_miRNAs$miRNA) & expression > 0)
mirexpression_keep_FI = subset(highest_expressed_organ, acc %in% paste0("hsa-",significant_testFI_miRNAs$miRNA) & expression > 0)
mirexpression_keep_mort = subset(highest_expressed_organ, acc %in% paste0("hsa-",significant_testmort_miRNAs$miRNA) & expression > 0)

#Tissue-specificity
tsi = data.table::fread("../../../Data/BiologicalInterpretation/sncRNA_tsi_group.csv", data.table=F) #Tissue Atlas v2
tsi = subset(tsi, sncRNA %in% paste0("hsa-",all591mirnas$miRNA)) #575 in tsi set
tsi = tsi %>% mutate(., acc = sncRNA) %>%
  left_join(., highest_expressed_organ) %>%
  mutate(organ = case_when(TSI <= 0.8 ~ "not specific",
                           TSI > 0.8 ~ organ),
         organ = factor(organ))

in_all <- tsi %>%
  filter(acc %in% paste0("hsa-", common_miRNAs))


# Define a function for performing Fisher's exact tests on organ counts
perform_fisher_tests <- function(tsi, significant_miRNAs, dataset_name1, dataset_name2, p_threshold = 0.05) {
  # Create a subset of the dataset for the significant miRNAs
  specific_miRNAs <- subset(tsi, acc %in% paste0("hsa-", significant_miRNAs$miRNA2))
  
  # Create a contingency table of organ counts in each dataset
  organ_counts <- table(
    organ = tsi$organ,
    dataset = ifelse(tsi$acc %in% specific_miRNAs$acc, dataset_name1, dataset_name2)
  )
  
  # Total number of rows in the dataset
  total_miRNAs <- nrow(tsi)
  sig_miRNAs <- nrow(significant_miRNAs)
  not_sig = total_miRNAs - sig_miRNAs
  
  # Initialize lists to store all results and significant results
  pairwise_results <- list()
  significant_results <- list()
  
  # Initialize a data frame to store organ counts and expected counts
  expected_results <- data.frame(organ = rownames(organ_counts), 
                                 count_set1 = numeric(nrow(organ_counts)), 
                                 count_set2 = numeric(nrow(organ_counts)),
                                 expected_count1 = numeric(nrow(organ_counts)), 
                                 expected_count2 = numeric(nrow(organ_counts)),
                                 p_value = numeric(nrow(organ_counts)),
                                 stringsAsFactors = FALSE)
  
  # Loop through each organ level and test separately
  for (org in rownames(organ_counts)) {
    # Calculate expected counts based on the proportion of total miRNAs
    expected_count1 <- sum(organ_counts[org, ]) / total_miRNAs * sig_miRNAs
    expected_count2 <- sum(organ_counts[org, ]) / total_miRNAs * not_sig
    
    # Create a contingency table for the current organ
    organ_table <- matrix(c(
      organ_counts[org, dataset_name2], sum(organ_counts[org, ]) - organ_counts[org, dataset_name2],
      organ_counts[org, dataset_name1], sum(organ_counts[org, ]) - organ_counts[org, dataset_name1]
    ), nrow = 2)
    
    # Check if any counts in the organ_table are zero
    if (any(organ_table == 0)) {
      expected_results[expected_results$organ == org, ] <- c(org, 
                                                             organ_counts[org, dataset_name1], 
                                                             organ_counts[org, dataset_name2], 
                                                             expected_count1, 
                                                             expected_count2, 
                                                             NA)
    }
    
    # Perform Fisher's exact test for this organ
    test_result <- fisher.test(organ_table, workspace = 2e8, conf.int = T, conf.level = 1-p_threshold)  # Increase workspace to 200 million
    pairwise_results[[org]] <- test_result
    
    # Store counts and expected counts in the results data frame
    expected_results[expected_results$organ == org, ] <- c(org, 
                                                           organ_counts[org, dataset_name1], 
                                                           organ_counts[org, dataset_name2], 
                                                           expected_count1, 
                                                           expected_count2, 
                                                           test_result$p.value)
    
    # Check if the p-value is significant and if the count of significant miRNAs is greater than expected
    if (test_result$p.value < p_threshold) {
      significant_results[[org]] <- test_result  # Save significant result
    }
  }
  expected_results$pFDR = p.adjust(expected_results$p_value, method = "fdr")
  expected_results$outcome = dataset_name1
  sig_expected_results = subset(expected_results, pFDR < p_threshold & pFDR !=0)
  
  
  # Return results and expected counts
  list(all_results = pairwise_results, 
       significant_results = significant_results, 
       expected_counts = expected_results,
       sig_expected = sig_expected_results)
}

# Test results
results_testage <- perform_fisher_tests(tsi, significant_testage_miRNAs, "tsi_specific_age", "tsi")
results_testPA <- perform_fisher_tests(tsi, significant_testPA_miRNAs, "tsi_specific_PA", "tsi")
results_testFI <- perform_fisher_tests(tsi, significant_testFI_miRNAs, "tsi_specific_FI", "tsi")
results_testmort <- perform_fisher_tests(tsi, significant_testmort_miRNAs, "tsi_specific_mort", "tsi")

#Combine in one file
tissue_allpheno = rbind(results_testage$expected_counts, results_testPA$expected_counts, results_testFI$expected_counts, results_testmort$expected_counts)
tissue_allpheno$pFDR = p.adjust(tissue_allpheno$p_value, method = "fdr")
# Convert the relevant columns to numeric
tissue_allpheno$count_set1 <- as.numeric(tissue_allpheno$count_set1)
tissue_allpheno$expected_count1 <- as.numeric(tissue_allpheno$expected_count1)
tissue_allpheno$count_set2 <- as.numeric(tissue_allpheno$count_set2)
tissue_allpheno$expected_count2 <- as.numeric(tissue_allpheno$expected_count2)

additionalfile8 = tissue_allpheno %>% rename(., 
                                             "N miRNAs differentially expressed"="count_set1",
                                             "N miRNAs not differentially expressed"="count_set2",
                                             "Expected N miRNAs differentially expressed" = "expected_count1", 
                                             "Expected N miRNAs not differentially expressed" = "expected_count2") %>%
  mutate(., Outcome = case_when(grepl("age", outcome) ~ "Chronological age",
                                grepl("PA", outcome) ~ "PhenoAge",
                                grepl("FI", outcome) ~ "Frailty Index", 
                                grepl("mort", outcome) ~ "Ten-Year Mortality"),
         Organ = stringr::str_to_title(gsub("_", " ", organ)),
         pFDR = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                          pFDR>= 0.01 ~ sprintf(pFDR, fmt = "%.2f"))) %>%
  dplyr::select(Outcome, 
                Organ,
                "N miRNAs differentially expressed", 
                "N miRNAs not differentially expressed", 
                "Expected N miRNAs differentially expressed", 
                "Expected N miRNAs not differentially expressed", 
                pFDR)

# Save the table to an Excel file
openxlsx::write.xlsx(additionalfile8, "../../AdditionalFile8.xlsx")

sig_tissue_allpheno = subset(tissue_allpheno, pFDR < 0.05)
overrep = subset(sig_tissue_allpheno, count_set1 > expected_count1)
underrep = subset(sig_tissue_allpheno, count_set1 < expected_count1)

### Comparison to other datasets:
#Chronological age
Wuage = openxlsx::read.xlsx("../../../WuEstimates.xlsx", sheet = "Age")
names(Wuage)[names(Wuage) == "miRNA"] <- "miRNA2"
Wu_age_not_included = anti_join(Wuage, significant_testage_miRNAs) #All associated
Wu_age_included = inner_join(Wuage, significant_testage_miRNAs)
sum(Wu_age_included$`Coefficient.(Log2)` < 0 & Wu_age_included$log2FoldChange > 0) #0
sum(Wu_age_included$`Coefficient.(Log2)` > 0 & Wu_age_included$log2FoldChange < 0) #0

#PhenoAge
WuPA = openxlsx::read.xlsx("../../../WuEstimates.xlsx", sheet = "PhenoAge")
names(WuPA)[names(WuPA) == "miRNA"] <- "miRNA2"
Wu_PA_not_included = anti_join(WuPA, significant_testPA_miRNAs) #All associated
Wu_PA_included = inner_join(WuPA, significant_testPA_miRNAs)
20 - nrow(Wu_PA_included) #4
head(Wu_PA_not_included$miRNA2)
sum(Wu_PA_included$`Coefficient.(Log2)` < 0 & Wu_PA_included$log2FoldChange > 0) #0
sum(Wu_PA_included$`Coefficient.(Log2)` > 0 & Wu_PA_included$log2FoldChange < 0) #0


## Additional File 5
create_tables5 <- function(outcome, set) {
  # Dynamically reference the data frame
  df <- get(paste0("DE_", outcome, "_", set))
  df = df %>% mutate(miRNA = rownames(df),
                     pFDR = case_when(padj < 0.01 ~ sub("e","x10^",sprintf(padj, fmt="%.2e")),
                                      padj>= 0.01 ~ sprintf(padj, fmt = "%.2f")))
  
  # Set outcome label based on outcome name
  outcome_label <- switch(outcome,
                          "AGE_Standard50" = "Chronological Age",
                          "PhenoAge_Standard50" = "PhenoAge",
                          "FI_Standard50" = "the Frailty Index",
                          "event_Standard50" = "Ten-Year Mortality",
                          "CHECK")
  
  # Select columns based on outcome
  cols <- if (outcome == "AGE_Standard50") c("miRNA", "baseMean", "log2FoldChange", "pFDR", "padj") else c("miRNA", "log2FoldChange", "pFDR", "padj")
  df <- df[, cols, drop = FALSE]
  
  # Define column names based on outcome and set
  names(df) <- c(
    "miRNA",
    if ("baseMean" %in% cols) paste("Average expression level across all samples normalised by sequencing depth in the", set),
    paste("log2 Fold Change estimate for", outcome_label, "in the", set, "set"),
    paste("pFDR for", outcome_label, "in the", set, "set"),
    paste("padj for", outcome_label, "in the", set, "set")
  )
  
  return(df)
}

combine_tables <- function(outcomes, sets) {
  # Initialize an empty list to collect tables
  all_tables <- list()
  
  # Iterate through sets first, then outcomes
  for (set in sets) {
    set_tables <- list()
    
    # Loop through outcomes
    for (outcome in outcomes) {
      # Skip "Validation" for outcomes that don't have it
      if (set == "Validation" && outcome != "AGE_Standard50") next
      
      # Create the table for the current outcome and set
      new_table <- create_tables5(outcome, set)
      
      # Add to set_tables list
      set_tables[[outcome]] <- new_table
    }
    
    # Combine all tables for the current set
    set_combined <- purrr::reduce(set_tables, full_join, by = "miRNA")
    all_tables[[set]] <- set_combined
  }
  
  # Combine all set tables into a single table
  final_table1 <- purrr::reduce(all_tables, full_join, by = "miRNA")
  final_table = inner_join(mibase20_df, final_table1)
  return(final_table)
}

# Define outcomes and sets
outcomes <- c("AGE_Standard50", "PhenoAge_Standard50", "FI_Standard50", "event_Standard50")
sets <- c("Training", "Test", "Validation")

# Load the sequence of mirBase v20
mibase20 = read.fasta("../../../Data/mature.fa")
# Convert to dataframe
mibase20_df <- data.frame(
  Name = names(mibase20),                # Extract sequence names
  Sequence = sapply(mibase20, function(seq) paste(seq, collapse = "")),  # Combine sequence elements into a single string
  Annotation = sapply(mibase20, function(seq) attr(seq, "Annot")),       # Extract annotation
  stringsAsFactors = FALSE              # Prevent factors
) %>%
  mutate(
    ID = gsub("^>(\\S+).*", "\\1", Annotation),           # Extract ID (e.g., "hsa-let-7a-5p")
    Accession = gsub("^>\\S+\\s(\\S+).*", "\\1", Annotation),  # Extract accession (e.g., "MIMAT0000062")
    Species = gsub("^>\\S+\\s\\S+\\s([A-Za-z]+\\s[A-Za-z]+).*", "\\1", Annotation), # Extract species (e.g., "Homo sapiens")
    miRNA = gsub("^>\\S+\\s\\S+\\s[A-Za-z]+\\s[A-Za-z]+\\s(.*)", "\\1", Annotation) # Extract description (e.g., "let-7a-5p")
  ) %>%
  subset(., Species == "Homo sapiens") %>%
  dplyr::select(., c(miRNA, ID, Accession, Sequence))  

# Create Additional File 5 dynamically
tables5 <- combine_tables(outcomes, sets)
top20 = tables5[order(tables5$`padj for Chronological Age in the Test set`), ]
top20_keep = top20[c(1:20),]
tables5 = tables5[order(tables5$`padj for Chronological Age in the Validation set`), ]
top20_keep_val = tables5[c(1:20),]
top20_keep[top20_keep$miRNA %in% top20_keep_val$miRNA,]$miRNA
tables5 = dplyr::select(tables5, -starts_with("padj")) #only included to 
openxlsx::write.xlsx(tables5, "../../TableS5.xlsx")

#Additional File 7
create_tables7 <- function(set_used) {
  # Dynamically reference the data frame
  df <- get(paste0("DE_", set_used, "_Training"))
  df = df %>% mutate(miRNA = rownames(df),
                     pFDR = case_when(padj < 0.01 ~ sub("e","x10^",sprintf(padj, fmt="%.2e")),
                                      padj>= 0.01 ~ sprintf(padj, fmt = "%.2f")))
  
  # Set outcome label based on outcome name
  set_used_label <- switch(set_used,
                          "AGE_Standard50" = "regular Trainin Set, with cut-off 50%",
                          "AGE_Sensitivity25" = "sensitivity analysis with cut-off 25%",
                          "AGE_NoCutoff" = "sensitivity analysis using all 2083 miRNAs",
                          "CHECK")
  
  # Select columns based on outcome
  cols <- c("miRNA", "baseMean", "log2FoldChange", "pFDR", "padj")
  df <- df[, cols, drop = FALSE]
  
  # Define column names based on outcome and set
  names(df) <- c(
    "miRNA",
    if ("baseMean" %in% cols) paste("Average expression level across all samples normalised by sequencing depth in the", set_used_label),
    paste("log2 Fold Change estimate for Chronological Age in the", set_used_label),
    paste("pFDR for Chronological Age in the", set_used_label),
    paste("padj for Chronological Age in the", set_used_label)
  )
  
  return(df)
}

combine_tables2 <- function(outcomes) {
  # Initialize an empty list to collect tables
  all_tables <- list()
  
  # Loop through outcomes
  for (outcome in outcomes) {
    # Create the table for the current outcome
    new_table <- create_tables7(outcome)
    
    # Add the table to the list
    all_tables[[outcome]] <- new_table
  }
  
  # Combine all tables by miRNA
  final_table1 <- purrr::reduce(all_tables, full_join, by = "miRNA")
  final_table = inner_join(mibase20_df, final_table1)
  return(final_table)
}

# Define outcomes
outcomes <- c("AGE_Standard50", "AGE_Sensitivity25", "AGE_NoCutoff")

# Create Additional File 7 dynamically
tables7 <- combine_tables2(outcomes)

# Sort the tables based on the pFDR column for Chronological Age
tables7 <- tables7[order(tables7$`padj for Chronological Age in the regular Trainin Set, with cut-off 50%`), ]

# Exclude 'padj' columns
tables7 <- dplyr::select(tables7, -starts_with("padj"))

# Save the table to an Excel file
openxlsx::write.xlsx(tables7, "../../TableS7.xlsx")

#Answer question Reviewer 3
old_associated = read_excel("../Old Results For Reviewer 3/Additional_File4_Table_S4.xlsx")
old_associated2 = read_excel("../Old Results For Reviewer 3/Additional_File6_Table_S5.xlsx")
old_associated2[, -1] <- lapply(old_associated2[, -1], as.numeric)

old_associated2 = subset(old_associated2, `pFDR Age in train set with LLOQ 50%` < 0.05)
old_associated[, -1] <- lapply(old_associated[, -1], as.numeric)
old_age = subset(old_associated, `pFDR Age in test set` < 0.05)

sum(old_age$miRNA %in% significant_testage_miRNAs$miRNA2) / nrow(old_age)
sum(significant_testage_miRNAs$miRNA2 %in% old_age$miRNA) / nrow(significant_testage_miRNAs)

sig_val_age = subset(DE_AGE_Standard50_Validation, padj < 0.05)
old_age_val = subset(old_associated, `pFDR Age in validation set` < 0.05)
sum(old_age_val$miRNA %in% rownames(sig_val_age)) / nrow(old_age_val)
sum(rownames(sig_val_age) %in% old_age_val$miRNA) / nrow(sig_val_age)

old_PA = subset(old_associated, `pFDR PhenoAge in test set` < 0.05)
sum(old_PA$miRNA %in% significant_testPA_miRNAs$miRNA2) / nrow(old_PA)
sum(significant_testPA_miRNAs$miRNA2 %in% old_PA$miRNA) / nrow(significant_testPA_miRNAs)

old_FI = subset(old_associated, `pFDR Frailty Index in test set` < 0.05)
sum(old_FI$miRNA %in% significant_testFI_miRNAs$miRNA2) / nrow(old_FI)
sum(significant_testFI_miRNAs$miRNA2 %in% old_FI$miRNA) / nrow(significant_testFI_miRNAs)

old_mort = subset(old_associated, `pFDR 10-years-mortality in test set` < 0.05)
sum(old_mort$miRNA %in% significant_testmort_miRNAs$miRNA2) / nrow(old_mort)
sum(significant_testmort_miRNAs$miRNA2 %in% old_mort$miRNA) / nrow(significant_testmort_miRNAs)


