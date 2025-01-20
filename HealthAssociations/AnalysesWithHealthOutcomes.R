rm(list=ls())
setwd("/Volumes/Biological_Age/miRNA")
library("haven")
library("dplyr")
library("ggplot2")
library("survival")
library("nnet")
library("tidyr")
library("ggh4x")
library("ggpubr")

source("AssociationOutcomeFunctions.R")

setwd("/Volumes/Biological_Age/miRNA")

#### DATA ####
#miRNA biomarkers
test4 = read.csv("./Data/e4_mirna_test.csv", row.names=1) 
young = read.csv("./Data/ex_mirna_vali.csv", row.names=1)
test4 = select(test4, c(ergoid, age, sex, rs_cohort, lympho, mono, RBC, PlateNR, border, starts_with("miRNA_")))
young = select(young, c(ergoid, age, sex, rs_cohort, lympho, mono, RBC, PlateNR, border, starts_with("miRNA_")))

#Multimorbidity
multim4 = read.csv("./Data/miRNA_phenotypes_RSI_II.csv", row.names=1)

#Mortality
mort=read_sav("Data/updated clinical endpoints/fp_VitalStatus_2023-47.sav")
mort <- mort %>%
  group_by(ergoid) %>%
  arrange(desc(fp_censordate)) %>%
  distinct(ergoid, .keep_all = TRUE)
bld12=read_spss("Data/e4_(4)_BLDAFNAM_(10-jul-2017).sav") #blood date
bld12=select(bld12, c(ergoid, e4_2686))
mort4=inner_join(mort, bld12)
mortx=subset(mort, rs_cohort == "RS-IV")
mort4=mort4%>%mutate(
  blooddate=e4_2686,
  enddate_mortality=fp_censordate,
  inc_mortality=case_when(!is.na(enddate_mortality) ~ fp_vitalstatus,
                          TRUE ~ NA_real_)) %>%
  select(., ergoid, enddate_mortality, inc_mortality)

mortx=mortx%>%mutate(
  event=fp_vitalstatus,
  enddate_mortality=fp_censordate,
  blooddate=fp_startdate, #blooddate not available
  inc_mortality=case_when(!is.na(enddate_mortality) ~ fp_vitalstatus,
                          TRUE ~ NA_real_))%>%
  select(., ergoid, blooddate, enddate_mortality, inc_mortality)

#ERGO-EXTRA Questions
ergox = read.csv("./Data/miRNA_phenotypes_RSIV.csv", row.names = 1)
ergox = ergox %>% dplyr::mutate(
  selfhealth = case_when(ex_EXTI_37_GEZLFT > 3 ~ NA_character_,
                         ex_EXTI_37_GEZLFT == 1 ~ "better",
                         ex_EXTI_37_GEZLFT == 2 ~ "same",
                         ex_EXTI_37_GEZLFT == 3 ~ "worse"),
  selfhealth = as.factor(selfhealth),
  hospital = case_when(ex_EXTI_04_HOSP == 0 ~ 0,
                       ex_EXTI_04_HOSP == 1 ~ 1),
  IADL_Telephone = case_when(ex_EXTI_02_28 == 7 ~ 3, #different annotation, do it before changing the rest
                             ex_EXTI_02_28 > 3 & ex_EXTI_02_28 < 7 ~ 2,
                             ex_EXTI_02_28 <= 3 ~ ex_EXTI_02_28 -1),
  across(starts_with("ex_EXTI_02_"), ~ ifelse(. > 4, NA, . - 1)),
  ADL_Dressing = pmax(ex_EXTI_02_03, ex_EXTI_02_04, ex_EXTI_02_05, na.rm = TRUE),
  ADL_Rising = pmax(ex_EXTI_02_06, ex_EXTI_02_07, na.rm = TRUE),
  ADL_Eating = pmax(ex_EXTI_02_08, ex_EXTI_02_09, na.rm = TRUE),
  ADL_Walking = pmax(ex_EXTI_02_10, ex_EXTI_02_11, na.rm = TRUE),
  ADL_Hygiene = pmax(ex_EXTI_02_12, ex_EXTI_02_13, ex_EXTI_02_15, na.rm = TRUE),
  ADL_Grip = pmax(ex_EXTI_02_14, ex_EXTI_02_19, ex_EXTI_02_20, na.rm = TRUE),
  ADL_Reach = pmax(ex_EXTI_02_17, ex_EXTI_02_18, na.rm = TRUE),
  ADL_Activities = pmax(ex_EXTI_02_22, ex_EXTI_02_23, ex_EXTI_02_25, na.rm = TRUE),
  IADL_Travel = ex_EXTI_02_24,
  IADL_Mealprep = ex_EXTI_02_29, 
  IADL_Medication = ex_EXTI_02_32, 
  IADL_Finances = ex_EXTI_02_33,
  IADL_Groceries = ex_EXTI_02_22,
  IADL_Laundry = ex_EXTI_02_31,
  IADL_Housekeeping = ex_EXTI_02_26) %>%
  rowwise() %>% 
  dplyr::mutate(BADL =  sum(c_across(starts_with("ADL_")), na.rm = F), #HAQ ADL Score
                IADL = sum(c_across(starts_with("IADL_")), na.rm = F) ) %>% #IADL score
  ungroup() %>% 
  select(., ergoid, selfhealth, hospital, BADL, IADL)

#Make feeling as healthy as peers reference
ergox$selfhealth <- relevel(ergox$selfhealth, ref = "same")

#Make one big dataframe
pheno4 = list(test4, mort4, multim4)
phenox = list(young, mortx, ergox)

testset = Reduce(function(x,y)left_join(x,y), pheno4, accumulate =F)
valiset = Reduce(function(x,y)left_join(x,y), phenox, accumulate =F)

biomarkers <- c("miRNA_age", "miRNA_PhenoAge", "miRNA_FI", "miRNA_mortality")
linear_outcomes <- c("FI_delta", "FI_E4", "BADL", "IADL", "PhenoAge")
multi_outcomes <- c("selfhealth")
cox_outcomes <- c("inc_mortality", "inc_multimorbid", "inc_morbidity")

models <- c("age + sex + lympho + mono + RBC + border + PlateNR", 
            "sex + lympho + mono + RBC + border + PlateNR", 
            "FI_E4 + age + sex + lympho + mono + RBC + border + PlateNR")
testvali = c("test", "vali")
# Create a data frame of selected combinations using expand.grid
params <- expand.grid(
  clocks = biomarkers,
  outcomes = cox_outcomes,
  modeluse = c(models[2],paste(models[2], "+ prev_cancer + prev_CHD"), paste(models[2], "+ PhenoAge"), paste(models[2], "+ FI_E4")),
  dataset = testvali
)
params = subset(params, !(dataset == "vali" & outcomes != "inc_mortality") & !(dataset == "vali" & grepl("prev_|PhenoAge|FI_E4", modeluse)) & !(grepl("prev_", modeluse) & (!outcomes %in% c("inc_mortality"))))

data = function(input){
  if(input == "test"){
    df = testset
  }else if(input == "vali"){
    df = valiset
  }
  return(df)
}

pb <- txtProgressBar(min = 0, max = nrow(params), style = 3)

# Loop through all combinations of frailty, i, models, and dets
for (j in 1:nrow(params)){
  # Update de progress bar every iteration
  setTxtProgressBar(pb, j)
  
  df = as.character(params[j, "dataset"]) #save datafile
  model = as.character(params[j,"modeluse"]) #save models
  clock = as.character(params[j,"clocks"])
  outs = as.character(params[j,"outcomes"])
  
  datafile = preprocess_data(data(df), clock, model, TRUE, outs) #create datafile
  fitcox = fit_cox_model(datafile, "biological_age", outs, model) #perform Cox Proportional Hazard model
  print(fitcox)
  #Saving output in Params next to input variables
  params[j,"N"]            =fitcox$n
  params[j,"Nevent"]       =fitcox$nevent
  params[j,c("BETA","HR","SE","Z","P")]=summary(fitcox)$coefficients[1,]
  params[j,"LL"]  = summary(fitcox)$conf.int[1,3]
  params[j,"UL"]  = summary(fitcox)$conf.int[1,4]
  params[j,"Concordance"]  = summary(fitcox)$concordance[1]
  # params[j,"interaction"] = summary(fitcox)$coefficients["biological_age:sex",5]
  
  close(pb)
}

mort <- mort_keep   <- dplyr::select(params  , -Z)
mort <- mort_keep <- mort %>% mutate(P_handig = case_when(P < 0.01 ~ sub("e","x10^",sprintf(P, fmt="%.2e")),
                                                          P>= 0.01 ~ sprintf(P, fmt = "%.2f")),
                                     #pFDR = p.adjust(P, method = "fdr"),
                                     HR_CI = paste0(sprintf(HR, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), ";", sprintf(UL, fmt = "%.2f"), ")"),
                                     
                                     C = case_when(Concordance < 0.01 ~ sub("e","x10^",sprintf(Concordance, fmt="%.2e")),
                                                   Concordance>= 0.01 ~ sprintf(Concordance, fmt = "%.2f")),
                                     #interactionFDR = p.adjust(interaction, method = "fdr")
)
mort <- subset(mort, outcomes %in% c("inc_mortality", "inc_multimorbid", "inc_morbidity") & !grepl("PhenoAge|FI", modeluse))

lmparams <- expand.grid(
  clocks = biomarkers,
  outcomes = linear_outcomes,
  modeluse = models,
  dataset = testvali
)

lmparams = subset(lmparams, 
                  (dataset == "test" & outcomes %in% c("FI_delta", "FI_E5") & modeluse == models[3]) | 
                    (dataset == "test" & outcomes%in% c("FI_E4", "PhenoAge") & modeluse == models[1]) |
                    (modeluse == models[1] & !grepl("FI|PhenoAge", outcomes)))


pb <- txtProgressBar(min = 0, max = nrow(lmparams), style = 3)

for (j in 1:nrow(lmparams)){
  # Update de progress bar every iteration
  setTxtProgressBar(pb, j)
  df = as.character(lmparams[j, "dataset"]) #save datafile
  fr= as.character(lmparams[j,"outcomes"]) #save frailty
  model = as.character(lmparams[j,"modeluse"]) #save models
  clock = as.character(lmparams[j,"clocks"])
  
  
  datafile = preprocess_data(data(df), clock, model, FALSE, fr) #create datafile
  
  fitlm   = fit_lm_model(datafile, fr, "biological_age", model) #perform linear regresion
  conf = confint(fitlm)
  
  #Saving output in Params next to input variables
  lmparams[j,"N"]            =length(fitlm$fitted.values)
  lmparams[j,c("BETA","SE")] =summary(fitlm)$coefficients["biological_age",c("Estimate","Std. Error")]
  lmparams[j,"P"]            =summary(fitlm)$coefficients["biological_age","Pr(>|t|)"]
  lmparams[j,"LL"]           =conf[2,1]
  lmparams[j,"UL"]           =conf[2,2]
  lmparams[j,"R2"]           =summary(fitlm)$adj.r.squared
  #lmparams[j,"interaction"]  =summary(fitlm)$coefficients["biological_age:sex", "Pr(>|t|)"]
  close(pb)
  
}

#### Multi-nomial ####

# Create a data frame of all combinations of datasets, models, and dets
multinom <- expand.grid("vali", multi_outcomes, models[1], biomarkers)
names(multinom) = c("dataset", "multiouts","model", "clocks")
pValue_extract <- function(x){
  z <- summary(x)$coefficients/summary(x)$standard.errors
  # 2-tailed Wald z tests to test significance of coefficients
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  p
}

# Loop through all combinations of binary outcomes, i, models, and dets
for (j in 1:nrow(multinom)){
  # Update de progress bar every iteration
  setTxtProgressBar(pb, j)
  
  df = as.character(multinom[j, "dataset"]) #save datafile
  model = as.character(multinom[j,"model"]) #save models
  clock = as.character(multinom[j,"clocks"])
  multi = as.character(multinom[j,"multiouts"])
  
  
  datafile = preprocess_data(data(df), clock, model, FALSE, bin) #create datafile
  fitmnom = fit_multinom_model(datafile, multi,"biological_age", model) #perform GLM
  # print(fitmnom)
  
  bci = as.data.frame(exp(confint(fitmnom)))
  
  a =  levels(data(df)[,multi])[2]
  b =  levels(data(df)[,multi])[3]
  if(length(levels(data(df)[,multi])) > 3){
    c =  levels(data(df)[,multi])[4]
  }
  pValue_extract(fitmnom)
  #Saving output in Params next to input variables
  multinom[j,"Totaal N"]             =nrow(fitmnom$fitted.values)
  multinom[j,"Level L1"]       =a
  multinom[j,"N L1"]     =sum(data(df)[,multi] ==a, na.rm =T)
  multinom[j,"OR L1"]           =exp(summary(fitmnom)$coefficients[a,"biological_age"])
  multinom[j, "LL L1"]= bci["biological_age", paste0("2.5 %.",a)]
  multinom[j, "UL L1"]= bci["biological_age", paste0("97.5 %.",a)]
  multinom[j, "P L1"]= pValue_extract(fitmnom)[a,"biological_age"]
  multinom[j,"Level L2"]       =b
  multinom[j,"N L2"]     =sum(data(df)[,multi] ==b, na.rm =T)
  multinom[j,"OR L2"]           =exp(summary(fitmnom)$coefficients[b,"biological_age"])
  multinom[j, "LL L2"]= bci["biological_age", paste0("2.5 %.",b)]
  multinom[j, "UL L2"]= bci["biological_age", paste0("97.5 %.",b)]
  multinom[j, "P L2"]= pValue_extract(fitmnom)[b,"biological_age"]
  if(length(levels(data(df)[,multi])) > 3){
    multinom[j,"Level 3"]       =c
    multinom[j,"N Level 3"]     =sum(data(df)[,multi] ==c, na.rm =T)
    multinom[j,"OR L3"]           =exp(summary(fitmnom)$coefficients[c,"biological_age"])
    multinom[j, "LL L3"]= bci["biological_age", paste0("2.5 %.",c)]
    multinom[j, "UL L3"]= bci["biological_age", paste0("97.5 %.",c)]
    multinom[j, "P L3"]= pValue_extract(fitmnom)[c,"biological_age"]
  }
  close(pb)
}


# Reshape the dataframe from wide to long format
multinom_long <- pivot_longer(
  multinom,
  cols = starts_with("Level") | starts_with("OR ") | starts_with("N") |starts_with("LL ") | starts_with("UL ") | starts_with("P "),
  names_to = c(".value", "group"),
  names_pattern = "(Level|OR|LL|UL|N|P)_?(.*)"
)

multinom_long = multinom_long %>% mutate(
  OR_handig = paste0(sprintf(OR, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), ";", sprintf(UL, fmt = "%.2f"), ")")
)

outcome_names <- list(
  'FI_E4'="Frailty index",
  'FI_delta'="Delta frailty",
  'PhenoAge'="PhenoAge",
  'IADL'="IADL",
  'BADL'="BADL"
)

# Create the ggplot2 graph
ggplot(lmparams, aes(x = BETA, y = clocks, color = as.factor(clocks))) +
  geom_point() +
  geom_errorbarh(aes(xmin = LL, xmax = UL)) +
  
  facet_grid(outcomes ~ dataset, scales = "fixed", space = "free_x", switch = "y") +  # Updated facet_grid call
  labs(x = "BETA Value", y = NULL) +
  scale_color_discrete(name = "Clocks") +
  theme_minimal() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(0.8, "lines"),
        legend.position = "bottom") +
  guides(color = guide_legend(nrow = 1))

Tol_muted <- c('#88CCEE', '#44AA99', '#117733', '#332288', '#DDCC77', '#999933','#CC6677', '#882255', '#AA4499', '#DDDDDD')

# Create custom order and labels for outcomes
lmparams = lmparams %>% mutate(
  B_CI = paste0(sprintf(BETA, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), ";", sprintf(UL, fmt = "%.2f"), ")"),
  outcome_name_use = case_when(outcomes == "FI_E4" ~ "Frailty index",
                               outcomes == "FI_delta" ~ "Delta frailty",
                               outcomes == "PhenoAge" ~ "PhenoAge",
                               outcomes == "BADL" ~ "BADL",
                               outcomes == "IADL" ~ "IADL"),
  datasetname = case_when(dataset == "test" ~ "Test set (n=772)",
                          dataset == "vali" ~ "Younger validation set (n=754)"),
  Tolmutedouts = case_when(outcomes == "FI_E4" ~ Tol_muted[1],
                           outcomes == "FI_delta" ~ Tol_muted[2],
                           outcomes == "PhenoAge" ~ Tol_muted[3],
                           outcomes == "BADL" ~ Tol_muted[4],
                           outcomes == "IADL" ~ Tol_muted[5]))
lmparams$outcome_name_use = factor(lmparams$outcome_name_use, levels = c("Frailty index",
                                                                         "Delta frailty",
                                                                         "PhenoAge",
                                                                         "BADL",
                                                                         "IADL"))
lmparams$clocks = factor(lmparams$clocks, levels = c("miRNA_mortality",
                                                     "miRNA_FI",
                                                     "miRNA_PhenoAge",
                                                     "miRNA_age"))
# Color blind friendly colours
colorblind=c(
  "#361ae5",
  "#dc267f",
  "#fe6100",
  "#ffb000")

Tol_light2 <- c( '#B57DC2', '#AAAA00', '#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF', '#44BB99', '#DDDDDD', '#BFE39A', '#A6B5E3', "#D4A4CC")


extracolorblind= c(
  "#f46d43",
  "#fdae61",
  "#fee090",
  "#ffffbf",
  "#e0f3f8",
  "#abd9e9",
  "#74add1",
  "#4575b4")

# Make plots
strip_make <- strip_themed(background_x = elem_list_rect(fill = Tol_light2[c(8,6)], linetype = "solid",
                                                         color = "black", linewidth = c(1,1)), 
                           background_y = elem_list_rect(fill = Tol_light2[c(1,3,4,5,7)], linetype = "solid",
                                                         color = "black", linewidth = c(1,1,1,1)))
# Create the ggplot2 graph
cross = ggplot(lmparams, aes(x = BETA, y = clocks, xmin = LL, xmax = UL, color = clocks)) +
  geom_point(size = 3) +
  geom_errorbarh(height=.1) +
  facet_grid2(outcome_name_use ~ datasetname, scales = "fixed", space = "free_x", switch = "y", strip = strip_make) +  # Updated facet_grid call with labeller
  labs(x = "Adjusted Beta", y = NULL) +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_minimal() +
  scale_color_manual(values = colorblind, "miRNA biomarker", labels = rev(c("miRNA Age", "miRNA PhenoAge", "miRNA FI", "miRNA Mortality")),na.translate = F) +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        legend.position = "bottom"
  ) +
  guides(colour = guide_legend(order = 1, nrow = 2, byrow=T, reverse = T))+
  theme(text = element_text(size = 18))
ggsave(cross, filename = "./Results/Cross_sectional_Biomarkers.png", width = 15, height =  12, units = "in", limitsize = F)


# Create custom order and labels for outcomes
mort = mort %>% mutate(
  outcome_name_use = case_when(outcomes == "inc_mortality" ~ "Mortality",
                               outcomes == "inc_multimorbid" ~ "Multimorbidity",
                               outcomes == "inc_morbidity" ~ "First morbidity"), 
  datasetname = case_when(dataset == "test" ~ "Test set (n=772)",
                          dataset == "vali" ~ "Younger validation set (n=754)"))

mort$outcome_name_use = factor(mort$outcome_name_use, levels = c("Mortality",
                                                                 "Multimorbidity",
                                                                 "First morbidity"))
mort$clocks = factor(mort$clocks, levels = c("miRNA_mortality",
                                             "miRNA_FI",
                                             "miRNA_PhenoAge",
                                             "miRNA_age"))

strip_make_cox <- strip_themed(background_x = elem_list_rect(fill = Tol_light2[c(8,6)], linetype = "solid",
                                                             color = "black", linewidth = c(1,1)), 
                               background_y = elem_list_rect(fill = Tol_light2[c(11, 2, 12)], linetype = "solid",
                                                             color = "black", linewidth = c(1,1,1,1)))
mort2 = subset(mort, (outcomes == "inc_mortality" | outcomes ==  "inc_morbidity"  | outcomes ==  "inc_multimorbid") & modeluse == models[2])

cox = ggplot(mort2, aes(x = HR, y = clocks, xmin = LL, xmax = UL, color = clocks)) +
  geom_point(size = 3) +
  geom_errorbarh( height=.1) +
  facet_grid2(outcome_name_use ~ datasetname, scales = "fixed", space = "free_x", switch = "y", strip = strip_make_cox) +  # Updated facet_grid call with labeller
  labs(x = "Adjusted Hazard Ratio", y = NULL) +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal() +
  scale_color_manual(values = colorblind, "miRNA biomarker", labels = rev(c("miRNA Age", "miRNA PhenoAge", "miRNA FI", "miRNA Mortality")),na.translate = F) +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        legend.position = "bottom"
  ) +
  guides(colour = guide_legend(order = 1, nrow = 2, byrow=T, reverse = T))+
  theme(text = element_text(size = 18))
ggsave(cox, filename = "./Results/Cox_sectional_Biomarkers.png", width = 14, height =  15, units = "in", limitsize = F)
ggsave(cox, filename = "./Results/Cox_sectional_Biomarkers_small.png", width = 14, height =  4, units = "in", limitsize = F)

strip_make_selfhealth <- strip_themed(background_x = elem_list_rect(fill = Tol_light2[c(8,6)], linetype = "solid",
                                                                    color = "black", linewidth = c(1,1)), 
                                      background_y = elem_list_rect(fill = Tol_light2[c(10,9)], linetype = "solid",
                                                                    color = "black", linewidth = c(1,1,1,1)))

# Extract unique combinations of multiouts and clocks
unique_combinations <- multinom_long %>% 
  distinct(multiouts, clocks) 


# Create a dataframe with the same number of rows as unique combinations
new_rows <- data.frame(
  dataset = rep("test", nrow(unique_combinations)),
  multiouts = unique_combinations$multiouts,
  model = NA,
  clocks = unique_combinations$clocks,
  `Totaal N` = NA,
  group = NA,
  Level =  rep(c("better", "worse"), nrow(unique_combinations)/2),
  OR = NA,
  N = NA,
  LL = NA,
  UL = NA,
  OR_handig = NA
)

# Combine the original dataframe and the new rows
multinom_long2 <- bind_rows(multinom_long, new_rows)

multinom_long2 = multinom_long2 %>%
  mutate(OR_CI = paste0(sprintf(OR, fmt = "%.2f"), " (", sprintf(LL, fmt = "%.2f"), ";", sprintf(UL, fmt = "%.2f"), ")"),
         datasetname = case_when(dataset == "test" ~ "Test set (n=772)",
                                 dataset == "vali" ~ "Younger validation set (n=754)"),
         outsname = case_when(Level == "better" ~ "Better",
                              Level == "worse" ~ "Worse"))

selfhealthplot = ggplot(multinom_long2, aes(x = OR, y = clocks, xmin = LL, xmax = UL, color = clocks)) +
  geom_point(size = 3) +
  geom_errorbarh( height=.1) +
  facet_grid2(outsname ~ datasetname, scales = "fixed", space = "free_x", switch = "y", strip = strip_make_selfhealth) +  # Updated facet_grid call with labeller
  labs(x = "Adjusted Odds Ratio On Self-Rated Health Compared to Peers", y = NULL) +
  geom_vline(xintercept=1, color='black', linetype='dashed', alpha=.5) +
  theme_minimal() +
  scale_color_manual(values = colorblind, "miRNA biomarker", labels = rev(c("miRNA Age", "miRNA PhenoAge", "miRNA FI", "miRNA Mortality")),na.translate = F) +
  theme(axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank(),  #remove y axis ticks
        legend.position = "bottom"
  ) +
  guides(colour = guide_legend(order = 1, nrow = 2, byrow=T, reverse = T))+
  theme(text = element_text(size = 18))


fig5 = ggarrange(cross, selfhealthplot, cox, nrow = 3, common.legend = T, legend = "bottom", heights = c(4,2,3), labels = c("a", "b", "c"), font.label = list(size = 18, face = "bold"))+ bgcolor("white")
ggsave(fig5, filename = "./Results/All_outs_Biomarkers2.png", width = 13, height =  17, units = "in", limitsize = F)


##### To answer questions from reviewer ####
#Direction betas miRNA-based aging biomarkers
biomarkers = readxl::read_excel("./Results/Results_miRNA_29Nov2024_EN.xlsx", sheet = "All") %>% subset(., !grepl("Intercept", miRNA))
# Convert all but miRNA to numeric
biomarkers[, -1] <- lapply(biomarkers[, -1], as.numeric)

# Replace 0 with NA
biomarkers[, -1] <- lapply(biomarkers[, -1], function(x) ifelse(x == 0, NA, x))

# Filter rows with conflicting directions
differentdirect <- biomarkers[apply(biomarkers[, -1], 1, function(row) {
  any(row > 0, na.rm = TRUE) && any(row < 0, na.rm = TRUE)
}), ]

#Selected miRNAs
count <- sum(rowSums(!is.na(biomarkers[, -1])) > 0)
print(count)
biomarkers$selected_in <- rowSums(!is.na(biomarkers[, -1]))
table(biomarkers$selected_in)

#Prevalent disease
mort_prev = subset(mort_keep, outcomes == "inc_mortality" & dataset == "test" & grepl("prev_", modeluse))

#Circularity of biomarkers
mort_circ = subset(mort_keep, outcomes == "inc_mortality" & dataset == "test" & !grepl("prev_", modeluse))
mort_circu = mort_circ %>% mutate(
  outcomes = case_when(grepl("PhenoAge", modeluse) ~ paste0(outcomes, "_PA"),
                       grepl("FI", modeluse) ~ paste0(outcomes, "_FI"))
) %>% subset(., !outcomes == "inc_mortality")

# Save results and calculate FDR-corrected p-value
multip = select(multinom_long, c(clocks, Level, P))
mortp  = select(mort2, c(clocks, outcomes, P))
lmparp = select(lmparams, c(clocks, outcomes, P))
mort_circup = select(mort_circu, c(clocks, outcomes, P))
names(multip) = names(mortp)
alloutsp = list(multip, mortp, lmparp, mort_circup)
pfdrcalc = Reduce(function(x,y)full_join(x,y), alloutsp)
pfdrcalc = pfdrcalc %>%
  mutate(pFDR = p.adjust(pfdrcalc$P, method = "fdr"),
         pFDR_handig = case_when(pFDR < 0.01 ~ sub("e","x10^",sprintf(pFDR, fmt="%.2e")),
                                 pFDR>= 0.01 ~ sprintf(pFDR, fmt = "%.2f")))
mort2 = left_join(mort2, pfdrcalc)
mort_circu = left_join(mort_circu, pfdrcalc)
mort2 = full_join(mort2, mort_circu)
lmparams = left_join(lmparams, pfdrcalc)
multinom_long = left_join(multinom_long, pfdrcalc, by = c("Level" = "outcomes", "clocks", "P"))

write.csv(mort2, "./Results/Incident.csv", row.names = F)
write.csv(lmparams, "./Results/LinearTestVali.csv", row.names = F)
write.csv(multinom_long, "./Results/SelfratedHealthVali.csv", row.names = F)

#Back to question from reviewer
circularity_table = mort2 %>% subset(., grepl("mortality", outcomes) & dataset != "vali") %>% 
  mutate(Nevent_N = paste0(Nevent,"/", N)) %>%
  select(., c(clocks, outcomes, Nevent_N, HR_CI, pFDR_handig)) 

#Make wide table
circularity_table_wide <- circularity_table %>%
  pivot_wider(
    names_from = outcomes,
    values_from = c(Nevent_N, HR_CI, pFDR_handig),
    names_sep = "_"
  )

# Define the order of metrics and suffixes
metrics <- c("Nevent_N", "HR_CI", "pFDR_handig")
suffixes <- c("inc_mortality", "inc_mortality_PA", "inc_mortality_FI")

# Generate a vector of all column names by adding suffixes to the metrics
ordered_columns <- c("clocks")  # Start with 'clocks'
for (suffix in suffixes) {
  # Add all columns for each metric with the current suffix
  ordered_columns <- c(ordered_columns, paste0(metrics, "_", suffix))
}

# Now apply the new column order to the dataframe
circularity_table_wide <- circularity_table_wide[, ordered_columns]

#Make aging biomarker names better:
circularity_table_wide = circularity_table_wide %>% mutate(
  clocks = case_when(clocks == "miRNA_age" ~ "mirAge",
                     clocks == "miRNA_PhenoAge" ~ "mirPA",
                     clocks == "miRNA_FI" ~ "mirFI",
                     clocks == "miRNA_mortality" ~ "mirMort")
)
#openxlsx::write.xlsx(circularity_table_wide, "Results/CircularityMortality.xlsx")
