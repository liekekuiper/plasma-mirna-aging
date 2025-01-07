#Creating data for DESEQ2 analyses

###################################################################################################
####################
# Folder locations
####################
rm(list=ls())
setwd("/Volumes/Biological_Age/miRNA")
library("dplyr")
library("openxlsx")
library("haven")

#Load data
#miRNAs RS-I and RS-II
rna12 = read.csv("Data/cpm_rs12.csv", row.names = 1)

rnatr=data.frame(read.csv("Data/rnaallageJulia.csv"))
rnatr=select(rnatr, c(ergoid, train, PhenoAge))

rna12=inner_join(rna12, rnatr)
demo12=read_spss("Data/RoterdamStudy_Basics2014.sav")
bld12=read_spss("Data/e4_(4)_BLDAFNAM_(10-jul-2017).sav")
demo12=inner_join(demo12, bld12)
demo12$AGE= as.numeric(demo12$e4_2686 - demo12$date_of_birth)/365.25
demo12=select(demo12, c(ergoid, rs_cohort, sex, AGE))

#miRNAs RS-IV
#Only select well-expressed miRNAs
rna4 = read.csv("Data/cpm_rs4.csv", row.names = 1)

#Demographics
demo4=read.xlsx("Data/miRNA_755Participants_RS4.xlsx", sheet = "Demographic info 755 poeple")
names(rna4)[names(rna4) == "ErgoID"]<- "ergoid"
names(demo4)=c("ergoid", "AGE", "sex")

## Cell count for technical adjustments
cell12=data.frame(read_sav("Data/e4_(4)_LAB_(10-mar-2010)b.sav"))
cell12=select(cell12, c(ergoid, e4_12794, e4_12801, e4_12841)) #keep percentage Lymphocytes and monocytes
cell4=data.frame(read_sav("Data/ex_(1)_LAB_(30-Mar-2023).sav"))
bl4 = read_spss("Data/ex_(1)_BLDAFNAM_(02-Nov-2023).sav")
bl4 = bl4 %>%
  mutate(blooddate = ex_2686_bld01) %>%
  select(ergoid, blooddate)

cell4=select(cell4, c(ergoid, ex_18850_lab01, ex_18852_lab01, ex_18812_lab01)) #keep percentage Lymphocytes and monocytes & RBC
names(cell12) = names(cell4) = c("ergoid", "lympho", "mono", "RBC")

#Plate count for technical adjustments
plate12=read.xlsx("Data/PlateNumber_2000SamplesERGO4_miRNAs.xlsx", sheet = 1)
names(plate12)=c("ergoid", "PlateNR", "Well")
plate4=read.xlsx("Data/PlateNumber_754SamplesERGOX_miRNAs.xlsx", sheet = 1)
names(plate4)=c("ergoid", "PlateNR", "Well")

plate12 = plate12 %>% mutate(
  border = case_when(grepl("A|G|H", Well) | grepl("1$|11$|12$", Well)~ "Yes",
                     TRUE ~ "No")
)

plate4 = plate4 %>% mutate(
  border = case_when(grepl("A|G|H", Well) | grepl("1$|11$|12$", Well)~ "Yes",
                     TRUE ~ "No")
)

e4=list(rna12, demo12, cell12, plate12)
e4rna=Reduce(function(x,y) inner_join(x,y), e4, accumulate = F)

ex=list(rna4, demo4, cell4, plate4)
exrna=Reduce(function(x,y) inner_join(x,y), ex, accumulate = F)
exrna$rs_cohort = 4
exrna$PhenoAge = NA

train <- subset(e4rna, train == 1)
test  <- subset(e4rna, train == 0)
train = select(train, -train)
test = select(test, -train)

trainclean = select(train, ergoid:miR99b5p)
standardadjustments = c("rs_cohort", "sex", "lympho", "mono", "RBC", "PlateNR", "Well", "border", "AGE")


sens = read.csv("Data/cpm_rs12_sens.csv", row.names = 1)
senstrain = subset(sens, ergoid %in% train$ergoid)
senstest  = subset(sens, ergoid %in% test$ergoid)
senstrain = left_join(senstrain, train[,c("ergoid", all_of(standardadjustments))])
senstest  = left_join(senstest, test[,c("ergoid",all_of(standardadjustments))])

sens25 = read.csv("Data/cpm_rs12_sens25.csv", row.names = 1)
sens25train = subset(sens25, ergoid %in% train$ergoid)
sens25test  = subset(sens25, ergoid %in% test$ergoid)
sens25train = left_join(sens25train, train[, c("ergoid", all_of(standardadjustments))])
sens25test  = left_join(sens25test, test[, c("ergoid", all_of(standardadjustments))])

#### Other outcomes
FI=data.frame(read.csv("Data/FIE4_imputed.csv"), row.names =1) #Frailty index in ERGO-4
FI <- FI %>%
  #mutate(across(starts_with("FI"), ~ . * 100)) %>%
  select(., ergoid, FI)

mort_all =read_sav("Data/updated clinical endpoints/fp_VitalStatus_2023-47.sav")
#Remove duplicated ergoids
mort_unique <- mort_all %>%
  group_by(ergoid) %>%
  arrange(desc(fp_censordate)) %>%
  distinct(ergoid, .keep_all = TRUE)

mort4 = mort_unique %>% subset(., ergoid %in% exrna$ergoid) %>%
  left_join(., bl4) %>%
  mutate(event = fp_vitalstatus) %>%
  select(ergoid, event, blooddate, fp_mortdat, fp_censordate)

mort <- mort_unique %>%
  # Join with selected columns from bld12
  inner_join(bld12[, c("ergoid", "e4_2686")], by = "ergoid") %>%
  
  # Add new variables
  mutate(
    blooddate = e4_2686, 
    event = case_when(
      fp_vitalstatus == 0 ~ 0,  # No event if vital status is 0
      fp_vitalstatus == 1 & fp_mortdat <= blooddate + (10 * 365.25) ~ 1,  # Event occurs within 10 years
      fp_vitalstatus == 1 & fp_mortdat > blooddate + (10 * 365.25) ~ 1    # Event does not occur within 10 years
    )
  ) %>%
  # Select relevant columns
  select(ergoid, event, fp_mortdat, fp_censordate, blooddate) %>%
  # Filter rows where censor date is after or equal to blooddate
  filter(fp_censordate >= blooddate)


listtrain = list(train, FI, mort)
listtest = list(test, FI, mort)
listvali = list(exrna, mort4)
train = Reduce(function(x,y) left_join(x,y), listtrain, accumulate = F)
test = Reduce(function(x,y) left_join(x,y), listtest, accumulate = F)
vali = Reduce(function(x,y) left_join(x,y), listvali, accumulate = F)

data.table::fwrite(train, "Data/train_standardadjustments.csv")
data.table::fwrite(test, "Data/test_standardadjustments.csv")
data.table::fwrite(vali, "Data/ex_standardadjustments.csv")
