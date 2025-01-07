############ Creating FI in ERGO-4 and ERGO-5 with same variables
rm(list=ls())

library("haven")
library("dplyr")
setwd("/Volumes/Biological_Age/miRNA/")
##Frailty index ERGO-5
FIE5        <- read_sav("/Volumes/Biological_Age/DATA_RAW/FrailtyIndex_ERGO5(eline).sav")
#Mobility is missing but available in the data
mobilityE5 <- read_sav("/Volumes/Nutrition_Lifestyle_Frailty/Raw_data_FI/ERGO-5/e5_intvw_JOINT_(29-apr-2016).sav")
mobilityE5 = mobilityE5 %>%
  mutate(FI_Mobility = case_when(e5_EI3_50 == 1 ~ 0,
                                 e5_EI3_50 == 7 ~ 1,
                                 e5_EI3_50 > 1 & e5_EI3_50 < 7 ~ 0.5)) %>%
  select(., c(ergoid, FI_Mobility))
#Add mobility
FIE5 = inner_join(FIE5, mobilityE5)

attr(FIE5$FI_FOLLOWUP, "label") #"FRAILTY INDEX FOLLOW-UP ERGO5 based on 38 variables (IMPUTATIE BASED)"
FIE5 = select(FIE5, -c(FI_IADL_MEDICATIE, FI_FOLLOWUP, TOTAL_KNOWN_FU, KNOWN_20_FU, KNOWN_30_FU)) #Medication not measured in ERGO4
FIE5 = subset(FIE5,  rowSums(!is.na(select(FIE5, starts_with("FI_"))))>=20)
ncolE5= length(grep(x = colnames(FIE5), pattern = "^FI_"))

FIE5_2 <- FIE5 %>%
  mutate(across(starts_with("FI_"), ~as.numeric(as.character(.)), .names = "{col}_num")) %>%
  mutate(
    FIcount = rowSums(select(., ends_with("_num"))),
    FI = FIcount / ncolE5
  )

FIE5_2 = select(FIE5_2, -ends_with("_num"))

write.csv(FIE5_2, file="./Data/FIE5_revised.csv")

##Frailty index ERGO-4
FIE4=read.csv("./Data/FI_E4_not_imputed.csv")
FIE4 = select(FIE4, -Def_CRP) #CRP Not measured in ERGO5
FIE4=subset(FIE4,  rowSums(!is.na(select(FIE4, starts_with("Def_"))))>=20)

require(mice)
require(lattice)
set.seed(123)

md.pattern(head(FIE4))
FIE4[, !names(FIE4) %in% c("ergoid", "age")] = lapply(FIE4[,!names(FIE4) %in% c("ergoid", "age")],  as.factor);
##### Imputation ####
init = mice(FIE4, maxit=0) #halt imputation
meth = init$method
predM = init$predictorMatrix
predM[, c("ergoid")]=0 #Do not take these into account for imputation
#meth[c("fp_csmt15_mortdat")]=""  #Do not impute mortdat
meth[c("Def_Grip")]="polyreg" #Def_Grip was seen as non-missing while missings were present

imputed = mice(FIE4, method=meth, predictorMatrix=predM, m=10, maxit=10)
FI_imputed <- complete(imputed)

sapply(FI_imputed, function(x) sum(is.na(x))) #Looking at missing values

ncolFI= length(grep(x = colnames(FI_imputed), pattern = "^Def_"))

library(dplyr)

FI_imputed2 <- FI_imputed %>%
  mutate(across(starts_with("Def_"), ~as.numeric(as.character(.)), .names = "{col}_num")) %>%
  mutate(
    FIcount = rowSums(select(., ends_with("_num"))),
    FI = FIcount / ncolFI
  )

FI_imputed2 = select(FI_imputed2, -ends_with("_num"))

write.csv(FI_imputed2, file="./Data/FIE4_imputed.csv")

selectscoreFIE5 = select(FIE5_2, c(ergoid, FI))
names(selectscoreFIE5) = c("ergoid", "FI_E5")
selectscoreFIE4 = select(FI_imputed2, c(ergoid, FI))
names(selectscoreFIE4) = c("ergoid", "FI_E4")
# Inner join the two data frames on 'ergoid'
delta_data = inner_join(selectscoreFIE5, selectscoreFIE4)

# Calculate the difference between FI (originating from selectscoreFIE5) and FI (originating from selectscoreFIE4)
delta_data = delta_data %>%
  mutate(FI_delta = FI_E5 - FI_E4) %>%
  select(., ergoid, FI_E4, FI_E5, FI_delta)

write.csv(delta_data, "./Data/FI_delta.csv")

