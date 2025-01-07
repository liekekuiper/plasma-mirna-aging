######## File creation
setwd("/Volumes/Biological_Age/miRNA")
rm(list=ls())

library("dplyr")
library("tidyr")
library("haven")
library("openxlsx")

#miRNAs RS-I and RS-II
rnatr=data.frame(read.csv("Data/rnaallageJulia.csv")) #Getting data split and PhenoAge
rnatr=select(rnatr, c(ergoid, train, PhenoAge)) #Keep calculated PhenoAge and whether someone is in the train set
#rna12 = read.csv("Data/cpm_rs12.csv", row.names = 1) #newly normalized data
#rna12 = inner_join(rnatr, rna12) 
rna12 = rnatr
  
#miRNAs RS-IV
rna4 = read.csv("Data/cpm_rs4.csv", row.names = 1)
rna4 = select(rna4, ergoid)

#Basic characteristics to calculate age
#RS-I and RS-II
demo12=read_spss("Data/RoterdamStudy_Basics2014.sav") #basic characteristics
bld12=read_spss("Data/e4_(4)_BLDAFNAM_(10-jul-2017).sav") #blood date
bld12$blooddate = bld12$e4_2686
bld12 = select(bld12, c(ergoid, blooddate))
demo12=inner_join(demo12, bld12)
demo12$age= as.numeric(demo12$blooddate - demo12$date_of_birth)/365.25
demo12=select(demo12, c(ergoid, rs_cohort, sex, age))
#RS-IV
demo4=read.xlsx("Data/miRNA_755Participants_RS4.xlsx", sheet = "Demographic info 755 poeple")
names(rna4)[names(rna4) == "ErgoID"]<- "ergoid"
names(demo4)=c("ergoid", "age", "sex")

## Cell count for technical adjustments
cell12=data.frame(read_sav("Data/e4_(4)_LAB_(10-mar-2010)b.sav"))
cell12=select(cell12, c(ergoid, e4_12794, e4_12801, e4_12841)) #keep percentage Lymphocytes and monocytes

#Loading ERGO-EXTRA
cell4=data.frame(read_sav("Data/ex_(1)_LAB_(30-Mar-2023).sav"))
cell4=select(cell4, c(ergoid, ex_18850_lab01, ex_18852_lab01, ex_18812_lab01)) #keep percentage Lymphocytes and monocytes & RBC
names(cell12) = names(cell4) = c("ergoid", "lympho", "mono", "RBC")
cell12=cell12[complete.cases(cell12),]
cell4=cell4[complete.cases(cell4),]

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
                     TRUE ~ "No"))


#### RS IV Endpoints ####
adl_ex = read_spss("./Data/RSIV_endpoints/ex_intvw_ADL_(04-mar-2021).sav")
adl_ex = select(adl_ex, c(ergoid, starts_with("ex_EX")))
med_ex = read_spss("./Data/RSIV_endpoints/ex_intvw_MEDHIST_(04-mar-2021).sav")
med_ex = select(med_ex, c(ergoid, starts_with("ex_EX")))
pain_ex= read_spss("./Data/RSIV_endpoints/ex_intvw_PAINPERC_(04-mar-2021).sav")
pain_ex = select(pain_ex, c(ergoid, starts_with("ex_EX")))

#Join what we have
rnatechlist12 = list(demo12, cell12, plate12, rna12)
rnatechlistex = list(demo4 , adl_ex, med_ex, pain_ex, cell4 , plate4 , rna4)
mirna12 = Reduce(function(x,y)inner_join(x,y), rnatechlist12)
mirnaex = Reduce(function(x,y)inner_join(x,y), rnatechlistex)


#Use 60% of dataset as training set and remaining 40% as testing set 
rnatr=data.frame(read.csv("Data/rnaallageJulia.csv")) #Getting data split and PhenoAge
rnatr=select(rnatr, c(ergoid, train))
trainyes = subset(rnatr, train == 1)

trainset <- mirna12 %>% subset(., ergoid %in% trainyes$ergoid)
testset  <- dplyr::anti_join(mirna12, trainset, by = 'ergoid')

mirna12 <- mirna12 %>% mutate(
  train = case_when(ergoid %in% trainset$ergoid ~ 1,
                    ergoid %in% testset$ergoid ~ 0)
)

#info
afspraak = read_spss("./Data/updated clinical endpoints/e4_(4)_RESPONS_(12-mar-2018)_excerpt.sav")
#mortality
mort  = read_spss("./Data/updated clinical endpoints/fp_VitalStatus_2023-47.sav")
mort <- mort %>%
  group_by(ergoid) %>%
  arrange(desc(fp_censordate)) %>%
  distinct(ergoid, .keep_all = TRUE)
mort  = select(mort, c(ergoid, fp_vitalstatus, fp_censordate))
mortex= subset(mort, ergoid %in% mirnaex$ergoid)
mortE4= subset(mort, ergoid %in% mirna12$ergoid)
rm(mort)

#### morbidity ####
#COPD
COPDE4 = read.csv("./Data/updated clinical endpoints/COPDstatus_byXB.csv") 
COPDE4 = inner_join(COPDE4, bld12)
COPDE4 = COPDE4 %>% mutate(
  prev_COPD = case_when(COPDstatus_byXB == 1 & incCOPDdate16_byXB <= blooddate ~ 1,
                        COPDstatus_byXB == 0 | (COPDstatus_byXB == 1 & incCOPDdate16_byXB > blooddate) ~ 0),
  inc_COPD = case_when(prev_COPD == 1 ~ NA_real_,
                       COPDstatus_byXB == 0 ~ 0,
                       COPDstatus_byXB == 1 & incCOPDdate16_byXB > blooddate ~ 1),
  enddate_COPD = case_when(is.na(incCOPDdate16_byXB) ~ as.Date("2023-04-01"),
                           TRUE ~ as.Date(incCOPDdate16_byXB)))
COPDE4 = select(COPDE4, c(ergoid, prev_COPD, inc_COPD, enddate_COPD))

#Dementia
dementiaE4 = read_spss("./Data/updated clinical endpoints/Dementia_RSI_II_III_2021_12_03.sav")
dementiaE4 = inner_join(dementiaE4, bld12)
dementiaE4 = dementiaE4 %>% mutate(
  prev_dementia = case_when(dementia_incident == 1 & censor_date <= blooddate ~ 1,
                  dementia_prevalent == 9 ~ NA_real_,
                  TRUE ~ dementia_prevalent),
  inc_dementia = case_when(prev_dementia == 1 | (dementia_incident == 1 & censor_date <= blooddate) ~ NA_real_,
                                TRUE ~ dementia_incident),
  enddate_dementia = censor_date
    )
dementiaE4 = select(dementiaE4, c(ergoid, prev_dementia, inc_dementia, enddate_dementia))

#High blood pressure
hbpE4 = read_spss("./Data/updated clinical endpoints/HT2018_analysisfile_(15-may-2018).sav")
visitdateE5 = read_spss("./Data/updated clinical endpoints/e5_(5)_RESPONS_(07-oct-2020)_excerpt.sav")
visitdateE5 = select(visitdateE5, c(ergoid, e5_3494))
names(visitdateE5) = c("ergoid", "e5_visit")
visitdateE6 = read_spss("./Data/updated clinical endpoints/e6_(6)_RESPONS_(10-feb-2017)_EXCERPT.sav")
visitdateE6 = select(visitdateE6, c(ergoid, e6_3494))
names(visitdateE6) = c("ergoid", "e6_visit")
hbplist = list(hbpE4, visitdateE5, visitdateE6)
hbpE4 = Reduce(function(x,y)full_join(x,y), hbplist)
hbpE4 = hbpE4 %>% mutate (
  preval_hbp = case_when(e4_HT2018 == 1 ~ 1,
                            TRUE ~ 0),
  incid_hbp = case_when(preval_hbp == 1 ~ NA_real_,
                           e5_HT2018 == 1 | e6_bpldrug == 1 | e6_diastolicBP > 90 | e6_systolicBP > 140 ~ 1,
                           TRUE ~0 ),
  enddate_hbp = case_when(incid_hbp == 1 & e5_HT2018 == 1 ~ e5_visit,
                          TRUE ~ e6_visit)
)
hbpE4 = select(hbpE4, c(ergoid, preval_hbp, incid_hbp, enddate_hbp))

#CHD
MI = read_spss("./Data/updated clinical endpoints/fp_MI.censor2020_(2024-03-08).sav")
AF = read_spss("./Data/updated clinical endpoints/Atrial Fibrillation - truncated 01.01.2020 - JAE van Oortmerssen.sav")
chdincE4 = read_spss("./Data/updated clinical endpoints/Incident Total Coronary heart disease - MJG Leening - updated.sav")
chdprevE4= read_spss("./Data/updated clinical endpoints/Prevalent CHD at RSI-4 and RSII-2.sav")
chdincE4[,"inc_CHD"][chdincE4[,"inc_CHD"] >= 7]<-NA #Make NA
names(chdprevE4)=c("ergoid", "prev_CHD")
chdprevE4[,"prev_CHD"][chdprevE4[,"prev_CHD"] >= 7]<-NA #Make NA
chdincE4 = select(chdincE4, -c(rs_cohort, fp_startdate))
chdlist = list(chdprevE4, chdincE4, bld12)
chdE4 = Reduce(function(x,y) full_join(x,y), chdlist)
chdE4 = chdE4 %>% mutate(
  inc_CHD = case_when(prev_CHD == 1 | (inc_CHD == 1 & enddat_CHD <= blooddate) ~ NA_real_,
                      TRUE ~ inc_CHD)
)
names(chdE4)[names(chdE4) == "enddat_CHD"] <- "enddate_CHD"
rm(chdincE4)
rm(chdprevE4)

#Oncology
oncologyE4 = read_spss("./Data/updated clinical endpoints/ONCOLOGY_prevalenceANDincidence_morbidityANDmortality_29.02.2020_First_event_no_NMNSC_2015FUP.sav")
oncologyE4 = select(oncologyE4, c(ergoid, Cancer, eventdat))
names(oncologyE4) = c("ergoid", "cancer", "enddate_cancer")

#Diabetes
#diabetesE4 = read_spss("./Data/updated clinical endpoints/rs123_diabetes_lifetimerisk_nieuwe_censordate_age.sav")
diabetesE4 = read_spss("./Data/updated clinical endpoints/fp_IFG. DM. set.sav") %>% subset(., ergoid %in% rnatr$ergoid)
diabetesE4 = select(diabetesE4, c(ergoid, Inci_DM_2015_Frank, Incidentdate_DM, Prevalent_DM))

names(diabetesE4) = c("ergoid", "inc_DM", "enddate_DM", "prev_DM")
diabetesE4 = inner_join(diabetesE4, bld12)
diabetesE4 = diabetesE4 %>% mutate(
  prev_DM = case_when(prev_DM == 1 | (inc_DM == 1 & enddate_DM <= blooddate) ~ 1,
                      TRUE ~ prev_DM),
  inc_DM = case_when(prev_DM == 1 ~ NA_real_,
                     inc_DM == 8 ~ NA_real_,
                     TRUE ~ inc_DM),
  enddate_DM = case_when(is.na(enddate_DM) ~ as.Date("2015-01-01"),
                         TRUE ~ enddate_DM)
)

strokeE4 = read_spss("./Data/updated clinical endpoints/Stroke_2020_detailed (14-02-2023).sav")
strokeE4 = strokeE4 %>% subset(., IC == 1) %>%
  select(., c(ergoid, incident_stroke, observed_end, prevalent_stroke))
names(strokeE4) = c("ergoid", "inc_stroke", "enddate_stroke", "prev_stroke")
strokeE4 = inner_join(strokeE4, bld12)
strokeE4 = strokeE4 %>% mutate(
  prev_stroke = case_when(prev_stroke == 1 | (prev_stroke == 1 & enddate_stroke <= blooddate) ~ 1,
                      TRUE ~ prev_stroke),
  inc_stroke = case_when(prev_stroke == 1 ~ NA_real_,
                     TRUE ~ inc_stroke)
)

library(purrr)

# identify data frames that end on "E4"
E4_names <- ls(pattern = "E4$")
outcomes = Reduce(function(x,y)full_join(x,y), mget(E4_names), accumulate=F)
outcomes = subset(outcomes, ergoid %in% c(mirna12$ergoid, mirnaex$ergoid))
outcomes <- outcomes %>% 
  mutate(
    prev_cancer = case_when(cancer == 1 & enddate_cancer <= blooddate ~ 1, TRUE ~ 0),
    inc_cancer = case_when(is.na(cancer) ~ 0,
                           prev_cancer == 1 | (cancer == 1 & enddate_cancer <= blooddate) ~ NA_real_,
                           cancer == 1 & enddate_cancer > blooddate ~ 1),
    eventdate_cancer = enddate_cancer,
    enddate_cancer = case_when(is.na(cancer) ~ as.Date("2020-02-29"),
                               TRUE ~ enddate_cancer),
    eventdate_CHD = case_when(inc_CHD == 0 ~ as.Date(NA_character_),
                              inc_CHD == 1 ~ enddate_CHD),
    eventdate_COPD= case_when(inc_COPD== 0 ~ as.Date(NA_character_),
                              inc_COPD== 1 ~ enddate_COPD),
    eventdate_stroke = case_when(inc_stroke == 0 ~ as.Date(NA_character_),
                                 inc_stroke == 1 ~ enddate_stroke),
    eventdate_dementia = case_when(inc_dementia == 0 ~ as.Date(NA_character_),
                                   inc_dementia == 1 ~ enddate_dementia),
    eventdate_DM  = case_when(inc_DM == 0 ~ as.Date(NA_character_),
                              inc_DM == 1 ~ enddate_DM),
    prev_multimorbid = case_when(rowSums(select(., contains("prev_")), na.rm = TRUE) > 1 ~ 1, TRUE ~ 0),
    inc_multimorbid = case_when(prev_multimorbid == 1 ~ NA_real_,
                                rowSums(select(., contains(c("prev_", "inc_"))), na.rm = TRUE) > 1 ~ 1,
                                TRUE ~ 0),
    #mindate = apply(select(., starts_with("eventdate_")), 1, function(x) if(all(is.na(x))) NA else min(x, na.rm = TRUE)) %>% as.Date(), #HIER GAAT IETS MIS
    mindate = as.Date(apply(across(starts_with("eventdate_")), 1, function(x) {
      x <- na.omit(x)  # Remove NA values before sorting
      if (length(x) > 0) sort(x, na.last = NA)[1] else NA
    })),
    secondsmall = as.Date(apply(across(starts_with("eventdate_")), 1, function(x) {
      x <- na.omit(x)  # Remove NA values before sorting
      if (length(x) > 1) sort(x, na.last = NA)[2] else NA
    })),
    enddate_multimorbid = case_when(
      inc_multimorbid == 1 & rowSums(select(., contains("prev_")), na.rm = TRUE) == 1 ~ mindate,
      inc_multimorbid == 1 & rowSums(select(., contains("inc_")), na.rm = TRUE) > 1 ~ secondsmall,
      inc_multimorbid == 0 ~ as.Date("2015-01-01"),
      is.na(inc_multimorbid) ~ as.Date(NA_real_)
    ),
    prev_morbidity = case_when(rowSums(select(., contains("prev_")), na.rm = TRUE) != 0 ~ 1, TRUE ~ 0),
    inc_morbidity = case_when(prev_morbidity == 1 ~ NA_real_,
                              rowSums(select(., contains("inc_")), na.rm = TRUE) >= 1 ~ 1,
                              TRUE ~ 0),
    enddate_morbidity = case_when(
      inc_morbidity == 1 ~ mindate,
      inc_morbidity == 0 ~ as.Date("2015-01-01"),
      is.na(inc_morbidity) ~ as.Date(NA_real_))
  ) %>%
  mutate(
    enddate_multimorbid = case_when(
      enddate_multimorbid > as.Date("2015-01-01") ~ as.Date("2015-01-01"), #As diabetes is only available until then
      TRUE ~ enddate_multimorbid
    ),
    enddate_morbidity = case_when(
      enddate_morbidity > as.Date("2015-01-01") ~ as.Date("2015-01-01"), #As diabetes is only available until then
      TRUE ~ enddate_morbidity
    ),
    inc_multimorbid = case_when(
      enddate_multimorbid > as.Date("2015-01-01") ~ 0,
      TRUE ~ inc_multimorbid
    ),
    inc_morbidity = case_when(
      enddate_morbidity > as.Date("2015-01-01") ~ 0,
      TRUE ~ inc_morbidity
    )
  )

                              
outcomes = select(outcomes, -c(cancer, starts_with("eventdate_")))

#weirdly now have missingness for age
stadjtrain = read.csv('Data/train_standardadjustments.csv')
stadjtest = read.csv('Data/test_standardadjustments.csv')
stadjex = read.csv('Data/ex_standardadjustments.csv')
stadjtrain = stadjtrain %>% select(., -c(starts_with("miR"), starts_with("let"),event, fp_mortdat, fp_censordate, blooddate))
stadjtest = stadjtest %>% select(., -c(starts_with("miR"), starts_with("let"),event, fp_mortdat, fp_censordate, blooddate))
stadjex = stadjex %>% select(., -c(starts_with("miR"), starts_with("let"),event, fp_mortdat, fp_censordate, blooddate))

#mirna12 = outcomes %>% inner_join(., mirna12) %>% select(., -c(age, sex, PhenoAge)) %>% left_join(full_join(stadjtrain[,c("ergoid", "AGE", "sex", "PhenoAge")], stadjtest[,c("ergoid", "AGE", "sex","PhenoAge")]))%>% rename(., age = AGE)
#trainset = trainset %>% left_join(., mirna12) %>% select(., -c(age, sex, PhenoAge)) %>% left_join(stadjtrain[,c("ergoid", "AGE", "sex", "PhenoAge")]) %>% rename(., age = AGE)
#testset = testset %>% left_join(., mirna12) %>% select(., -c(age, sex, PhenoAge)) %>% left_join(stadjtest[,c("ergoid", "AGE", "sex", "PhenoAge")]) %>% rename(., age = AGE)

trainset = left_join(stadjtrain, outcomes)
testset = left_join(stadjtest, outcomes)
mirna12 = left_join(full_join(stadjtrain, stadjtest), outcomes)


ADL_E4 = haven::read_spss("/Volumes/Nutrition_Lifestyle_Frailty/Raw_data_FI/ERGO-4/e4_intvw_ADL_(07-nov-2011).sav")
ADL_E4 = ADL_E4 %>% mutate(
  IADL_Telephone = case_when(e4_di1_28 == 7 ~ 3, #different annotation, do it before changing the rest
                             e4_di1_28 > 3 & e4_di1_28 < 7 ~ 2,
                             e4_di1_28 <= 3 ~ e4_di1_28 -1),
  across(starts_with("e4_di1_"), ~ ifelse(. > 4, NA, . - 1)),
  ADL_Dressing = pmax(e4_di1_03, e4_di1_04, e4_di1_05, na.rm = TRUE),
  ADL_Rising = pmax(e4_di1_06, e4_di1_07, na.rm = TRUE),
  ADL_Eating = pmax(e4_di1_08, e4_di1_09, na.rm = TRUE),
  ADL_Walking = pmax(e4_di1_10, e4_di1_11, na.rm = TRUE),
  ADL_Hygiene = pmax(e4_di1_12, e4_di1_13, e4_di1_15, na.rm = TRUE),
  ADL_Grip = pmax(e4_di1_14, e4_di1_19, e4_di1_20, na.rm = TRUE),
  ADL_Reach = pmax(e4_di1_17, e4_di1_18, na.rm = TRUE),
  ADL_Activities = pmax(e4_di1_22, e4_di1_23, e4_di1_25, na.rm = TRUE),
  IADL_Travel = e4_di1_24,
  IADL_Mealprep = e4_di1_29, 
  IADL_Medication = e4_di1_32, 
  IADL_Finances = e4_di1_33,
  IADL_Groceries = e4_di1_22,
  IADL_Laundry = e4_di1_31,
  IADL_Housekeeping = e4_di1_26) %>%
  rowwise() %>% 
  mutate(BADL =  sum(c_across(starts_with("ADL_")), na.rm = F), #HAQ ADL Score
         IADL = sum(c_across(starts_with("IADL_")), na.rm = F) ) %>% #IADL score
  ungroup() %>% 
  select(., ergoid, BADL, IADL)

# Frailty index
FI_delta=data.frame(read.csv("Data/FI_delta.csv", row.names = 1)) %>% select(., ergoid, FI_delta)
FI=read.csv("Data/FIE4_imputed.csv") %>% rename(., FI_E4 = FI) %>% select(., c(ergoid, FI_E4)) %>% left_join(., FI_delta)
FI = subset(FI, ergoid %in% rna12$ergoid)

mirna12 = left_join(mirna12, ADL_E4)
mirna12 = left_join(mirna12, FI)
mirna12 = left_join(mirna12, rna12[, c("ergoid", "train")])

write.csv(mirna12, "./Data/miRNA_phenotypes_RSI_II.csv")
write.csv(mirnaex, "./Data/miRNA_phenotypes_RSIV.csv")
write.csv(trainset, "./Data/trainset.csv")
write.csv(testset, "./Data/testset.csv")

mirna12$train = factor(mirna12$train, levels = c(1,0))
#Demographics
ancestry_selfreported = read_spss("./Data/Ancestry/ERGO_ethnic_backgrounddef_(2024-10-09).sav")
ancestryRSI = read_spss("./Data/Ancestry/RS_I_GeneticEtnicity_(11-JAN-2016).sav")

ancestryRSII = read_spss("./Data/Ancestry/RS_II_GeneticEtnicity_(11-JAN-2016).sav")
ancestryRSIV = read_spss("./Data/Ancestry/RS_IV_GeneticEtnicity_(20-MAY-2021).sav")
ancestry_selfreported = ancestry_selfreported %>%
  mutate(self_ances = case_when(ASCCEG %in% c(0, 1, 2, 12, 14) ~ "CEU",
                                ASCCEG %in% c(3, 4) ~ "AFR",
                                ASCCEG %in% c(5, 6, 7) ~ "ASN",
                                ASCCEG %in% c(13) ~ "ADMIXED",
                                TRUE ~ NA_character_))

E4ancestry = list(ancestryRSI, ancestryRSII, ancestry_selfreported[,c("ergoid", "self_ances")])
ancestryRSI_II = Reduce(function(x,y) full_join(x,y), E4ancestry)
ancestryRSI_II = ancestryRSI_II %>% mutate(
  RSI_ancestry = case_when(!RSI_ancestry %in% c("ASN", "CEU", "AFR", "ADMIXED") ~ NA_character_, TRUE ~ RSI_ancestry),
  RSII_ancestry = case_when(!RSII_ancestry %in% c("ASN", "CEU", "AFR", "ADMIXED") ~ NA_character_, TRUE ~ RSII_ancestry),
  ancestry_use = coalesce(RSI_ancestry, RSII_ancestry, self_ances)) 

ancestryRSIV = full_join(ancestryRSIV, ancestry_selfreported[,c("ergoid", "self_ances")])
ancestryRSIV = ancestryRSIV %>% mutate(
  RSIV_ancestry = case_when(!RSIV_ancestry %in% c("ASN", "CEU", "AFR", "ADMIXED") ~ NA_character_, TRUE ~ RSIV_ancestry),
  ancestry_use = coalesce(RSIV_ancestry, self_ances)) 
mirna12 = left_join(mirna12, ancestryRSI_II)

table1::table1(~ factor(sex) + AGE + factor(fp_vitalstatus) + PhenoAge + FI_E4 + BADL + IADL + ancestry_use + as.factor(prev_morbidity) + as.factor(prev_multimorbid) + FI_delta| train, data=mirna12)

valiset = mirnaex
valiset = left_join(valiset, ancestryRSIV)
t.test(mirna12$age ~ mirna12$train)$p.value
t.test(mirna12$PhenoAge ~ mirna12$train)$p.value
t.test(mirna12$FI_E4 ~ mirna12$train)$p.value
t.test(mirna12$FI_delta ~ mirna12$train)$p.value
t.test(mirna12$BADL ~ mirna12$train)$p.value
t.test(mirna12$IADL ~ mirna12$train)$p.value

chisq.test(mirna12$sex, mirna12$train)
chisq.test(mirna12$fp_vitalstatus, mirna12$train)
chisq.test(mirna12$ancestry_use, mirna12$train)
trainset = left_join(trainset, ancestryRSIV[,c("ergoid", "ancestry_use")])

table1::table1(~ AGE + factor(sex) + ancestry_use, data = valiset)
#table1::table1(~ AGE + factor(sex) + BADL + IADL + factor(inc_mortality), data = valiset)
#vali = select(valiset, c(age, BADL, IADL, sex, inc_mortality)) #get valiset from Analyses.R
vali = select(valiset, c(ergoid, age, sex, ancestry_use)) #get valiset from Analyses.R

vali = vali %>% mutate(#fp_vitalstatus = inc_mortality,
              train = 0) %>%
  full_join(.,trainset[,c("ergoid", "age","sex", "ancestry_use", "train")])

vali$train <- as.factor(vali$train)
t.test(vali$age ~ vali$train)$p.value
#t.test(vali$BADL ~ vali$train)$p.value
#t.test(vali$IADL ~ vali$train)$p.value




chisq.test(vali$sex, vali$train)
chisq.test(vali$fp_vitalstatus, vali$train)
chisq.test(vali$ancestry_use, vali$train)
