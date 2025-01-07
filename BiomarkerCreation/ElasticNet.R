rm(list=ls())
setwd("/Volumes/Biological_Age/miRNA")
library("caret")
library("haven")
library("dplyr")
library("glmnet")
library("survival")
library("openxlsx")
library("ggfortify") 
library("paletteer")

#miRNAs RS-I and RS-II
rnatr=data.frame(read.csv("Data/rnaallageJulia.csv")) #Getting data split and PhenoAge
rnatr=select(rnatr, c(ergoid, train, PhenoAge))
rna12 = read.csv("Data/cpm_rs12.csv", row.names = 1) #newly normalized data
rna12 = inner_join(rna12, rnatr)
demo12=read_spss("Data/RoterdamStudy_Basics2014.sav") #basic characteristics
bld12=read_spss("Data/e4_(4)_BLDAFNAM_(10-jul-2017).sav") #blood date
demo12=inner_join(demo12, bld12)
demo12$age= as.numeric(demo12$e4_2686 - demo12$date_of_birth)/365.25
demo12=select(demo12, c(ergoid, rs_cohort, sex, age))

#miRNAs RI-IV
rna4 = read.csv("Data/cpm_rs4.csv", row.names = 1)
#Demographics
demo4=read.xlsx("Data/miRNA_755Participants_RS4.xlsx", sheet = "Demographic info 755 poeple")
names(rna4)[names(rna4) == "ErgoID"]<- "ergoid"
names(demo4)=c("ergoid", "age", "sex")

## Cell count for technical adjustments
cell12=data.frame(read_sav("Data/e4_(4)_LAB_(10-mar-2010)b.sav"))
cell12=select(cell12, c(ergoid, e4_12794, e4_12801, e4_12841)) #keep percentage Lymphocytes and monocytes
cell4=data.frame(read_sav("Data/ex_(1)_LAB_(30-Mar-2023).sav"))
cell4=select(cell4, c(ergoid, ex_18850_lab01, ex_18852_lab01, ex_18812_lab01)) #keep percentage Lymphocytes and monocytes & RBC
names(cell12) = names(cell4) = c("ergoid", "lympho", "mono", "RBC")

#Plate count for technical adjustments
plate12=read.xlsx("Data/PlateNumber_2000SamplesERGO4_miRNAs.xlsx", sheet = 1)
names(plate12)=c("ergoid", "PlateNR", "Well")
plate4=read.xlsx("Data/PlateNumber_754SamplesERGOX_miRNAs.xlsx", sheet = 1)
names(plate4)=c("ergoid", "PlateNR", "Well")

#Create inner/outer-well
plate12 = plate12 %>% mutate(
  border = case_when(grepl("A|G|H", Well) | grepl("1$|11$|12$", Well)~ "Yes",
                     TRUE ~ "No")
)

plate4 = plate4 %>% mutate(
  border = case_when(grepl("A|G|H", Well) | grepl("1$|11$|12$", Well)~ "Yes",
                     TRUE ~ "No"))

e4=list(rna12, demo12, cell12, plate12)
e4rna=Reduce(function(x,y) inner_join(x,y), e4, accumulate = F)

ex=list(rna4, demo4, cell4, plate4)
exrna=Reduce(function(x,y) inner_join(x,y), ex, accumulate = F)
exrna$rs_cohort = 4
exrna$ergoid <- as.character(exrna$ergoid)

#Read in end-points
FI=data.frame(read.csv("Data/FIE4_imputed.csv"), row.names =1) #Frailty index in ERGO-4
FI <- FI %>%
  #mutate(across(starts_with("FI"), ~ . * 100)) %>%
  select(., ergoid, FI)

mort=read_sav("Data/updated clinical endpoints/fp_VitalStatus_2023-47.sav")
#Remove duplicated ergoids
mort_unique <- mort %>%
  group_by(ergoid) %>%
  arrange(desc(fp_censordate)) %>%
  distinct(ergoid, .keep_all = TRUE)

mort_unique=select(mort_unique, -c(rs_cohort, ond_cod, qes_cod))
mort_unique=inner_join(mort_unique, bld12)
mort_unique=mort_unique%>%mutate(
  event=fp_vitalstatus,
  time=fp_censordate,
  blooddate=e4_2686,
  followup=case_when(time > blooddate ~ as.numeric(time - blooddate)/365.25,
                     TRUE ~ NA_real_),
  event=case_when(!is.na(followup) ~ event,
                  TRUE ~ NA_real_))
morttrain = subset(mort_unique, !is.na(event) & !is.na(followup))
morttrain$mortality = Surv(morttrain$followup, morttrain$event)
morttrain=dplyr::select(morttrain, c(ergoid, mortality))

listouts = list(e4rna, FI, morttrain)
e4rna = Reduce(function(x,y) left_join(x,y), listouts, accumulate = F)
names(e4rna) <- gsub("[^[:alnum:]]", "_", names(e4rna)) 
e4rna$ergoid = as.character(e4rna$ergoid)
test_data = subset(e4rna, train == "0")

eid = match("ergoid",names(e4rna))
mrn=591

displaynames = read.csv("displaynames.csv", row.names = 1)

set.seed(123)

# create a new dataframe to store the residuals
data_elastic = function(trainortest, outcome){
  if(trainortest == "train"){
    df = subset(e4rna, train == "1" & !is.na(mono) & !is.na(lympho) & !is.na(e4rna[,outcome]))
  }else if (trainortest == "test"){
    df = subset(e4rna, train == "0" & !is.na(mono) & !is.na(lympho))
  }else if (trainortest == "all"){
    df = subset(e4rna, !is.na(mono) & !is.na(lympho))
  }else if (trainortest == "vali"){
    df = subset(exrna, !is.na(mono) & !is.na(lympho))
    if(!outcome %in% names(df)){
      df[,outcome] = NA
    }
  }
  
  if(outcome == "age"){
  for (i in (eid + 1):(eid + mrn)) {
    lm_model <- lm(formula = paste(names(df)[i], "~ sex + rs_cohort + mono + lympho + RBC + border + PlateNR"), data = df)}}  
 else{
   for (i in (eid + 1):(eid + mrn)) {
   lm_model <- lm(formula = paste(names(df)[i], "~ age + sex + rs_cohort + mono + lympho + RBC + border + PlateNR"), data = df)
 }}
    df[, i] <- rstandard(lm_model)
  outcomenr = match(outcome,names(df))
    residuals_df = df[,c(outcomenr, (eid + 1):(eid + mrn))]
  rownames(residuals_df)=df$ergoid
  return(residuals_df)
}

### PCA plots to check removal of technical variation ####
pca_check = data_elastic("all", "age")
pca_check$ergoid = rownames(pca_check)
pcaplate = plate12
pcaplate$ergoid = as.character(pcaplate$ergoid)
pcaplate$PlateNR = as.factor(pcaplate$PlateNR)
pcaplate$border = as.factor(pcaplate$border)
pca_check = left_join(pca_check, pcaplate)

check.pca <- prcomp(pca_check[,c(2:592)], 
                   center = TRUE, 
                   scale. = TRUE) 

# Generate a colorblind-friendly palette with 26 colors
colorblind_palette <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", 
  "#0072B2", "#D55E00", "#CC79A7", "#999999", 
  "#66CCEE", "#228833", "#CCBB44", "#EE6677", 
  "#AA3377", "#BBBBBB", "#332288", "#88CCEE", 
  "#44AA99", "#117733", "#999933", "#882255", 
  "#661100", "#6699CC", "#AA4466", "#DDCC77", 
  "#117777", "#888888"
)

# Preview the palette
scales::show_col(colorblind_palette)


check.pca.plot <- autoplot(check.pca, 
                           data = pca_check, 
                           colour = 'PlateNR') + 
  theme_classic() + 
  scale_colour_manual(values = colorblind_palette) + 
  labs(colour = "Plate number")

check.pca.plot2 <- autoplot(check.pca, 
                            data = pca_check, 
                            colour = 'border') + 
  theme_classic() + 
  labs(colour = "Placed on\nouter well")
library("ggpubr")
save_pca_plots = ggarrange(check.pca.plot, check.pca.plot2, labels = c("a", "b"))

ggsave(plot = save_pca_plots, filename = "Results/Rebuttal/PCA_plot_technical.png", width = 12, height = 5)


wb <- createWorkbook() #create workbook to save Excel and avoid xlsx package
y=c("age", "PhenoAge", "FI", "mortality")#Outcomes to train on
coefs <- data.frame(matrix(NA, nrow=eid + mrn, ncol=length(y)))
colnames(coefs) <- y
for (outcome in y) {
  
  set.seed(123)
  # Remove rows with missing outcome values
  train.data <- data_elastic("train", outcome)
  test.data  <- data_elastic("test", outcome)
  vali.data  <- data_elastic("vali", outcome)
  x.train = as.matrix(train.data[,2:ncol(train.data)])
  x.test  = as.matrix( test.data[,2:ncol( test.data)])
  x.vali  = as.matrix( vali.data[,2:ncol( vali.data)])
  test.data$ergoid = rownames(test.data)
  vali.data$ergoid = rownames(vali.data)
  
  if (outcome == "mortality"){
    set.seed(123)
    
    
    #Tune for alpha 
    alphalist <- seq(0,1,by=0.1)
    # Fit glmnet model with automatic lambda selection
    cv_alpha <- lapply(alphalist, function(a){
      cv.glmnet(x.train, y=train.data[[outcome]], nfolds = 10, family = "cox",alpha = a)
    })
    # Initialize an empty list to store the minimum values
    min_values_list <- list()
    
    #Store minimum value per alpha
    for (i in 1:11) {
      min_values_list[[i]] <- min(cv_alpha[[i]]$cvm)
    }
    
    # Find the overall minimum value after the loop
    lowest_index <- which.min(unlist(min_values_list))

    #two was lowest
    best_alpha <- alphalist[lowest_index] 
    model3 <- cv_alpha[[lowest_index]]

    # Extract the optimal value of alpha and lambda
    best_lambda <- model3$lambda.min
    print(paste("best alpha:", best_alpha, "best lambda:", best_lambda))
    
    #Find coefficients
    coeff3 <- coef(model3, s = best_lambda)
    # Predict scores on train and test data
    predtrain = predict(model3, newx = x.train, s = best_lambda, type="link")
    train.data[,paste0("miRNA_", outcome)] <- predict(model3, newx = x.train, s = best_lambda, type="link")
    predtest = predict(model3, newx = x.test, s = best_lambda, type="link")
    test.data[,paste0("miRNA_", outcome)] <- predict(model3, newx = x.test, s = best_lambda, type="link")
    predvali = predict(model3, newx = x.vali, s = best_lambda, type="link")
    vali.data[,paste0("miRNA_", outcome)] <- predict(model3, newx = x.vali, s = best_lambda, type="link")
    coefs[1,outcome]=NA #No intercept, thus empty row
    coefs[2:nrow(coefs),outcome] <- as.matrix(coeff3)[,1] 
    coeff3 <- data.frame(rownames(coeff3), as.matrix(coeff3))
    names(coeff3) = c("miRNA", "Weight")
    coeff3 %>%
      left_join(., displaynames) %>%
      select(., c(miRNA2, Weight)) %>%
      dplyr::rename(miRNA = miRNA2)
      

    
    # Calculate the evaluation metrics
    print(Cindex(predtrain, train.data[[outcome]]))
    print(coxnet.deviance(predtrain, train.data[[outcome]]))
    print(Cindex(predtest, test.data[[outcome]]))
    test2 = subset(test.data, !is.na(test.data[[outcome]]))
    print(coxnet.deviance(predtest, test2[[outcome]]))
    
  }
  else{
    set.seed(123)
    model3 <- train(x.train, train.data[[outcome]], method="glmnet", trControl= trainControl("cv", number=10), tuneLength = 10)#Fit elastic net
    
    #Best tuning parameter 
    print(model3$bestTune)
    train.data[,paste0("miRNA_", outcome)] <- predict(model3$finalModel, model3$bestTune$lambda, newx = as.matrix(x.train))
     test.data[,paste0("miRNA_", outcome)] <- predict(model3$finalModel, model3$bestTune$lambda, newx = as.matrix(x.test))
     vali.data[,paste0("miRNA_", outcome)] <- predict(model3$finalModel, model3$bestTune$lambda, newx = as.matrix(x.vali))
    
     #Coefficient of the final model. You need # to specify the best lambda 
    coeff3 <-coef(model3$finalModel, model3$bestTune$lambda)
    coeff3 <- data.frame(rownames(coeff3), as.matrix(coeff3))
    names(coeff3) = c("miRNA", "Weight")
    coeff3 %>%
      left_join(., displaynames) %>%
      select(., c(miRNA2, Weight)) %>%
      dplyr::rename(miRNA = miRNA2)
    
    coefs[,outcome] <- coeff3$Weight
    rownames(coefs) <- rownames(coeff3)
    
    # Calculate the evaluation metrics
    msetrain <- mean((train.data[[outcome]] - train.data[[paste0("miRNA_", outcome)]])^2, na.rm = TRUE)
    msetest  <- mean((test.data[[outcome]]  -  test.data[[paste0("miRNA_", outcome)]])^2, na.rm = TRUE)
    rmsetr <- sqrt(msetrain)
    rmsete <- sqrt(msetest)
    rsqtr <- cor(train.data[[outcome]], train.data[[paste0("miRNA_", outcome)]], use="na.or.complete")^2
    rsqte <- cor(test.data[[outcome]], test.data[[paste0("miRNA_", outcome)]], use="na.or.complete")^2
    
    # Print the evaluation metrics
    cat(paste0("MSE train: ", msetrain, " MSE test: ", msetest, "\n"))
    cat(paste0("RMSE train: ", rmsetr, " RMSE test: ", rmsete, "\n"))
    cat(paste0("R-squared train: ", rsqtr, " R-squared test: ", rsqte, "\n"))
    
    
  }
  save = test.data[, c("ergoid", paste0("miRNA_", outcome))]
  save2= vali.data[, c("ergoid", paste0("miRNA_", outcome))]
  test_data = left_join(test_data, save)
  exrna = left_join(exrna, save2)
  # Create a new worksheet with the outcome name
  addWorksheet(wb, sheetName = outcome)
  writeData(wb, outcome, as.matrix(coeff3))

  if (outcome == y[length(y)]) {
    # Create a new worksheet for all coefficients
    addWorksheet(wb, sheetName = "All")
    coefsave = data.frame(rownames(coefs), coefs)
    names(coefsave)[names(coefsave) == "rownames.coefs."] <- "miRNA"
    coefsave %>%
      left_join(., displaynames) %>%
      select(., c(miRNA2, all_of(outcome))) %>%
      dplyr::rename(miRNA = miRNA2)
    
    writeData(wb, "All", as.matrix(coefsave))
   }
}
saveWorkbook(wb, "Results/Results_miRNA_29Nov2024_EN.xlsx")

sum(coefs$age != 0) -1
sum(coefs$PhenoAge != 0) -1
sum(coefs$FI != 0) -1
sum(coefs$mortality != 0, na.rm=T)

write.csv(test_data, "./Data/e4_mirna_test.csv")
write.csv(exrna, "./Data/ex_mirna_vali.csv")


#### Correlation ####
correlatie=select(test_data, c(age, starts_with("miRNA_")))
correlatie = correlatie[complete.cases(correlatie),]
correlatiex=select(exrna, c(age, starts_with("miRNA_")))
correlatiex = correlatiex[complete.cases(correlatiex),]
miRNAclocks=c("miRNA_age", "miRNA_PhenoAge", "miRNA_FI", "miRNA_mortality")
for(c in miRNAclocks){
  fitage = lm(formula=paste(c,"~age"), data=correlatie)
  correlatie[,paste0("AgeAccel_", c)]=residuals(fitage)
  fitagex = lm(formula=paste(c,"~age"), data=correlatiex)
  correlatiex[,paste0("AgeAccel_", c)]=residuals(fitagex)
}

onlyaa_correlatie = select(correlatie, starts_with("AgeAccel_"))
onlyaa_correlatiex= select(correlatiex, starts_with("AgeAccel_"))

notaa_correlatie = select(correlatie, !starts_with("AgeAccel_"))
notaa_correlatiex= select(correlatiex,!starts_with("AgeAccel_"))

clrs=c(
  "black",
  "#361ae5",
  "#dc267f",
  "#fe6100",
  "#ffb000",
  "#361ae5",
  "#dc267f",
  "#fe6100",
  "#ffb000")
names(correlatie) = c("Age", "mAge", "mPA", "mFI", "mMort", "mAge", "mPA", "mFI", "mMort")

names(onlyaa_correlatie) = names(onlyaa_correlatiex) = c("mAge", "mPA", "mFI", "mMort")
names(notaa_correlatie) = names(notaa_correlatiex) = c("Age", "mAge", "mPA", "mFI", "mMort")

clrs=c(
  "black",
  "#361ae5",
  "#dc267f",
  "#fe6100",
  "#ffb000")

panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  c <<- c +1 
  rect(breaks[-nB], 0, breaks[-1], y, col = clrs[c] , ...) #DNAm and metabolomics in different colours
}


library("RColorBrewer")
cols = brewer.pal(11, "RdYlBu")   # goes from red to white to blue
pal = colorRampPalette(cols)
cor_colors = data.frame(correlation = seq(-1,1,0.01), 
                        correlation_color = pal(201)[1:201])  # assigns a color for each r correlation value
cor_colors$correlation_color = as.character(cor_colors$correlation_color)


panel.cor <- function(x, y, digits = 2, prefixr = "r=", prefixp= "p=",cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  u <- par('usr') 
  names(u) <- c("xleft", "xright", "ybottom", "ytop")
  r <- cor.test(x, y, method = "spearman",exact=F)
  bgcolor = cor_colors[2+(-r$estimate+1)*100,2]    # converts correlation into a specific color
  do.call(rect, c(col = bgcolor, as.list(u))) # colors the correlation box
  txt <- sprintf(r$estimate, fmt = paste0("%.",digits,"f"))
  pval<- formatC(r$p.value, format = "e", digits = 1)
  txtr <- paste0(prefixr, txt)
  txtr <- substitute(paste(italic(prefixr), txt))
  txtp <- substitute(paste(italic(prefixp), pval))
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 12)
  #text(0.5, 0.7, txtr, cex = 2.5)
  #text(0.5, 0.2, txtp, cex = 2.5)
  #abline(h = 0.5, lty = 2) # draws a line between correlatoin coefficient and p value
  
}

png(filename="./Results/CorrBiomarkersTest.png", res=400, width = 17.1, height = 10, units = 'in') 
c <- 0
pairs(notaa_correlatie,
      na.action = na.omit,
      upper.panel = panel.cor,        
      diag.panel = panel.hist,# Adding the histograms
      text.panel = function(x, y, labels, cex, font, ...) 
        text(x, y, labels, cex = 5.5, font=c(1,1,1,1,1,3,3,3,3)[c]),
      oma=c(0.7,0.7,0.7,0.7),
      xaxt="n", yaxt="n")   
dev.off()

png(filename="./Results/CorrBiomarkersValidation.png", res=400, width = 17.1, height = 10, units = 'in') 
c <- 0
pairs(notaa_correlatiex,
      na.action = na.omit,
      upper.panel = panel.cor,        
      diag.panel = panel.hist,# Adding the histograms
      text.panel = function(x, y, labels, cex, font, ...) 
        text(x, y, labels, cex = 5.5, font=c(1,1,1,1,1,3,3,3,3)[c]),
      oma=c(0.7,0.7,0.7,0.7),
      xaxt="n", yaxt="n")   
dev.off()
