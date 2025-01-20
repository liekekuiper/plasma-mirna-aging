##### Define helper functions ####

#To preprocess data and create age accelerated biomarkers
preprocess_data <- function(data, clock_var, covariates, coxyn, outcome) {
  
  if(coxyn == TRUE){
    endvar = gsub("inc_", "enddate_", outcome)  
    data$followup = as.numeric(as.Date(data[,endvar]) - as.Date(data$blooddate))/365.25
    # Subset data
    sub_data <- subset(data, followup > 0 )
    # Add study time variable
    sub_data$Studytime <- sub_data$age + sub_data$followup
  }else if(coxyn == FALSE){
    sub_data = subset(data, !is.na(data[,outcome]))
  }
  
  var_names <- strsplit(covariates, "\\+")[[1]]  # split by + sign
  var_names <- trimws(var_names)  # remove leading/trailing white space
  
  sub_data <- sub_data[complete.cases(sub_data[, var_names]), ] #we want participants to have information on all variables of interest
  
  # Compute residuals for outcome variable
  fit_age <- lm(paste0(clock_var, " ~ age"), data = sub_data)
  residuals <- scale(residuals(fit_age), center = TRUE, scale = TRUE)#Z-scoring
  # Add residuals to data
  sub_data$biological_age <- residuals
  
  # Convert Visit variable to factor
  sub_data$rs_cohort <- as.factor(sub_data$rs_cohort)
  return(sub_data)
}

# Run Cox-PH model
fit_cox_model <- function(datafile, clock_var, outcome_var, covariates) {
  if(length(levels(datafile$rs_cohort)) >1 ){
    covariates2 = paste(covariates, "+ rs_cohort")
  }
  else{
    covariates2 = covariates
  }
  
  formula <- as.formula(paste("Surv(Studytime,", outcome_var,") ~ ", clock_var, "+",paste(covariates2, collapse = " + ")))
  fit <- coxph(formula, data = datafile)
  return(fit)}

# Run linear regression model
fit_lm_model <- function(datafile, outcome, clock_var, covariates) {
  if(length(levels(datafile$rs_cohort)) >1 ){
    covariates2 = paste(covariates, "+ rs_cohort")
  }
  else{covariates2 = covariates
  }
  
  datafile[,"oucome_use"] <- scale(datafile[[outcome]])
  formula <- as.formula(paste("oucome_use", "~",clock_var,"+", paste(covariates2, collapse = " + ")))
  fit <- lm(formula, data = datafile)
  return(fit)
}

# Multinomial model
fit_multinom_model <- function(datafile, multioutcome, clock_var, covariates) {
  if(length(levels(datafile$rs_cohort)) >1 ){
    covariates2 = paste(covariates, "+ rs_cohort")
  }
  else{covariates2 = covariates
  }
  formula <- as.formula(paste(multioutcome, "~",clock_var,"+", paste(covariates2, collapse = " + ")))
  fit <- multinom(formula, data = datafile)
  
  return(fit)}


outcome_labeller <- function(variable,value){
  return(outcome_names[value])
}
