# This script reads in input data and then produces 2-year incident TB risk estimates for each individual in the dataset, using the PERISKOPE-TB model
# Please see README for instructions on software requirements and how to prepare the input data
# Further details on the model can be found at http://periskope.org
# For advice, please contact Dr Rishi Gupta (r.gupta@ucl.ac.uk)

# Load required R packages

library(tidyverse)
library(rstpm2)
library(rms)
library(lubridate)

# Load required objects

load("qfn_lookup") # QFT percentile look-up table
load("tspot_lookup") # T=SPOT percentile look-up table
load("tst_lookup") # TST percentile look-up table
load("fit_final_github_2020-08-14") # This is the model object

# Load WHO TB incidence by country data and preprocess for merging later

country_tb_inc <- read.csv("TB_burden_countries_2020-08-14.csv")
country_tb_inc <- country_tb_inc %>% select(country, year, e_inc_100k) %>% 
  rename(country_of_birth=country, year_of_entry=year)

#Load input data (named as input_data.csv). A dummy CSV of 10 simulated patients is provided for testing

input_data <- read.csv("input_data.csv")

# Preprocess the input_data

## Test results ##

### Create QFN and T-SPOT 'difference' variables
input_data$qfn_tbag_nil <- input_data$qfn_tbag_max - input_data$qfn_negative_control
input_data$tspot_max_diff <- input_data$tspot_tbag_max - input_data$tspot_negative_control

### Replace inderminate IGRA results as missing, since these are not valid for entry into the model
input_data$qfn_tbag_nil[input_data$qfn_result=="Indeterminate"] <- NA
input_data$tspot_max_diff[input_data$tspot_result=="Indeterminate"] <- NA

### Replace values <0 as zero
input_data$qfn_tbag_nil[input_data$qfn_tbag_nil<0]<-0
input_data$tspot_max_diff[input_data$tspot_max_diff<0]<-0
input_data$mantoux_result[input_data$mantoux_result<0]<-0

### Look up normalised percentile for each test result
input_data$pct_qfn <- qfn.lookup$pct_qfn[findInterval(input_data$qfn_tbag_nil, qfn.lookup$qfn.min)]
input_data$pct_tspot <- tspot.lookup$pct_tspot[findInterval(input_data$tspot_max_diff, tspot.lookup$tspot.min)]
input_data$pct_tst <- tst.lookup$pct_tst[findInterval(input_data$mantoux_result, tst.lookup$tst.min)]

### Compostite percentile variable (uses QFT > TSPOT > TST)
### If quantitative IGRA result is missing, this is imputed as the median for LTBI positive and negative groups, respectively
input_data <- input_data %>% mutate(pct_testspl1 = case_when(!is.na(pct_qfn) ~ as.integer(pct_qfn),
                                                         !is.na(pct_tspot) ~ as.integer(pct_tspot),
                                                         qfn_result=="Positive" ~ as.integer(87),
                                                         qfn_result=="Negative" ~ as.integer(1),
                                                         tspot_result=="Positive" ~ as.integer(87),
                                                         tspot_result=="Borderline positive" ~ as.integer(79),
                                                         tspot_result=="Borderline negative" ~ as.integer(76),
                                                         tspot_result=="Negative" ~ as.integer(1),
                                                         !is.na(pct_tst) ~ as.integer(pct_tst)))

### Add test result splines (5 knots at fixed positions)
pct_test_spline5 <- as.data.frame(rcspline.eval(input_data$pct_testspl1, knots = c(5, 27.5, 50, 72.5, 95)))
colnames(pct_test_spline5) <- c("pct_testspl2", "pct_testspl3", "pct_testspl4")
input_data <- data.frame(input_data,pct_test_spline5)

## Age splines (5 knots at fixed positions)
input_data <- input_data %>% rename(agespl1=age)
age_spline5 <- as.data.frame(rcspline.eval(input_data$agespl1, knots = c(8, 25, 33.07, 45, 64)))
colnames(age_spline5) <- c("agespl2", "agespl3", "agespl4")
input_data <- data.frame(input_data,age_spline5)

## Composite 'exposure' variable with 4 levels
## If the person was tested through contact screening, this will either be coded as "Household, smear+" or "Other contacts"
## For non-contacts, it is coded as "No contact, migrant" if the participant was born in a country with annual TB incidence >100/100,000
## For non-contacts who are not born in high TB burden countries, it is coded as "No contact, non-migrant"

### Contact levels
input_data$exposure_cat4b[input_data$indexcase_proximity=="Household" & input_data$indexcase_sputumsmear=="Positive"] <- "Household, smear+"
input_data$exposure_cat4b[input_data$indexcase_proximity=="Non-Household" | input_data$indexcase_sputumsmear=="Negative"] <- "Other contacts"

### Migrant level

#### First ensure entry date is a date, and create new 'year' of entry variable (to merge with TB incidence by country data). 
#### 2018 is the year with most recent data available for TB incidence by country. Therefore, years later than 2018 are replaced as 2018
#### The model is developed & validated for a range of 0-15 years since entry among migrants. Therefore, entry pre-2005 is replaced as 2015
input_data$date_of_entry <- ymd(input_data$date_of_entry)
input_data$year_of_entry <- year(input_data$date_of_entry)
input_data$year_of_entry[input_data$year_of_entry<2005] <- 2005
input_data$year_of_entry[input_data$year_of_entry>2018] <- 2018

#### Get TB incidence in country of birth at year of migration (earliest 2005, latest 2017) and replace variable if appropriate
input_data <- left_join(input_data, country_tb_inc, all.x=T)
input_data$exposure_cat4b[input_data$e_inc_100k>=100 & input_data$contact=="No"] <- "No contact, migrant"

### No exposure level
input_data$exposure_cat4b[(input_data$migrant=="No" | input_data$e_inc_100k<100) & input_data$contact=="No"]<- "No contact, non-migrant"

#### Create months since migration variable for migrant non-contacts (maximum of 180 months = 15 years allowed)
#### This is set to 0 for all participants who are not in the "No contact, migrant" exposure category
input_data$test_date <- ymd(input_data$test_date)
input_data$months_migrant <- as.numeric((input_data$test_date-input_data$date_of_entry)/30)
input_data$months_migrant[input_data$exposure_cat4b!="No contact, migrant"] <- 0
input_data$months_migrant[input_data$months_migrant>12*15] <- 12*15
input_data$months_migrant[input_data$months_migrant<0] <- 0

## Concise input dataframe - removes redundant variables
input_data <- input_data %>% select(agespl1, agespl2, agespl3, agespl4, 
                                    pct_testspl1, pct_testspl2, pct_testspl3, pct_testspl4, 
                                    exposure_cat4b, months_migrant, 
                                    hivpos, transplant, 
                                    ltbi_treatment, qfn_tbag_nil, qfn_result, tspot_max_diff,tspot_result, mantoux_result) %>%
  rename(transplant_assumed=transplant)

# Generate predictions 
# (for 2-year interval = 730-42 days as indicated by 'studytime' variable, since the first 42 days after testing is excluded as prevalent TB interval)
input_data$studytime <- 730-42
predictions <- as.data.frame(100*(predict(fit.final.github, input_data, type="fail", se.fit=T)))
input_data <- bind_cols(input_data, predictions)

## The columns #Estimate, 'lower' and 'upper' indicate 2-year incident TB risk % with lower and upper 95% confidence intervals for each individual
## Predictions will be missing for any participants with missing input values (as shown for participant #2 in the input data, where HIV status was missing)