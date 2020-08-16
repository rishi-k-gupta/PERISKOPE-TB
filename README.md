# PERISKOPE-TB

# Overview

- This repository provides code required to reconstruct the PERISKOPE-TB prognostic model for incident tuberculosis in settings with low transmission. 
- Please refer to http://periskope.org and our peer-reviewed publication for full details of the model, including the methods used during development.
- Please [get in touch](mailto:r.gupta@ucl.ac.uk) with us if you would like to collaborate on PERISKOPE-TB. 
- If you use code from this repository, please cite our paper as a condition of use.

# Requirements

The code is run in R and requires the following packages to be installed:

```
install.packages("rstpm2")
install.packages("tidyverse")
install.packages("rms")
install.packages("lubridate")
```

# Files

The repository includes:
- `PERISKOPE-TB-replication.R` - the script to run the model
- `input_data.csv` - simulated dataset including 10 example patients
- `input_data_dictionary.csv` - data dictionary for `input_data.csv`
- `qfn_lookup` - look-up table to normalise quantitative Quantiferon results to a percentile scale
- `tspot_lookup` - look-up table to normalise quantitative T-SPOT.TB results to a percentile scale
- `tst_lookup` - look-up table to normalise quantitative tuberculin skin test results to a percentile scale
- `TB_burden_countries_2020-08-14.csv` - [WHO TB burden estimates](https://www.who.int/tb/country/data/download/en/)

# Input data

The input dataset must include the following variables (as shown in the simulated `input.data.csv`):

| Variable               | Description                                                                    | Type                 | Levels / range                                                                | Required?                                                     |
|------------------------|--------------------------------------------------------------------------------|----------------------|-------------------------------------------------------------------------------|---------------------------------------------------------------|
| age                    | Age (years)                                                                    | Numeric              | 0-80                                                                          | Yes                                                           |
| qfn_tbag_max           | Maximum Quantiferon TB antigen response (higher of TB1 and TB2 tubes; IU/mL) | Numeric              | 0-20                                                                          | One valid quantitative or binary test result must be provided |
| qfn_negative_control   | Quantiferon negative control (IU/mL)                                        | Numeric              | 0-20                                                                          | One valid quantitative or binary test result must be provided |
| qfn_result             | Quantiferon result                                                             | Factor with 3 levels | Positive, Negative, Indeterminate                                             | One valid quantitative or binary test result must be provided |
| tspot_tbag_max         | Maximum T-SPOT.TB antigen spot count (higher of Panel A and Panel B)         | Integer              | 0-400                                                                         | One valid quantitative or binary test result must be provided |
| tspot_negative_control | T.SPOT.TB negative control spot count                                          | Integer              | 0-400                                                                         | One valid quantitative or binary test result must be provided |
| tspot_result           | T-SPOT.TB result                                                               | Factor with 3 levels | Positive, Borderline positive, Borderline negative, Negative,   Indeterminate | One valid quantitative or binary test result must be provided |
| mantoux_result         | Tuberculin skin test induration (mm)                                           | Integer              | 0-100                                                                         | One valid quantitative or binary test result must be provided |
| contact                | Was the person tested through contact tracing?                                 | Binary               | Yes, No                                                                       | Yes                                                           |
| indexcase_proximity    | Proximity of index case                                                        | Factor               | Household, Non-Household                                                      | For contacts                                                  |
| indexcase_sputumsmear  | Sputum smear status of index case                                              | Binary               | Positive, Negative                                                            | For contacts                                                  |
| migrant                | Was the person tested born abroad?                                             | Binary               | Yes, No                                                                       | For non-contacts                                              |
| country_of_birth       | Country of birth                                                               | Factor               | Formatted as per `TB_burden_countries_2020-08-14.csv`                                        | For migrants                                                  |
| date_of_entry          | Approximate date of migration                                                  | Date (yyyy-mm-dd)    | 0-15 years prior to test_date                                                 | For migrants                                                  |
| hivpos                 | Is the person tested known to be living with HIV?                              | Binary               | Positive, Negative                                                            | Yes                                                           |
| transplant             | Has the person tested receieved a solid organ or haematological transplant?  | Binary               | Yes, No                                                                       | Yes                                                           |
| ltbi_treatment         | Preventative treatment commenced?                                              | Binary               | Yes, No                                                                       | Yes                                                           |
| test_date              | Date of latent TB test                                                         | Date (yyyy-mm-dd)    | N/A                                                                           | Yes                                                           |

# Pre-processing

After reading in the input dataset, the code will perform a number of pre-processing steps, as annotated in the script. These steps include:

- Transforming age to restricted cubic spline variables (using 5 knots at fixed positions).
- Transforming latent TB test results to a normalised percentile scale and then to restricted cubic spline variables (using 5 knots at fixed positions). Note: If >1 latent TB test has been done, the QuantiFERON result will enter the model in preference to T-SPOT.TB, while the tuberculin skin test will only be used if no valid QuantiFERON of T-SPOT.TB result is available.
- Generating a composite TB exposure variable (`exposure_cat4b`) with 4 levels (household contact of smear positive index case (`Household, smear+`); other contact (`Other contacts`); migrant from high TB burden country with no contact (`No contact, migrant`); or no exposure (`No contact, non-migrant`)).
- Generating a 'months since migration' (`months_migrant`) variable (months from migration date to latent TB test date) for people in the `No contact, migrant` exposure category. This is coded as 0 people in other exposure categories.

# Predictions

The code will generate 2-year estimates for predicted risk of incident TB. The interval for the predictions can be altered in the code by amending the `studytime` variable (in days).  
Note: predictions will be missing for any participants with missing input values (as shown for participant #2 in the input data, where HIV status was missing). For external validation of the PERISKOPE-TB model, we would recommend consideration of multiple imputation to deal with missing data. Please get in touch](mailto:r.gupta@ucl.ac.uk) with us if you would like further advice. 
