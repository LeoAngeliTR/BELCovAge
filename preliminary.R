#___________________________________________________________________________
# This file uses parts of the Social Contact Rates (SOCRATES) modelling project: please visit https://socialcontactdata.org/ 
# 
# Purpose: To perform preliminary operations for the longitudinal sensitivity analysis of the COVID-19 impact in Belgium.
#
# Copyright 2020, SIMID, UNIVERSITY OF ANTWERP & HASSELT UNIVERSITY
#___________________________________________________________________________

# Required Libraries
# Loading necessary R libraries for data manipulation, plotting, and statistical analysis
library(rstudioapi)
library(matrixcalc)
library(popdemo)
library(demography)
library(smoothAPC)
library(plot.matrix)
library(plot3D)
library(blockmatrix)
library(lubridate)
library(readxl)
library(geometry)
library(dplyr)
library(ggplot2)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(reshape2)
library(ie2misc)
library(ggpubr)
library(expm)
library(jsonlite)
library(xtable)
library(BSDA)
library(readr)
library(scales)
library(data.table)
library(tidyverse)
library(qs)

# Preliminary Operations
# Suppress package startup messages and source main project scripts for analysis
suppressPackageStartupMessages(source('./R/socrates_main.R'))  # Load main Socrates functions
source('./add_functions.R')  # Load additional custom functions used in the analysis

# Survey Configuration
# Set parameters for the type of contact data to retrieve from the survey

country <- "Belgium 2020 CoMix (Coletti 2020)"  # Specify the country and survey year

# Define the type of contacts to be included in the analysis
daytype = "All contacts"   # Options: "Saturday-Sunday", "Monday-Friday", "All contacts"
touch = "Physical contacts"  # Type of contact
duration = "All"   # Duration of contact
gender = "All"  # Gender consideration in contacts
cnt_location = opt_location  # Location of contacts, use 'opt_location' to specify particular locations
bool_reciprocal = TRUE  # Whether the data should generate a symmetrical matrix
bool_suppl_professional_cnt = FALSE  # Include imputed supplementary working contacts
bool_hhmatrix_selection = FALSE  # Consider only household member contacts

# Susceptibles Data Loading
# Load simulation data of susceptible individuals, including both naive susceptibles and those with waning immunity from vaccine and infection
sm_susceptible <- readRDS('./data/sm_susc_last_calibration_2024.rds')
# sm_susceptible structure explanation:
#  - sm_susceptible[,n,k] -> 100 simulated observations of susceptibles in age group k on day n
#  - sm_susceptible[n,,k] -> Days of observation from 1/03/2020
#  - sm_susceptible[n,k,] -> nth observation on day k across all age-groups
# Note: [n,k,1] is 0 for all n,k


# Load Survey Data
# Retrieve and store survey data based on the specified parameters
survey_object_CoMix <- get_survey_object(country = "Belgium 2020 CoMix (Coletti 2020)",
                                         daytype = daytype,
                                         touch = touch,
                                         duration = duration,
                                         gender = gender,
                                         cnt_location = cnt_location,
                                         bool_reciprocal = bool_reciprocal,  # Generate symmetrical matrix
                                         bool_suppl_professional_cnt = TRUE,  # Include supplementary working contacts
                                         bool_hhmatrix_selection = FALSE,
                                         wave = 'All',  # Include all waves
                                         quiet = TRUE)

# Retrieve sorted list of unique wave identifiers from the survey data
id_waves <- gtools::mixedsort(unique(survey_object_CoMix$participants$wave))




# WAVES DATES------

# Retrieve and process the dates of survey waves for Belgium from the CoMix survey
waves_df <- gen_w_dates()  # Returns a data_frame with columns: DATE, Participants (Freq), and Id

# Initialize simulated susceptible dates starting from 1st March 2020
sm_dates <- as.Date('2020-03-01') + sm_susceptible[1,,1]  # Vector of dates corresponding to simulation days

# Initialize arrays to store start, end, and weighted average dates of each wave
w_dates <- rep('', length(id_waves))
w_dates_start <- rep('', length(id_waves))
w_dates_end <- rep('', length(id_waves))
w_dates_weights <- rep(as.Date('2020-03-08'), length(id_waves))

# Calculate weighted average dates for each wave where days with more participants have greater influence
for (i in 1:length(id_waves)) {
  # Identify the first and last dates for each wave
  w_dates_start[i] <- waves_df$DATE[which(waves_df$Id == i)][1]
  w_dates_end[i] <- tail(waves_df$DATE[which(waves_df$Id == i)], 1)
  n_days <- length(waves_df$DATE[which(waves_df$Id == i)])
  weights <- rep(0, n_days)
  
  # Calculate weights based on the frequency of participants and the number of days since the start of the wave
  for (j in 2:n_days) {
    d_delta <- as.numeric(difftime(dmy(waves_df$DATE[which(waves_df$Id == i)][j]), dmy(w_dates_start[i]), units = "days"))
    weights[j] <- waves_df$Freq[which(waves_df$Id == i)][j] * d_delta
  }
  weig_delta <- round(sum(weights) / sum(waves_df$Freq[which(waves_df$Id == i)]))
  w_dates_weights[i] <- as.Date(dmy(w_dates_start[i]) + days(weig_delta))
}

# Parse dates to ensure they are in the correct format
w_dates_start <- dmy(w_dates_start)
w_dates_end <- dmy(w_dates_end)
w_dates <- as.character(w_dates_weights)

# Define the horizon of the observation to February 2022 (last_date represents the day count from the start)
last_date <- 42

# Adapt the wave dates to match the observation time interval, filtering to include only relevant dates
w_dates <- as.Date(w_dates)
w_dates_start <- w_dates_start[1:last_date]
w_dates_end <- w_dates_end[1:last_date]
w_dates <- w_dates[1:last_date]
w_dates <- w_dates[w_dates <= sm_dates[length(sm_dates)]]  # Filter to include dates not beyond the last simulation date

# EPIDEMIOLOGICAL ASSUMPTIONS-----
# Set epidemiological parameters using posterior mean values from Abrams et al.
# This configuration specifically mimics the age distribution system of schools in Belgium (case2).

params <- set_epi_params("case2")  # Retrieve parameters for the specified epidemiological setting
age_breaks <- params$Age          # Age intervals for epidemiological data
Asusc_vec <- params$A             # Age-specific susceptibility to infection
Hinf_vec <- params$H              # Age-specific infectiousness
p <- params$p_asym                # Age-specific proportion of asymptomatic infections
phi_0 <- params$phi_0             # Age-specific proportion of mildly symptomatic infections
omega_1 <- params$omega_1         # Removal rate for severe symptoms
p_mask <- params$p_mask           # Age-specific mask usage rates
prop_ratio <- 0.51                # Infectiousness ratio (Symptomatic to Asymptomatic)
gamma <- 0.729                    # Removal rate for exposed individuals
theta <- 0.475                    # Removal rate for pre-symptomatic individuals
sigma1 <- 0.24                    # Removal rate for asymptomatic infections
sigma2 <- 0.756                   # Recovery rate for mildly symptomatic infections

# Construct age-specific vectors and matrices for disease progression and transmission
Sigma <- sigma2 * phi_0           # Recovery rate matrix for mildly symptomatic individuals
Psi <- sigma2 * (1 - phi_0)       # Progression rate matrix from symptomatic to more severe stages

## Q-SUSCEPTIBILITY
# Construct a diagonal matrix for age-specific susceptibility
A <- diag(Asusc_vec)

## Q-INFECTIVITY
# Construct a diagonal matrix for age-specific infectivity
H <- diag(Hinf_vec)

# Setup an identity matrix for model operations
nbreaks <- length(age_breaks)
I <- diag(rep(1, nbreaks))  # Identity matrix used in matrix operations



# SUSCEPTIBLE MATRIX----
# Construct matrices containing estimated susceptibles S(s_ij) for each age group i during wave j.
# Average across all iterations to compute the mean value s_ij.

# Initialize matrices for susceptibles and their lower and upper confidence bounds for age groups (0,10,...,80,90+)
SW_90 <- matrix(0, 10, length(w_dates))  # Susceptible counts
SW_90_low <- matrix(0, 10, length(w_dates))  # Lower bounds
SW_90_up <- matrix(0, 10, length(w_dates))  # Upper bounds

# Fill the matrices with simulated susceptible counts and confidence intervals
for (i in 1:10) {
  for (j in 1:length(w_dates)) {
    index <- which(sm_dates == w_dates[j])
    SW_90[i, j] <- mean(sm_susceptible[, index, i + 1], trim = 0.1)  # Mean with 10% trimmed mean
    SW_90_low[i, j] <- quantile(sm_susceptible[, index, i + 1], .025)  # 2.5th percentile
    SW_90_up[i, j] <- quantile(sm_susceptible[, index, i + 1], .975)  # 97.5th percentile
  }
}

# Adjust susceptibles to aggregated age classes, summing older age groups for broader categories
# Aggregating age groups [70,80), [80,90), 90+ into a single class 70+
SW_70 <- apply(SW_90[8:10,], 2, sum)  # Sum values across selected age groups for each wave
SW_70 <- as.matrix(rbind(SW_90[1:7,], SW_70))  # Combine with other age groups
SW_70_low <- apply(SW_90_low[8:10,], 2, sum)
SW_70_low <- as.matrix(rbind(SW_90_low[1:7,], SW_70_low))
SW_70_up <- apply(SW_90_up[8:10,], 2, sum)
SW_70_up <- as.matrix(rbind(SW_90_up[1:7,], SW_70_up))


# Extract demographic data to match the structure of the susceptible simulations for age brackets 0-10, 10, 20..., 70-80, 90+
set.seed(2345)  # Ensure reproducibility
matrixCoMix <- contact_matrix(
  survey = survey_object_CoMix,  # Primary survey data
  age.limits = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90),  # Define age breaks
  n = 1,  # Number of bootstrap samples
  missing.contact.age = "sample",  # Handling missing age data by sampling
  estimated.contact.age = "sample",  # Same for estimated ages
  weigh.age = TRUE,  # Apply weighting based on age
  weigh.dayofweek = TRUE,  # Apply weighting based on the day of the week
  return.demography = TRUE,  # Return demographic data
  weight.threshold = 3  # Maximum weight for a single participant
)
N_0_90 <- matrixCoMix$demography[['population']]  # Extract demographic population vector

# Preparing for contact matrix adaptation by matching demographic data to simulation data
set.seed(2345)
survey_object_CoMix <- get_survey_object(
  country = country,
  daytype = "All contacts",
  touch = touch,
  duration = duration,
  gender = gender,
  cnt_location = cnt_location,
  bool_reciprocal = bool_reciprocal,
  bool_suppl_professional_cnt = bool_suppl_professional_cnt,
  bool_hhmatrix_selection = bool_hhmatrix_selection,
  wave = id_waves[9],  # Specific wave selection
  quiet = TRUE
)

matrixCoMix <- contact_matrix(
  survey = survey_object_CoMix,
  age.limits = age_breaks,  # Updated age breaks
  n = 1,
  missing.contact.age = "sample",
  estimated.contact.age = "sample",
  weigh.age = TRUE,
  weigh.dayofweek = TRUE,
  return.demography = TRUE,
  weight.threshold = 3
)
S_vec <- matrixCoMix$demography[['population']]  # Extract population data
names(S_vec) <- matrixCoMix$participants$age.group  # Label demographic data with age groups
N_vec <- S_vec  # Copy demographic data

# Adjust numerically simulated susceptible numbers to demographic age group distributions
# By imposing that the age-distribution over non-coinciding age intervals matches the demographic distribution
p1 <- N_vec[1] / N_0_90[1]
p2 <- 1 - p1
n1 <- N_0_90[1] - N_vec[1]
p3 <- (N_vec[2] - n1) / N_0_90[2]
p4 <- N_vec[3] / N_0_90[2]
n2 <- N_vec[4] - N_0_90[3]
p5 <- n2 / N_0_90[2]

props <- c(p1, p2, p3, p4, p5)


# Adjust susceptible estimates to the demographic age groups for analysis
SW <- adapt_susc(M = SW_70, age_groups = age_breaks, dates = w_dates, proportions = props)
SW_low <- adapt_susc(M = SW_70_low, age_groups = age_breaks, dates = w_dates, proportions = props)
SW_up <- adapt_susc(M = SW_70_up, age_groups = age_breaks, dates = w_dates, proportions = props)


# PRE-PANDEMIC CONTACT DATA
# Utilize social contact data from 2011 as described by Van Hoang et al. for Belgium.

# Set a fixed seed for reproducibility
set.seed(2345)

# Fetch the pre-pandemic survey object from Belgium 2010 survey by Van Hoang, configuring various parameters for contact data extraction
survey_object_CoMix <- get_survey_object(
  country = "Belgium 2010* (Van Hoang 2020)",
  daytype = daytype,
  touch = touch,
  duration = duration,
  gender = gender,
  cnt_location = cnt_location,
  bool_reciprocal = bool_reciprocal,              # Should the data generate a symmetrical matrix?
  bool_suppl_professional_cnt = bool_suppl_professional_cnt,   # Include supplementary working contacts?
  bool_hhmatrix_selection = bool_hhmatrix_selection,      # Select just household members' contact?
  wave = "All",                      # Include all waves 
  quiet = TRUE
)

# Generate contact matrices using the survey data
matrixCoMix_start <- contact_matrix(
  survey = survey_object_CoMix,
  age.limits = age_breaks,
  n = 1,  # Number of bootstrap samples
  missing.contact.age = "sample",  # Handle missing contact ages by sampling from population
  estimated.contact.age = "sample",
  weigh.age = TRUE,  # Apply weighting based on age
  weigh.dayofweek = TRUE,  # Apply weighting based on the day of the week
  return.demography = TRUE,
  weight.threshold = 3  # Maximum weight for one participant
)

# Assign age group names for demography
age_names <- matrixCoMix_start$demography$age.group
di_names <- list(age_names, age_names)

# Initialize matrices for asymptomatic and pre-symptomatic contact configurations
B_mean_asym <- matrix(0, nbreaks, nbreaks)
B_asym <- matrix(0, nbreaks, nbreaks)

# Process and normalize the contact data based on demographic data
if ('demography' %in% names(matrixCoMix_start) && !any(is.na(matrixCoMix_start$matrix))) {
  num_age_groups <- nrow(matrixCoMix_start$demography)
  pop_matrix <- matrix(rep(S_vec, num_age_groups), ncol = num_age_groups, byrow = TRUE)
}

# Symmetrize the contact matrix if specified
if (bool_reciprocal == TRUE && all(cnt_location == opt_location)) {
  B_mean_asym <- symmetrize(matrixCoMix_start$matrix, N = S_vec)
  matrixCoMix_start$matrix_per_capita <- B_mean_asym / pop_matrix
  B_asym <- matrixCoMix_start$matrix_per_capita
} else if (bool_reciprocal == TRUE && !all(cnt_location == opt_location)) {
  B_mean_asym <- symmetrize(matrixCoMix_start$matrix, N = S_vec)
  matrixCoMix_start$matrix_per_capita <- B_mean_asym / pop_matrix
  B_asym <- matrixCoMix_start$matrix_per_capita
  warning('The contact matrices have been symmetrized, but some reported contacts may not be reciprocal, e.g., customer(work) - client(leisure)')
} else {
  B_mean_asym <- matrixCoMix_start$matrix
  matrixCoMix_start$matrix_per_capita <- B_mean_asym / pop_matrix
  B_asym <- matrixCoMix_start$matrix_per_capita
}

# Handle location-specific behavioral changes and adjust the symptomatic contact matrices accordingly
B_mean_sym <- matrix(0, nbreaks, nbreaks)
B_sym <- matrix(0, nbreaks, nbreaks)
for (loc in opt_location) {
  set.seed(2345)
  survey_object_CoMix <- get_survey_object(
    country = "Belgium 2010* (Van Hoang 2020)",
    daytype = daytype,
    touch = touch,
    duration = duration,
    gender = gender,
    cnt_location = loc,
    bool_reciprocal = FALSE,  # Should the data generate a symmetrical matrix?
    bool_suppl_professional_cnt = bool_suppl_professional_cnt,
    bool_hhmatrix_selection = bool_hhmatrix_selection,
    wave = "All",
    quiet = TRUE
  )
  
  matrixCoMix_start <- contact_matrix(
    survey = survey_object_CoMix,
    age.limits = age_breaks,
    n = 1,
    missing.contact.age = "sample",
    estimated.contact.age = "sample",
    weigh.age = TRUE,
    weigh.dayofweek = TRUE,
    return.demography = TRUE,
    weight.threshold = 3
  )
  
  if ('demography' %in% names(matrixCoMix_start) && !any(is.na(matrixCoMix_start$matrix))) {
    num_age_groups <- nrow(matrixCoMix_start$demography)
    pop_matrix <- matrix(rep(S_vec, num_age_groups), ncol = num_age_groups, byrow = TRUE)
  }
  
  C_mean <- matrixCoMix_start$matrix
  # Weighted sum of contact matrices across locations with specific factors
  B_mean_sym <- B_mean_sym + (1 * (loc == "Home") + 0.09 * (loc == "Work") + 0.09 * (loc == "School") + 0.13 * (loc == "Transport") + 0.06 * (loc == "Leisure") + 0.25 * (loc == "Otherplace")) * C_mean
}

# Symmetrize and normalize the symptomatic contact matrices if required
if (bool_reciprocal == TRUE && all(cnt_location == opt_location)) {
  B_mean_sym <- symmetrize(B_mean_sym, N = S_vec)
  B_sym <- B_mean_sym / pop_matrix
} else if (bool_reciprocal == TRUE && !all(cnt_location == opt_location)) {
  B_mean_sym <- symmetrize(B_mean_sym, N = S_vec)
  B_sym <- B_mean_sym / pop_matrix
  warning('The contact matrices have been symmetrized, but some reported contacts may not be reciprocal: e.g., customer(work) - client(leisure)')
}


# NGM COMPOSITION-----

### Calibrate q-factor (residual) to match an initial R_0 = 3.4
# This step adjusts the proportionality factor, non-age-specific, so that the resulting R_0 matches 3.4 as reported in Coletti et al. 2021.

R_0 <- 0  # Initialize R_0
q_sym <- 0.1  # Initial guess for symptomatic infectiousness scaling factor
q_asym <- prop_ratio * q_sym  # Compute asymtomatic infectiousness scaling factor based on symptomatic

# Define matrices for asymptomatic and symptomatic transmission coefficients
Beta_asym <- matrix()
Beta_sym <- matrix()

# Define matrices for the various stages of the disease progression
G = diag(gamma, nbreaks, nbreaks)  # Exposed removal rate
O = diag(theta, nbreaks, nbreaks)  # Presymptomatic removal rate
S1 = diag(sigma1, nbreaks, nbreaks)  # Asymptomatic removal rate
S2 = diag(Sigma)  # Recovery rate for mildly symptomatic
R = diag(Psi)  # Progression from symptomatic to more severe stages
Omega = diag(omega_1)  # Removal rate for severe symptoms

# Calculate factors for disease progression matrix calculations
factor1 <- (theta * sigma1) ** (-1)  # Inverse of product of theta and sigma1
factor2 <- solve(Omega * (R + S2))  # Solve linear equations for progression factors

# Compose matrices for different disease states
D_1_asym <- (S1 + (O * p)) %*% diag(factor1, nbreaks, nbreaks)  # Asymptomatic disease matrix
D_2_sym <- factor2 * ((Omega + R) * (1 - p))  # Symptomatic disease matrix

# Compose transmission terms
Beta_sym <- A %*% B_sym %*% (q_sym * H)  # Symptomatic transmission matrix
Beta_asym <- A %*% B_asym %*% (q_asym * H)  # Asymptomatic transmission matrix
Beta_asym <- diag(S_vec) %*% Beta_asym  # Adjust by susceptible population
Beta_sym <- diag(S_vec) %*% Beta_sym
K_1 <- Beta_asym %*% D_1_asym + Beta_sym %*% D_2_sym  # Total NGM
R_0 <- max(Re(eigen(K_1)$values))  # Compute the dominant eigenvalue (R_0)

# Adjust q factors based on observed R_0
q = 3.4 / R_0 
q_asym <- q_asym * q
q_sym <- q_sym * q
Beta_asym <- q * Beta_asym
Beta_sym <- q * Beta_sym
K_1 <- K_1 * q  # Scale NGM by new q

# Compute dominant eigenstructure for R_0
R_0_1 <- max(Re(eigen(K_1)$values))
w_1 <- abs(eigen(K_1)$vectors[,1])
v_1 <- abs(eigen(t(K_1))$vectors[,1])
w_1 <- w_1 / sum(w_1)  # Normalize leading right eigenvector
v_1 <- v_1 / sum(v_1 * w_1)  # Normalize leading left eigenvector


