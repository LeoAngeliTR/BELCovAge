#___________________________________________________________________________
# This script executes the longitudinal sensitivity analysis as described in the study by Angeli et al., 2024.
# It accounts for the impact of reporting fatigue and recalculates contact matrices for each wave.
# Then the NGM is evaluated for each survey wave as weel as the sensitivity and elasticity matrices
# (LONGITUDINAL SENSITIVITY ANALYSIS (BELGIUM))
# Copyright 2020, SIMID, UNIVERSITY OF ANTWERP & HASSELT UNIVERSITY
#___________________________________________________________________________

# PRELIMINARY OPERATIONS------
# Load required preliminary configurations and data setup scripts.
source('./preliminary.R')

## Fatigue effect correction----
# Adjust contact matrices up to wave 43 to account for underreporting due to participant fatigue as per Loedy et al. 2023, Belgium.
folder1 <- './fatigue_correction/SCM_reci_predict_f.csv'  # File containing fatigue-corrected data
folder2 <- './fatigue_correction/SCM_reci_predict.csv'  # File containing original survey data
df_corr <- read.csv(folder1)  # Load corrected contact data
df_orig <- read.csv(folder2)  # Load original contact data

# Prepare a list to store corrections for each wave from the 9th to the last studied date.
correction_list_reci <- vector(mode = "list", length = length(9:last_date))

# Process each wave's data to calculate corrected contact matrices.
row_names <- unique(df_orig$age.group)  # Get unique age groups
col_names <- row_names  # Column names same as row names for matrix symmetry
cicle <- 1  # Initialize cycle counter for tracking progress

for (i in 9:length(w_dates)) {  # Start from wave 9 to the end
  temp_df_corr <- df_corr[df_corr$wave == paste0('Wave ', i),] 
  temp_df <- df_orig[df_orig$wave == paste0('Wave ', i),]
  M_wave <- matrix(nrow = length(row_names), ncol = length(col_names))
  M_wave_up <- matrix(nrow = length(row_names), ncol = length(col_names))
  M_wave_low <- matrix(nrow = length(row_names), ncol = length(col_names))
  
  # Calculate the matrix of correction factors for each age group interaction
  for (j in seq_len(length(row_names))) {
    for (k in seq_len(length(col_names))) {
      numerator <- temp_df_corr[temp_df_corr$age.group == row_names[j] & temp_df_corr$contact.age.group == col_names[k],]
      denominator <- temp_df[temp_df$age.group == row_names[j] & temp_df$contact.age.group == col_names[k],]
      M_wave[j, k] <- numerator$mean / denominator$mean
      M_wave_up[j, k] <- numerator$upper / denominator$upper
      M_wave_low[j, k] <- numerator$lower / denominator$lower
    }
  }
  correction_list_reci[[i]] <- list(Wave = i, M_perc_f = M_wave, M_perc_f_up = M_wave_up, M_perc_f_low = M_wave_low)
  cicle <- cicle + 1
}

# Select the range of waves for analysis
wrange <- 9:length(w_dates)  # Analyze from November 2020 to February 2022 
                             #(please select a wave range starting from wave 9 onwards)

### Initialize the list containing wave-specific contacts and matrices (NGM)
ngms <- vector(mode = "list", length = length(wrange))
c_asym <- vector(mode = "list", length = length(wrange))
c_sym <- vector(mode = "list", length = length(wrange))
w1_list <- vector(mode = "list", length = length(wrange))
v1_list <- vector(mode = "list", length = length(wrange))
S_tot_index <- vector(mode = "list", length = length(wrange))
S_tot_index_mid <- vector(mode = "list", length = length(wrange))

## Store Belgian demography vector----
# Load demographic data to adjust contact matrices according to population distribution
N_vec <- matrixCoMix$demography[['population']]
num_age_groups <- nrow(matrixCoMix$demography)
pop_matrix <- matrix(rep(N_vec, num_age_groups), ncol = num_age_groups, byrow = TRUE)
kl_weight = TRUE  # Weighting factor for adjusting contact matrices

# Load effective reproduction number data (R_t) from PCR analysis by Gressani et al. 2022
Rt_list <- readRDS("./data/list_avgRt_1_42.rds")
q_factor_sym <- list()
R_t = TRUE  # Set TRUE to evaluate effective reproduction number during the analysis
# Note: The cycle counter 'cicle' and timestamp 'old' are used for performance tracking and will be initialized before looping through waves.
cicle <- 1
old <- Sys.time()


#Example: evaluate wave-specific contact matrix-----

wave <- 11
survey_object_CoMix_w<-get_survey_object(country=country,
                                         daytype=daytype,
                                         touch=touch,
                                         duration=duration,
                                         gender=gender,
                                         cnt_location=cnt_location,
                                         bool_reciprocal=bool_reciprocal,              ### should the data generate symmetrical matrix?
                                         bool_suppl_professional_cnt=bool_suppl_professional_cnt,   ### should supplementary working contacts be included?
                                         bool_hhmatrix_selection=bool_hhmatrix_selection,      ###select just household members contact
                                         wave=id_waves[wave],                      ### wave to be included (optional)
                                         quiet = TRUE)

matrixCoMix_w<-contact_matrix(survey=survey_object_CoMix, ## First and most important parameter: the survey object
                              age.limits = age_breaks,
                              n=1,                             ## number of bootstraps
                              missing.contact.age="sample",    ## What to do with contacts with missing contact age (possibilities: sample (sample from population), remove (remove them))
                              estimated.contact.age="sample",
                              weigh.age=T,                  ## TRUE if weight on the age is performed
                              weigh.dayofweek=T,            ## TRUE if weight on the day of the week is performed
                              return.demography   = TRUE,
                              weight.threshold = 3             ## maximum weight for one participant
)
#Visualize contact matrix
matrix_plot(matrixCoMix_w$matrix,main = "Daily Contacts")


# This loop processes each wave specified in 'wrange', adjusting for changes in contact patterns,
# susceptibles shifts, and epidemiological parameters to update the next-generation matrix (NGM) for each wave.

for(wave in wrange){
  
  # Store previous wave's data for comparison and reference in current calculations
  B_mean_ref_asym <- B_mean_asym
  B_mean_ref_sym <- B_mean_sym
  B_ref_asym <- B_asym
  B_ref_sym <- B_sym
  R_0_ref <- R_0_1
  K_1_ref <- K_1
  
  # Initialize contact matrices for the current wave
  B_mean_asym <- matrix(0, nbreaks, nbreaks)
  B_asym <- matrix(0, nbreaks, nbreaks)
  B_mean_sym <- matrix(0, nbreaks, nbreaks)
  B_sym <- matrix(0, nbreaks, nbreaks)
  # Adjust susceptible populations based on current wave data and reporting fatigue effects
  Susc <- SW[,wave]  # Extract susceptible counts for the current wave
  if(length(age_breaks) != 9){  # Adjust for different age group breakdowns if necessary
    Susc <- c(mean(SW[,wave][1:3]), mean(SW[,wave][4:5]), mean(SW[,wave][6:7]), mean(SW[,wave][8:9])) 
  }
  S_vec <- S_vec * (!R_t) + Susc * R_t  # Update susceptible vector dynamically based on R_t flag
  
  # Load previously evaluated contact matrices ensuring reciprocity
  fname <- "./data/raw_contact_list.rds"
  matrices_epi_list <- readRDS(fname)
  matrices_epi_list <- lapply(matrices_epi_list, symmetrize, N = N_vec)
  B_mean_asym <- matrices_epi_list[[wave]]
  B_asym <- B_mean_asym / pop_matrix
  
  # Compute symptomatic contact matrices considering specific locations and settings
  B_mean_sym_loc <- matrix(0, nbreaks, nbreaks)
  for(loc in opt_location){
    set.seed(2345)
    survey_object_CoMix_loc <- get_survey_object(
      country = country, daytype = daytype, touch = "All", duration = duration,
      gender = gender, cnt_location = loc, bool_reciprocal = F,
      bool_suppl_professional_cnt = bool_suppl_professional_cnt,
      bool_hhmatrix_selection = bool_hhmatrix_selection, wave = id_waves[wave], quiet = TRUE)
    
    matrixCoMix_loc <- contact_matrix(
      survey = survey_object_CoMix_loc, age.limits = age_breaks, n = 1,
      missing.contact.age = "sample", estimated.contact.age = "sample",
      weigh.age = TRUE, weigh.dayofweek = TRUE, return.demography = TRUE, weight.threshold = 3)
    
    C_mean <- matrixCoMix_loc$matrix
    B_mean_sym_loc <-  B_mean_sym_loc + (1 * (loc == "Home") + 0.09 * (loc == "Work") + 0.09 * (loc == "School") +
                         0.13 * (loc == "Transport") + 0.06 * (loc == "Leisure") + 0.25 * (loc == "Otherplace")) * C_mean
  }
  B_mean_sym <- B_mean_sym_loc
  
  # Impose reciprocity and ensure non-negative values
  if(bool_reciprocal && all(cnt_location == opt_location)){
    B_mean_sym <- symmetrize(B_mean_sym, N = N_vec)
  } else if(bool_reciprocal && !all(cnt_location == opt_location)){
    B_mean_sym <- symmetrize(B_mean_sym, N = N_vec)
    warning('The contact matrices have been symmetrized, but some reported contacts may not be reciprocal: e.g., customer(work) - client(leisure)')
  }
  B_sym <- non_neg(B_mean_sym / pop_matrix)
  
  ## Apply fatigue correction to the contact matrices -----
  ## (from wave 9 onwards) See Loedy et al. 2022
    M_corr_f <- matrix(0,nrow(B_asym),ncol(B_asym))
    M_corr_f_up<- matrix(0,nrow(B_asym),ncol(B_asym))
    M_corr_f_low<- matrix(0,nrow(B_asym),ncol(B_asym))
    for (ir in seq_len(nrow(B_asym))) {
      for (jc in seq_len(ncol(B_asym))) {
        i_row <- 1*(ir<=3)+2*(ir>3 & ir<9)+3*(ir>8)
        j_col <- 1*(jc<=3)+2*(jc>3 & jc<9)+3*(jc>8)
        M_corr_f[ir,jc] <- correction_list_reci[[wave]]$M_perc_f[i_row,j_col] 
        M_corr_f_up[ir,jc] <- correction_list_reci[[wave]]$M_perc_f_up[i_row,j_col] 
        M_corr_f_low[ir,jc] <- correction_list_reci[[wave]]$M_perc_f_low[i_row,j_col] 
      }
    }
    B_asym <- B_asym*M_corr_f
    B_sym <- B_sym*M_corr_f
  
  
  
  #Calculate wave-specific NGM ----
  
  G=diag(gamma, nbreaks,nbreaks)
  O=diag(theta, nbreaks,nbreaks)
  S1=diag(sigma1, nbreaks,nbreaks)
  S2=diag(Sigma)
  R=diag(Psi)
  Omega=diag(omega_1)
  factor1 <- (theta*sigma1)**(-1)
  factor2 <- solve(Omega*(R+S2))
  D_1_asym <- (S1+(O*p))%*%diag(factor1,nbreaks,nbreaks)
  D_2_sym<- factor2*((Omega+R)*(1-p))
  
  
  ###Calculate q-factor (residual)####
  R_0 <- 0
  q_sym <- 0.1  
  q_asym <- prop_ratio*q_sym
  Beta_asym <- matrix()
  Beta_sym <- matrix()
  K_1 <- matrix()
  #composing transmission terms
  Beta_sym <-A%*%B_sym%*%(q_sym*H)
  Beta_asym <- A%*%B_asym%*%(q_asym*H)
  Beta_asym <- diag(S_vec)%*%Beta_asym
  Beta_sym <- diag(S_vec)%*%Beta_sym
  K_1 <- Beta_asym%*%D_1_asym +Beta_sym%*%D_2_sym
  R_0 <-max(Re(eigen(K_1)$values))
  q=Rt_list$Mean[wave]/R_0 
  q_asym <- q_asym*q
  q_sym <- q_sym*q
  Beta_asym <- q*Beta_asym
  Beta_sym <- q*Beta_sym
  K_1 <- K_1*q
  q_factor_sym[[cicle]] <- q_sym
  
  ##for plots
  wave_title <- toString(waves)
  
  ##Beta-Matrices composed by transmission terms
  Su <-diag(S_vec)
  Beta_sym <-Su%*%A%*%B_sym%*%(q_sym*H)
  Beta_asym <- Su%*%A%*%B_asym%*%(q_asym*H)
  
 
  ###NGM & eigenstructure----
  K_1 <- Beta_asym%*%D_1_asym +Beta_sym%*%D_2_sym
  R_0_1 <- max(Re(eigen(K_1)$values))
  w_1 <- abs(eigen(K_1)$vectors[,1])
  v_1 <- abs(eigen(t(K_1))$vector[,1])
  w_1 <- w_1/sum(w_1)  
  v_1 <- v_1/sum(v_1*w_1)
  age_names <- matrixCoMix_start$demography$age.group
  di_names <- list(age_names,age_names)
  names(v_1)<- age_names
  names(w_1)<- age_names

  
  #Mid-point approx NGM
  
  if(kl_weight==T){
    weight1 <- sum(B_mean_ref_asym%*%diag(p) + B_mean_ref_sym%*%diag(1-p))
    weight2 <- sum(B_mean_sym%*%diag(p)+B_mean_asym%*%diag(1-p))
    mid_point_KL <- (K_1_ref*weight1+K_1*weight2)/(weight1+weight2)
  }else{
    mid_point_KL <- (K_1_ref+K_1)/2
  }
  
  #Overall sensitivity Index-----
  #upper bound for the magnitude of change in Rt,
  #produced by K perturbation
  S_tot_index[cicle]<-sqrt(sum(hadamard.prod(S_built_1,S_built_1))) 
  S_tot_index_mid[cicle]<-sqrt(sum(calc.sens(mid_point_KL)*calc.sens(mid_point_KL)))
  
 
  ##Store wave-specific contact matrices----
  c_asym[[cicle]] <- B_mean_asym
  c_sym[[cicle]] <- B_mean_sym
  
  
  ##Store wave-specific NGM----
  ngms[[cicle]] <- K_1
  
  ##Store wave-specific eigenvectors
  w1_list[[cicle]] <- w_1
  v1_list[[cicle]]<-v_1
   
  cicle <- cicle+1
  print(paste0('Wave:',wave," - Rt_w: ",R_0_1, "- Rt_pcr: ",Rt_list$Mean[wave]))
  
} 
# print elapsed time
new <- Sys.time() - old # calculate difference
print(new) # print in nice format

#Sensitivity & Elasticity matrices (lists) - set to identity matrix the first 8 elements
sens_list<- lapply(ngms, calc.sens)
el_list <- list()
el_list<- lapply(ngms, calc.el)




