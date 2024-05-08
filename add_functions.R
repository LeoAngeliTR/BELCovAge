#___________________________________________________________________________
# This file contains the functions necessary to run the analysis described in Angeli et al., 2024
# (LONGITUDINAL SENSITIVITY ANALYSIS (BELGIUM))
# 
# Copyright 2020, SIMID, UNIVERSITY OF ANTWERP & HASSELT UNIVERSITY
#___________________________________________________________________________

# Define necessary classes
setClass("MatrixProp", slots=list(
  determinant="numeric", 
  primitive="logical", 
  irreducible="logical",
  ev="eigen",           # Eigenvalues
  ev_left="eigen",      # Left eigenvalues
  R_0='numeric'))       # Basic reproduction number

# Function to calculate per capita matrices
get_percapita_matrix <- function (cnt_matrix, pop_vec, byrow=FALSE) {
  res <- matrix(0, nrow(cnt_matrix), ncol(cnt_matrix))
  # byrow=TRUE assumes the contact is from the perspective of the infected group
  for(i in 1:nrow(cnt_matrix)) {
    for (j in 1:ncol(cnt_matrix)) {
      res[i,j] <- if(byrow) cnt_matrix[i,j] / pop_vec[i] else cnt_matrix[i,j] / pop_vec[j]
    }
  }
  return(res)
}

# Function to build the next-generation matrix (NGM)
built_ngm <- function(S, qsusc, qinf, C_asym, C_sym, D1, D2, q, N_vec, prop_ratio) {
  S <- diag(S)
  C_asym <- get_percapita_matrix(C_asym, N_vec)
  C_sym <- get_percapita_matrix(C_sym, N_vec)
  q_a <- prop_ratio * q
  Asym <- S %*% qsusc %*% C_asym %*% (q_a * qinf)
  Sym <- S %*% qsusc %*% C_sym %*% (q * qinf)
  Kbuilt <- Asym %*% D1 + Sym %*% D2
  return(list(Asymp=Asym, Symp=Sym, NGM=Kbuilt))
}

# Function to ensure matrix symmetry or reciprocity
symmetrize <- function(A, N=NULL) {
  res <- matrix(0, nrow(A), ncol(A))
  if(is.null(N)) {
    res <- if(is.symmetric.matrix(A)) A else (A + t(A)) / 2
  } else {
    for (i in 1:nrow(A)) {
      for (j in 1:ncol(A)) {
        res[i, j] <- (A[i, j] * N[i] + A[j, i] * N[j]) / (2 * N[i])
      }
    }
  }
  return(res)
}

# Function to create an object with matrix properties
mat_obj <- function(A) {
  obj <- new("MatrixProp")
  obj@determinant <- det(A)
  obj@ev <- eigen(A)
  obj@ev_left <- eigen(t(A))
  return(obj)
}

# Function to set epidemiological parameters based on age groups
set_epi_params <- function(case) {
  
  # Define epidemiological parameters based on the specified case
  age_breaks <- switch(case,
                       "default" = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90),
                       "case1" = c(0, 10, 20, 30, 40, 50, 60, 70),
                       "case2" = c(0, 6, 12, 18, 30, 40, 50, 60, 70),
                       "case3" = c(0, 18, 40, 60))
  
  Asusc_vec <- switch(case, # age-specific susceptibility
                      "default" = c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74, 0.74),
                      "case1" = c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74),
                      "case2" = c(0.4, 0.39, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74),
                      "case3" = c(0.4, 0.83, 0.81, 0.81))
  
  Hinf_vec <- switch(case, # age-specific infectiousness
                     "default" = c(0.55, 0.55, 0.59, 0.7, 0.76, 0.9, 0.9, 0.9, 0.9, 0.9),
                     "case1" = c(0.55, 0.55, 0.59, 0.7, 0.76, 0.9, 0.9, 0.9),
                     "case2" = c(0.54, 0.55, 0.56, 0.59, 0.7, 0.76, 0.9, 0.99, 0.99),
                     "case3" = c(0.55, 0.65, 0.83, 0.99))
  
  p <- switch(case, # age-specific probability of asymptomatic infection
              "default" = c(0.94, 0.90, 0.84, 0.61, 0.49, 0.21, 0.02, 0.02, 0.2, 0.2),
              "case1" = c(0.94, 0.90, 0.84, 0.61, 0.49, 0.21, 0.02, 0.02),
              "case2" = c(0.94, 0.92, 0.90, 0.84, 0.61, 0.49, 0.21, 0.02, 0.02),
              "case3" = c(0.92, 0.75, 0.35, 0.2))
  
  phi_0 <- switch(case, # Age-specific probability of developing only mild symptoms
                  "default" = c(0.972, 0.992, 0.984, 0.987, 0.977, 0.971, 0.958, 0.936, 0.956, 0.926),
                  "case1" = c(0.972, 0.992, 0.984, 0.987, 0.977, 0.971, 0.958, 0.936),
                  "case2" = c(0.972, 0.982, 0.992, 0.984, 0.987, 0.977, 0.971, 0.958, 0.936),
                  "case3" = c(0.972, 0.982, 0.992, 0.984))
  
  omega_1 <- switch(case, # Severe symptomatic removal rate
                    "default" = c(0.167, 0.095, 0.099, 0.162, 0.338, 0.275, 0.343, 0.338, 0.334, 0.302),
                    "case1" = c(0.167, 0.095, 0.099, 0.162, 0.338, 0.275, 0.343, 0.338),
                    "case2" = c(0.167, 0.131, 0.095, 0.099, 0.162, 0.338, 0.275, 0.343, 0.338),
                    "case3" = c(0.12, 0.13, 0.3, 0.34))
  
  p_mask <- switch(case, # Mask compliance rates
                   "default" = c(0.37, 0.37, 0.37, 0.41, 0.41, 0.41, 0.57, 0.57, 0.57, 0.57),
                   "case1" = c(0.37, 0.37, 0.37, 0.41, 0.41, 0.41, 0.57, 0.57),
                   "case2" = c(0.37, 0.37, 0.37, 0.37, 0.41, 0.41, 0.41, 0.57, 0.57),
                   "case3" = c(0.37, 0.39, 0.41, 0.57))
  
  # Compile all epidemiological parameters into a list
  params <- list(Age=age_breaks, A=Asusc_vec, H=Hinf_vec, p_asym=p, phi_0=phi_0, omega_1=omega_1, p_mask=p_mask)
  
  # Error handling to ensure all parameters are properly initialized
  if (any(sapply(params, is.null)) || any(sapply(params, function(x) any(is.na(x))))) {
    stop("An error occurred: check the values of 'case'. Ensure all parameters are defined.")
  } else {
    print("Epidemiological parameters initialized successfully.")
  }
  
  return(params)
}

#This function adjusts susceptibility data based on demographic proportions.

adapt_susc <- function(M, age_groups, dates, proportions) {
  # Initialize a matrix to store adapted susceptibility values
  R <- matrix(NA, nrow = length(age_groups), ncol = length(dates))
  for (i in 1:ncol(M)) {
    # Adjust the susceptibility across age groups based on predefined proportions
    R[1, i] <- M[1, i] * proportions[1]
    R[2, i] <- M[1, i] * proportions[2] + M[2, i] * proportions[3]
    R[3, i] <- M[2, i] * proportions[4]
    R[4, i] <- M[2, i] * proportions[5] + M[3, i]
    R[5:9, i] <- M[4:nrow(M), i]
  }
  return(round(R))
}


#This function plots a matrix with a heatmap, potentially used for visualizing contact matrices.

plot_matrix <- function(mij, title='', axisnames=NULL, xlab=NULL, ylab=NULL) {
  # Set column names if provided
  if (!is.null(axisnames)) {
    colnames(mij) <- axisnames
  }
  
  # Return NA if the matrix is empty
  if (all(is.na(mij))) {
    return(NA)
  }
  
  # Define the color palette for the heatmap
  redc <- rev(heat.colors(100))
  
  # Set margins and plot parameters
  par(mar = c(5, 6, 2, 2), mgp = c(3, 0.5, 0))
  p <- image(mij, 
             xlab = ifelse(is.null(xlab), "Age of participant (year)", xlab),
             ylab = ifelse(is.null(ylab), "Age of contact (year)", ylab),
             col = redc, 
             main = title)
  axis(1, labels = axisnames, cex.axis = 0.9)
  axis(2, labels = axisnames, cex.axis = 0.9, las = 1)
  return(p)
}

#This function ensures that all entries in a matrix are non-negative, replacing negative values with 0.1.
non_neg <- function(M = matrix()) {
  # Replace negative values in the matrix with 0.1
  if (any(M < 0)) {
    M[M < 0] <- 0.1
  }
  return(M)
}
#This function deals with NA values in a matrix by using another matrix to estimate missing entries.

deal_nas <- function(A, B, perc = 1) {
  # Handle missing participant data in matrix A using matrix B
  missing_pop <- NULL
  for (j in 1:nrow(A$matrix)) {
    if (all(is.na(A$matrix[j,]))) {
      A$matrix[j,] <- B$matrix[j,] * perc
      missing_pop <- c(missing_pop, B$demography$population[j])
    }
  }
  return(list(contact_matrix = A$matrix, population = c(missing_pop, A$demography$population)))
}


#his function creates a bar plot for visualizing data arrays.
plot_bars <- function(mij, title = '', extratitle = NULL, lab_y = '', col = "lightcyan", adj_ylim = 1,
                      adj_labels = 0, setylim = FALSE, yaxis = NULL, barnames = NULL, change_xlab = FALSE, x_lab = NULL) {
  if (all(is.na(mij))) {
    return(NA)  # Return NA if all data points are NA
  }
  
  # Determine the y-axis limits
  if (all(mij >= 0)) {
    lim_y = c(0, max(mij) * adj_ylim)
  } else {
    # Adjust limits for negative values in the data
    lim_y = c(min(mij), max(mij)) * adj_ylim
  }
  
  y_lim = if (setylim) yaxis else lim_y  # Set custom or automatic y-axis limits
  
  # Bar plot setup
  if (!is.null(barnames)) {
    names(mij) <- barnames  # Assign names to bars if provided
  }
  
  # Set x-axis labels based on condition
  x_lab_final = if (change_xlab && !is.null(x_lab)) x_lab else 'Age groups-participants'
  
  x <- barplot(mij, main = paste(title, extratitle, sep = ' '), ylim = y_lim, xlab = x_lab_final, ylab = lab_y, col = col)
  
  # Adjust labels on bars for better visibility
  lab_pos = adjustTextPosition(mij, adj_labels)
  
  # Place text labels on the bars
  text(x, lab_pos, labels = formatNumbers(mij, digits = 2), pos = 3)
  
  return(x)
}

# Helper function to adjust the position of labels on the bars
adjustTextPosition <- function(values, adjustment) {
  ifelse(values >= 0, values + adjustment, values - adjustment)
}

# Helper function to format numbers, showing significant digits
formatNumbers <- function(numbers, digits = 2) {
  signif(numbers, digits = digits)
}



#This function creates perturbation matrices based on specific vector and matrix inputs.

build_perturbation_matrices <- function(nrows, ncols=NULL, w, v, dimnames, KL, R_0) {
  # Set ncols to nrows if it is not specified
  if (is.null(ncols)) { ncols <- nrows }
  
  # Initialize sensitivity (S) and elasticity (E) matrices
  S <- matrix(0, nrow = nrows, ncol = ncols)
  E <- matrix(0, nrow = nrows, ncol = ncols)
  
  # Populate the matrices with calculated values
  for (i in 1:nrows) {
    for (j in 1:ncols) {
      S[i, j] <- v[i] * w[j]
      E[i, j] <- S[i, j] * KL[i, j] / R_0
    }
  }
  return(list(S = S, E = E))
}

#This function calculates a sensitivity matrix based on the given matrix G, using eigenvalues and eigenvectors.

calc.sens <- function(G, normalize = TRUE) {
  ev <- eigen(G)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  U <- ev$vectors
  u <- abs(Re(U[, lmax]))
  ev_left <- eigen(t(G))
  lmax_l <- which(Re(ev_left$values) == max(Re(ev_left$values)))
  V <- ev_left$vectors
  v <- abs(Re(V[, lmax_l]))
  
  # Normalize the vectors if specified
  if (normalize) {
    u <- u / sum(u)
    v <- v / (as.numeric(v %*% u))
  }
  
  s <- as.matrix(v %o% u)  # Compute outer product of v and u
  return(s)
}

#This function calculates an elasticity matrix for a given matrix G, based on the method in calc.sens.

calc.el <- function(G, normalize = TRUE) {
  M <- calc.sens(G, normalize)  # Calculate sensitivity matrix
  lambda <- Re(eigen(G)$values[1])  # Compute the leading eigenvalue
  E <- M * G / lambda  # Adjust the matrix by the leading eigenvalue
  return(E)
}

#This function extracts the dominant right and left eigenvectors of a matrix G, with optional scaling.

top_eigenvectors <- function(G, scale = TRUE) {
  ev <- eigen(G)
  lmax <- which(Re(ev$values) == max(Re(ev$values)))
  U <- ev$vectors
  u <- abs(Re(U[, lmax]))
  ev_left <- eigen(t(G))
  lmax_l <- which(Re(ev_left$values) == max(Re(ev_left$values)))
  V <- ev_left$vectors
  v <- abs(Re(V[, lmax_l]))
  
  # Scale the vectors if specified
  if (scale) {
    u <- u / sum(u)
    v <- v / sum(v * u)
  }
  
  return(list(right = u, left = v))
}

gen_w_dates <- function() {
  # Purpose: Generates a data frame containing dates and IDs for survey waves from the CoMix survey data.
  
  # Initialize an empty data frame to store wave date and id information
  waves <- data.frame(DATE = character(0),
                      id = integer(0),
                      stringsAsFactors = FALSE)
  
  # Retrieve survey data for Belgium 2020 from the CoMix dataset
  survey_object_CoMix <- get_survey_object(country = "Belgium 2020 CoMix (Coletti 2020)",
                                           daytype = "All contacts",
                                           touch = "All contacts",
                                           duration = "All contacts",
                                           gender = "All",
                                           cnt_location = opt_location,
                                           bool_reciprocal = FALSE,            # Exclude symmetrical matrix generation
                                           bool_suppl_professional_cnt = TRUE, # Include supplementary working contacts
                                           bool_hhmatrix_selection = FALSE,
                                           wave = 'All',                      # Include all waves
                                           quiet = TRUE)
  
  # Extract and sort unique wave identifiers
  id_waves <- gtools::mixedsort(unique(survey_object_CoMix$participants$wave))
  
  # Iterate over each wave to process data
  for (w in id_waves) {
    # Retrieve survey data for each specific wave
    survey_object_CoMix <- get_survey_object(country = "Belgium 2020 CoMix (Coletti 2020)",
                                             daytype = "All contacts",
                                             touch = "All contacts",
                                             duration = "All contacts",
                                             gender = "All",
                                             cnt_location = opt_location,
                                             bool_reciprocal = FALSE,
                                             bool_suppl_professional_cnt = TRUE,
                                             bool_hhmatrix_selection = FALSE,
                                             wave = w,                          # Process the current wave
                                             quiet = TRUE)
    
    # Process the survey data to extract date information and ensure non-zero frequencies
    newrows <- data.frame(table(survey_object_CoMix[["participants"]][["sday_id"]])) %>%
      dplyr::filter(Freq != 0) %>%
      dplyr::mutate(DATE = format(as.Date(as.character(Var1), format = "%Y.%m.%d"), "%d/%m/%Y")) %>%
      dplyr::select(DATE, Freq)
    
    # Extract wave ID from the wave string and add to the dataset
    newrows$Id <- as.integer(strsplit(w, ':')[[1]][1])
    
    # Append the new rows to the main dataset
    waves <- rbind(waves, newrows)
  }
  
  return(waves)
}



# Define a function to generate the matrix for a given location
generate_location_matrix <- function(loca, N=1,
                                     ages,
                                     profession=NULL,
                                     country_ref,
                                     wave_id="All",
                                     corr_list=NULL) {
  
  survey_object <- get_survey_object(country=country_ref,    
                                     daytype=daytype,
                                     touch=touch,
                                     duration=duration,
                                     gender=gender,
                                     cnt_location=loca,
                                     cnt_profession = profession,
                                     bool_reciprocal=bool_reciprocal,
                                     bool_suppl_professional_cnt=bool_suppl_professional_cnt,
                                     bool_hhmatrix_selection=bool_hhmatrix_selection,
                                     wave=wave_id,
                                     quiet=TRUE)
  N_boot <- N
  if(wave_id!="All"){
    wave <-  which(wave_ids==wave_id)
    N_boot <- N_boot*(wave>8)+1*(wave<9)
    
  }
  
  matrix <- contact_matrix(survey=survey_object,
                           age.limits=ages,
                           n= N_boot,
                           missing.contact.age="sample",
                           estimated.contact.age="sample",
                           weigh.age=TRUE,
                           weigh.dayofweek=TRUE,
                           return.demography=TRUE,
                           weight.threshold=3)
  
  if (N==1) {
    ret <- matrix$matrix
    if(wave_id!="All"){
      
      if (wave<9 && country_ref== "Belgium 2020 CoMix (Coletti 2020)"){
        output_folder="./output_matrices_imputed_test/"
        fname=paste0(output_folder,"boot_",loca,"_wave_",wave,".rds")
        matrices_loc_list <- readRDS(fname)
        ret<-apply(simplify2array(matrices_loc_list), 1:2, mean)
        
      }else if(wave>8 && country_ref== "Belgium 2020 CoMix (Coletti 2020)" && !is.null(corr_list)){
        #adapt correction matrices to my age structure
        M_corr_f <- matrix(0,length(ages),length(ages))
        M_corr_f_up<- matrix(0,length(ages),length(ages))
        M_corr_f_low<- matrix(0,length(ages),length(ages))
        for (ir in seq_len(length(ages))) {
          for (jc in seq_len(length(ages))) {
            i_row <- 1*(ir<=3)+2*(ir>3 & ir<9)+3*(ir>8)
            j_col <- 1*(jc<=3)+2*(jc>3 & jc<9)+3*(jc>8)
            M_corr_f[ir,jc] <- corr_list[[wave]]$M_perc_f[i_row,j_col] 
            M_corr_f_up[ir,jc] <- corr_list[[wave]]$M_perc_f_up[i_row,j_col] 
            M_corr_f_low[ir,jc] <- corr_list[[wave]]$M_perc_f_low[i_row,j_col] 
          }
        }
        
        ret<- ret*M_corr_f
        
      }else{
        warning("WARNINNG NO CONTACT SELECTED: set country_ref = 'Belgium 2020 CoMix (Coletti 2020)' ")
        ret <- NULL
        }
      
      
    }
  }else{
    
    if(wave_id!="All"){
      wave <-  which(wave_ids==wave_id)
      if (wave<9 && country_ref== "Belgium 2020 CoMix (Coletti 2020)"){
        output_folder="./output_matrices_imputed_test/"
        fname=paste0(output_folder,"boot_",loca,"_wave_",wave,".rds")
        matrices_loc_list <- readRDS(fname)
        ret<-matrices_loc_list
        
      }else if(wave>8 && country_ref== "Belgium 2020 CoMix (Coletti 2020)" && !is.null(corr_list)){
        ret <- matrix$matrices
        corrected_mat <- correct_matrix_f(ages,corr_list,wave)
        ret<- lapply(ret, function(m) multiply_matrix(matrix = m$matrix, corrected_mat$M_corr_f ))
      }else{
        warning("WARNINNG NO CONTACT SELECTED: set country_ref = 'Belgium 2020 CoMix (Coletti 2020)' ")
        ret <- NULL
      }
      
      
    }else{
      ret <- lapply(matrix$matrices, function(m) multiply_matrix(matrix = m$matrix, 1 ))
    }
  }
  
  
  
  return(ret)
  #returns a single matrix if no bootstrapping, else a list of matrices 
}

correct_matrix_f <- function(ages, corr_list, wave) {
  # Purpose: Corrects a matrix for specified age groups using correction factors provided for a specific wave.
  
  # Initialize matrices for corrected values and their upper and lower bounds
  M_corr_f <- matrix(0, length(ages), length(ages))
  M_corr_f_up <- matrix(0, length(ages), length(ages))
  M_corr_f_low <- matrix(0, length(ages), length(ages))
  
  # Iterate over the matrix indices to apply corrections
  for (ir in seq_len(length(ages))) {
    for (jc in seq_len(length(ages))) {
      # Determine the group indices based on age bands
      i_row <- 1 * (ir <= 3) + 2 * (ir > 3 & ir < 9) + 3 * (ir >= 9)
      j_col <- 1 * (jc <= 3) + 2 * (jc > 3 & jc < 9) + 3 * (jc >= 9)
      
      # Assign corrected values from the correction list for the given wave
      M_corr_f[ir, jc] <- corr_list[[wave]]$M_perc_f[i_row, j_col]
      M_corr_f_up[ir, jc] <- corr_list[[wave]]$M_perc_f_up[i_row, j_col]
      M_corr_f_low[ir, jc] <- corr_list[[wave]]$M_perc_f_low[i_row, j_col]
    }
  }
  
  # Return the matrices containing the corrected values and their uncertainty bounds
  return(list(M_corr_f = M_corr_f, M_corr_f_up = M_corr_f_up, M_corr_f_low = M_corr_f_low))
}


multiply_matrix <- function(matrix, multiplier) {
  # Purpose: Multiplies each element of a matrix by a specified multiplier.
  matrix * multiplier
}
