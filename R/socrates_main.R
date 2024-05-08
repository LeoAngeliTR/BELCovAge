#___________________________________________________________________________
# This file is part of the SOcial Contact RATES (SOCRATES) modelling project
# 
# => LOAD AND SELECT SOCIAL CONTACT SURVEY DATA
#
#  Copyright 2020, SIMID, UNIVERSITY OF ANTWERP & HASSELT UNIVERSITY
#___________________________________________________________________________

# load packages and help functions
source('R/load_config_base.R')
source('R/load_config_comix.R')
source('R/contact_matrix_fix.R')



## MAIN FUNCTION ####
get_contact_matrix <- function(country,daytype,touch,duration,gender,
                               cnt_location,cnt_matrix_features,age_breaks_text,
                               weight_threshold,wave){
  
  # parse age intervals
  age_breaks_num <- parse_age_values(age_breaks_text,bool_unique = TRUE)
  
  # if no breaks specified, group all participants
  if(length(age_breaks_num)==0){
    age_breaks_num <- 0
  }
  
  bool_reciprocal      <- opt_matrix_features[[1]]  %in% cnt_matrix_features
  bool_weigh_age       <- opt_matrix_features[[2]]  %in% cnt_matrix_features
  bool_weigh_dayofweek <- opt_matrix_features[[3]]  %in% cnt_matrix_features
  bool_age_range       <- opt_matrix_features[[4]]  %in% cnt_matrix_features
  bool_age_missing     <- opt_matrix_features[[5]]  %in% cnt_matrix_features
  bool_suppl_professional_cnt <- opt_matrix_features[[6]]  %in% cnt_matrix_features
  bool_hhmatrix_selection    <- opt_matrix_features[[7]]  %in% cnt_matrix_features
  
  # get specific social_mixr survey object
  survey_object <- get_survey_object(country      = country,
                                     daytype      = daytype,
                                     touch        = touch,
                                     duration     = duration,
                                     gender       = gender,
                                     cnt_location = cnt_location,
                                     bool_reciprocal   = bool_reciprocal,
                                     bool_suppl_professional_cnt =  bool_suppl_professional_cnt,
                                     bool_hhmatrix_selection = bool_hhmatrix_selection,
                                     wave         = wave)
  
  if(nrow(survey_object$participants)==0){
    return(list(matrix=NA,
                participants = NA,
                warning="Participant selection too strict... no data left!")
    )
  }
  
  if(nrow(survey_object$contacts) == 0){
    return(list(matrix=NA,
                participants = NA,
                warning="Contact selection too strict... no data left!")
    )
  }
  
  # (re)set rng seed (if ages are sampled from the reported range)
  set.seed(rng_seed)
  
  # run social_mixr function
  matrix_out <- contact_matrix(survey          = survey_object, 
                               age.limits      = age_breaks_num,
                               symmetric       = bool_reciprocal,
                               weigh.age       = bool_weigh_age,
                               weigh.dayofweek = bool_weigh_dayofweek,
                               weight.threshold = weight_threshold,
                               estimated.contact.age = ifelse(bool_age_range,'sample','mean'),
                               missing.contact.age = ifelse(bool_age_missing,'remove','ignore'),
                               return.part.weights = TRUE,
                               #return.demography   = TRUE,
                               quiet           = TRUE)
  
  
  # add per capita contact rate (if demography data)
  if('demography' %in% names(matrix_out) && !any(is.na(matrix_out$matrix))){
    num_age_groups <- nrow(matrix_out$demography)
    pop_matrix     <- matrix(rep(matrix_out$demography$population,num_age_groups),ncol=num_age_groups,byrow = T)
    matrix_out$matrix_per_capita <- matrix_out$matrix / pop_matrix
  }
  
  ## TMP: remove weights from output
  if('weights' %in% names(matrix_out)){
    tmp <- matrix_out$weights
    matrix_out$weights <- NULL
    matrix_out$weights <- tmp
    
  }
  
  # ## add date
  dates_str <- paste0(survey_object$participants$year,'-',
                      survey_object$participants$month,'-',
                      survey_object$participants$day)
  dates_str <- dates_str[!grepl('NA',dates_str)]
  
  if(length(dates_str)==0){
    dates_all <- unique(survey_object$participants$year)
    matrix_out$survey_period <- paste(c('Survey period: ', paste(dates_all, collapse=', ')),collapse=' ')
  } else{
    dates_all <- as.Date(dates_str)
    matrix_out$survey_period <- paste(c('From', paste(range(dates_all), collapse=' to ')),collapse=' ')
  }
  

  # return
  matrix_out
}

## GET SURVEY DATA ####
get_survey_object <- function(country,
                              daytype,
                              touch,
                              duration,
                              gender,
                              cnt_location,
                              cnt_profession=NULL,
                              bool_reciprocal,
                              bool_suppl_professional_cnt,
                              bool_hhmatrix_selection,
                              missing.contact.age = "remove",  # adopted from socialmixr package
                              #missing.contact.age = "keep",  # adopted from socialmixr package
                              wave,
                              quiet = FALSE){
  
  # select dataset filename and load #####
  sel_dataset <- opt_country_admin[opt_country_admin$name == country,]
  
  # get original data
  survey_data <- readRDS(sel_dataset$dataset)
  data_part   <- survey_data$participants
  data_cnt    <- survey_data$contacts
  
  # include wave id notation  
  data_part <- add_wave_id(data_part)

  # option to select country-specific participant and contact data
  if(nchar(sel_dataset$country)>0){
    bool_country <- (data_part$country == sel_dataset$country)
    data_part    <- data_part[bool_country,]
    data_cnt     <- data_cnt[data_cnt$part_id %in% data_part$part_id,]
  }
  
  # select type of day ####
  if(!daytype %in% names(opt_day_type[c(1,6)])){
    bool_dayofweek <- data_part$dayofweek >= 0 # all
    if(daytype == opt_day_type[[3]]){ # weekend
      bool_dayofweek <- data_part$dayofweek %in% c(0,6)
      data_part      <- data_part[bool_dayofweek,]
    } else{
      bool_dayofweek <- data_part$dayofweek %in% 1:5
      data_part      <- data_part[bool_dayofweek,]
    }
  }
  
  # select period ####
  if(daytype %in% names(opt_day_type[4:6])){
    if(!any(data_part$holiday)){
      load('data/holiday_all.RData')
      country_iso3 <- countrycode(unlist(country), 'country.name', 'iso3c')
      country_holiday_data <- holiday_all[holiday_all$iso3 == country_iso3,]
      data_part$date <- as.Date(paste(data_part$day,
                                      data_part$month,
                                      data_part$year,
                                      sep='/')
                                ,'%d/%m/%Y')
      data_part$holiday <- data_part$date %in% country_holiday_data$date
    }
    
    # check if holiday is a boolean... if not, try to convert and throw warning
    if(typeof(data_part$holiday) != 'logical'){
      warning("holiday variable is not a 'logical', try to convert binary 0/1 to FALSE/TRUE")
      data_part$holiday <- data_part$holiday == 1
    }
    
    if(daytype == opt_day_type[[4]]){
      # if(!any(data_part$is_holiday)){ # if no holiday period data
      #   print("NO HOLIDAY DATA... USE REGULAR PERIOD DATA")
      # } else{ # select holiday period data
      data_part <- data_part[data_part$holiday,]
      # }
    } else{ # select regular period data
      data_part <- data_part[!data_part$holiday,]
    }
    # select cnt data of remaining participants
    data_cnt <- data_cnt[data_cnt$part_id %in% data_part$part_id]
  }
  
  # select contact duration ####
  if(duration != opt_duration[[1]]){
    #print(duration)
    duration_code <- which(opt_duration == duration)-1
    
    if(duration %in% names(opt_duration[2:3]) ){
      bool_duration <- !is.na(data_cnt$duration_multi)  & data_cnt$duration_multi <= duration_code
      data_cnt      <- data_cnt[bool_duration,]
    } else{
      bool_duration <- !is.na(data_cnt$duration_multi)  & data_cnt$duration_multi >= duration_code
      data_cnt      <- data_cnt[bool_duration,]
    }
  }
  
  # select contact intensity ####
  if(touch != opt_touch[[1]]){
    touch_code    <- which(opt_touch == touch)-1
    bool_touching <- !is.na(data_cnt$phys_contact) & data_cnt$phys_contact == touch_code
    data_cnt      <- data_cnt[bool_touching,]
    #print(touch)
  }
  
  # select gender ####
  if(gender != opt_gender[[1]]){
    
    # first select cnt data of remaining participants
    data_cnt <- data_cnt[data_cnt$part_id %in% data_part$part_id]
    
    # make sure the gender variable is used as character
    data_cnt$cnt_gender   <- as.character(data_cnt$cnt_gender)
    data_part$part_gender <- as.character(data_part$part_gender)
    
    # set gender-specific booleans
    bool_cnt_female  <- data_cnt$cnt_gender   == 'F'
    bool_part_female <- data_part$part_gender == 'F'
    bool_cnt_male    <- data_cnt$cnt_gender   == 'M'
    bool_part_male   <- data_part$part_gender == 'M'
    
    # merge dataset to compare participant and contact gender
    data_cnt_gender  <- merge(data_cnt[,c('part_id','cnt_gender')],data_part[,c('part_id','part_gender')],by='part_id')
    data_cnt_gender$cnt_gender[!data_cnt_gender$cnt_gender %in% c('M','F')] <- NA
    data_cnt_gender$part_gender[!data_cnt_gender$part_gender %in% c('M','F')] <- NA
    bool_gender_diff <- data_cnt_gender$cnt_gender != data_cnt_gender$part_gender
    bool_gender_diff <- bool_gender_diff & !is.na(bool_gender_diff)
    
    if(gender == opt_gender[[2]]){                  # female-female
      data_cnt       <- data_cnt[bool_cnt_female,]
      data_part      <- data_part[bool_part_female,]
    } else if(gender == opt_gender[[5]]){           # male-male
      data_cnt       <- data_cnt[bool_cnt_male,]
      data_part      <- data_part[bool_part_male,]
    } else if(bool_reciprocal){
      data_cnt       <- data_cnt[bool_gender_diff,]
    } else {
      if(gender == opt_gender[[3]]){                # female-male
        data_cnt       <- data_cnt[bool_cnt_male,]
        data_part      <- data_part[bool_part_female,]
      } else if(gender == opt_gender[[4]]){         # male-female
        data_cnt       <- data_cnt[bool_cnt_female,]
        data_part      <- data_part[bool_part_male,]
      }
    }
  }
  
  ## select wave (optional) ----
  if(wave != opt_waves[[1]]){
    if(!is.null(data_part$wave) & wave %in% data_part$wave){
      # print(table(data_part$wave))
      # print(table(data_part$wave == wave))
      
      bool_part_wave <- data_part$wave == wave
      data_part <- data_part[bool_part_wave,]
      
      bool_cnt_wave <- data_cnt$part_id %in%  data_part$part_id
      data_cnt      <- data_cnt[bool_cnt_wave,]
      # print(paste('select wave', wave))
      # print(table(data_part$wave))
    }
  }
  
  #__________________________________________________________________________
  # adjust location data: missing and multiple locations ####
  if(nrow(data_cnt)>0){
    # set data.table to data.frame
    data_cnt_tmp <- data.frame(data_cnt)
    
    # select all location-specific columns
    cnt_location_colnames <- c(paste0('cnt_',tolower(opt_location)))
    data_cnt_tmp <- data_cnt_tmp[,cnt_location_colnames]
    dim(data_cnt_tmp)
    
    # replace value 'NA' for a location to 'false' (=not-present)
    data_cnt_tmp[is.na(data_cnt_tmp)] <- 0
    
    # add missing location to "other"
    # note: missing could also be "other locations than specified in opt_location"
    cnt_loc_missing <- rowSums(data_cnt_tmp,na.rm=T) == 0
    data_cnt_tmp$cnt_otherplace  <- as.numeric(data_cnt_tmp$cnt_otherplace | cnt_loc_missing)

    # 1. calculate cumulative sum (from left to right)
    tmp_loc_cumsum <- t(apply(data_cnt_tmp,1,cumsum))
    
    # 2. set locations with cummulative sum >1 (== not unique and not the "main location") to 0
    data_cnt_tmp[tmp_loc_cumsum>1] <- 0
    
    # 3. copy adjusted location data back
    data_cnt[,cnt_location_colnames] <- data_cnt_tmp
  }

  #__________________________________________________________________________
  
  # household members: get matrix with only household members? ####
  if('is_hh_member' %in% names(data_cnt) && !is.na(bool_hhmatrix_selection) && 
     bool_hhmatrix_selection == TRUE){
    flag_cnt_adapt                       <- data_cnt$cnt_home == 1 & data_cnt$is_hh_member == FALSE
    data_cnt$cnt_home[flag_cnt_adapt]    <- 0
    data_cnt$cnt_leisure[flag_cnt_adapt] <- 1
  }
  
  #select location ####
  if(length(cnt_location)==0){
    warning("WARNING: NO LOCATIONS SPECIFIED...")
    data_cnt <- data_cnt[0,]
  } else if(!identical(as.character(cnt_location),as.character(opt_location))
            && nrow(data_cnt)>0){
    
    # set data.table to data.frame
    data_cnt_tmp <- data.frame(data_cnt)
    
    # select requested location-specific columns
    cnt_location_colnames <- c(paste0('cnt_',tolower(cnt_location)))
    
    # select columns
    if(length(cnt_location)>1){
      is_present <- rowSums(data_cnt_tmp[,cnt_location_colnames] == 1,na.rm=T)
    } else{
      is_present <- data_cnt_tmp[,cnt_location_colnames]
    }
    
    # select contact at specified location(s) 
    bool_location <- is_present > 0
    
    # select
    data_cnt <- data_cnt[bool_location,]
    
    # add warning
    if(!any(bool_location)){
      warning("WARNING: NO CONTACTS LEFT AFTER LOCATION SELECTION...")
    } 
  }
  
  
  # remove temporal data from contact data.frame
  if('dayofweek' %in% names(data_cnt)){
    data_cnt$dayofweek <- NULL
  }
 
  # suppl. professional contact? ####
  # ==>> remove imputed supplementary professional contacts?
  if('is_imputed' %in% names(data_cnt) && !is.na(bool_suppl_professional_cnt) && 
     bool_suppl_professional_cnt == FALSE){
    data_cnt <- data_cnt[data_cnt$is_imputed == 0,]
  }
  
  
  
  
  
  # #new column: contact intensity BE####
  # Collapsing redundant working categories
  if(country=="Belgium 2020 CoMix (Coletti 2020)"){
    data_part <- data_part %>%
      mutate(part_occupation = case_when(
        part_occupation %in% c("member of the general management, senior executive responsible for 5 employees or less",
                               "member of the general management senior executive responsible for 5 employees or less") ~
          "General management, senior executive (responsible <=5 employees)",
        
        part_occupation %in% c("other employee who does not do office work (eg teacher, nurses ...)",
                               "other employee who does not do office work (eg teacher nurses ...)") ~
          "Other non-office worker (e.g., teacher, nurse)",
        
        part_occupation %in% c("craftsman, trader with 5 employees or less",
                               "craftsman trader with 5 employees or less") ~
          "Craftsman, trader (<=5 employees)",
        
        part_occupation %in% c("member of the general management, senior executive responsible for 6 to 10 employees",
                               "member of the general management senior executive responsible for 6 to 10 employees") ~
          "General management, senior executive (6-10 employees)",
        
        part_occupation %in% c("industrial, wholesaler with 6 employees or more",
                               "industrial wholesaler with 6 employees or more") ~
          "Industrial, wholesaler (>6 employees)",
        
        part_occupation %in% c("member of the general management, senior executive responsible for 11 employees or more",
                               "member of the general management senior executive responsible for 11 employees or more") ~
          "General management, senior executive (>10 employees)",
        part_occupation %in% c("middle management that is not part of the general management responsible for 5 employees or less",
                               "middle management, that is not part of the general management, responsible for 5 employees or less") ~
          "Middle management, not part of the general management (<=5 employees)",
        part_occupation %in% c("middle management, that is not part of the general management, responsible for 6 employees or more",
                               "middle management that is not part of the general management responsible for 6 employees or more") ~
          "Middle management, not part of the general management (>6 employees)",
        
        TRUE ~ part_occupation
      ))
    
    data_part$part_occupation <- as.factor(data_part$part_occupation)
    levels(data_part$part_occupation) <- c(levels(data_part$part_occupation),"I do not know")
    occupations <- c(
      "Craftsman, trader (<=5 employees)",                                    
      "farmer",                                                               
      "General management, senior executive (>10 employees)" ,                
      "General management, senior executive (6-10 employees)" ,               
      "General management, senior executive (responsible <=5 employees)" ,    
      "house man or housewife",                                               
      "I do not know"  ,                                                      
      "Industrial, wholesaler (>6 employees)",                                
      "liberal profession or profession for which qualification is required" ,
      "Middle management, not part of the general management (<=5 employees)",
      "Middle management, not part of the general management (>6 employees)" ,
      "never worked",                                                         
      "non-skilled worker" ,                                                  
      "other" ,                                                               
      "other employee who mainly performs office work" ,                      
      "Other non-office worker (e.g., teacher, nurse)" ,                      
      "pre-retired",                                                          
      "retired" ,                                                             
      "skilled worker" ,                                                      
      "student" ,                                                             
      "unable for work" ,                                                     
      "unemployed")
    # occupations <- unique(data_part$part_occupation)
    data_part <- data_part %>%
      mutate(part_occupation_contact = case_when(
        part_occupation %in% occupations[c(13,9,16)]~ "High", #conact intensity classification
        
        part_occupation %in% occupations[c(1,2,10,11,15)] ~ "Medium",
        
        part_occupation %in% occupations[c(3,4,5,8,19)]~ "Low",
        
        #part_occupation %in% occupations[c(6,7,12,14,17,18,20,21,22)] | is.na(part_occupation)~ "Other",
        part_occupation %in% occupations[c(6,7,12,14,17,18,20,21,22)] ~ "Other",
        
        TRUE ~ NA_character_
      ))
    
    data_part$part_occupation_contact <- as.factor(data_part$part_occupation_contact)
    
    
    if(!is.null(cnt_profession) && length(cnt_location)==1){
      if(!cnt_profession %in% opt_profession ){
        warning("WARNING NO CONTACT SELECTED: select correct cnt_profession from opt_profession elements")
        data_cnt <- data_cnt[0,]
      }else if(cnt_location!='Work'){
        warning("WARNING SELECTION cnt_profession IGNORED:\n
      the selection of work contact intensity is only valid when cnt_location== 'Work' ")
      }else{
        # set data.table to data.frame
        data_part_tmp <- data.frame(data_part)
        
        #select ID partecipant related to professions "cnt_profession"
        id_profession <- filter(data_part_tmp, part_occupation_contact==cnt_profession)$part_id
        
        #filter contact data
        data_cnt <- data_cnt %>%
          filter(part_id  %in% id_profession)  
        
        #fill the column with selected work-contact intensity
        data_cnt <- data_cnt %>%
          mutate(cnt_high = 0,               
                 cnt_medium = 0,               
                 cnt_low = 0,
                 cnt_other=0) 
        # select requested profession-specific columns
        cnt_profession_colname <- c(paste0('cnt_',tolower(cnt_profession)))
        data_cnt[,cnt_profession_colname] <- data_cnt$cnt_work
      }
    }else if (!is.null(cnt_profession) && length(cnt_location)>1){
      warning("WARNING selection cnt_profession ignored:\n
    the selection of work contact intensity is only valid when cnt_location== 'Work' ")
      
    }
    
  }
  
  # create new survey object
  mixr_survey <- survey(data_part, data_cnt)
  
  # return
  return(mixr_survey)
  
}


#mija <- contact_matrix(polymod, countries = "Belgium", age.limits = c(0, 18, 45,65))$matrix*c(1,0.5,0.6,1)
#mijb <- contact_matrix(polymod, countries = "Belgium", age.limits = c(0, 18, 45, 65))$matrix
compare_contact_matrices <- function(mija,mijb,
                                     bool_transmission_param,age_susceptibility_text,age_infectiousness_text){
  
  # mij ratio
  mij_ratio     <- mija/mijb
  
  # adjust for age-specific transmission?
  if(bool_transmission_param){
    mija <- adjust_mij_transmission(mija,age_susceptibility_text,age_infectiousness_text)
    mijb <- adjust_mij_transmission(mijb,age_susceptibility_text,age_infectiousness_text)
  }
  
  if(any(is.na(mija))|any(is.na(mijb))){
    warning('Social contact matrix contains NA... no comparison possible!')
    out <- list(notes='Social contact matrix contains NA... no comparison possible!')
  } else{
    R0_ratio      <- max(Re(eigen(mija)$values))/max(Re(eigen(mijb)$values))
    
    # relative incidence 
    RIa             <- standardize_RI(eigen(mija)$vectors[,1])
    RIb             <- standardize_RI(eigen(mijb)$vectors[,1])
    RI_ratio        <- RIa/RIb
    names(RI_ratio) <- colnames(mija)
    
    #output 
    out <- list(R0_ratio=R0_ratio,mij_ratio=mij_ratio,RI_ratio=RI_ratio,
                notes="ratio = with intervention / without intervention")
    
    # fix NA-results
    if(identical(mija,mijb)){ # set 1 if mija == mijb
      for(i in 1:length(out)) { 
        out[[i]][] <- 1 
      }
    } else if(sum(mija) == 0){ # set 0 if mija[] == 0
      for(i in 1:length(out)) { 
        out[[i]][] <- 0 
      }
    }
  }
  
  return(out)
}



get_location_matrices <- function(country,daytype,touch,duration,gender,
                                  cnt_location,
                                  cnt_matrix_features,
                                  age_breaks_text,
                                  weight_threshold,
                                  wave){
  
  
  # location specific ==> NOT reciprocal
  sel_cnt_matrix_features <- cnt_matrix_features[!grepl('recipocal',cnt_matrix_features,ignore.case = T)]
  
  # initialise list
  matrix_list <- list()
  
  for(i_loc in 1:length(cnt_location)){
    matrix_list[i_loc] <- list(get_contact_matrix(country,daytype,touch,duration,gender,
                                                  cnt_location = cnt_location[i_loc],
                                                  sel_cnt_matrix_features,
                                                  age_breaks_text,
                                                  weight_threshold,
                                                  wave = wave))
  }
  
  # add location names
  names(matrix_list) <- cnt_location
  
  return(matrix_list)
}

