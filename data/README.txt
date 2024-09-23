Data explaination:

###Rt (Reproduction number)###

#
File name:"list_avgRt_1_42.rds"
Type: list [3 x 42]
Content: 
each atom [[i]] of the list contains a vector of the estimated Rt (pcr- Gressani et al). The three atoms are named as "Mean", "Low", "Up" in corrispondence with the vector they contain of mean and the extrema of the 95% CI, respectively.
N.B. The Rt with its confidence interval is referring to the mean value of the effective reproduction number over a 7 days interval centered at the survey wave date. 

###Susceptibles###

#
File name:"sm_susc_last_calibration_2024.rds"
Type: array [400 x 730  x  11]
Content: 
It is an array containing [400,,] different simulated number of susceptible individuals, for every day [,1:730,],  for each of the age classes [,,2:11]. Here the considered classes are from 0 to 90+ years of age, with an age interval of 10yrs. (Abrams et al. 2021)


###Contact Matrices####

#
File name:"raw_contact_list.rds"
Type: list , matrix [[9x9], 43]
Content: 
each atom [[k]] of the list contains a 9x9 matrix, whose entry [i,j] contains the average number of daily contacts reported by a participant to the survey of age in i, with an individual of age in j, during Wave k. These matrices are not reciprocal.

N.B. The matrices list is considering already the fatigue effect compensation, exploiting the work of Loedy et al.

#

###Vaccinations & VOCs###

#File name: "df_vaccines_1stDose_1_43.rds"
            "df_vaccines_2ndDose_1_43.rds"
            
Type: list , dataframe 
Content: 
Per-day cumulative vaccinations (1st,2nd dose) for the time interval between waves 1 and 43 of CoMix (Belgium)- (Sciensano).
The dataframe contains both absolute numbers and proportions at each day of the above interval.

#File name: "df_vocs.rds"
Type: list , dataframe 
Content: 
Generic information about main VOCs. 
Initial VOIs(wild type & other variants)
Alpha,Delta, Omicron, Omicron BA.2

Date_first: date of first detection/sequencing (set to 2020-02-04 for Initial VOIs)
Date_start: date of switching from <50% of cases to >50% of cases  (set to 2020-02-04 for Initial VOIs)
Date_end: date of switching from > 50% of cases to <50% of cases
Date80_start: date of prevalence >=80%
Date80_end: date of   50%< prevalence <80% 
(Sciensano, Oxford, KU Leuven)


#File name: "df_measures.rds"
Type: list , dataframe 
Content: 
Evolution of the main interventions took by the Belgian government to counteract the virus' spread. The records report the date of the main changes in policies (Date_start), the date of end of that stringency level (Date_end, it generally coincides with the beginnning of the subsequent period) and a brief description (Description).

#File name: "survey_****.rds"
These are some of the available data on social contact data both from the CoMix survey and from other studies around the word. Check also https://socialcontactdata.org/socrates/ for more details. 