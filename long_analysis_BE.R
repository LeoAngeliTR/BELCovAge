#___________________________________________________________________________
# This script generates the graphical outputs for the paper by Angeli et al. 2024,
# focusing on the longitudinal sensitivity analysis in Belgium.
# Copyright 2020, SIMID, UNIVERSITY OF ANTWERP & HASSELT UNIVERSITY
#___________________________________________________________________________

library(rstudioapi)
library(matrixcalc)
library(popdemo)
library(demography)
library(smoothAPC)
library('plot.matrix')
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
library (readr)
library(scales)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(qs)
library(janitor)
library(showtext)
library(MetBrewer)
library(scico)
library(ggtext)
library(patchwork)
library(gghighlight)
library(ggrepel)
library(sysfonts)


# Load preliminary operations
suppressPackageStartupMessages(source('./R/socrates_main.R'))  # Main function source for modeling
source('./add_functions.R')  # Additional function definitions
source('./wave_specific_an.R')  # Wave-specific analysis functions

# Data Loading
data_eu <- read.csv('https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_nat_latest.csv')
stringency <- data_eu[data_eu$CountryName == 'Belgium',] %>% 
  filter(!is.na(StringencyIndex_Average))
stringency$date <- as.Date(as.character(stringency$Date), format = "%Y%m%d")

#Load NPIs
output='./data/'
df_measures <- readRDS(paste0(output,"df_measures.rds"))

# Load genomic surveillance data and vaccination data from local sources
df_vocs <- readRDS("./data/df_vocs.rds")  # Genomic surveillance data
df_vaccines_1st = readRDS("./data/df_vaccines_1stDose_1_43.rds")  # First dose vaccine data
df_vaccines_2nd = readRDS("./data/df_vaccines_2ndDose_1_43.rds")  # Second dose vaccine data

#### Styling for plots ####
fam <- font <- "Gudea"  # Define the font for plots
font_add_google(family = font, font, db_cache = TRUE)  # Add Google font
fa_path <- systemfonts::font_info(family = "Font Awesome 6 Brands")[["path"]]
font_add(family = "fa-brands", regular = fa_path)
theme_set(theme_minimal(base_family = font, base_size = 20))  # Set the default theme for ggplot2
bg <- "white"  # Background color
txt_col <- "black"  # Text color
showtext_auto(enable = TRUE)  # Enable auto-text rendering
col_blind_pal <- c("#FF8C00", "#D62728", "#F0E442", "#117733", "#44AA99", "#0072B2", "#56B4E9", "#CC6677", "#882255")  # Color-blind-friendly palette
pal <- c("#FFFFF0", "#F0EAD6", "#FFFDD0", "#F3E5AB")  # Custom color palette


# Prepare data for plotting and analysis by calculating matrix differences between consecutive NGMs
wave8_NGM <- readRDS("./data/wave8_NGM.rds")
matrix_diff <- list(ngms[[1]] - wave8_NGM)  # Initialize with difference from wave 8
matrix_diff <- c(matrix_diff, Map(`-`, ngms[-1], ngms[-length(ngms)]))  # Calculate differences between successive NGMs


pert_df <- data.frame(Wave=integer(),
                      Date=Date(),
                      Age=character(),
                      K_col=double(),
                      K_rows=double(),
                      Sens=double(),
                      Sens_standard=double(),
                      Sens_rows=double(),
                      El=double(),
                      Rep_val=double(),
                      Stable_stage=double(),
                      Cont=double(),
                      Cont_perc=double(),
                      Cont_cols=double(),
                      Tot_S=double(),
                      Tot_S_delta=double(),
                      Tot_S_delta_mid=double(),
                      
                      stringsAsFactors=FALSE)
for (i in 1:length(ngms)) {
  for (k in 1:length(age_breaks)) {
    
    new_record <- list(as.integer(8+i),
                       w_dates[8+i],
                       age_names[k],
                       colSums(ngms[[i]])[k],
                       rowSums(ngms[[i]])[k],
                       colSums(calc.sens(ngms[[i]]))[k],
                       colSums(calc.sens(ngms[[i]]))[k]/sum(calc.sens(ngms[[i]])),
                       rowSums(calc.sens(ngms[[i]]))[k],
                       colSums(calc.el(ngms[[i]]))[k],
                       v1_list[[i]][k],
                       w1_list[[i]][k],
                       rowSums(c_asym[[i]])[k],
                       rowSums(c_asym[[i]])[k]/sum(c_asym[[i]]),
                       colSums(c_asym[[i]])[k],
                       S_tot_index[[i]],
                       S_tot_index[[i]]*norm(matrix_diff[[i]],type = 'F'),
                       S_tot_index_mid[[i]]*norm(matrix_diff[[i]],type = 'F')
                       
    )
    pert_df[nrow(pert_df) + 1,] = new_record 
      }
}
# Add a column to classify age groups based on their names
pert_df <- pert_df %>%
  mutate(AgeGroup = case_when(
    Age %in% age_names[1:3] ~ 'Minors',
    Age %in% age_names[4:7] ~ 'Adults',
    Age %in% age_names[8:9] ~ 'Elderly',
    TRUE ~ 'Other'  # For any age not covered in the above categories
  ))


#R_t dataframe

rep_df<- data.frame(Wave=integer(),
                    Date=Date(),
                    R0= double(),
                    R0_up= double(),
                    R0_low= double(),
                    stringsAsFactors=FALSE)
for (i in 9:length(w_dates)) {
  new_record <- list(as.integer(i),
                     w_dates[i],
                     Rt_list$Mean[i],
                     Rt_list$Up[i],
                     Rt_list$Low[i]
                     
  )
  
  rep_df[nrow(rep_df) + 1,] = new_record 
}



mean_df<- data.frame(Wave=integer(),
                     Date=Date(),
                     Cont=double(),
                     El=double(),
                     El_w=double(),
                     Sens_cols=double(),
                     Sens_cols_w=double(),
                     Sens_rows=double(),
                     Sens_rows_w=double(),
                     Rep_val_w=double(),
                     Stable_w=double(),
                     
                     stringsAsFactors=FALSE)
for (i in 1:length(ngms)) {
  new_record <- list(as.integer(8+i),
                     w_dates[8+i],
                     weighted.mean(filter(pert_df,Wave==8+i)$Cont,N_vec),
                     mean(filter(pert_df,Wave==8+i)$El),
                     weighted.mean(filter(pert_df,Wave==8+i)$El, N_vec),
                     mean(filter(pert_df,Wave==8+i)$Sens),
                     weighted.mean(filter(pert_df,Wave==8+i)$Sens,N_vec),
                     mean(filter(pert_df,Wave==8+i)$Sens_rows),
                     weighted.mean(filter(pert_df,Wave==8+i)$Sens_rows,N_vec),
                     weighted.mean(filter(pert_df,Wave==8+i)$Rep_val,N_vec),
                     weighted.mean(filter(pert_df,Wave==8+i)$Stable_stage,N_vec)
  )
  mean_df[nrow(mean_df) + 1,] = new_record 
  
}

## Data frame for susceptibility data across age groups and waves
susc_df<- data.frame(Wave=integer(),
                     Age=character(),
                     Susc=double(),
                     Susc_prop=double(),
                     Susc_up=double(),
                     Susc_prop_up=double(),
                     Susc_low=double(),
                     Susc_prop_low=double(),
                     Date=Date(),
                     stringsAsFactors=FALSE)
for (i in 1:length(w_dates)) {
  for (k in 1:length(age_breaks)) {
    new_record <- list(as.integer(i),
                       age_names[k],
                       SW[k,i],
                       SW[k,i]/N_vec[k],
                       SW_up[k,i],
                       SW_up[k,i]/N_vec[k],
                       SW_low[k,i],
                       SW_low[k,i]/N_vec[k],
                       w_dates[i])
    susc_df[nrow(susc_df) + 1,] = new_record
  }
}





#ELASTICITY EVOLUTION----

#BACKGROUND PLOT
xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
annotations_size <- 31.5
size_text <- 98
plot.margin = margin(t = 0.3, r = 0.5, b = 0.3, l = 0.5, unit = "cm")
base_size = 90
legend.key.size = unit(2,'cm')
legend.key.height = unit(1.5,'cm')
legend.key.width = unit(4,'cm')
y_max <-0.91
gb <-ggplot()
base<- gb+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  #annotate("text", x=w_dates[14], y=max(rep_df$R0*scaling_Rt), label="Schools & non-ess. shops reopen",size=annotations_size, angle=0)+
  annotate("text", x=w_dates[13]+5, y=y_max, vjust=2, label="Schools & non-ess. shops reopen",size=annotations_size, angle=0)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[8]+15,y=y_max, vjust=2, label=  df_measures$Descr_short[8],size=annotations_size, angle=0)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=w_dates[29],y=y_max, vjust=2, label= "School re-opening & Lifting",size=annotations_size, angle=0)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[10]+18,y=y_max, vjust=2, label= df_measures$Descr_short[10] ,size=annotations_size, angle=0)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[11]+20,y=y_max, vjust=2, label= "Lifting" ,size=annotations_size, angle=0)


#Stringency index
#annotate("text", x=w_dates[24], y=0.53, label="Stringency index(%)",size=6, angle=0)+

#VOCs
y_rect <- -0.02
y_min <- -0.06
y_text <- -0.04
base <- base+
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = y_min, ymax =y_rect)+
  annotate("text", x=w_dates[12], y=y_text, label= "Wild Type",size=annotations_size)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=w_dates[20], y=y_text, label= "Alpha",size=annotations_size)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=y_text, label= "Delta",size=annotations_size)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=y_text, label= "Omicron",size=annotations_size)+
  theme(
    plot.margin=plot.margin
  )







scaling_Rt <- max(pert_df$El)*1.5
pert_df$Age <- factor(pert_df$Age,levels =age_names)

##Minors----
#size_text <- 100
minor_plot <- base+geom_line(data =pert_df %>% filter(AgeGroup == 'Minors'),aes(x=Date,y=El,color=Age), lwd=3.4)+
  geom_line(data = rep_df,aes(x=Date,y=R0*scaling_Rt,linetype="Rt"), lwd=1.8) +
  geom_line(data =stringency %>% filter(date>=w_dates[9] & date<=w_dates[length(w_dates)]),aes(x=date,y=StringencyIndex_Average/102, linetype="S.I.(%)"), lwd=1.8, color='black')+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 0.11, lwd=1.4, linetype="longdash")+
  
  scale_x_date( breaks=w_dates[9:length(w_dates)],
                labels = rep('', length(w_dates[9:length(w_dates)])),
                #labels = xaxis_el, 
                limits=c(w_dates[9],w_dates[length(w_dates)]),
                expand = c(0, 0),
                sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave")
  )+
  scale_linetype_manual(values=c('solid', "twodash"))+
  scale_y_continuous(expand = c(0, 0),sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names)+
  theme_bw(base_size = base_size) +
  theme(plot.margin=plot.margin,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.spacing.y = unit(1, 'cm'),
        axis.text.x.bottom = element_text(size=size_text, angle=90, vjust=0.5, family =fam),
        axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
        axis.text.y  = element_text(size=size_text, family =fam),
        legend.title = element_text(size=size_text, family =fam),
        legend.text  = element_text(size=size_text, family =fam),
        legend.key.size = legend.key.size,
        legend.key.height = legend.key.height,
        legend.key.width = legend.key.width,
        axis.title.y.left = element_text(size=size_text*1.1, family =fam),
        axis.title.x.top = element_text(size=size_text, family =fam),
        axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  guides(
    linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
    color = guide_legend(order = 2)     # 'color' represents the Age group legend
  ) +
  labs(y=expression(tilde(e)[j]),
       x="",linetype='',
       color='Age group')

ggsave("Minors.png", plot = minor_plot, width = 40, height = 25, dpi = 200)

##Adults----

adults_plot <- base+geom_line(data =pert_df %>% filter(AgeGroup == 'Adults'),aes(x=Date,y=El,color=Age), lwd=3.4)+
  geom_line(data = rep_df,aes(x=Date,y=R0*scaling_Rt), lwd=1.8)+
  geom_line(data =stringency %>% filter(date>=w_dates[9] & date<=w_dates[length(w_dates)]),aes(x=date,y=StringencyIndex_Average/102),linetype="twodash", lwd=1.8, color='black')+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 0.11, lwd=1.4, linetype="longdash")+
  
  scale_x_date(breaks=w_dates[9:length(w_dates)],
               labels = rep('', length(w_dates[9:length(w_dates)])),
               limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0,0),
               sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates), name = "")
  )+
  scale_y_continuous( expand = c(0,0), sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal[4:7],
                     labels=age_names[4:7])+
  theme_bw(base_size = base_size) +
  theme(
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x.bottom = element_text(size=size_text, angle=90,vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text, family =fam),
    legend.text  = element_text(size=size_text, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.1, family =fam),
    axis.title.x.top = element_text(size=size_text, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  labs(y=expression(tilde(e)[j]),
       x="",linetype='',
       color='')

ggsave("Adults.png", plot = adults_plot, width = 40, height = 25, dpi = 200)

##Elderly----

elderly_plot <- base+geom_line(data =pert_df %>% filter(Wave > 8 & AgeGroup == 'Elderly'),aes(x=Date,y=El,color=Age), lwd=3.4)+
  
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt), lwd=1.8)+
  geom_line(data =stringency %>% filter(date>=w_dates[9] & date<=w_dates[length(w_dates)]),aes(x=date,y=StringencyIndex_Average/102), linetype="twodash", lwd=1.8, color='black')+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 0.11, lwd=1.4, linetype="longdash")+
  geom_segment(x = w_dates[13], xend = w_dates[13], y = -Inf, yend = 0.5, linetype = 2, col = 'grey20',lwd=2.7) +
  
  scale_x_date(breaks=w_dates[9:length(w_dates)],
               labels = xaxis_el,
               limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0,0),
               sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "")
  )+
  scale_y_continuous( expand = c(0,0), sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal[8:9],
                     labels=age_names[8:9])+
  theme_bw(base_size = base_size) +
  theme(
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.x.bottom = element_text(size=size_text, angle=90, vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text, family =fam),
    axis.title.y.left = element_text(size=size_text*1.1, family =fam),
    legend.text  = element_text(size=size_text, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.x.top = element_text(size=size_text, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  labs(y=expression(tilde(e)[j]),
       x="",linetype='',
       color='')

ggsave("Elderly.png", plot = elderly_plot, width = 40, height = 25, dpi = 200)

# Combine the plots
combined_plot <- ggarrange(minor_plot, adults_plot, elderly_plot, ncol = 1, nrow = 3, legend = "right")

# Save the combined plot: this generates Figure 4 in the main text of the manuscript
ggsave( "AgeGroupPlot.png", plot = combined_plot, width = 58, height = 58, dpi =178,limitsize = FALSE)

###To generate Figures S4 to S8 of the SM of the manuscript, it is sufficient to select the correspondent attributes of the data frame "pert_df" in the code above.
###namely: K_col,K_rows, Cont, Sens, Rep_val


#SUSCEPTIBLES evolution----

#fam <- 'sans'
xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
susc_df$Age <- factor(susc_df$Age,levels =age_names)
y_ann <- 1
ann_size <- 20
gb <-ggplot()+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=w_dates[14]-7, y=y_ann, label="Schools & non-ess. shops reopen",size=ann_size, angle=0)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[8]+18, y=y_ann, label=  df_measures$Descr_short[8],size=ann_size, angle=0)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=w_dates[29], y=y_ann, label= "School re-opening + gen. lifting",size=ann_size, angle=0)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[10]+20, y=y_ann, label= df_measures$Descr_short[10] ,size=ann_size, angle=0)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[11]+25, y=y_ann, label= "Lifting" ,size=ann_size, angle=0)+
  
  
  #VOCs
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=w_dates[12], y=-0.06, label= "Wild Type",size=ann_size)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=w_dates[20], y=-0.06, label= "Alpha",size=ann_size)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=-0.06, label= "Delta",size=ann_size)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=-0.06, label= "Omicron",size=ann_size)
  
  
 #Susc proportion:FIG 6 bottom----
text_size <- 55
 #Susceptibility
  susc_prop <- gb+
  geom_line(data =susc_df %>% filter(Wave>8),aes(x=Date,y=Susc_prop,color=Age), lwd=1.6)+
  geom_line(data=df_vaccines_1st,aes(x=Day,y=Percnt,linetype="1st dose"),lwd=0.6)+
  geom_line(data=df_vaccines_2nd,aes(x=Day,y=Percnt, linetype="2nd dose"),lwd=0.6)+
  
  
  scale_x_date(breaks=w_dates[9:length(w_dates)],labels = xaxis_el,   limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0.0, 0),sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "Wave"))+
  scale_color_manual("Age",values=col_blind_pal, 
                     labels=age_names)+
  scale_y_continuous(sec.axis=sec_axis(~./1,breaks = seq(0,1, by=0.25)))+
  
  theme_bw(base_size = 20) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x.bottom = element_text(size=text_size, angle=90, vjust=0.5, family =fam), 
          axis.text.y  = element_text(size=text_size, family =fam),
          axis.text.y.right = element_text(color = 'darkred'),
          legend.text  = element_text(size=text_size, family =fam),
          legend.key.height = unit(1.5, "cm"),
          legend.key.size = unit(0.7,'cm'),
          legend.key.width = unit(2,'cm'),
          legend.spacing.y = unit(1,'cm'),
          legend.title = element_text(size=text_size, family =fam),
          axis.title = element_text(size=text_size, family =fam),
          axis.text.x.top = element_text(size=text_size, vjust=0.5,angle=30, family =fam),
          axis.text.y.left = element_text(size=text_size, vjust=0.5,angle=90, family =fam)
    ) +
    guides(
      linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
      color = guide_legend(order = 2)     # 'color' represents the Age group legend
    ) +
  labs(y=bquote("Susceptibles proportion"),
       x="",
       color='Age', linetype="")

  # Save the combined plot
  ggsave("Susc_prop.png", plot = susc_prop, width = 25, height = 13, dpi =200,limitsize = FALSE)
  

#Susc Number: FIG 6 top----
text_size <- 55
ann_size <- 20
y_text <- max(filter(susc_df, Wave >8)$Susc)*1.05
vocs_y <- 5*(10^4)
sec_axis_scale <- 1.8*10^6
gb <-ggplot()
susc_numb <- gb+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=w_dates[14]-7, y=y_text, label="Schools & non-ess. shops reopen",size=ann_size, angle=0)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[8]+18, y=y_text, label=  df_measures$Descr_short[8],size=ann_size, angle=0)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=w_dates[29], y=y_text, label= "School re-opening + gen. lifting",size=ann_size, angle=0)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[10]+20, y=y_text, label= df_measures$Descr_short[10] ,size=ann_size, angle=0)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[11]+25, y=y_text, label= "Lifting" ,size=ann_size, angle=0)+
  
  
  #VOCs
  
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=w_dates[12], y=-vocs_y, label= "Wild Type",size=ann_size)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=w_dates[20], y=-vocs_y, label= "Alpha",size=ann_size)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=-vocs_y, label= "Delta",size=ann_size)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=-vocs_y, label= "Omicron",size=ann_size)+
  
  
  #Susceptibility
  geom_line(data =susc_df %>% filter(Wave>8),aes(x=Date,y=Susc,color=Age), lwd=1.6)+
  geom_line(data=df_vaccines_1st %>% filter(Day>=w_dates[9]),aes(x=Day,y=Percnt*(sec_axis_scale),linetype="1st dose"),lwd=0.6)+
  geom_line(data=df_vaccines_2nd %>% filter(Day>=w_dates[9]),aes(x=Day,y=Percnt*(sec_axis_scale), linetype="2nd dose"),lwd=0.6)+
  
  scale_color_manual("Age Group",values=col_blind_pal, 
                     labels=age_names)+
  scale_x_date(breaks=w_dates[9:length(w_dates)],labels = xaxis_el,   limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0.0, 0),sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave"))+
  scale_y_continuous(sec.axis=sec_axis(~./sec_axis_scale,breaks = seq(0,1, by=0.25)),
  breaks = c(500000, 10^6, 1.5*10^6), # Choose appropriate breaks for your data
  labels = c("0.5*10^6", "10^6", "1.5*10^6"))+
  theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x.bottom = element_text(size=text_size, angle=90, vjust=0.5, family =fam), 
        axis.text.y  = element_text(size=text_size, family =fam),
        axis.text.y.right = element_text(color = 'darkred'),
        legend.text  = element_text(size=text_size, family =fam),
        legend.key.height = unit(1.5, "cm"),
        legend.key.size = unit(0.7,'cm'),
        legend.key.width = unit(2,'cm'),
        legend.spacing.y = unit(1,'cm'),
        legend.title = element_text(size=text_size, family =fam),
        axis.title = element_text(size=text_size, family =fam),
        axis.text.x.top = element_text(size=text_size, vjust=0.5,angle=30, family =fam),
        axis.text.y.left = element_text(size=text_size, vjust=0.5,angle=90, family =fam)
  ) +
  labs(y=bquote("# susceptible individuals"),
       x="",
       color='Age Group', linetype="")



# Save the combined plot
ggsave( "Susc_number.png", plot = susc_numb, width = 25, height = 13, dpi =200,limitsize = FALSE)


# Combine the plots
combined_plot <- ggarrange(susc_numb, susc_prop,  ncol = 1, nrow = 2, legend = "right", common.legend = T)

# Save the combined plot
ggsave("AgeGroupPlotSusc.png", plot = combined_plot, width = 25, height = 20, dpi =250,limitsize = FALSE)


#SCATTERPLOTS FIGs 1,2,3 -----

plot_scatter_el <- list()

# Set a fixed size for the plots
plot_size <- theme(plot.title = element_text(size = 40),
                   axis.title = element_text(size = 40),
                   axis.text = element_text(size = 36),
                   legend.title = element_text(size=40),
                   legend.text = element_text(size=40)
)

#Set the margins for the axis labels and titles
plot_margin <- theme(plot.margin = unit(c(1.2, 0.85, 0.85,  0.85), "lines"),
                     axis.title.y = element_text(margin = margin(r = 10)),
                     axis.title.x = element_text(margin = margin(t = 10)))


##Elast vs Sens (col)----
wave_set <- c(9,14,17) #FIGURE 1
for (i in wave_set) {
  df_plot <-filter(pert_df,Wave==i)
  df_plot$Age <- factor(df_plot$Age, levels=age_names)
  pos <- position_jitter(height = -0.6,width =-0.05, seed = 28)
  scat_el_sens <-  ggplot(df_plot , aes(x=Sens, y=El, color=Age) )+
    plot_size +
    plot_margin+
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth=1) +
    geom_hline(yintercept = 0.11, linetype = "dashed", color = "black", linewidth=1)+
    
    # Annotate each quadrant
    annotate("text", x = 0.17, y = -0.08, label = "MARGINAL", size=35) + 
    annotate("text", x = 2.2, y = -0.08, label = "EFFECTIVE",size=35) + 
    annotate("text", x = 0.17, y = 0.42, label = "INCIDENTAL",size=35) + 
    annotate("text", x = 2.2, y = 0.42, label = "MAIN",size=35,fontface = 'bold',color='darkred') +
    geom_point(aes(size=Rep_val), position=pos)+
    geom_label_repel(aes(label = Age),
                     position=pos,
                     size          = 26,
                     max.overlaps = Inf,
                     box.padding   = 0.5,
                     point.padding = 0.5,
                     force = 70,
                     force_pull = 14,
                     segment.size  = 0,
                     segment.color = NA)+
    
    # labs(y=bquote('Elasticity'~tilde(e)[j]),x= bquote('Sensitivity'~tilde(s)[j]),color='',size=bquote(v[j]), 
    #      title = paste0("Wave ",wrange[i],",  ", w_dates[i], '                              ',expression(~R[t])," =",round(R_series_mean[i],2)))+
    labs(y = bquote('Elasticity' ~ tilde(e)[j]), 
         x = bquote('Sensitivity' ~ tilde(s)[j]), 
         color = '', 
         size = bquote(V[j]), 
         #size = expression("v"[j] ~ "(infective" * "\n" * "value)"), 
         title = bquote("Wave " ~ .(wrange[i]) ~ ",  " ~ .(format(as.Date(w_dates[i], origin = "1970-01-01"), "%d %b %Y")) ~ 
                          '                              ' ~ R[t] ~ " =" ~ .(round(R_series_mean[i], 2))))+
  
    theme(axis.title.x = element_text(margin = margin(t = 10, unit = "pt"), size=100, face = 'bold'),
          axis.title.y = element_text(margin = margin(r = 10, unit = "pt"), size=100, face = 'bold'),
          axis.text = element_text(size=90),
          legend.text = element_text(size=90),
          legend.title = element_text(size = 100),
          plot.title =  element_text(size = 100)
    )+
    scale_x_continuous(limits = c(-0.1,2.5))+
    scale_y_continuous(limits = c(-0.1,0.45))+
    scale_color_manual(values=col_blind_pal)+
    
    guides(color = "none")
  
  # Set the desired range for size scaling
  desired_range <- c(10, 18)
  # Update the size scale
  scat_el_sens <- scat_el_sens + scale_size_continuous(range = desired_range)
  
  ggsave(paste0("El_vs_Sens_",i, ".png"), plot = scat_el_sens, width = 10, height = 7.5, dpi = 500)
}




  #INDEX S_j^t------


sens_index_df <- data.frame(Wave=integer(),
                            Age=character(),
                            Sens_index=double(),
                            Date=Date(),
                            stringsAsFactors=FALSE)
count <- 0
for (i in 1:length(wrange)) {
  for (k in 1:length(age_breaks)) {
    new_record <- list(as.integer(8+i),
                       age_names[k],
                       (colSums(el_list[[i]])[k]>0.11) + 
                         (v1_list[[i]][k]>mean(unlist(v1_list)) | colSums(sens_list[[i]])[k]>1) +
                         (v1_list[[i]][k]>mean(unlist(v1_list)) & colSums(sens_list[[i]])[k]>1),
                       w_dates[8+i])
    sens_index_df[nrow(sens_index_df) + 1,] = new_record 
    
  }
}



  
time_hor <- length(w_dates)
xaxis_el <-format(as.Date(w_dates[9:time_hor]),  "%d %b %y")
scaling_Rt <- 5
sens_index_df$Age <- factor(sens_index_df$Age,levels =age_names)


# 
# gb <-ggplot()
# gb+
#     #NDICES
#   geom_point(data =sens_index_df %>% filter(Wave>8),aes(x=Date,y=Sens_index,color=Age), lwd=1.6, position=position_jitter(h=0.15,w=0.15))+
#   geom_hline(yintercept = 1*scaling_Rt,linetype="dashed",color='darkred')+
#   geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt,linetype="Rt"), lwd=1, color='darkred')+
#   
#   
#   scale_x_date(breaks=w_dates[9:time_hor],labels = xaxis_el,   limits=c(w_dates[9],w_dates[time_hor]),
#                expand = c(0.02, 0),sec.axis = sec_axis(~ . ,  breaks = w_dates[9:time_hor], labels = 9:time_hor , name = "Wave"))+
#   scale_linetype_manual(values=c('solid'))+
#   scale_y_continuous(sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
#   
#   theme_bw(base_size = 20) +
#   theme(axis.line = element_line(colour ="black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.text.x.bottom = element_text(size=18, angle=75, vjust=0.5, family =fam), 
#         axis.text.y  = element_text(size=18, family =fam),
#         axis.text.y.right = element_text(color = 'darkred'),
#         legend.text  = element_text(size=18, family =fam),
#         legend.key.size = unit(0.7,'cm'),
#         legend.key.width = unit(1,'cm'),
#         axis.text.x.top = element_text(size=18, vjust=0.5,angle=25, family =fam)
#   ) +
#   labs(y=bquote("Composite Index S_I"),
#        x="",linetype='',
#        color='Age') 
# 
# 
# 




##TABLE 1-----
row_names <- c(age_names,"Stringency index", "School & Work" , "Travel ban", 
                "Avg. contacts/day", "Avg. contacts/day (-18)" ,"Avg. contacts/day(+18)","Rt (mean)")

table_df <- data.frame(rep(0,length(row_names)), row.names =row_names )

cicle <- 1
for (i in 9:length(w_dates)) {
  index_df <- filter(stringency, date==w_dates[i])
  col <- c(round(as.integer(filter(sens_index_df,Wave==i)$Sens_index),0),  
           index_df$StringencyIndex_Average,
           (index_df$C1M_School.closing+index_df$C2M_Workplace.closing)/(max(stringency$C1M_School.closing)+max(stringency$C2M_Workplace.closing))*100,
           (index_df$C8EV_International.travel.controls+ index_df$C5M_Close.public.transport)/(max(stringency$C8EV_International.travel.controls)+max(stringency$C5M_Close.public.transport))*100,
           weighted.mean(filter(pert_df, Wave==i)$Cont,N_vec),
           weighted.mean(filter(pert_df, Wave ==i)$Cont[1:3],N_vec[1:3]),
           weighted.mean(filter(pert_df, Wave ==i)$Cont[4:9],N_vec[4:9]),
           Rt_list$Mean[i]
  )
  table_df <- cbind(table_df,col)
  colnames(table_df)[cicle+1] <- paste0('Wave ',i)
  cicle <- cicle+1
}
table_df <- table_df[,-1] 
xtable(table_df[,1:12], caption = "Sensitivity and Stringency indices Nov-Dec 2020", label = "tab:ss_ND2020")
xtable(table_df[,12:21])
xtable(table_df[,22:34])
#ComiX Wave dates
df_dates <- data.frame(format(w_dates[27:43], "%d/%b/%y"))
t(df_dates)
xtable(t(df_dates))



#TABLE S1----

df_tab <- data.frame(Wave = 9:length(w_dates), Rmax = unique(filter(pert_df, Wave %in% 9:length(w_dates))$Tot_S_delta), Date=unique(pert_df$Date)[9:42])
df_tab$Date <- format(as.Date(df_tab$Date),  "%d %b %y")
n_1 <- round(length(w_dates)/2)
n_2 <- length(w_dates)-n_1
df1 <- head(df_tab, n_1) 
df2 <- tail(df_tab, n_2) 
df_combined <- cbind(df1, df2)
names(df_combined) <- c("Wave", "R_max","Date", "Wave", "R_max","Date")
table_combined <- xtable(df_combined)
print(table_combined, include.rownames=FALSE)


#FIGURE S3----


get_max_indices <- function(mat) {
  which(mat == max(mat, na.rm = TRUE), arr.ind = TRUE)
}

inf_val_list <- lapply(v1_list, function(l) outer(l,l,'/'))

max_cont_ind <-lapply(c_asym, get_max_indices)
max_k_ind <-lapply(ngms, get_max_indices)
max_sens_ind <- lapply(sens_list, get_max_indices)
max_el_ind<- lapply(el_list, get_max_indices)
max_v_ind<- lapply(inf_val_list, get_max_indices)

x_coords_c <- sapply(max_cont_ind, function(ind) ind[1,2]) #Columns
y_coords_c <- sapply(max_cont_ind, function(ind) ind[1,1]) #Rows
x_coords_k <- sapply(max_k_ind, function(ind) ind[1,2]) #Columns
y_coords_k <- sapply(max_k_ind, function(ind) ind[1,1]) #Rows
x_coords_s <- sapply(max_sens_ind, function(ind) ind[1,2]) #Columns
y_coords_s <- sapply(max_sens_ind, function(ind) ind[1,1]) #Rows
x_coords_e <- sapply(max_el_ind, function(ind) ind[1,2]) #Columns
y_coords_e <- sapply(max_el_ind, function(ind) ind[1,1]) #Rows
x_coords_v <- sapply(max_v_ind, function(ind) ind[1,2]) #Columns
y_coords_v <- sapply(max_v_ind, function(ind) ind[1,1]) #Rows

df_c <- data.frame(x = x_coords_c, y = y_coords_c)
df_k <- data.frame(x = x_coords_k, y = y_coords_k)
df_s <- data.frame(x = x_coords_s, y = y_coords_s)
df_e <- data.frame(x = x_coords_e, y = y_coords_e)
df_v <- data.frame(x = x_coords_v, y = y_coords_v)


# CONTACTS
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_c)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)

p_c <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = 8) +
  
  scale_fill_gradient(low="white", high="red") +
  labs(title="Contacts", x="Age (contactee)", y="Age (contacted)", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         axis.text = element_text(size=20),
                         legend.text = element_text(size=18),
                         legend.title = element_text(size = 30),
                         plot.title =  element_text(size = 35, face = 'bold')  
  )




# NGM
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_k)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)

p_k <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = 8) +
  
  scale_fill_gradient(low="white", high="red") +
  labs(title="NGM ", x="Age (Infectious)", y="Age (Susceptible)", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         axis.text = element_text(size=20),
                         legend.text = element_text(size=18),
                         legend.title = element_text(size = 30),
                         plot.title =  element_text(size = 35, face = 'bold')
  )




# Sensitivity 
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_s)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)


p_s <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = 8) +
  
  scale_fill_gradient(low="white", high="red") +
  labs(title="Sensitivity", x="Age (columns)", y="Age (rows)", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         axis.text = element_text(size=20),
                         legend.text = element_text(size=18),
                         legend.title = element_text(size = 30),
                         plot.title =  element_text(size = 35, face = 'bold')
                         )


# Elasticity
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_e)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)


p_e <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = 8) +
  
  scale_fill_gradient(low="white", high="red") +
  labs(title="Elasticity", x="Age (columns)", y="Age (rows)", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         axis.text = element_text(size=20),
                         legend.text = element_text(size=18),
                         legend.title = element_text(size = 30),
                         plot.title =  element_text(size = 35, face = 'bold')
  )

# Inv val
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_v)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)


p_v <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = 8) +
  
  scale_fill_gradient(low="white", high="red") +
  labs(title=bquote("Inf. value "~v[i]/v[j]), x="Age (j)", y="Age (i)", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=20),
                         axis.title.y = element_text(size=20),
                         axis.text = element_text(size=20),
                         legend.text = element_text(size=18),
                         legend.title = element_text(size = 30),
                         plot.title =  element_text(size = 35, face = 'bold')
  )

#FIGURE S3-----
combined_plot <- ggarrange(p_c, p_k, p_s, p_e, p_v, 
                           ncol=3, nrow=2)



