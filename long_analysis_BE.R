# ==============================================================================
# Title: Longitudinal Sensitivity Analysis: Graphical Outputs for Belgium (Angeli et al., 2024)
# Author: Leonardo Angeli  
# Institution: Univesrsity of Hasselt
# Date: 09/2024
# ==============================================================================

# Description:
# This script generates the graphical outputs for the main text of the study by Angeli et al. (2024),
# focusing on the longitudinal sensitivity analysis in Belgium.
# It includes:
# - Importing required libraries and setting up the environment.
# - Loading and preparing data, including NPIs, genomic surveillance, and vaccine data.
# - Calculating matrix differences between Next Generation Matrices (NGMs) across waves.
# - Generating several figures (scatter plots, contact matrices) that visualize 
#   the sensitivity, elasticity, and contact data across age groups and waves.

# Dependencies:
# This script relies on pre-calculated variables (such as `ngms`, `sens_list`, `el_list`)
# and external data sources, including genomic surveillance data, vaccine data, 
# and NPIs. Ensure these datasets are available before running the script.

# Output:
# The script outputs several graphical files (scatter plots, contact matrices)
# to visualize the sensitivity analysis results, as described in the paper.

# ==============================================================================
old <- Sys.time()
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
library(gridExtra)
library(cowplot)
library(ggtext)
library(patchwork)
library(gghighlight)
library(ggrepel)
library(sysfonts)

rm(list=ls(all=TRUE))
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}
# Load preliminary operations
suppressPackageStartupMessages(source('./R/socrates_main.R'))  # Main function source for modeling
source('./add_functions.R')  # Additional function definitions
source('./helpers.R')  # Additional function definitions (supplementary results)
source('./wave_specific_an.R')  # Call Wave-specific analysis 

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
#matrix_diff <- list(ngms[[1]] - wave8_NGM)  # Initialize with difference from wave 8
#for now we duplicate the first ngm of the list, but we will have to evaluate the imputed matrix for wave 8
matrix_diff<- list(ngms[[1]])
matrix_diff <- c(matrix_diff, Map(`-`, ngms[-1], ngms[-length(ngms)]))  # Calculate differences between successive NGMs


pert_df <- data.frame(Wave=integer(),
                      Date=Date(),
                      Age=character(),
                      K_col=double(),
                      K_rows=double(),
                      Sens=double(),
                      Sens_standard=double(),
                      Sens_rows=double(),
                      Sens_susc=double(),
                      Sens_susc_stand=double(),
                      El=double(),
                      Inf_val=double(),
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
                       colSums(sens_list[[i]])[k],
                       colSums(sens_list[[i]])[k]/sum(sens_list[[i]]),
                       rowSums(sens_list[[i]])[k],
                       rowSums(lower_level_sens[[i]]$S_Susc)[k],
                       rowSums(lower_level_sens[[i]]$S_Susc)[k]/sum(lower_level_sens[[i]]$S_Susc),
                       colSums(el_list[[i]])[k],
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
switch (case,
        "case2" = {
          pert_df <- pert_df %>%
            mutate(AgeGroup = case_when(
              Age %in% age_names[1:3] ~ 'Minors',
              Age %in% age_names[4:7] ~ 'Adults',
              Age %in% age_names[8:9] ~ 'Elderly',
              TRUE ~ 'Other'  # For any age not covered in the above categories
            ))
        },
        "case3" = {
          pert_df <- pert_df %>%
            mutate(AgeGroup = case_when(
              Age %in% age_names[1] ~ 'Minors',
              Age %in% age_names[2:3] ~ 'Adults',
              Age %in% age_names[4:5] ~ 'Elderly',
              TRUE ~ 'Other'  # For any age not covered in the above categories
            ))
        },
        # Default case to handle unexpected 'case' values
        {
          warning("Unexpected case value: ", case)
        }
)

pert_df$Age <- factor(pert_df$Age,levels =age_names)

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
                     round(Rt_list$Mean[i],2),
                     round(Rt_list$Up[i]),
                     round(Rt_list$Low[i])
                     
  )
  
  rep_df[nrow(rep_df) + 1,] = new_record 
}



mean_df<- data.frame(Wave=integer(),
                     Date=Date(),
                     Cont=double(),
                     El=double(),
                     El_w=double(),
                     Sens=double(),
                     Sens_cols_w=double(),
                     Sens_rows=double(),
                     Sens_rows_w=double(),
                     Sens_susc=double(),
                     Inf_val=double(),
                     Inf_val_w=double(),
                     Stable_stage=double(),
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
                     mean(filter(pert_df,Wave==8+i)$Sens_susc),
                     mean(filter(pert_df,Wave==8+i)$Inf_val),
                     weighted.mean(filter(pert_df,Wave==8+i)$Inf_val,N_vec),
                     mean(filter(pert_df,Wave==8+i)$Stable_stage),
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
susc_df$Age <- factor(susc_df$Age,levels =age_names)




#FIGs 1,2,3 -----

plot_scatter_el <- list()
size_text_plot <- 135  # Increase the text size for better readability in papers
plot_margin <- theme(plot.margin = unit(c(1.5, 1.2, 1.2,  1.2), "lines"),  # Adjust margins slightly for better spacing
                     axis.title.y = element_text(margin = margin(r = 12)),
                     axis.title.x = element_text(margin = margin(t = 12)))
##Elast vs Sens (col)
wave_set <- c(31,32,34) #FIGURE 1
for (i in wave_set) {
  df_plot <-filter(pert_df,Wave==i)
  df_plot$Age <- factor(df_plot$Age, levels=age_names)
  pos <- position_jitter(height = -0.6,width =-0.05, seed = 28)
  scat_el_sens <-  ggplot(df_plot , aes(x=Sens, y=El, color=Age) )+
    # plot_size +
    plot_margin+
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth=1) +
    geom_hline(yintercept = 0.11, linetype = "dashed", color = "black", linewidth=1)+
    
    # Annotate each quadrant
    annotate("text", x = 0.25, y = -0.08, label = "MARGINAL", size=size_text_plot*0.3) + 
    annotate("text", x = 2.2, y = -0.08, label = "EFFECTIVE",size=size_text_plot*0.3) + 
    annotate("text", x = 0.28, y = 0.42, label = "INCIDENTAL",size=size_text_plot*0.30) + 
    annotate("text", x = 2.2, y = 0.42, label = "MAIN",size=size_text_plot*0.3,fontface = 'bold',color='darkred') +
    geom_point(aes(size=Inf_val), position=pos)+
    geom_label_repel(aes(label = Age),
                     size= size_text_plot*0.3,
                     max.overlaps = Inf,
                     box.padding   = 1.1,  # Increased padding around the label
                     point.padding = 0.8,  # Increased distance from the point to the label
                     force = 110,          # Stronger repelling force
                     force_pull = 40,      # Stronger pull back force
                     fill = alpha("white",0.5),
                     segment.size = 1.2,   # Thicker segment lines
                     segment.color =  col_blind_pal,
                     nudge_x = ifelse(df_plot$Sens > 1, 0.4, -0.4)) +
    
    
    # labs(y=bquote('Elasticity'~tilde(e)[j]),x= bquote('Sensitivity'~tilde(s)[j]),color='',size=bquote(v[j]), 
    #      title = paste0("Wave ",wrange[i],",  ", w_dates[i], '                              ',expression(~R[t])," =",round(R_series_mean[i],2)))+
    labs(y = bquote('Elasticity' ~ tilde(e)[j]), 
         x = bquote('Sensitivity' ~ tilde(s)[j]), 
         color = '', 
         size = bquote(V[j]), 
         #size = expression("v"[j] ~ "(infective" * "\n" * "value)"), 
         title = bquote("Wave " ~ .(i) ~ ",  " ~ .(format(as.Date(w_dates[i], origin = "1970-01-01"), "%d %b %Y")) ~ 
                          '                 ' ~ R[t] ~ " =" ~ .(round(Rt_list$Mean[i], 2))))+
    
    theme(axis.title.x = element_text(margin = margin(t = 12, unit = "pt"), size=size_text_plot, face = 'bold'),
          axis.title.y = element_text(margin = margin(r = 12, unit = "pt"), size=size_text_plot, face = 'bold'),
          axis.text = element_text(size=size_text_plot),
          legend.text = element_text(size=size_text_plot),
          legend.title = element_text(size = size_text_plot),
          plot.title =  element_text(size = size_text_plot)
    )+
    scale_x_continuous(limits = c(-0.1,2.5))+
    scale_y_continuous(limits = c(-0.1,0.45))+
    scale_color_manual(values=col_blind_pal)+
    
    guides(color = "none")
  
  # Set the desired range for size scaling
  desired_range <- c(10, 20)
  # Update the size scale
  scat_el_sens <- scat_el_sens + scale_size_continuous(range = desired_range)
  
  ggsave(paste0("El_vs_Sens_",i, ".png"), plot = scat_el_sens, width = 12, height = 9.5, dpi = 400)
}

######Alternative FIG 1,2,3####
#Combine 3 scatterplots together for better readability

plot_list <- list()
size_text_plot <- 190  # Text size for better readability in papers
plot_margin <- theme(plot.margin = unit(c(1.5, 1.2, 1.2,  1.2), "lines"),  # Adjust margins slightly for better spacing
                     axis.title.y = element_text(margin = margin(r = 12)),
                     axis.title.x = element_text(margin = margin(t = 12)))

## Elast vs Sens (col)
wave_set <- c(31,32,34) # Select the waves to plot
panel_labels <- c("a)", "b)", "c)")
# Determine the global range for the size scale across all waves
global_size_range <- range(filter(pert_df, Wave %in% wave_set)$Inf_val, na.rm = TRUE)

for (i in wave_set) {
  df_plot <- filter(pert_df, Wave == i)
  df_plot$Age <- factor(df_plot$Age, levels = age_names)
  pos <- position_jitter(height = -0.6, width = -0.05, seed = 28)
  
  scat_el_sens <- ggplot(df_plot, aes(x = Sens, y = El, color = Age, size = Inf_val)) +
    plot_margin +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
    geom_hline(yintercept = 0.11, linetype = "dashed", color = "black", linewidth = 1) +
    
    # Annotate each quadrant
    annotate("text", x = 0.25, y = -0.08, label = "MARGINAL", size = size_text_plot * 0.3) + 
    annotate("text", x = 2, y = -0.08, label = "EFFECTIVE", size = size_text_plot * 0.3) + 
    annotate("text", x = 0.28, y = 0.42, label = "INCIDENTAL", size = size_text_plot * 0.3) + 
    annotate("text", x = 2.2, y = 0.42, label = "MAIN", size = size_text_plot * 0.3, fontface = 'bold', color = 'darkred') +
    
    geom_point(position = pos) +
    annotate("text", x = -0.05, y = -0.1, label = panel_labels[i], size = size_text_plot * 0.35, hjust = 0, vjust = 0, fontface = "bold") +
    
    labs(y = bquote('Elasticity' ~ tilde(e)[j]), 
         x = bquote('Sensitivity' ~ tilde(s)[j]), 
         color = 'Age', 
         size = bquote(V[j]), 
         title = bquote("Wave " ~ .(i) ~ ",  " ~ .(format(as.Date(w_dates[i], origin = "1970-01-01"), "%d %b %Y")) ~ 
                          ', ' ~ R[t] ~ " =" ~ .(round(Rt_list$Mean[i], 2)))) +
    
    theme(axis.title.x = element_text(margin = margin(t = 12, unit = "pt"), size = size_text_plot, face = 'bold'),
          axis.title.y = element_text(margin = margin(r = 12, unit = "pt"), size = size_text_plot, face = 'bold'),
          axis.text = element_text(size = size_text_plot * 0.9),
          legend.text = element_text(size = size_text_plot * 0.85),
          legend.title = element_text(size = size_text_plot),
          plot.title = element_text(size = size_text_plot),
          legend.spacing.x = unit(1, 'cm'),  # Increase horizontal space between legend keys
          legend.spacing.y = unit(0.7, 'cm')) +  # Increase vertical space between legend rows
    scale_x_continuous(limits = c(-0.1, 2.5)) +
    scale_y_continuous(limits = c(-0.1, 0.45)) +
    scale_color_manual(values = col_blind_pal) +
    scale_size_continuous(range = c(20, 30), limits = global_size_range) +  # Apply global range for the size scale
    
    guides(color = guide_legend(override.aes = list(size = 15)), size = guide_legend(legend.position = 'right'))  # Increase point size in Age legend
  
  # Store the plot in the list
  plot_list[[as.character(i)]] <- scat_el_sens
}

# Combine the plots using patchwork
combined_plot1 <- plot_list[[1]] + theme(legend.position = "none")
combined_plot2 <- plot_list[[2]] + theme(legend.position = "bottom") & guides(size = "none")
# Now, separate the size legend and move it to the left
size_legend_plot <- plot_list[[3]] +
  guides(color = "none")  # Remove the color legend, keeping only the size legend



# Combine the size legend and the three plots
final_plot <- (combined_plot1 + plot_spacer() + combined_plot2 + plot_spacer() + size_legend_plot) +
  plot_layout(widths = c(1, 0.085, 1, 0.085, 1)) +  # Adjust width ratios with spacers
  plot_annotation(tag_levels = "a") +  # Add labels (a), (b), (c)
  theme(plot.margin = margin(5, 10, 5, 10))

# Save the final combined plot
name_plot <- paste0("Combined_El_vs_Sens",wave_set[1],"_",wave_set[3],".png" )

ggsave(name_plot, plot = final_plot, width = 40, height = 19, dpi = 400)







#######FIG 2 (d-f)----
#Plotting social contact matrix
w_set <- c(31, 32, 34)
mat_plot <- vector(mode = "list", length = 3)
cicle <- 1

for (wave_obs in w_set) {
  melted_mij <- reshape2::melt(c_asym[[wave_obs - 8]])
  
  # Round the values to 2 decimal places for display
  melted_mij$value_rounded <- round(melted_mij$value, 2)
  
  # Define text size variables
  text_size_var <- 150  # Base text size
  
  c_mat <- ggplot(melted_mij, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    geom_text(aes(label = value_rounded), size = text_size_var * 0.31, family = fam, color = "black") +  # Text labels
    scale_fill_gradient(low = "#FFFFF0", high = "red") +
    labs(
      title = paste0("Social contacts     Wave ", wave_obs),
      x = "",
      y = "",
      fill = ""
    ) +
    theme_minimal() +
    
    # Adjustments for axis labels, titles, and spacing
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = text_size_var*1.05 , family = fam, 
                                 margin = margin(t = 20),color = "black"),  # Add space above x-axis labels
      axis.text.y = element_text(size = text_size_var*1.05 , family = fam, 
                                 margin = margin(r = 20),color = "black"),  # Add space to the right of y-axis labels
      plot.title = element_text(size = text_size_var * 1.266, family = fam, 
                                margin = margin(b =40)),  # Add more space below the title
      plot.margin = unit(c(1, 1, 1, 1), "cm")  # Adjust the overall plot margin
    ) +
    
    scale_x_discrete(expand = expansion(add = 0.5)) +
    scale_y_discrete(expand = expansion(add = 0.5)) +
    
    # Remove the legend bar for the gradient
    guides(fill = "none")
  
  mat_plot[[cicle]] <- c_mat
  cicle <- cicle + 1
}

# Combine the plots and ensure the legend is positioned horizontally at the bottom
plotC <- ggarrange(
  mat_plot[[1]], plot_spacer(), mat_plot[[2]], plot_spacer(), mat_plot[[3]],
  ncol = 5, widths = c(1, 0.1, 1, 0.1, 1),  # Adjust widths to add space between plots
  common.legend = TRUE, legend = "bottom"
)

# Save the final combined plot
name_plot_c <- paste0("c_mat_",w_set[1],"_",w_set[3],".png" )

ggsave(name_plot_c, plot = plotC, width = 40, height = 19, dpi = 400)




#FIG4: ELASTICITY EVOLUTION----

#BACKGROUND PLOT
xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
annotations_size <- 35
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
  annotate("text", x=w_dates[13]+5, y=y_max, vjust=2, label="Schools & non-ess. shops reopen",size=annotations_size,family=fam)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[8]+15,y=y_max, vjust=2, label=  df_measures$Descr_short[8],size=annotations_size, family=fam)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=w_dates[29],y=y_max, vjust=2, label= "School re-opening & Lifting",size=annotations_size, family=fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[10]+18,y=y_max, vjust=2, label= df_measures$Descr_short[10] ,size=annotations_size, family=fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[11]+28,y=y_max, vjust=2, label= "Lifting" ,size=annotations_size, family=fam)


#Stringency index
#annotate("text", x=w_dates[24], y=0.53, label="Stringency index(%)",size=6, angle=0)+

#VOCs
y_rect <- -0.02
y_min <- -0.12
y_text <- -0.07
base <- base+
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = y_min, ymax =y_rect)+
  annotate("text", x=w_dates[12], y=y_text, label= "Wild Type",size=annotations_size, family=fam)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=w_dates[20], y=y_text, label= "Alpha",size=annotations_size, family=fam)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=y_text, label= "Delta",size=annotations_size, family=fam)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=y_text, label= "Omicron",size=annotations_size, family=fam)+
  theme(
    plot.margin=plot.margin
  )







scaling_Rt <- max(pert_df$El)*1.5
pert_df$Age <- factor(pert_df$Age,levels =age_names)

#####Minors----
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
        legend.title = element_text(size=size_text*1.3, family =fam),
        legend.text  = element_text(size=size_text, family =fam),
        legend.key.size = legend.key.size,
        legend.key.height = legend.key.height,
        legend.key.width = legend.key.width,
        axis.title.y.left = element_text(size=size_text*1.5, family =fam),
        axis.title.x.top = element_text(size=size_text*1.3, family =fam),
        axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  guides(
    linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
    color = guide_legend(order = 2)     # 'color' represents the Age group legend
  ) +
  labs(y=expression(tilde(e)[j]),
       x="",linetype='',
       color='Age group')

#ggsave("Minors.png", plot = minor_plot, width = 40, height = 25, dpi = 200)

#####Adults----

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
    legend.title = element_text(size=size_text*1.3, family =fam),
    legend.text  = element_text(size=size_text, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.5, family =fam),
    axis.title.x.top = element_text(size=size_text*1.3, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  labs(y=expression(tilde(e)[j]),
       x="",linetype='',
       color='')

#ggsave("Adults.png", plot = adults_plot, width = 40, height = 25, dpi = 200)

#####Elderly----

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
    legend.title = element_text(size=size_text*1.3, family =fam),
    axis.title.y.left = element_text(size=size_text*1.5, family =fam),
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

#ggsave("Elderly.png", plot = elderly_plot, width = 40, height = 25, dpi = 200)

# Combine the plots
combined_plot <- ggarrange(minor_plot, adults_plot, elderly_plot, ncol = 1, nrow = 3, legend = "right")

# Save the combined plot: this generates Figure 4 in the main text of the manuscript
ggsave( "AgeGroupPlot.png", plot = combined_plot, width = 58, height = 58, dpi =178,limitsize = FALSE)


#FIG 6 (susceptibles evolution)----
text_size <-90 
ann_size <- 30
y_text <- max(filter(susc_df, Wave >8)$Susc)*1.08
vocs_y <- 5*(10^4)
sec_axis_scale <- 1.8*10^6
gb <-ggplot()
susc_numb <- gb+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=w_dates[14]-7, y=y_text, label="Schools & non-ess. shops reopen",size=ann_size, angle=0, family=fam)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[8]+18, y=y_text, label=  df_measures$Descr_short[8],size=ann_size, angle=0, family=fam)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=w_dates[29], y=y_text, label= "School re-opening + gen. lifting",size=ann_size, angle=0, family=fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[10]+20, y=y_text, label= df_measures$Descr_short[10] ,size=ann_size, angle=0, family=fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = Inf)+
  annotate("text", x=df_measures$Date_start[11]+25, y=y_text, label= "Lifting" ,size=ann_size, angle=0, family=fam)+
  
  
  #VOCs
  
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=w_dates[12], y=-vocs_y, label= "Wild Type",size=ann_size, family=fam)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=w_dates[20], y=-vocs_y, label= "Alpha",size=ann_size, family=fam)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=-vocs_y, label= "Delta",size=ann_size, family=fam)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = 0)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=-vocs_y, label= "Omicron",size=ann_size, family=fam)+
  
  
  #Susceptibility
  geom_ribbon(data = susc_df, aes(x = Date, ymin = Susc_low, ymax = Susc_up, fill = Age), alpha = 0.2,show.legend = FALSE) +
  geom_line(data = susc_df, aes(x = Date, y = Susc, color = Age), linewidth = 0.8, linetype = "solid") +
  geom_point(data = susc_df, aes(x = Date, y = Susc, color = Age), size = 3.8) +
  
  #geom_line(data =susc_df %>% filter(Wave>8),aes(x=Date,y=Susc,color=Age), lwd=1.6)+
  geom_line(data=df_vaccines_1st %>% filter(Day>=w_dates[9]),aes(x=Day,y=Percnt*(sec_axis_scale),linetype="1st dose"),lwd=0.6)+
  geom_line(data=df_vaccines_2nd %>% filter(Day>=w_dates[9]),aes(x=Day,y=Percnt*(sec_axis_scale), linetype="2nd dose"),lwd=0.6)+
  
  scale_color_manual("Age Group",values=col_blind_pal, 
                     labels=age_names)+
  scale_x_date(breaks=w_dates[9:length(w_dates)],labels = xaxis_el,   limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0.0, 0),sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave"))+
  scale_y_continuous(sec.axis=sec_axis(~./sec_axis_scale,breaks = seq(0,1, by=0.25)),
                     breaks = c(500000, 10^6, 1.5*10^6), # Choose appropriate breaks for your data
                     labels = c("0.5*10^6", "10^6", "1.5*10^6"))+
  #theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x.bottom = element_text(size=text_size*0.9, angle=90, vjust=0.5, family =fam), 
        axis.text.y  = element_text(size=text_size, family =fam),
        axis.text.y.right = element_text(color = 'darkred'),
        legend.text  = element_text(size=text_size, family =fam),
        legend.key.height = unit(2, "cm"),
        legend.key.size = unit(4,'cm'),
        legend.key.width = unit(3.5,'cm'),
        legend.spacing.y = unit(1.5,'cm'),
        legend.title = element_text(size=text_size, family =fam),
        axis.title = element_text(size=text_size, family =fam),
        axis.text.x.top = element_text(size=text_size*0.9, vjust=0.5,angle=30, family =fam),
        axis.text.y.left = element_text(size=text_size, vjust=0.5,angle=90, family =fam)
  ) +
  labs(y=bquote("# susceptible individuals"),
       x="",
       color='Age Group', linetype="")+
  guides(color = guide_legend(override.aes = list(size = 7, linetype = 1, linewidth = 2)),
         fill = guide_legend(override.aes = list(alpha = 0.2)),
         linetype= guide_legend(override.aes = list(size = 5, linetype = c(1,3), linewidth = 1.5)))



# Save the combined plot
ggsave( "Susc_number.png", plot = susc_numb, width = 35, height = 17, dpi =230,limitsize = FALSE)


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



#FIGs S3,5,6: counterfactual scenarios-----
source('./helpers.R')  # Additional helper functions 
# Waves 22-25
range_1 <- 25:30 #select wave range
variants <- c("Wild Type","Alpha","Delta","Omicron")
VOC <- variants[2]
bas_index_1 <- range_1[1]-8
age_obs <- age_names[c(3:5)]

## Scenario 1: Evolving contact pattern with static susceptibles (for age group "age_obs")

## Scenario 2: Evolving susceptibility with static contact pattern (for age group "age")

## Scenario 3: Evolving contact pattern with static susceptibles (for all age groups)

## Scenario 4: Evolving susceptibility with static contact pattern (for all age groups)

for(s in 1:2){
  # Set the scenario flag
  scenario<- s  # Change this value to 1, 2, 3, or 4 to select the scenario (see helper.R function "select_scenario")
  
  # Evaluate the selected scenario
  scenario_result <- select_scenario(scenario_flag=scenario, range=range_1, bas_index=bas_index_1,
                                     susceptibles=Susc_list,
                                     contact_sym_pc=c_sym_per_c,
                                     contact_asym_pc= c_asym_per_c,
                                     age = age_obs, age_ints=age_names, contact_dir='both')
  
  # Prepare the data for plotting
  c_analysis <- prepare_plotting(scenario_result, age_breaks, w_dates, age_names, range_1)
  
  
  # Generate and save the plots for different attributes
  xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
  annotations_size <- 70
  size_text <- 180
  plot.margin = margin(t = 0.3, r = 0.5, b = 0.3, l = 0.5, unit = "cm")
  base_size =150
  legend.key.size = unit(8,'cm')
  legend.key.height = unit(10,'cm')
  legend.key.width = unit(6,'cm')
  
  
  age_obs_str <- if (length(age_obs) > 1) paste(age_obs, collapse = "_") else age_obs
  width <- 55
  height <- 40
  dpi <- 150
  #attributes <- c("El", "Sens", "Inf_val")
  attributes <- "El"
  for (attribute in attributes) {
    file_name <- paste0("./scenarios/",age_obs_str,"_",range_1[1],"_",range_1[length(range_1)] ,"/", attribute, "_S", scenario,"_",age_obs_str,"_",range_1[1],"_",range_1[length(range_1)],".png")
    generate_plot(c_analysis, attribute = attribute,range =  range_1, w_dates = w_dates,age = age_obs,age_names =  age_names,file_name =  file_name,variant = VOC,extra_variant = T,
                  title_lab=scenario, highlight = T, width =width, height = height,dpi = dpi, size_text = size_text )
    # file_name <- paste0("./scenarios/", attribute, "_baseline_",range_1[1],"_",range_1[length(range_1)],".png")
    # generate_plot(c_analysis = pert_df,attribute =  attribute, range = range_1,w_dates =  w_dates,age_names =  age_names,age = age_obs,file_name =  file_name,variant = VOC,
    #               title_lab='OBSERVED', width =width, height = height,dpi = dpi,size_text = size_text )
  }
  
}

# FIGs S4,7: Elasticity gradients-----

# Define necessary variables
wave_obs <- 13 
age_group <- 3  # Position age group (j): the sensitivity indicators of this groups will be derived for the columns and rows of the NGM
K <- ngms[[wave_obs - 8]]
iter <- seq_len(ncol(K))

# Preallocate lists to store results
second_order_sens <- vector("list", length(iter))
second_order_sens_row <- vector("list", length(iter))
second_order_elast <- vector("list", length(iter))
second_order_elast_row <- vector("list", length(iter))

# Loop through columns of K and compute sensitivity and elasticity
for (m in iter) {
  sens_result <- gradient_sens(A = K, j = age_group, m = m)
  sens_result_row <- gradient_sens(A = K, j = age_group, m = m, col = FALSE)
  
  atom_sens <- sens_result$sens
  atom_sens_row <- sens_result_row$sens
  atom_el <- sens_result$el
  atom_el_row <- sens_result_row$el
  
  names(atom_sens) <- names(atom_sens_row) <- names(atom_el) <- names(atom_el_row) <- age_names
  
  second_order_sens[[m]] <- atom_sens
  second_order_sens_row[[m]] <- atom_sens_row
  second_order_elast[[m]] <- atom_el
  second_order_elast_row[[m]] <- atom_el_row
}

# Summarize data
data <- data.frame(
  Age = factor(age_names, levels = age_names),
  Sens_col = simplify2array(lapply(second_order_sens, sum)),
  Sens_row = simplify2array(lapply(second_order_sens_row, sum)),
  El_col = simplify2array(lapply(second_order_elast, sum)),
  El_row = simplify2array(lapply(second_order_elast_row, sum))
)
data$Age <- factor(data$Age,levels =age_names)

# Select sensitivity indicator to plot
index_2_plot <- 3
indices <- c('Sens_col', 'Sens_row', 'El_col', 'El_row')
plot_title <- switch(index_2_plot,
                     "1" = paste0(age_names[age_group], ": Sensitivity gradient (columns)     wave ", wave_obs),
                     "2" = paste0(age_names[age_group], ": Sensitivity gradient (rows)     wave ", wave_obs),
                     "3" = paste0(age_names[age_group], ": Elasticity gradient (columns)     wave ", wave_obs),
                     "4" = paste0(age_names[age_group], ": Elasticity gradient (rows)     wave ", wave_obs))
direction <- switch(index_2_plot, "1" = "_columns_", "2" = "_rows_", "3" = "_columns_","4"= "_rows_")
measure <- switch(index_2_plot, "1" = "sens", "2" = "sens", "3" = "el", "4"= "el")
# Text size variables
title_size <-80
axis_title_size <- title_size*0.9
axis_text_size <- title_size*0.8
# Generate the barplot
sens_2_plot <- ggplot(data, aes_string(x = "Age", y = indices[index_2_plot], fill = "Age")) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = plot_title,
    x = "",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = title_size, family = fam),
    axis.title.x = element_text(size = axis_title_size, family = fam),
    axis.title.y = element_text(size = axis_title_size),
    axis.text.x = element_text(angle = 90, size = axis_text_size, family = fam),
    axis.text.y = element_text(size = axis_text_size),
    legend.position = 'none'
  ) + scale_fill_manual(values = col_blind_pal)

ggsave(paste0("second_order/", measure,"_sec_order_", age_names[age_group], direction, wave_obs, ".png"), plot = sens_2_plot, width = 20, height = 15, dpi = 150)


#FIGs S8,9. Exercise:vaccinating children-----
col_blind_pal <- c("#FF8C00", "#D62728", "#F0E442", "#117733", "#44AA99", "#0072B2", "#56B4E9", "#CC6677", "#882255")

# Generate lighter version of the palette
col_blind_pal_lighter <- alpha(col_blind_pal, 0.5) # Adjust alpha to lighten the colors


wave_obs <- 31 #
K <- ngms[[wave_obs - 8]]
#Modify susceptibles
susc <- Susc_list[[wave_obs-8]]
susc[1] <- susc[1]-10^5
susc[2] <- susc[2]-10^5
#susc[3] <- Susc_list[[26-8]][3]*0.9 #Assume [12,18) had a higher susceptible proportion ( uncomment for Fig 9b)
El_mod <- NGM_analysis(wave_obs-8,Susc_n = susc,C_sym = c_sym_per_c[[wave_obs-8]],C_asym =c_asym_per_c[[wave_obs-8]] )$El
K_mod <- NGM_analysis(wave_obs-8,Susc_n = susc,C_sym = c_sym_per_c[[wave_obs-8]],C_asym =c_asym_per_c[[wave_obs-8]] )$NGM

# Prepare the data
data <- data.frame(
  Age = age_names,
  El = colSums(calc.el(K)),
  El_mod = colSums(calc.el(K_mod))
)
data$Age <- factor(data$Age, levels = age_names)

data_long <- melt(data, id.vars = "Age")


# Combine the palettes into a named vector
color_map <- setNames(
  c(rep(col_blind_pal, length.out = length(unique(data_long$Age))), 
    rep(col_blind_pal_lighter, length.out = length(unique(data_long$Age)))),
  interaction(rep(unique(data_long$variable), each = length(unique(data_long$Age))), unique(data_long$Age))
)

text_size <- 90

ex_plot <- ggplot(data_long, aes(x = Age, y = value, fill = interaction(variable, Age))) +
  geom_bar(stat = "identity", position = "dodge", color = "white") +
  scale_fill_manual(values = color_map) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(family = fam, size = text_size), 
    axis.title.y = element_text(family = fam, size = text_size*1.2,angle = 90),   
    axis.text.x= element_text(family = fam, size = text_size*1.2,angle = 90),
    axis.text.y= element_text(family = fam, size = text_size*1.2) 
  ) +
  labs(
    title = bquote("Wave " * .(wave_obs) * "             " * R[t] * " = " * .(round(eigen(K)$value[1], 2)) * ",   " * tilde(R)[t] * " = " * .(round(eigen(K_mod)$value[1], 2))),
    x = "",
    y = bquote(tilde(e)[j])
  )

ggsave(paste0("Ex2_susc_wave_",wave_obs,".png"), plot = ex_plot, width = 25, height = 18, dpi = 150)





#FIGURE S10 (Overall trends)----


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

text_size <- 110
# CONTACTS
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_c)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)

p_c <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = text_size*0.3,fontface = 'plain') +
  
  scale_fill_gradient(low = "white", high = "red", guide = "none") +  
  
  labs(title="Contacts", x="", y="", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=text_size*0.7),
                         axis.title.y = element_text(size=text_size*0.7),
                         axis.text = element_text(size=text_size),
                         axis.text.x = element_text(angle=90),
                         legend.text = element_text(size=text_size*0.7),
                         plot.title =  element_text(size = text_size, face = 'plain'),
                         text = element_text(family=fam)
  )




# NGM
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_k)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)

p_k <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = text_size*0.3,fontface = 'plain') +
  
  scale_fill_gradient(low = "white", high = "red", guide = "none") +  
  
  labs(title="NGM ", x="", y="", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=text_size*0.7),
                         axis.title.y = element_text(size=text_size*0.7),
                         axis.text = element_text(size=text_size),
                         axis.text.x = element_text(angle=90),
                         legend.text = element_text(size=text_size*0.7),
                         plot.title =  element_text(size = text_size, face = 'plain'),
                         text = element_text(family=fam)
  )




# Sensitivity 
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_s)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)


p_s <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = text_size*0.3,fontface = 'plain') +
  
  scale_fill_gradient(low = "white", high = "red", guide = "none") +  
  
  labs(title="Sensitivity", x="", y="", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=text_size*0.7),
                         axis.title.y = element_text(size=text_size*0.7),
                         axis.text = element_text(size=text_size),
                         axis.text.x = element_text(angle=90),
                         legend.text = element_text(size=text_size*0.7),
                         plot.title =  element_text(size = text_size, face = 'plain'),
                         text = element_text(family=fam)
  )

# Elasticity
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_e)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)


p_e <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = text_size*0.3,fontface = 'plain') +
  
  scale_fill_gradient(low = "white", high = "red", guide = "none") +  
  
  labs(title="Elasticity", x="", y="", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=text_size*0.7),
                         axis.title.y = element_text(size=text_size*0.7),
                         axis.text = element_text(size=text_size),
                         axis.text.x = element_text(angle=90),
                         legend.text = element_text(size=text_size*0.7),
                         plot.title =  element_text(size = text_size, face = 'plain'),
                         text = element_text(family=fam)
  )

# Inv val
all_couples <- expand.grid(x = seq(1, 9), y = seq(1, 9))
df_freq <- merge(all_couples, as.table(table(df_v)), by = c("x", "y"), all.x = TRUE)
df_freq$Freq[is.na(df_freq$Freq)] <- 0
df_heatmap <- as.data.frame(df_freq)


p_v <- ggplot(df_heatmap, aes(x=as.numeric(as.character(x)), y=as.numeric(as.character(y)), fill=Freq) ) + 
  geom_tile() +
  geom_text(aes(label = ifelse(Freq > 0, Freq, "")), vjust = 0.5, color = "black", size = text_size*0.3,fontface = 'plain') +
  
  scale_fill_gradient(low = "white", high = "red", guide = "none") +  
  
  labs(title=bquote("Inf. value "~v[i]/v[j]), x="", y="", fill='') +
  scale_x_continuous(breaks = 1:9, labels = age_names)+
  scale_y_continuous(breaks = 1:9, labels = age_names)+
  theme_minimal()+ theme(panel.grid.major = element_blank(), 
                         panel.grid.minor = element_blank(),
                         axis.ticks = element_line(colour = "black"),
                         axis.title.x = element_text(size=text_size*0.7),
                         axis.title.y = element_text(size=text_size*0.7),
                         axis.text = element_text(size=text_size),
                         axis.text.x = element_text(angle=90),
                         legend.text = element_text(size=text_size*0.7),
                         plot.title =  element_text(size = text_size, face = 'plain'),
                         text = element_text(family=fam)
  )

combined_plot <- ggarrange(p_c, p_k, p_s, p_e, p_v, 
                           ncol=3, nrow=2)
ggsave( "Max_index.png", plot = combined_plot, width = 40, height = 25, dpi =230)































#FIG S11,12:NGM columns & rows----

#BACKGROUND PLOT
scaling_Rt <- max(pert_df$K_col)*1.3
xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
annotations_size <- 34
size_text <- 100
plot.margin = margin(t = 0.3, r = 0, b = 0.3, l = 0.5, unit = "cm")
#base_size = 45
legend.key.size = unit(2,'cm')
legend.key.height = unit(1.5,'cm')
legend.key.width = unit(4,'cm')
y_max <-3.3
gb <-ggplot()
base<- gb+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  #annotate("text", x=w_dates[14], y=max(rep_df$R0*scaling_Rt), label="Schools & non-ess. shops reopen",size=annotations_size, angle=0)+
  annotate("text", x=w_dates[13]+8, y=y_max, vjust=2, label="Schools & non-ess. shops reopen",size=annotations_size,family=fam)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[8]+22,y=y_max, vjust=2, label=  df_measures$Descr_short[8],size=annotations_size, family=fam)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=w_dates[29],y=y_max, vjust=2, label= "School re-opening & Lifting",size=annotations_size, family=fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[10]+18,y=y_max, vjust=2, label= df_measures$Descr_short[10] ,size=annotations_size, family=fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[11]+28,y=y_max, vjust=2, label= "Lifting" ,size=annotations_size, family=fam)


#Stringency index
#annotate("text", x=w_dates[24], y=0.53, label="Stringency index(%)",size=6, angle=0)+

#VOCs
y_rect <- 0.15
y_min <- -0.02
y_text <- 0.08
base <- base+
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = y_min, ymax =y_rect)+
  annotate("text", x=w_dates[12], y=y_text, label= "Wild Type",size=annotations_size, family=fam)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=w_dates[20], y=y_text, label= "Alpha",size=annotations_size, family=fam)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=y_text, label= "Delta",size=annotations_size, family=fam)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=y_text, label= "Omicron",size=annotations_size, family=fam)+
  theme(
    plot.margin=plot.margin
  )

#Rows or columns
display_sum <- c("rows","col")
plot_dir <- display_sum[2]

#####Minors----
minor_col <- base+
  
  geom_line(data =pert_df %>% filter(Wave>8 & AgeGroup == 'Minors'),aes_string(x="Date",y=paste0("K_",plot_dir),color="Age"), lwd=3.4)+
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt,linetype="Rt"), lwd=1.8)+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 1, lwd=1.4, linetype="longdash")+
  #geom_segment(x = w_dates[13], xend = w_dates[13], y = -Inf, yend = 0.5, linetype = 2, col = 'grey20',lwd=2.7) +
  
  scale_x_date( breaks=w_dates[9:length(w_dates)],
                #labels = rep('', length(w_dates[9:length(w_dates)])),
                labels = xaxis_el, 
                limits=c(w_dates[9],w_dates[length(w_dates)]),
                expand = c(0, 0),
                sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave")
  )+
  scale_linetype_manual(values=c('solid'))+
  scale_y_continuous(expand = c(0, 0),sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names,drop = FALSE)+
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
    legend.title = element_text(size=size_text*1.4, family =fam),
    legend.text  = element_text(size=size_text*1.25, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.2, family =fam),
    axis.title.x.top = element_text(size=size_text*1.2, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  guides(
    linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
    color = guide_legend(order = 2)     # 'color' represents the Age group legend
  ) +
  labs(y="",#paste0("NGM ",plot_dir," sum"),
       x="",linetype='',
       color='') 

#ggsave("Minors_col.png", plot = minor_col, width = 40, height = 15, dpi = 100)

#####Adults----
adults_col <- base+geom_line(data =pert_df %>% filter(Wave>8 & AgeGroup != 'Minors'),aes_string(x="Date",y=paste0("K_",plot_dir),color="Age"), lwd=3.4)+
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt,linetype="Rt"), lwd=1.8)+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 1, lwd=1.4, linetype="longdash")+
  
  scale_x_date(breaks=w_dates[9:length(w_dates)],
               #labels = xaxis_el,
               labels = rep('', length(w_dates[9:length(w_dates)])),
               limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0,0),
               sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates), name = "Comix Wave")
  )+
  scale_linetype_manual(values=c('solid'))+
  scale_y_continuous( expand = c(0,0), sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names,drop = FALSE)+
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
    legend.title = element_text(size=size_text*1.4, family =fam),
    legend.text  = element_text(size=size_text*1.25, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.2, family =fam),
    axis.title.x.top = element_text(size=size_text*1.2, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  labs(y="",#paste0("NGM ",plot_dir," sum"),
       x="",linetype='',
       color='') 

#ggsave(path = path_graphs ,"Adults_col.png", plot = adults_col, width = 40, height = 15, dpi = 100)

# Combine the plots
combined_col <- ggarrange(adults_col, minor_col, ncol = 1, nrow = 2,common.legend = TRUE, legend = "right" )
# Save the combined plot
ggsave( paste0("AgeGroup",plot_dir,".png"), plot = combined_col, width = 40, height =30, dpi =250,limitsize = FALSE)

#FIG S13:CONTACTS----

scaling_Rt <- max(pert_df$Cont)*1.1
pert_df$Age <- factor(pert_df$Age,levels =age_names)
xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
annotations_size <- 35
size_text <- 120
plot.margin = margin(t = 0.3, r = 0.5, b = 0.3, l = 0.5, unit = "cm")
base_size = 120
legend.key.size = unit(2,'cm')
legend.key.height = unit(2.5,'cm')
legend.key.width = unit(4,'cm')
y_max <- 31
trans_function <- function(x) x/max(pert_df$Cont)/3
gb <-ggplot()
base<- gb+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  #annotate("text", x=w_dates[14], y=max(rep_df$R0*scaling_Rt), label="Schools & non-ess. shops reopen",size=annotations_size, angle=0)+
  annotate("text", x=w_dates[13]+5, y=y_max, vjust=2, label="Schools & non-ess. shops reopen",size=annotations_size, angle=0,family =fam)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[8]+15,y=y_max, vjust=2, label=  df_measures$Descr_short[8],size=annotations_size, angle=0,family =fam)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=w_dates[29],y=y_max, vjust=2, label= "School re-opening & Lifting",size=annotations_size, angle=0,family =fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[10]+18,y=y_max, vjust=2, label= df_measures$Descr_short[10] ,size=annotations_size, angle=0,family =fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[11]+20,y=y_max, vjust=2, label= "Lifting" ,size=annotations_size, angle=0,family =fam)


#Stringency index
#annotate("text", x=w_dates[24], y=0.53, label="Stringency index(%)",size=6, angle=0)+

#VOCs
y_rect <- 1
y_min <- -1
y_text <- 0
base <- base+
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = y_min, ymax =y_rect)+
  annotate("text", x=w_dates[12], y=y_text, label= "Wild Type",size=annotations_size,family =fam)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=w_dates[20], y=y_text, label= "Alpha",size=annotations_size,family =fam)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=y_text, label= "Delta",size=annotations_size,family =fam)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=y_text, label= "Omicron",size=annotations_size,family =fam)+
  theme(
    plot.margin=plot.margin
  )


cont_plot <- base+
  geom_line(data =pert_df %>% filter(Wave>8),aes(x=Date,y=Cont,color=Age), lwd=3.4)+
  geom_line(data =stringency %>% filter(date>=w_dates[9] & date<=w_dates[length(w_dates)]),aes(x=date,y=StringencyIndex_Average/100*2*max(pert_df$Cont), linetype="S.I.(%)"), lwd=1.8, color='black')+
  scale_x_date(breaks=w_dates[9:length(w_dates)],labels = xaxis_el,   limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0, 0),sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave"))+
  scale_linetype_manual(values=c("twodash"))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names)+
  scale_y_continuous(expand = c(0, 0),sec.axis=sec_axis(trans_function,name = "",
                                                        breaks = seq(0, 1, by = 0.2)))+
  theme_bw(base_size = base_size) +
  theme(#axis.line = element_line(colour ="black"),
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.spacing.y = unit(2.5,'cm'),
    axis.text.x.bottom = element_text(size=size_text, angle=90, vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text, family =fam),
    legend.text  = element_text(size=size_text, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.1, family =fam),
    #axis.text.x = element_blank(), #hide bottom x-axis labels
    axis.title.x.top = element_text(size=size_text, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  guides(
    linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
    color = guide_legend(order = 2)     # 'color' represents the Age group legend
  ) +
  labs(y=bquote("Daily contacts"),
       x="",linetype='',
       color='Age')

ggsave("Contacts_evo.png", plot = cont_plot, width = 48, height = 27, dpi = 220)

#FIG S14: Cumulative sensitivity----
#BACKGROUND PLOT
scaling_Rt <-max(pert_df$Sens)*1.15
pert_df$Age <- factor(pert_df$Age,levels =age_names)

xaxis_el <-format(as.Date(w_dates[9:length(w_dates)]),  "%d %b %y")
annotations_size <- 34
size_text <- 100
plot.margin = margin(t = 0.3, r = 0.3, b = 0.3, l = 0.5, unit = "cm")
base_size =45
legend.key.size = unit(2,'cm')
legend.key.height = unit(1.5,'cm')
legend.key.width = unit(4,'cm')
y_max <-3.5
gb <-ggplot()
base<- gb+
  
  #School & non-essential shope reopenings
  annotate("rect", fill = "darkgrey", alpha = 0.3, 
           xmin = w_dates[9], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  #annotate("text", x=w_dates[14], y=max(rep_df$R0*scaling_Rt), label="Schools & non-ess. shops reopen",size=annotations_size, angle=0)+
  annotate("text", x=w_dates[13]+5, y=y_max, vjust=2, label="Schools & non-ess. shops reopen",size=annotations_size, angle=0,family =fam)+
  
  #Easter pause
  annotate("rect", fill = "darkgrey", alpha = 0.5, 
           xmin = df_measures$Date_start[8], xmax = df_measures$Date_end[8],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[8]+15,y=y_max, vjust=2, label=  df_measures$Descr_short[8],size=annotations_size, angle=0,family =fam)+
  
  #School re-openings and general relaxation measures
  annotate("rect", fill = "darkgrey", alpha = 0.20, 
           xmin = df_measures$Date_start[9], xmax = df_measures$Date_end[9],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=w_dates[29],y=y_max, vjust=2, label= "School re-opening & Lifting",size=annotations_size, angle=0,family =fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.10, 
           xmin = df_measures$Date_start[10], xmax = df_measures$Date_end[10],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[10]+18,y=y_max, vjust=2, label= df_measures$Descr_short[10] ,size=annotations_size, angle=0,family =fam)+
  annotate("rect", fill = "darkgrey", alpha = 0.05, 
           xmin = df_measures$Date_start[11], xmax = w_dates[length(w_dates)],
           ymin = -Inf, ymax = y_max)+
  annotate("text", x=df_measures$Date_start[11]+20,y=y_max, vjust=2, label= "Lifting" ,size=annotations_size, angle=0,family =fam)



#Stringency index
#annotate("text", x=w_dates[24], y=0.53, label="Stringency index(%)",size=6, angle=0)+

#VOCs
y_rect <- 0
y_min <- -0.3
y_text <- -0.15
base <- base+
  annotate("rect", fill = pal[1], alpha = 1,
           xmin = w_dates[9], xmax = df_vocs$Date_80_start[2],
           ymin = y_min, ymax =y_rect)+
  annotate("text", x=w_dates[12], y=y_text, label= "Wild Type",size=annotations_size,family =fam)+
  
  
  annotate("rect", fill = pal[2], alpha = 1,
           xmin = df_vocs$Date_start[2], xmax = df_vocs$Date_80_start[3],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=w_dates[20], y=y_text, label= "Alpha",size=annotations_size,family =fam)+
  
  
  annotate("rect", fill = pal[3], alpha = 1,
           xmin = df_vocs$Date_start[3], xmax = df_vocs$Date_80_start[4],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_80_start[3]+90, y=y_text, label= "Delta",size=annotations_size,family =fam)+
  
  
  annotate("rect", fill = pal[4], alpha = 1,
           xmin = df_vocs$Date_start[4], xmax = w_dates[length(w_dates)],
           ymin = y_min, ymax = y_rect)+
  annotate("text", x=df_vocs$Date_start[4]+28, y=y_text, label= "Omicron",size=annotations_size,family =fam)+
  theme(
    plot.margin=plot.margin
  )

##Minors----
minor_sens <- base+geom_line(data =pert_df %>% filter(Wave>8 & AgeGroup == 'Minors'),aes(x=Date,Sens,color=Age), lwd=3)+
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt,linetype="Rt"), lwd=1.2)+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 1, lwd=1.2, linetype="longdash")+
  scale_x_date( breaks=w_dates[9:length(w_dates)],
                #labels = rep('', length(w_dates[9:length(w_dates)])),
                labels = xaxis_el, 
                limits=c(w_dates[9],w_dates[length(w_dates)]),
                expand = c(0, 0),
                sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave")
  )+
  scale_linetype_manual(values=c('solid'))+
  scale_y_continuous(expand = c(0, 0),sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names,drop = FALSE)+
  theme_bw(base_size = base_size) +
  theme(#axis.line = element_line(colour ="black"),
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.spacing.y = unit(1, 'cm'),
    axis.text.x.bottom = element_text(size=size_text, angle=90, vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text*1.4, family =fam),
    legend.text  = element_text(size=size_text*1.25, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.3, family =fam),
    #axis.text.x = element_blank(), #hide bottom x-axis labels
    axis.title.x.top = element_text(size=size_text*1.2, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  guides(
    linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
    color = guide_legend(order = 2)     # 'color' represents the Age group legend
  ) +
  labs(y=bquote("Sensitivity"~tilde(s)[j]),
       x="",linetype='',
       color='Age') 

#ggsave("Sens_minors.png", plot = minor_sens, width = 40, height = 20, dpi = 200)


##Adults----

adults_sens <- base+geom_line(data =pert_df %>% filter(Wave>8 & AgeGroup != 'Minors'),aes(x=Date,y=Sens,color=Age), lwd=3)+
  
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt), lwd=1.2)+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 1, lwd=1.2, linetype="longdash")+
  scale_x_date(breaks=w_dates[9:length(w_dates)],
               labels = rep('', length(w_dates[9:length(w_dates)])),
               #labels = xaxis_el,
               limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0,0),
               sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates), name = "Comix Wave")
  )+
  #scale_linetype_manual(values=c('solid'))+
  scale_y_continuous( expand = c(0,0), sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names,drop = FALSE)+
  theme_bw(base_size = base_size) +
  theme(#axis.line = element_line(colour ="black"),
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.spacing.y = unit(1, 'cm'),
    axis.text.x.bottom = element_text(size=size_text, angle=45, vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text*1.4, family =fam),
    legend.text  = element_text(size=size_text*1.25, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.3, family =fam),
    #axis.text.x = element_blank(), #hide bottom x-axis labels
    axis.title.x.top = element_text(size=size_text*1.2, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  labs(y=bquote("Sensitivity"~tilde(s)[j]),
       x="",linetype='',
       color='') 

#ggsave(path = path_graphs ,"Adults_col.png", plot = adults_col, width = 40, height = 15, dpi = 100)

# Combine the plots
combined_sens <- ggarrange(adults_sens, minor_sens,  ncol = 1, nrow = 2,common.legend = T, legend = "right")

# Save the combined plot
ggsave("AgeGroupPlotSens.png", plot = combined_sens, width = 42, height =30, dpi =250,limitsize = FALSE)



#FIG S15: Inf.value----
#Use the same "base" of FIG s14
##Minors----
minor_inf <- base+geom_line(data =pert_df %>% filter(Wave>8 & AgeGroup == 'Minors'),aes(x=Date,Inf_val,color=Age), lwd=3)+
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt,linetype="Rt"), lwd=1.2)+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 1, lwd=1.2, linetype="longdash")+
  scale_x_date( breaks=w_dates[9:length(w_dates)],
                #labels = rep('', length(w_dates[9:length(w_dates)])),
                labels = xaxis_el, 
                limits=c(w_dates[9],w_dates[length(w_dates)]),
                expand = c(0, 0),
                sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates) , name = "CoMix Wave")
  )+
  scale_linetype_manual(values=c('solid'))+
  scale_y_continuous(expand = c(0, 0),sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names,drop = FALSE)+
  theme_bw(base_size = base_size) +
  theme(#axis.line = element_line(colour ="black"),
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.spacing.y = unit(1, 'cm'),
    axis.text.x.bottom = element_text(size=size_text, angle=90, vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text*1.4, family =fam),
    legend.text  = element_text(size=size_text*1.25, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.3, family =fam),
    #axis.text.x = element_blank(), #hide bottom x-axis labels
    axis.title.x.top = element_text(size=size_text*1.2, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  guides(
    linetype = guide_legend(order = 1),  # Adjust the order here, 'linetype' represents the Rt and S.I. legend
    color = guide_legend(order = 2)     # 'color' represents the Age group legend
  ) +
  labs(y=bquote("Infective value "~v[j]),
       x="",linetype='',
       color='Age') 

#ggsave("Inf_minors.png", plot = minor_inf, width = 40, height = 15, dpi = 100)


##Adults----

adults_inf <- base+geom_line(data =pert_df %>% filter(Wave>8 & AgeGroup != 'Minors'),aes(x=Date,y=Inf_val,color=Age), lwd=3)+
  
  geom_line(data = rep_df%>% filter(Wave>8),aes(x=Date,y=R0*scaling_Rt), lwd=1.2)+
  geom_hline(yintercept = 1*scaling_Rt, lwd=0.8)+
  geom_hline(yintercept = 1, lwd=1.2, linetype="longdash")+
  scale_x_date(breaks=w_dates[9:length(w_dates)],
               labels = rep('', length(w_dates[9:length(w_dates)])),
               #labels = xaxis_el,
               limits=c(w_dates[9],w_dates[length(w_dates)]),
               expand = c(0,0),
               sec.axis = sec_axis(~ . ,  breaks = w_dates[9:length(w_dates)], labels = 9:length(w_dates), name = "Comix Wave")
  )+
  #scale_linetype_manual(values=c('solid'))+
  scale_y_continuous( expand = c(0,0), sec.axis=sec_axis(~./scaling_Rt,breaks = seq(0,1.5, by=0.5)))+
  scale_color_manual("Age group",values=col_blind_pal,
                     labels=age_names,drop = FALSE)+
  theme_bw(base_size = base_size) +
  theme(#axis.line = element_line(colour ="black"),
    plot.margin=plot.margin,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.spacing.y = unit(1, 'cm'),
    axis.text.x.bottom = element_text(size=size_text, angle=45, vjust=0.5, family =fam),
    axis.text.x = element_text(size=size_text, angle=75, vjust=0.5, family =fam), 
    axis.text.y  = element_text(size=size_text, family =fam),
    legend.title = element_text(size=size_text*1.4, family =fam),
    legend.text  = element_text(size=size_text*1.25, family =fam),
    legend.key.size = legend.key.size,
    legend.key.height = legend.key.height,
    legend.key.width = legend.key.width,
    axis.title.y.left = element_text(size=size_text*1.3, family =fam),
    #axis.text.x = element_blank(), #hide bottom x-axis labels
    axis.title.x.top = element_text(size=size_text*1.2, family =fam),
    axis.text.x.top = element_text(size=size_text, vjust=0.5,angle=30, family = fam)
  ) +
  labs(y=bquote("Infective value "~v[j]),
       x="",linetype='',
       color='') 

# Combine the plots
combined_inf <- ggarrange(adults_inf, minor_inf,  ncol = 1, nrow = 2,common.legend = T, legend = "right")

# Save the combined plot
ggsave("AgeGroupPlotInf.png", plot = combined_inf, width = 42, height =30, dpi =250,limitsize = FALSE)

#

