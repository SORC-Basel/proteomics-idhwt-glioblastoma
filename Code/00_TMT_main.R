
# ++++++++++++++++++++++++
# Load the packages
# ++++++++++++++++++++++++
## Loading libraries
library("limma")
library("stringr")
library("dplyr")
library("ggrepel")
library("MSnbase")
library("RColorBrewer")
library("tidyr")
library("clusterProfiler")
library("org.Hs.eg.db")
library("flextable")
library("patchwork")
library("ggplot2")
library("cowplot")
library("stringr")
library("enrichplot")
library("ComplexHeatmap")
library("FactoMineR")
library("factoextra")
library("corrplot")
library("openxlsx") 


# Directory
path_directory <- here::here()

### Load datasets
## TMT 
dta_full_tmt     <- read.csv(paste0(path_directory,"Data/Total_proteome_normalized.CSV"), stringsAsFactors = F)

## Confounders
dta_confounder   <- readxl::read_xlsx(paste0(path_directory,"Data/confounders_data.xlsx"))



## --------------------------------------------------------------------------------------------

## Load source files 
source(paste0(path_directory,"02_TMT_functions.R"))
source(paste0(path_directory,"02_TMT_global_variables.R"))

# Data preparation
source(paste0(path_directory,"01_TMT_datapreparation.R"))    


# DEP
source(paste0(path_directory,"03_TMT_analysis_low_vs_rest.R"))    


## Figures 
# Barplots
source(paste0(path_directory,"04_TMT_barplot.R"))    



