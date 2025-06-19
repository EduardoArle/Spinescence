## Calculating Proportion of Spinescent Mimosoid Species per Grid Cell ##

# Description: This script compiles trait data on mimosoid species defenses      from Ringelberg et al. (2023) with updates from Bruneau et al. (2024). 

#It calculates the proportion of spinescent species in each spatial grid cell      by dividing the number of spinescent species by the total number of mimosoid     species present.

# Species distribution data are based on quality-controlled digitized herbarium specimens from Ringelberg et al. (2023).

## Authors: Rachel Souza Ferreira (rachel.souza_ferreira@idiv.de) and Eduardo    Arlé 


# ---------------------------- #
# 1. Load Required Libraries   #
# ---------------------------- #

library(terra)
#detach(dplyr) 
#library(tidyr)
#library(plyr)
#library(stringr)


# ---------------------------- #
# 2. List Working Directories #
# ---------------------------- #

wd_tables <- '/Users/carloseduardoaribeiro/Documents/Collaborations/Rachel/Data/Tables'

# ---------------------------- #
# 3. Load Data                 #
# ---------------------------- #

setwd(wd_tables)
plants <- read.csv("MimosoidsTrimmed.csv") 