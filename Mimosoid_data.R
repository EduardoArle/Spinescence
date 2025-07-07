## Calculating Proportion of Spinescent Mimosoid Species per Grid Cell ##

# Description: This script compiles trait data on mimosoid species defenses      from Ringelberg et al. (2023) with updates from Bruneau et al. (2024). 

#It calculates the proportion of spinescent species in each spatial grid cell      by dividing the number of spinescent species by the total number of mimosoid     species present.

# Species distribution data are based on quality-controlled digitized herbarium specimens from Ringelberg et al. (2023).

# The data used in this script are available onhttps://zenodo.org/records/15773619

## Authors: Rachel Souza Ferreira (rachel.souza_ferreira@idiv.de) and Eduardo    Arlé 


# ---------------------------- #
# 1. Load Required Libraries   #
# ---------------------------- #

library(terra)
library(data.table)


# ---------------------------- #
# 2. List Working Directories #
# ---------------------------- #

#set working directory of the downloaded data
wd_download <- '/Users/carloseduardoaribeiro/Documents/Collaborations/Rachel/Data/Download/'
wd_results <- '/Users/carloseduardoaribeiro/Documents/Collaborations/Rachel/Data/Results'

# ---------------------------- #
# 3. Load Data                 #
# ---------------------------- #

#load mimosoids occurrences
setwd(paste0(wd_download, 'Data/Occurance'))
plants <- read.csv("Mimotrimmed_new.csv") 

#load mimosoids traits
setwd(paste0(wd_download, 'Data/Trait_Spinescence'))
spines <- read.csv('Mimosoid_spinescence.csv')

#load grid
setwd(paste0(wd_download, 'ELEData'))
grids <- terra::vect('SpatialData')


# ---------------------------- #
# 3. Overlap occ grid          #
# ---------------------------- #

#create SpatVector object from csv
plants_geom <- vect(plants, geom=c("decimalLongitude", "decimalLatitude"),
                    crs = crs(grids))

#list species
sps_list <- unique(plants_geom$SpeciesEdited)

#identify species within cells
sps_grid <- list()
for(i in 1:length(sps_list))
{
  sps_geom <- plants_geom[plants_geom$SpeciesEdited == sps_list[i],]
  within_pl <- relate(sps_geom, grids, "within")
  
  GRID <- numeric()
  for(j in 1:nrow(sps_geom))
  {
    a <- which(within_pl[j,] == TRUE)
    if(length(a) > 0){
      GRID[length(GRID) + 1] <- a
    }
  }
  Taxon <- gsub(' ', '.', sps_list[i])
  tab <- data.frame(GRID = GRID, Taxon = Taxon)
  tab2 <- unique(tab, by = c('GRID', 'Taxon'))
  sps_grid[[i]] <- tab2 
  
  print(i)
}

#rbind all species
sps_grid_all <- rbindlist(sps_grid)


#identify species touching cells 
sps_grid2 <- list()
for(i in 1:length(sps_list))
{
  sps_geom <- plants_geom[plants_geom$SpeciesEdited == sps_list[i],]
  touches_pl <- relate(sps_geom, grids, "touches")
  
  GRID <- numeric()
  for(j in 1:nrow(sps_geom))
  {
    a <- which(touches_pl[j,] == TRUE)
    if(length(a) > 0){
      GRID[length(GRID) + 1] <- a
    }
  }
  if(length(GRID) > 0){
    Taxon <- gsub(' ', '.', sps_list[i])
    tab <- data.frame(GRID = GRID, Taxon = Taxon)
    tab2 <- unique(tab, by = c('GRID', 'Taxon'))
    sps_grid2[[i]] <- tab2 
  }
  print(i)
}

#rbind all species
sps_grid2_all <- rbindlist(sps_grid2)

#join matrices
plant_relate <- rbind(sps_grid_all, sps_grid2_all)

#select unique combinatios of species per grid
plant_relate <- unique(plant_relate)

#save
setwd(wd_results)
write.csv(plant_relate, "plant_join_within_touches.csv")


# ---------------------------- #
# 3. Merge trait and occ data  #
# ---------------------------- #

#harmonise Taxon names
spines$Taxon <- gsub(' ', '.', spines$Taxon)

#merge data
merge_data <- merge(plant_relate, spines, by = "Taxon")

#save
setwd(wd_results)
write.csv(merge_data, "plants_grid_2025.csv")

