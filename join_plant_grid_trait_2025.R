# _____________________________________________________________________________
# ## find mimosoid in grig
# Description: 
## Project: The complex interaction between megaherbivores, climate and fire has shaped the evolution and distribution of plant spinescence across biogeographical realms
## Author: Rachel Souza Ferreira (rachel.souza_ferreira@idiv.de) 
# _____________________________________________________________________________



setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Git/Data")
library(terra)
library(dplyr)

grids <- terra::vect("SpatialData")
plants <- read.csv("MimosoidsTrimmed.csv")
plants <- plants[,c(2:4)]#Keeping only coordinates + species id
plants_geom <- vect(plants, geom=c("decimalLongitude", "decimalLatitude"), crs="")# Creating SpatVector object from csv


######### Plants -> overlapping points (plants) with polygons (grids)
# Job 1 
#Overlaping data - within / plants ####
within_pl <- relate(plants_geom, grids, "within")
which(within_pl==1)
within_pl <- relate(plants_geom, grids[2390], "within")
colnames(within_pl) <- grids$WorldID
write.csv(within_pl,"pre_result_plants_grid_12September_within.csv")

#Overlaping data - touches / plants ####
touches_pl <- relate(plants_geom, grids, "touches")
which(touches_pl==1)
touches_pl <- relate(plants_geom, grids, "touches")
colnames(within_pl) <- grids$WorldID
write.csv(touches_pl,"pre_result_plant_grid_7October_touches.csv")


####Load matrix of plant within and touching grids after Terra relate function 

pw <- readRDS(file="final_plant_within_7November.rds")
pt <- readRDS("final_plant_touches_7November.rds")

###Joining matrices 

plant_relate <- rbind(pw,pt)
duplicated(plant_relate)
sum(duplicated(plant_relate))
plant_relate <- unique(plant_relate)
write.csv(plant_relate, "plant_join_within_touches.csv")

###Load matrix of Mimosoideae traits
trait_mat <- read.csv("Mimosoid_spinescence.csv", sep=",")
joined <- read.csv("plant_join_within_touches.csv", sep=",")

str(trait_mat)
str(joined)

joined2 <- joined %>%
  # replace every “.” with a space so it matches trait_mat$Taxon
  mutate(Taxon = gsub("\\.", " ", Taxon)) %>%
  # now join by the (space-separated) Taxon column
  left_join(trait_mat, by = "Taxon")

write.csv(joined2, "plants_grid_2025.csv")







