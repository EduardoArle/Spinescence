## Mapping Extant Mammal Occurrences to Grid Cells ##

# Description: This script processes IUCN mammal range data, overlaps it with spatial grid polygons

# filters species based on native and extant criteria, and outputs species-grid matrices.

## Project: The complex interaction between megaherbivores, climate and fire has shaped the evolution and distribution of plant spinescence across biogeographical realms

## Author: Rachel Souza Ferreira (rachel.souza_ferreira@idiv.de) and Eduardo Arlé

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

wd_tables <- '/Users/carloseduardoaribeiro/Documents/Collaborations/Rachel/Data/Tables' # 

# ---------------------------- #
# 3. Load Data                 #
# ---------------------------- #

mmls <- terra::vect("redlist_species_data") # use function to split IUCN ranges


grids <- terra::vect("SpatialData")


plants <- read.csv("MimosoidsTrimmed.csv")

plants <- plants[,c(2:4)]#Keeping only coordinates + species id

plants_geom <- vect(plants, geom=c("decimalLongitude", "decimalLatitude"), crs="")# Creating SpatVector object from csv


# Correcting geometries ####
# necessary? ->  to check

test_mmls <- makeValid(mmls)
###############################################

######### Mammals -> overlapping polygons (mammals) with polygons (grids)
# Job 1 

mmls_cover <- terra::relate(test_mmls, grids, relation = "coveredby")

### done in cluster###
#mmls_over <- terra::relate(test_mmls, grids, relation = "overlaps")

mmls_over <- readRDS("mmls_overl_subset.rds")

write.csv(mmls_over,"pre_result_mmls_overlap_4october.csv")
write.csv(mmls_cover,"pre_result_mmls_covered_4october.csv")

mmls_valid <- is.valid(mmls)
which(mmls_valid == F)
test <- mmls[mmls_valid==T]
length(test$BINOMIAL)


test <- mmls[mmls_valid==T]
mmls_overl <- terra::relate(mmls[mmls_valid==TRUE], grids, relation="overlaps")

test <- cover
test[] <- lapply(test[],as.numeric)
test$X <- table$BINOMIAL

#rename columm in order to join dataframes ###

colnames(mmls_over) <- grids$WorldID
colnames(mmls_cover) <- grids$WorldID

anti_join(cover,over)

###### handling dataset####

over <- read.csv("pre_result_mmls_overlap_4october.csv")
cover <-read.csv("pre_result_mmls_covered_4october.csv")

#### join matrix of mammals "within" and "overlaps" #####

mammals_join <- join(over, cover, by="X")
df_merge1 <- merge(over, cover,by="row.names")

# Handling dataset ####
df <- read.csv("pre_result_mmls_covered_4october.csv")
table <- read.csv("IUCN_mammals.csv")

df1 <- df
df1[] <- lapply(df1[],as.numeric)
df1$X <- table$BINOMIAL
colnames((table))
all_mammals <- cbind(df1, table[,4:6])

# Filtering presence = "extant" (code 1) and origin = "native" (code 1) ####
ext_mammals <- filter(all_mammals, all_mammals[,19563] == 1)
ext_nat_mammals <- filter(ext_mammals, ext_mammals[,19564] == 1)


df1 <- ext_nat_mammals
df1 <- df1[,-c(19563:19565)]
df1$X <- make.unique(df1$X)
df2 <- t(df1)
df3 <- df2[-1,]  
colnames(df3) <- df2[1,]  
#write.csv(df3, "pre_result_mmls_grid_12September_within.csv")

# Merging columns of duplicated species ####
df3 <- apply(df3, 2 ,as.numeric)
colnames(df3) <- gsub(" ",".",colnames(df3))
species <- colnames(df3)
spu <- lapply(strsplit(species, "\\."), function(x){
  paste0(x[1],".",x[2])
})
spu <- unique(unlist(spu))



# create final data
df4 <- df3[,colnames(df3)%in%spu]
for (i in 1:length(spu)){
  df4[,spu[i]] <- rowSums(df3[,grep(pattern = spu[i], species), drop=F])
}
gd <- colnames(df1)[-1]
rownames(df4) <-gd[1:19561]
final <- df4
final[final[] > 1] <- 1
write.csv(final,"mmls_final_occurence_matrix_12September_within.csv")


# Getting list of grid per species ####
df <- read.csv("mmls_final_occurence_matrix_12September_within.csv")
df1 <- df[,-1]
rownames(df1) <- df[,1]
remove(df2)

df2 <- t(df1)
df2 <- as.data.frame(df2)
df2$GRID <- NA

for(i in 1:nrow(df2)){
  index <- which(df2[i,]==1)
  gd <- colnames(df2)[index]
  df2$GRID[i] <- paste(gd, collapse=",")
} 

df3 <- df2
df3$Taxon <- rownames(df2)
df4 <- df3[,c(19562,19563)]
test <- t(df4)


df5 <- df4 %>% 
  mutate(GRID = strsplit(as.character(GRID), ",")) %>% 
  unnest(GRID)

#write.csv(df5, "mmls_final_bc_list_03May2022_within")

###############################################
colnames(over) <- grids$WorldID
over <-read.csv("pre_result_mmls_overlap_4october.csv")

# Handling dataset ####
df <- read.csv("pre_result_mmls_overlap_4october.csv")
table <- read.csv("IUCN_mammals.csv")
####delete rows of invalid geometry from IUCN

table2 <- table
iucn_valid <- table2[-c(3326, 8075, 8957, 10857, 11015), ]

df1 <- over
df1[] <- lapply(df1[],as.numeric)
df1$X <- iucn_valid$BINOMIAL
colnames((iucn_valid))
all_mammals <- cbind(df1, iucn_valid[,4:6])

# Filtering presence = "extant" (code 1) and origin = "native" (code 1) ####
ext_mammals <- filter(all_mammals, all_mammals[,19563] == 1)
ext_nat_mammals <- filter(ext_mammals, ext_mammals[,19564] == 1)


df1 <- ext_nat_mammals
df1 <- df1[,-c(19563:19565)]
df1$X <- make.unique(df1$X)
df2 <- t(df1)
df3 <- df2[-1,]  
colnames(df3) <- df2[1,]  
#write.csv(df3, "pre_result_mmls_grid_12September_within.csv")

# Merging columns of duplicated species ####
df3 <- apply(df3, 2 ,as.numeric)
colnames(df3) <- gsub(" ",".",colnames(df3))
species <- colnames(df3)
spu <- lapply(strsplit(species, "\\."), function(x){
  paste0(x[1],".",x[2])
})
spu <- unique(unlist(spu))



# create final data
df4 <- df3[,colnames(df3)%in%spu]
for (i in 1:length(spu)){
  df4[,spu[i]] <- rowSums(df3[,grep(pattern = spu[i], species), drop=F])
}
gd <- colnames(df1)[-1]
rownames(df4) <-gd[1:19561]
final <- df4
final[final[] > 1] <- 1
write.csv(final,"mmls_final_occurence_matrix_12September_over.csv")


# Getting list grids per species ####
df <- read.csv("mmls_final_occurence_matrix_12September_over.csv")
df1 <- df[,-1]
rownames(df1) <- df[,1]

df2 <- t(df1)
df2 <- as.data.frame(df2)
df2$GRID <- NA

for(i in 1:nrow(df2)){
  index <- which(df2[i,]==1)
  gd <- colnames(df2)[index]
  df2$GRID[i] <- paste(gd, collapse=",")
} 

df3 <- df2
df3$Taxon <- rownames(df2)
df4 <- df3[,c(19562,19563)]
test <- t(df4)


df5 <- df4 %>% 
  mutate(GRID = strsplit(as.character(GRID), ",")) %>% 
  unnest(GRID)

write.csv(df5, "mmls_final_grid_over.csv")
to <- read.csv ("mmls_final_grid_over.csv")
tc <- read.csv ("mmls_final_grid_list_04_october_cover.csv")

###fim chique################

final <- rbind(to,tc)
final <- final[,-1]
final2 <- final %>% distinct()
write.csv(final2,"mammals_final_grid_within&overlaps_07Oct2022")

df <- read.csv("mammals_final_grid_within&overlaps_01_November_2022.csv", sep = ";")

IUCN <- merge(df,diet1)
df$Binomial.1.2 <- str_replace('.','_',df$Binomial.1.2)
write.csv(IUCN, "IUCN_final_grid_1_November.csv")



