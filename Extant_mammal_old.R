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

setwd(wd_tables)
plants <- read.csv("MimosoidsTrimmed.csv") #[, 2:4] ???

mmls <- terra::vect("redlist_species_data")
grids <- terra::vect("SpatialData") #include function to split IUCN shp per species
plants_geom <- vect(plants, geom = c("decimalLongitude", "decimalLatitude"), crs = "")

# ---------------------------- #
# 4. Validate Mammal Geometry  #
# ---------------------------- #

test_mmls <- makeValid(mmls)
mmls_valid <- is.valid(mmls)
valid_mmls <- mmls[mmls_valid == TRUE]
length(valid_mmls$BINOMIAL)

# ---------------------------- #
# 5. Spatial Relationships     #
# ---------------------------- #
mmls_cover <- terra::relate(test_mmls, grids, relation = "coveredby")
mmls_over <- readRDS("mmls_overl_subset.rds")

# Save Pre-Results
write.csv(mmls_over, "pre_result_mmls_overlap_4october.csv")
write.csv(mmls_cover, "pre_result_mmls_covered_4october.csv")

# ---------------------------- #
# 6. Load and Join IUCN Traits #
# ---------------------------- #
cover <- read.csv("pre_result_mmls_covered_4october.csv")
over <- read.csv("pre_result_mmls_overlap_4october.csv")
traits <- read.csv("IUCN_mammals.csv")

# ---------------------------- #
# 7. Preprocess and Filter     #
# ---------------------------- #
cover[] <- lapply(cover[], as.numeric)
cover$X <- traits$BINOMIAL
all_mammals <- cbind(cover, traits[, 4:6])

# Keep extant + native species
ext_mammals <- filter(all_mammals, all_mammals[, ncol(all_mammals) - 2] == 1)
ext_nat_mammals <- filter(ext_mammals, ext_mammals[, ncol(all_mammals) - 1] == 1)

# ---------------------------- #
# 8. Create Binary Matrix      #
# ---------------------------- #
df1 <- ext_nat_mammals[ , -c(ncol(ext_nat_mammals)-2:ncol(ext_nat_mammals))]
df1$X <- make.unique(df1$X)

# Transpose and format
df2 <- t(df1)
df3 <- df2[-1, ]
colnames(df3) <- df2[1, ]
df3 <- apply(df3, 2, as.numeric)
colnames(df3) <- gsub(" ", ".", colnames(df3))

# Merge duplicate species
species <- colnames(df3)
spu <- lapply(strsplit(species, "\\."), function(x) paste0(x[1], ".", x[2]))
spu <- unique(unlist(spu))
df4 <- df3[, colnames(df3) %in% spu]

for (i in 1:length(spu)) {
  df4[, spu[i]] <- rowSums(df3[, grep(pattern = spu[i], species), drop = FALSE])
}

rownames(df4) <- colnames(df1)[-1][1:nrow(df4)]
df4[df4 > 1] <- 1
write.csv(df4, "mmls_final_occurence_matrix_12September_within.csv")

# ---------------------------- #
# 9. Generate Grid List by Species #
# ---------------------------- #
df <- read.csv("mmls_final_occurence_matrix_12September_within.csv")
df1 <- df[ , -1]
rownames(df1) <- df[ , 1]
df2 <- as.data.frame(t(df1))
df2$GRID <- NA

for (i in 1:nrow(df2)) {
  index <- which(df2[i, ] == 1)
  gd <- colnames(df2)[index]
  df2$GRID[i] <- paste(gd, collapse = ",")
}

df2$Taxon <- rownames(df2)
df4 <- df2[ , c("Taxon", "GRID")]
df5 <- df4 %>%
  mutate(GRID = strsplit(as.character(GRID), ",")) %>%
  unnest(GRID)
# write.csv(df5, "mmls_final_bc_list_03May2022_within")

# ---------------------------- #
# 10. Repeat for 'Overlap'     #
# ---------------------------- #
over <- read.csv("pre_result_mmls_overlap_4october.csv")
traits <- read.csv("IUCN_mammals.csv")
traits <- traits[-c(3326, 8075, 8957, 10857, 11015), ]  # Remove invalid geometries

over[] <- lapply(over[], as.numeric)
over$X <- traits$BINOMIAL
all_mammals <- cbind(over, traits[, 4:6])

ext_mammals <- filter(all_mammals, all_mammals[, ncol(all_mammals) - 2] == 1)
ext_nat_mammals <- filter(ext_mammals, ext_mammals[, ncol(all_mammals) - 1] == 1)

df1 <- ext_nat_mammals[ , -c(ncol(ext_nat_mammals)-2:ncol(ext_nat_mammals))]
df1$X <- make.unique(df1$X)
df2 <- t(df1)
df3 <- df2[-1, ]
colnames(df3) <- df2[1, ]
df3 <- apply(df3, 2, as.numeric)
colnames(df3) <- gsub(" ", ".", colnames(df3))

species <- colnames(df3)
spu <- lapply(strsplit(species, "\\."), function(x) paste0(x[1], ".", x[2]))
spu <- unique(unlist(spu))
df4 <- df3[, colnames(df3) %in% spu]

for (i in 1:length(spu)) {
  df4[, spu[i]] <- rowSums(df3[, grep(pattern = spu[i], species), drop = FALSE])
}

rownames(df4) <- colnames(df1)[-1][1:nrow(df4)]
df4[df4 > 1] <- 1
write.csv(df4, "mmls_final_occurence_matrix_12September_over.csv")

# Grid list by species (overlap)
df <- read.csv("mmls_final_occurence_matrix_12September_over.csv")
df1 <- df[ , -1]
rownames(df1) <- df[ , 1]
df2 <- as.data.frame(t(df1))
df2$GRID <- NA

for (i in 1:nrow(df2)) {
  index <- which(df2[i, ] == 1)
  gd <- colnames(df2)[index]
  df2$GRID[i] <- paste(gd, collapse = ",")
}

df2$Taxon <- rownames(df2)
df4 <- df2[ , c("Taxon", "GRID")]
df5 <- df4 %>%
  mutate(GRID = strsplit(as.character(GRID), ",")) %>%
  unnest(GRID)

write.csv(df5, "mmls_final_grid_over.csv")

# ---------------------------- #
# 11. Combine & Export Final Matrix #
# ---------------------------- #
to <- read.csv("mmls_final_grid_over.csv")
tc <- read.csv("mmls_final_grid_list_04_october_cover.csv")
final <- rbind(to, tc)[ , -1]
final2 <- final %>% distinct()
write.csv(final2, "mammals_final_grid_within&overlaps_07Oct2022.csv")

# ---------------------------- #
# 12. Merge with Traits        #
# ---------------------------- #
df <- read.csv("mammals_final_grid_within&overlaps_01_November_2022.csv", sep = ";")
df$Binomial.1.2 <- str_replace(df$Binomial.1.2, "\\.", "_")
IUCN <- merge(df, diet1)  # Assumes 'diet1' is preloaded
write.csv(IUCN, "IUCN_final_grid_1_November.csv")


