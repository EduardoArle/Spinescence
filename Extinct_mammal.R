## Phylacine Range Data Processing
## Description: Processes present and natural species ranges from the PHYLACINE dataset.
## Project: The complex interaction between megaherbivores, climate and fire has shaped the evolution and distribution of plant spinescence across biogeographical realms
## Author: Rachel Souza Ferreira (rachel.souza_ferreira@idiv.de)

### 1. Load Libraries -------------------------------------------------------
required_packages <- c("raster", "rgdal", "data.table", "terra", "rgeos", 
                       "raadtools", "exactextractr", "tidyverse")
lapply(required_packages, require, character.only = TRUE)
save.image(file='phylacine.RData')

setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Data")
files <- list.files(pattern="*.tif")
grids <- terra::vect("SpatialData")

#first i
#import all files in a single folder as a list 

present_natural1 <- list.files(path = "/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Data/Phylacine/Ranges/Present_natural", pattern='.tif$', all.files=TRUE, full.names=TRUE)
present_natural_stack <- stack(present_natural1)
wd_grid <- ("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Data/SpatialData")
shape <- readOGR("BehrmannMeterGrid_WGS84_land", dsn = wd_grid)

phyla_current <- list.files(path = "/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Data/Phylacine/Ranges/Current", pattern='.tif$', all.files=TRUE, full.names=TRUE)
present_current_stack <- stack(phyla_current)

plot(present_current_stack$Eudorcas_rufina)
plot(present_natural_stack$Eudorcas_rufina)

####convert shapefile into raster format

shape_cea <- spTransform(shape,CRS('+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'))


#######make loop of rasterstack to get list of presences in grid cell  #######

setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Outputs/current loop")


for(i in 1:nlayers(present_current_stack))
{
  a <- which(present_current_stack[[i]][] != 0)
  saveRDS(a,names(present_current_stack[[i]]))
  print(i)
}
#######make loop for each grid cell to identify all species in cell ###

setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Outputs/current loop")

sps <- lapply(list.files(), readRDS) 
names(sps) <- list.files()

setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Outputs/current in grid")

for(i in 1:length(present_current_stack[[1]]))
{
  a <- lapply(sps, function(x){i %in% x})
  b <- a[which(a == TRUE)]
  
  if(length(b) > 0){
    spsmm2 <- names(b)
    c <- data.frame(species = spsmm2)
    saveRDS(c, paste(i))
  }
  print(i) 
}

#matrix of all species 
setwd("/Users/rachelsouzaferreira/Dropbox/PhD/Chapter 1/Outputs/current in grid")

test <- lapply(list.files(), readRDS) 
names(test) <- list.files()

unlist_test <-  unlist(test)
phyllacine_grid <- as.matrix(unlist_test)

colnames(phyllacine_grid)[1] <- "Binomial.1.2"


##create loop to relate raster ID fall in same GRID ID/regio

###make ID raster without info, just cell number in shapefile region 

rasterID <- present_current_stack[[1]]
rasterID[] <- c(1:length(rasterID))
plot(rasterID)


region <- matrix(data=NA,ncol=2,nrow=nrow(shape_cea))

for (i in 1:nrow(shape_cea))
{
  a <- try(crop(rasterID,shape_cea[i,]),silent=T)
  b <- ifelse(class(a) != "try-error", values(a), NA)
  region[i,1] <- shape_cea[i,]$WorldID
  region[i,2] <- b
  print(i)
}

region2 <- as.data.frame(region)
names(region2) <- c("WorldID","RasterID")
tail(region2)


# Assuming phyllacine_grid is a matrix and region2 is a data frame
# Convert them into data.tables
phyllacine_dt <- as.data.table(phyllacine_grid, keep.rownames = "row_names")


# Assuming phyllacine_grid is your matrix and it's loaded into R
# Convert the matrix to a tibble
phyllacine_tb <- as_tibble(phyllacine_dt, .name_repair = "unique")

# Remove the '.species' part from the new 'ID' column using regular expression
phyllacine_tb$row_names <- gsub("\\.species.*", "", phyllacine_tb$row_names)

### change column name 

# Assuming phyllacine_tb is your tibble
# Let's say you want to change the name of the first column to 'ID' and the second to 'Species'
phyllacine_tb <- phyllacine_tb %>%
 rename_with(~ "RasterID", .cols = 1)

phyla <- merge(phyllacine_tb,region2)
phylacine_all_current = phyla[,-1]

write.csv(phylacine_all_current,"phylacine_all_Current.csv")


