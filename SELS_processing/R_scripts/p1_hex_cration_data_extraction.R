#-----------------------
library(tidyverse)
library(sp)
library(raster)
library(sf)
library(tmap)
library(stars)
library(terra)
library(dplyr)
library(exactextractr)
getwd()
input_folder <-setwd("E:/users/input_data") #run sometimes twice

list.files(input_folder)

#### FIRST SECTION: RASTER    ####
#Review list of files so recognize the ones that are "special".
#---------------------------------------------------------------#
#       IMPORTANT:
# Be familiar with your data: some raster attributes 
# can be automatically extracted, while others need more work 
# (ex: MapBiomas "Buff_reclass0602_v2022.tif").
#---------------------------------------------------------------#

#Selecting TIF (Raster) files that are "normal"
# List only TIFF files in the folder
tif_files <- list.files(input_folder, pattern = "\\.tif$", full.names = TRUE)

# Extract only the file names without the full path
tif_file_names <- basename(tif_files)

# Print the list of TIFF file names
print(tif_file_names)

# Define exclusion list
#Lnd Cover, MapBiomas: "MapBiomas_and_SoyMap_v2_clipped_reclass_100mt"
#Burned pixels from MODIS: "burnedIM.tif"
#Pop. density: "GHS_POP_E2020_merged.tif"
exclusion_list <- c("reclass", "burned","Cattle", "GHS_POP_E2020_merged",
                    "mapbiomas",  "size", "TAXN","access"
                    )  # Add more names as needed

# Exclude files that include names from the exclusion list
exclude_pattern <- paste(exclusion_list, collapse = "|")
tif_file_names_sel <- tif_file_names[!grepl(exclude_pattern, tif_file_names)]

print(tif_file_names_sel)

#### READ: Study Area ####
study_area <- st_read(dsn=".", layer="AltoParana_adm_bound_50km_50perc_v0802" )
sa_terra <- terra::vect(study_area) 
# Get the extent from the sa_terra object
sa_terra_extent <- ext(sa_terra)
sa_t_4326 <- project(sa_terra, "epsg:4326") 
st_crs(sa_t_4326)
# Reproject the shapefile to EPSG:102033
# Define the target CRS (EPSG:3857)
target_crs <- st_crs(3857)
sa_t_3857 <- st_transform(study_area, target_crs)
st_crs(sa_t_3857)

# Plot the polygon
plot(sa_t_3857["AREA"])

#### CREATE HEXAGONAL GRIDS #### 

hex_grid<-st_make_grid(sa_t_3857,
                       cellsize = 5000, #km
                       square = FALSE, what = "polygons")

# To sf and add grid ID
hex_grid_sf_v2 = st_sf(hex_grid)
#plot(hex_grid)
#summary(hex_grid_sf) #takes time when hex width = 1000km


# Clip polygon1 with polygon2
{start_time <- Sys.time() #set timer
  sa_hex_3857 <- st_intersection(hex_grid_sf_v2,sa_t_3857)
  #x11()
  #plot(sa_hex_3857["grid_id"])
  end_time <- Sys.time();end_time - start_time} 
#plot(sa_hex_3857["burned"])
# Add a unique grid ID column
sa_hex_3857$grid_id <- 1:nrow(sa_hex_3857)
sa_hex_3857$AREA <- st_area(sa_hex_3857)

#Clean possible points or multipolygons:
# Filter out non-point geometries
pol_hex_3857 <- sa_hex_3857[sf::st_geometry_type(sa_hex_3857) != "POINT", ]
# Print summary of the resulting layer
print(summary(pol_hex_3857))

#In case you have multipolygons
# Filter out multipolygons
multipolygons <- pol_hex_3857[sf::st_geometry_type(pol_hex_3857) == "MULTIPOLYGON", ]
# Convert multipolygons to individual polygons
individual_polygons <- st_cast(multipolygons, "POLYGON")
# Print summary of the resulting layer
print(summary(individual_polygons))
# Compute areas
individual_polygons$AREA <- st_area(individual_polygons)
# Print summary of areas
print(summary(individual_polygons$AREA))


#--------------------------------------------------------------------#
# Load raster data - Burned Pixels in the past ten years
# Get the list of TIFF files that start with "burned"
burned_tif_files <- list.files(pattern = "^burned.*\\.tif$", full.names = TRUE)
# Extract the first band from each raster
first_bands <- lapply(burned_tif_files, function(file) {
  raster_obj <- raster(file)
  return(raster_obj)
})
# Combine the extracted bands into a raster stack
raster_stack <- stack(first_bands)
# Write the raster stack to a GeoTIFF file
#writeRaster(raster_stack, "raster_stack_v2.tif", format = "GTiff", overwrite = TRUE)

{start_time <- Sys.time() #set timer
  sum_raster <- sum(raster_stack)
  end_time <- Sys.time();end_time - start_time}#1.269911 mins
# Plot the mean raster
plot(sum_raster)

mean_raster <- sum_raster/10

#We have decided to leave it as the sum of the average pixels burned in the last 10 years.
#This makes the spixels that burned every year a value of 1 and the ones that never burned are NA.
{start_time <- Sys.time() #set timer
  pol_hex_3857$burned<- exact_extract(mean_raster,  pol_hex_3857, 'sum')
  #pol_hex_3857$burned<- exact_extract(mean_raster, pol_hex_3857, 'sum') #number of pixels burned per hexagon
  end_time <- Sys.time();end_time - start_time} #4.486882 secs

#st_write(pol_hex_3857  , dsn =input_folder, 
#         layer="hex_test_burned_4.shp", driver="ESRI Shapefile", delete_layer = TRUE)
#
#-----------------------------------------------------------# 
# RASTER: 
#Note: it is assumed that your dataset are in the same projections or 
# that you have opened the datasets and these are all placed correctly.

#### Extract raster data ####
tif_file_names_sel

# Loop through each TIFF file name
{start_time <- Sys.time() #set timer
  for (file_name in tif_file_names_sel) {
    # Read the raster file
    my_tif <- raster(file_name)
    
    # Assign column name based on file_name
    column_name <- gsub(".tif", "", file_name)
    column_name <- gsub(".", "_", column_name, fixed = TRUE)  # Replace '.' with '_' to make it a valid column name
    column_name <- make.names(column_name, unique = TRUE)  # Ensure the column name is unique
    
    # Print the current file name
    print(file_name)
    
    # If the raster has multiple bands, extract mean from each band separately
    if (nlayers(my_tif) > 1) {
      # Extract mean from the first three bands
      means <- sapply(1:min(3, nlayers(my_tif)), function(i) {
        exact_extract(my_tif[[i]], pol_hex_3857, 'mean')
      })
      
      # Calculate the overall mean for all bands
      overall_mean <- rowMeans(means)
      
      # Assign the overall mean to the column
      pol_hex_3857[[column_name]] <- overall_mean
    } else {
      # If the raster has only one band, extract mean directly
      pol_hex_3857[[column_name]] <- exact_extract(my_tif, pol_hex_3857, 'mean')
    }
  }
  end_time <- Sys.time();end_time - start_time} #

#Making copy:
pol_hex_3857_copy1<-pol_hex_3857
#--------------------------------------------------------------------#
# Load raster data - time travel to cities
print(tif_file_names)
timetravel <- raster(tif_file_names[grep("access", tif_file_names)])
# Define the reclassification matrix
reclass_matrix <- matrix(c(-Inf, -9999, NA), ncol = 3, byrow = TRUE)
# Reclassify the raster
reclassified_raster <- reclassify(timetravel, reclass_matrix)
# Plot the reclassified raster to verify the changes
plot(reclassified_raster)
{start_time <- Sys.time() #set timer
  pol_hex_3857$timetravel<- exact_extract(reclassified_raster, pol_hex_3857, 'mean') 
  end_time <- Sys.time();end_time - start_time}
#--------------------------------------------------------------------#
# Load raster data - population
population <- raster(tif_file_names[grep("POP", tif_file_names)])
#plot(population)
{start_time <- Sys.time() #set timer
  pol_hex_3857$pop<- exact_extract(population, pol_hex_3857, 'sum') #takes some m
  end_time <- Sys.time();end_time - start_time} #

#--------------------------------------------------------------------#
# Load raster data - Field size from crowd sourcing data.
print(tif_file_names)
field_size <- raster(tif_file_names[grep("size", tif_file_names)])
plot(field_size)

{start_time <- Sys.time() #set timer
  pol_hex_3857$field_size<- exact_extract(field_size, pol_hex_3857, 'mode') #category of field size that repeats the most in the hexagon
  end_time <- Sys.time();end_time - start_time} #2.5 secs

#Reclass classes:
# 0 - no fields;
# 3502 -Very large fields with an area of greater than 100 ha;
# 3503 - Large fields with an area between 16 ha and 100 ha;
# 3504 - Medium fields with an area between 2.56 ha and 16 ha;
# 3505 - Small fields with an area between 0.64 ha and 2.56 ha; and
# 3506 - Very small fields with an area less than 0.64 ha.
# Reclassify values in the data frame: max, min and there rest are averages
#pol_hex_3857$field_size[pol_hex_3857$field_size == 0] <- 0
pol_hex_3857$field_size[pol_hex_3857$field_size == 3502] <- 100
pol_hex_3857$field_size[pol_hex_3857$field_size == 3503] <- 58
pol_hex_3857$field_size[pol_hex_3857$field_size == 3504] <- 9.28
pol_hex_3857$field_size[pol_hex_3857$field_size == 3505] <- 1.6
pol_hex_3857$field_size[pol_hex_3857$field_size == 3506] <- 0.64

#--------------------------------------------------------------------#
#Total cattle production by hexagon according to the gridded Livestock of
#the World 2.0 by Livestock Geowiki (Robinson et al. 2014)
print(tif_file_names)
cattle <- raster(tif_file_names[grep("Cattle", tif_file_names)])
plot(cattle)
{start_time <- Sys.time() #set timer
  pol_hex_3857$cattle<- exact_extract(cattle, pol_hex_3857, 'sum') #number of cattle heads per hexagon
  end_time <- Sys.time();end_time - start_time} #

#Making copy:
pol_hex_3857_copy<-pol_hex_3857
#st_write(pol_hex_3857, dsn =input_folder, 
#         layer="hex_cattle_v22032024_v1.shp", driver="ESRI Shapefile")
#--------------------------------------------------------------------#
#Plant Species diversity:
# Construct the full path to the shapefile
shapefile_path <- file.path(input_folder, "sr_Ensemble_Prediction_7774_Eckert-IV.shp" )
# Read the shapefile
shapefile_data <- st_read(shapefile_path)
print(shapefile_data)
plantsrich_3857 <- st_transform(shapefile_data, target_crs)
st_crs(plantsrich_3857)

# Perform spatial join to find intersections
joined_data <- st_join(pol_hex_3857, plantsrich_3857)

# Group by hexagon and calculate average value
average_values <- aggregate(joined_data$value, by = list(grid_id = joined_data$grid_id), FUN = mean)
# Merge the data frames by "grid_id"
pol_hex_3857 <- merge(pol_hex_3857, average_values, by = "grid_id", all = TRUE)
# Rename the "value" column to "new_value"
pol_hex_3857 <- pol_hex_3857 %>% 
  rename(plantsrich = x)

#Making copy:
pol_hex_3857_copy2<-pol_hex_3857
#--------------------------------------------------------------------#
#Indigenous Communities:
# Construct the full path to the shapefile
indigenous_path <- file.path(input_folder, "indigenous_villages_apaf_25032024.shp" )
# Read the shapefile
indigenous_data <- st_read(indigenous_path)
# Add a new column filled with the value 1
indigenous_3857 <- st_transform(indigenous_data, target_crs)
# Add a new column filled with the value 1
indigenous_3857 <- indigenous_3857 %>%
  mutate(presence = 1)
# Remove the first 46 columns
indigenous_3857 <- indigenous_3857[, -c(1:48)]
# Print the resulting data
print(indigenous_3857)
#st_crs(indigenous_3857)
#plot(indigenous_3857)

# Perform spatial join to find intersections
indigenous_joined_data <- st_join(pol_hex_3857, indigenous_3857)
# Aggregate the "presence" column by "grid_id" and sum the values
ind_by_hex <- aggregate(presence ~ grid_id, data = indigenous_joined_data, FUN = sum)
# Merge the data frames by "grid_id"
pol_hex_3857 <- merge(pol_hex_3857, ind_by_hex, by = "grid_id", all.x = TRUE)
#in indi and agri replace 0 for NAs
pol_hex_3857$presence[is.na(pol_hex_3857$presence)] <- 0
#Making copy:
pol_hex_3857_copy3<-pol_hex_3857

#--------------------------------------------------------------------#
#Protected areas:
protectarea_path <- file.path(input_folder,"Intersect_Areas_countries_dissolved_10042024.shp")
# Read the shapefile
protectarea_data <- st_read(protectarea_path)
protectarea_3857 <- st_transform(protectarea_data, target_crs)
#intersected_polygons <- st_intersection(shape, shp_template)
{start_time <- Sys.time() #set timer
  int <-as_tibble(st_intersection(protectarea_3857, pol_hex_3857))
  int$areaPA <- st_area(int$geometry)
  tb_APbyHEX <- int %>%
    group_by(grid_id) %>%
    summarise(areaPA = sum(areaPA))
  end_time <- Sys.time();end_time - start_time} #31.02669 secs

pol_hex_3857 <-left_join(pol_hex_3857, tb_APbyHEX, by = 'grid_id')
pol_hex_3857$areaPA[is.na(pol_hex_3857$areaPA)]<-0
pol_hex_3857<- pol_hex_3857%>%mutate(areaHEX=st_area(geometry),
                                     propPA=areaPA/areaHEX)
#Making copy:
pol_hex_3857_copy4<-pol_hex_3857

#--------------------------------------------------------------------#
#### Extract Land Cover data ####
#MapBiomas: In this step, MapBiomas is complemented with Soy maps for South America in 2021.
#Pre processing together with a resampling happened in QGIS
# Load raster data
print(tif_file_names)
land_cover <- raster(tif_file_names[grep("reclass", tif_file_names)])
# Print the resolution
print(land_cover) #0.0008983332 (100mt) - before 

# Plot the raster
# plot(land_cover, main = "Raster")
# Plot the shapefile on top of the raster
# plot(sa_t_4326, add = TRUE, border = "red")

{start_time <- Sys.time() #set timer
  extracted_values <- exact_extract(land_cover, pol_hex_3857, 
                                    function(df)df %>% group_by(grid_id, value),
                                    summarize_df = TRUE, include_cols = 'grid_id')
  end_time <- Sys.time();end_time - start_time} 

extracted_LC_summary <- extracted_values %>%
  group_by(grid_id) %>%
  summarise(
    forest = sum(value == 1, na.rm = TRUE),
    pasture = sum(value == 15, na.rm = TRUE),
    productive = sum(value == 19, na.rm = TRUE), #productive temporary and perennial crops/plantations
    noveg = sum(value == 22, na.rm = TRUE), #Urban Area/Non vegetated area
    water = sum(value == 33, na.rm = TRUE)) #River, Lake and Ocean ->not used

extracted_LC_summary <- extracted_LC_summary %>%
  mutate(grid_id = as.integer(grid_id))

#join data with main table:
merged_data <- left_join(pol_hex_3857, extracted_LC_summary, by = "grid_id")
# Calculate the size of pixels
pixel_size <- prod(res(land_cover))

# Convert columns 6:8 to numeric if needed
# Assuming merged_data is an "sf" object
modified_data <- merged_data %>%
  mutate(across(forest:tea, ~as.numeric(as.character(.))))
#merged_data %>% select(forest:soy) <- lapply(merged_data %>% select(forest:soy), function(x) as.numeric(as.character(x)))
# Check for missing values and handle them if necessary
modified_data <- modified_data %>%
  mutate(across(forest:tea, ~replace_na(., 0)))
#merged_data[, 4:15][is.na(merged_data[, 4:15])] <- 0
# Calculate new area for polygons
modified_data$area <- as.numeric(st_area(modified_data))
# Use left_join function to join the shapefile with the data frame
modified_data <- modified_data %>%
  mutate(across(forest:tea)* pixel_size / (modified_data$area))
#---------------------------------------------------------------------#

# Land Cover diversity
#Select the classes of Land Use
classes<- extracted_LC_summary %>% dplyr::select(forest:tea)#select the classes columns
# Load required library
library(vegan)
data <- extracted_LC_summary %>%
  mutate(Shannon_Div = diversity(extracted_LC_summary %>% dplyr::select(forest:tea),index = "shannon")) %>%
  mutate(Simpson_Div = diversity(extracted_LC_summary %>% dplyr::select(forest:tea),index = "simpson")) %>%
  mutate(Inv_Simpson_Div = diversity(extracted_LC_summary %>% dplyr::select(forest:tea),index = "invsimpson"))#Inverse_Simpson_Diversity

total_selected_columns <- ncol(classes)

data$Simpson_Eve <- data$Inv_Simpson_Div/total_selected_columns #Simpson_Evenness 
#Replace columns with the proportions of hexagon covered with the class 

#join data with main table:
modified_data <- modified_data %>%
  left_join(data %>% dplyr::select(grid_id, Shannon_Div,Simpson_Div,Inv_Simpson_Div,Simpson_Eve), by = "grid_id")

#Making copy:
modified_data_copy<-modified_data

modified_data<-modified_data_copy
#--------------------------------------------------------------------#
#RENAME COLUMNS:
colnames(modified_data)

modified_data <- modified_data %>%
  rename(    
    timetravel = accessibility_to_cities_2015_v1_0,
    bio16 = CHELSA_bio16_1981.2010_V_2_1,
    iucn = Combined_SR_2022,
    distrivers = proximity_to_rivers_25032024.,
    slope = Slope_sa_4326,
    indi = presence,
    clay = isric_clay,
    wind = wind_speed
  )

st_write(modified_data, dsn =input_folder, 
         layer="hex_originalvars_5km.shp", driver="ESRI Shapefile")
