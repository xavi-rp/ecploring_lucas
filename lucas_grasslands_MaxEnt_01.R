

library(data.table)
#library(tidyr)
#library(devtools)
#install_github("xavi-rp/PreSPickR", 
#               ref = "v2", 
#               INSTALL_opts = c("--no-multiarch"))  # https://github.com/rstudio/renv/issues/162
#library(PreSPickR)
library(raster)
library(giscoR)
library(sf)
library(tidyverse)



#sessionInfo()

if(Sys.info()[4] == "D01RI1700308") {
  wd <- "D:/xavi_rp/D5_FFGRCC_lucas_grasslands/"
}else if(Sys.info()[4] == "S-JRCIPRAP320P") {
  wd <- "D:/rotllxa/D5_FFGRCC_lucas_grasslands/"
}else if(Sys.info()[4] %in% c("jeodpp-terminal-jd001-03", "jeodpp-terminal-03", "jeodpp-terminal-dev-12" )) {
  if(!dir.exists("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/")) 
    dir.create("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/")
  wd <- "/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/"
  gbif_creds <- "/home/rotllxa/Documents/"
}else{
  wd <- "C:/Users/rotllxa/D5_FFGRCC_lucas_grasslands/"
  gbif_creds <- "C:/Users/rotllxa/Documents/"
}

setwd(wd)

# file.edit('~/.Renviron')


## Lucas grassland 2018 data set ####

## Harmonized by Momo
dir_LucasGrassland_Momo <- "/eos/jeodpp/data/projects/REFOCUS/data/flowerpower/data/LUCAS_grassland_data/scrap_output/estat_version_2/"
list.files(dir_LucasGrassland_Momo)

# Point geometries with all LUCAS Grassland attributes minus releve - 
point_geom_surveyors <- paste0(dir_LucasGrassland_Momo, "/estat_s_attr_point_allattr_new.csv")
# Line geometries with all LUCAS Grassland attributes minus releve - 
line_geom_surveyors <- paste0(dir_LucasGrassland_Momo, "/estatdb_allattr_s_transects.kml")
# Expert point geometries all LUCAS Grassland attributes minus releve - 
point_geom_experts <- paste0(dir_LucasGrassland_Momo, "/estat_e_attr_point_allattr_new.csv")


## Surveyors data set (points)
point_geom_surveyors <- fread(point_geom_surveyors)
point_geom_surveyors

nrow(point_geom_surveyors)  # 2622
ncol(point_geom_surveyors)  #  127
names(point_geom_surveyors)


sum(is.na(point_geom_surveyors$SURVEY_GRASS_GPS_LAT))  # no NAs
range(point_geom_surveyors$SURVEY_GRASS_GPS_LAT)       # 34.67436   66.26501
sum(is.na(point_geom_surveyors$SURVEY_GRASS_GPS_LON))  # no NAs
range(point_geom_surveyors$SURVEY_GRASS_GPS_LON)       # -9.746292 34.059950

sort(unique(point_geom_surveyors$SURVEY_GRASS_RICHNESS_SPEC27))
sum(is.na(point_geom_surveyors$SURVEY_GRASS_RICHNESS_SPEC27))
summary(point_geom_surveyors$SURVEY_GRASS_RICHNESS_SPEC27)
table(point_geom_surveyors$SURVEY_GRASS_RICHNESS_SPEC27)
# I assume that 0 means no individuals in the point; 1:10 means number of individuals; 11 means >10 individuals


sort(unique(point_geom_surveyors$NUMBER_KEY_SPECIES)) # 0-17; this should be species richness in the point

sort(unique(point_geom_surveyors$NUMBER_FLOWERS_KEY_SPECIES)) # 0-131; number of flowers of all key species (?)
summary(point_geom_surveyors$NUMBER_FLOWERS_KEY_SPECIES) 
#    Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
#   0.00    4.00    16.00    22.95   33.00   131.00

sort(unique(point_geom_surveyors$SURVEY_GRASS_FLOWER_DENSITY)) # 1-10; do not know what this is


## Plotting the points (surveyors) ####
library(sf)
library(ggplot2)
library(ggExtra)
library(viridis)  
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
# The “viridis” and “magma” scales do better - they cover a wide perceptual range in brightness in brightness and blue-yellow, 
# and do not rely as much on red-green contrast
library(ggforce)
library(ggpubr)
library(patchwork)
library(giscoR)
library(dplyr)


# Gisco maps
# https://ropengov.github.io/giscoR/

eur_gisco <- gisco_get_countries(region = "Europe")
eur_gisco

eur_gisco <- st_crop(eur_gisco, xmin = -10.5, xmax = 32, ymin = 33, ymax = 70)

p <- ggplot() +
  geom_sf(data = eur_gisco) +
  geom_point(
    data = point_geom_surveyors, 
    aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT),
    size = 0.1,
    color = "darkred"
  ) +
  
  theme_light() +
  #scale_color_viridis(option = "viridis", discrete = TRUE) +
  labs(title = "LUCAS grasslands 2018 (surveyors)") + #, x = "TY [°C]", y = "Txxx") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 16))# +
  #guides(color = guide_legend("Species", override.aes = list(size = 2)))

# https://jtr13.github.io/cc21fall2/tutorial-for-scatter-plot-with-marginal-distribution.html
p1 <- ggMarginal(p,
                 #aes(colour = species),
                 #type = "density", 
                 type = "histogram", 
                 #type = "densigram", 
                 color = "darkred",
                 groupColour = FALSE, groupFill = FALSE)
p1

ggsave("point_geom_surveyors.png", p1)#, width = 20, height = 20, units = "cm")





## Plotting by key species richness 

#point_geom_surveyors <- paste0(dir_LucasGrassland_Momo, "/estat_s_attr_point_allattr_new.csv")
#point_geom_surveyors <- fread(point_geom_surveyors)

#point_geom_surveyors$NUMBER_KEY_SPECIES <- as.numeric(point_geom_surveyors$NUMBER_KEY_SPECIES)
sort(unique(point_geom_surveyors$NUMBER_KEY_SPECIES))

point_geom_surveyors[, sp_richness_class := cut(point_geom_surveyors$NUMBER_KEY_SPECIES,
                                                #breaks = c(0, 1, 6, 17),
                                                breaks = c(0, 1, 8, 17),
                                                include.lowest = TRUE,
                                                right = FALSE)]
sort(unique(point_geom_surveyors$sp_richness_class))

point_geom_surveyors$NUMBER_KEY_SPECIES <- as.factor(point_geom_surveyors$NUMBER_KEY_SPECIES)


p <- ggplot() +
  geom_sf(data = eur_gisco) +
  geom_point(
    data = point_geom_surveyors, 
    aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT,
        color = sp_richness_class),
    #size = 0.1
    size = 0.4
  ) +
  
  theme_light() +
  scale_color_viridis(option = "viridis", discrete = TRUE) +
  labs(title = "LUCAS grasslands 2018 (surveyors): Key species richness") + #, x = "TY [°C]", y = "Txxx") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 16)) +
  guides(color = guide_legend("Key species richness classes", override.aes = list(size = 2)))

# https://jtr13.github.io/cc21fall2/tutorial-for-scatter-plot-with-marginal-distribution.html
p1 <- ggMarginal(p,
                 aes(colour = sp_richness_class),
                 type = "density", 
                 #type = "histogram", 
                 #type = "densigram", 
                 #color = "darkred",
                 groupColour = TRUE, groupFill = TRUE)
p1

ggsave("point_geom_surveyors_SpRichness.png", p1)#, width = 20, height = 20, units = "cm")





## Joining species name ####
# 






## Experts; normal survey ####

point_geom_experts #<- paste0(dir_LucasGrassland_Momo, "/estat_e_attr_point_allattr.csv")


## Experts data set (points)
point_geom_experts <- fread(point_geom_experts)

point_geom_experts

nrow(point_geom_experts)   # 605
names(point_geom_experts)  # 124+2

# EUNIS habitat classification
sort(unique(point_geom_experts$SURVEY_GRASS_EUNIS_HABITAT))

table(point_geom_experts$SURVEY_GRASS_RICHNESS_SPEC27)
# I assume that 0 means no individuals in the point; 1:10 means number of individuals; 11 means >10 individuals


sort(unique(point_geom_experts$NUMBER_KEY_SPECIES)) # 0-17; this should be species richness in the point
summary(point_geom_experts$NUMBER_FLOWERS_KEY_SPECIES) 
#    Min. 1st Qu.  Median    Mean   3rd Qu.    Max. 
#   0.00  17.00    37.50    41.31   59.00  155.00


## Plotting Experts points
# Gisco maps
# https://ropengov.github.io/giscoR/

eur_gisco <- gisco_get_countries(region = "Europe")
eur_gisco

eur_gisco <- st_crop(eur_gisco, xmin = -10.5, xmax = 40, ymin = 33, ymax = 70)

p <- ggplot() +
  geom_sf(data = eur_gisco) +
  geom_point(
    data = point_geom_experts, 
    aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT),
    size = 0.1,
    color = "darkblue"
  ) +
  
  theme_light() +
  #scale_color_viridis(option = "viridis", discrete = TRUE) +
  labs(title = "LUCAS grasslands 2018 (experts)") + #, x = "TY [°C]", y = "Txxx") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 16))# +
#guides(color = guide_legend("Species", override.aes = list(size = 2)))

# https://jtr13.github.io/cc21fall2/tutorial-for-scatter-plot-with-marginal-distribution.html
p1 <- ggMarginal(p,
                 #aes(colour = species),
                 #type = "density", 
                 type = "histogram", 
                 #type = "densigram", 
                 color = "darkblue",
                 groupColour = FALSE, groupFill = FALSE)
p1

ggsave("point_geom_experts.png", p1, width = 20, height = 20, units = "cm")





## Plotting by key species richness 

sort(unique(point_geom_experts$NUMBER_KEY_SPECIES))

point_geom_experts[, sp_richness_class := cut(point_geom_experts$NUMBER_KEY_SPECIES,
                                                #breaks = c(0, 1, 6, 17),
                                                breaks = c(0, 1, 8, 17),
                                                include.lowest = TRUE,
                                                right = FALSE)]
sort(unique(point_geom_experts$sp_richness_class))

point_geom_experts$NUMBER_KEY_SPECIES <- as.factor(point_geom_experts$NUMBER_KEY_SPECIES)


p <- ggplot() +
  geom_sf(data = eur_gisco) +
  geom_point(
    data = point_geom_experts, 
    aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT,
        color = sp_richness_class),
    #size = 0.1
    size = 0.4
  ) +
  
  theme_light() +
  scale_color_viridis(option = "viridis", discrete = TRUE) +
  labs(title = "LUCAS grasslands 2018 (experts): Key species richness") + #, x = "TY [°C]", y = "Txxx") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 16)) +
  guides(color = guide_legend("Key species richness classes", override.aes = list(size = 2)))

# https://jtr13.github.io/cc21fall2/tutorial-for-scatter-plot-with-marginal-distribution.html
p1 <- ggMarginal(p,
                 aes(colour = sp_richness_class),
                 type = "density", 
                 #type = "histogram", 
                 #type = "densigram", 
                 #color = "darkred",
                 groupColour = TRUE, groupFill = TRUE)
p1

ggsave("point_geom_experts_SpRichness.png", p1, width = 20, height = 20, units = "cm")






## Releve' ####

dir_LucasGrassland <- "/eos/jeodpp/data/projects/REFOCUS/data/flowerpower/data/LUCAS_grassland_data/scrap_output/estat_version/"
list.files(dir_LucasGrassland)

estatdb_allattr_e <- fread(paste0(dir_LucasGrassland, "/estatdb_allattr_e.csv"))   # This seems the experts normal survey
estatdb_allattr_e


estat_attr_s_start_loc_EW_fromB <- fread(paste0(dir_LucasGrassland, "/estat_attr_s_start_loc_EW_fromB.csv"))   # This seems surveyors coords start transect
estat_attr_s_start_loc_EW_fromB


estatdb_allattr_s_GPSplus_plusGPSEXIF <- fread(paste0(dir_LucasGrassland, "/estatdb_allattr_s_GPSplus_plusGPSEXIF.csv"))   
estatdb_allattr_s_GPSplus_plusGPSEXIF


estatdb_allattr_th <- fread(paste0(dir_LucasGrassland, "/estatdb_allattr_s_GPSplus_plusGPSEXIF.csv"))   # This seems surveyors coords start transect
estatdb_allattr_th




#
dir_LucasGrassland_2 <- "/eos/jeodpp/data/projects/REFOCUS/data/flowerpower/data/LUCAS_grassland_data/databases/"
list.files(dir_LucasGrassland_2)


LUCAS2018_DMT_EXP_v5.5 <- readxl::read_excel(paste0(dir_LucasGrassland_2, "/LUCAS2018_DMT_EXP_v5.5.xlsx"))
LUCAS2018_DMT_EXP_v5.5
names(LUCAS2018_DMT_EXP_v5.5)
names(LUCAS2018_DMT_EXP_v5.5)[235:275]   # key species


Annex2_DLV2.1_Database_v2 <- readxl::read_excel(paste0(dir_LucasGrassland_2, "/Annex2_DLV2.1_Database_v2.xlsx")) # all together
Annex2_DLV2.1_Database_v2
names(Annex2_DLV2.1_Database_v2)
#[235:275]   # key species


DLV5_lot9_surveydatabase_LUCAS2018_v3 <- readxl::read_excel(paste0(dir_LucasGrassland_2, "/DLV5_lot9_surveydatabase_LUCAS2018_v3.xlsx"))   # experts releve
DLV5_lot9_surveydatabase_LUCAS2018_v3
names(DLV5_lot9_surveydatabase_LUCAS2018_v3)
nrow(DLV5_lot9_surveydatabase_LUCAS2018_v3)
head(DLV5_lot9_surveydatabase_LUCAS2018_v3[, 1:3])

DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point <- readxl::read_excel(paste0(dir_LucasGrassland_2, "/DLV5_lot9_surveydatabase_LUCAS2018_v3.xlsx"),
                                                                         sheet = "LUCAS relevé data w point",
                                                                         col_names = TRUE)   # experts releve
DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point
DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point
nrow(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point) # one row per species, plus botanist name, species richness in the point and an empty row (3)
ncol(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point) # one col per point, plus species name, type (subspecies or genus) and number of points present
colnames(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point) # 
colnames(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point)[1:5] # 
tail(colnames(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point))

table(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point$`number of points present`)
DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point[DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point$`number of points present` == 138, 1:4]

## remove empty row
if(sum(is.na(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point[3, ])) == ncol(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point)) 
  releve_point <- DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point[-3, ]


## is "species richness" consistent with column reported? 
releve_point
head(names(releve_point))
head((releve_point))
releve_point[2, 1]

releve_point_data <- releve_point[3:nrow(releve_point), 4:ncol(releve_point)] # keeping only plants occurrences/abundances 
releve_point_data[, 1:4]
sort(unique(unlist(apply(releve_point_data, 2, function(x) sort(unique(x))))))

all(as.vector(apply(releve_point_data, 2, function(x) sum(!is.na(x)))) == as.numeric(unlist(releve_point[2, 4:ncol(releve_point)])))

kk <- data.frame(names(releve_point_data), as.vector(apply(releve_point_data, 2, function(x) sum(!is.na(x)))), as.numeric(unlist(releve_point[2, 4:ncol(releve_point)])))
sp_richn_inconsitent_points <- kk[kk$as.vector.apply.releve_point_data..2..function.x..sum..is.na.x.... != kk$as.numeric.unlist.releve_point.2..4.ncol.releve_point...., ]
sp_richn_inconsitent_points <- sp_richn_inconsitent_points$names.releve_point_data.
# This points are different: 31122070, 54902758. But only by 1

as.vector(releve_point[, colnames(releve_point) %in% sp_richn_inconsitent_points[1]])
as.vector(releve_point[, colnames(releve_point) %in% sp_richn_inconsitent_points[2]])


## is reported abundance consistent?
# The releve' reports % of coverage of each species within the transect, instead of number of plants
# How can this be used to calculate diversity indexes (e,g. Shannon)??
kk2 <- apply(releve_point_data, 2, function(x) sum(as.numeric(x), na.rm = TRUE))
kk2

unlist(releve_point_data[, 1])
kk3 <- as.numeric(unlist(releve_point_data[, 1]))
sum(as.numeric(unlist(releve_point_data[, 1])), na.rm = TRUE)
sort(kk3[!is.na(kk3)])

summary(kk2)   # much more than half of the points report more than 100% coberture (adding all species)!!
#   Min.  1st Qu.  Median    Mean    3rd Qu.    Max. 
#  2.60    93.22   115.80   127.13   156.60    422.00 

sum(kk2 == 100)
sum(kk2 < 100)
sum(kk2 > 100)


### list of species ####
releve_point

species_list_releve <- as.vector(unlist(releve_point[, 1]))
species_list_releve <- species_list_releve[!species_list_releve %in% c("Botanist", "Species richness")]
species_list_releve

length(species_list_releve) # 2672

species_list_releve <- gsub("\\.", " ", species_list_releve)
species_list_releve

species_list_releve <- gsub("  ", " ", species_list_releve)
species_list_releve

species_list_releve <- as.data.table(species_list_releve)
species_list_releve

write.csv(species_list_releve, "species_list_releve.csv", row.names = FALSE, quote = FALSE)
species_list_releve <- fread("species_list_releve.csv", sep = ",")



#apply(releve_point_data, 1, function(x) sum(!is.na(x)))
#sum(is.na(apply(releve_point_data, 1, function(x) sum(!is.na(x)))))
#
#all(as.vector(apply(releve_point_data, 1, function(x) sum(!is.na(x)))) == as.numeric(unlist(releve_point[3:nrow(releve_point), 3])))
#
#kk4 <- data.frame(releve_point$...1[-c(1:2)], as.vector(apply(releve_point_data, 1, function(x) sum(!is.na(x)))), as.numeric(unlist(releve_point[3:nrow(releve_point), 3])))
#head(kk4)
#kk4[kk4$as.vector.apply.releve_point_data..1..function.x..sum..is.na.x.... != kk4$as.numeric.unlist.releve_point.3.nrow.releve_point...3..., ]
#kk4[is.na(kk4$as.numeric.unlist.releve_point.3.nrow.releve_point...3...), ]



### Reported Species Richness ####
head(names(releve_point))
tail(names(releve_point))

releve_point[2, 1]
releve_point[2, ]

releve_point_SpRichness <- releve_point[2, -c(1:3)]
releve_point_SpRichness <- as.data.table(releve_point_SpRichness)
ncol(releve_point_SpRichness)  # 728
releve_point_SpRichness

## Keeping only validated points (Momo's harmo)

head(point_geom_experts[, 1])

releve_valid_points <- point_geom_experts$POINT_ID
length(releve_valid_points)  # 605
releve_valid_points

sum(releve_valid_points %in% names(releve_point_SpRichness))  # 599 (6 points missed, why?)
sum(names(releve_point_SpRichness) %in% releve_valid_points)  # 599

# Only valid points in both data sets
releve_valid_points <- releve_valid_points[releve_valid_points %in% names(releve_point_SpRichness)]

releve_point_SpRichness_valid <- releve_point_SpRichness[, .SD, .SDcols = as.character(releve_valid_points)] 
releve_point_SpRichness_valid
ncol(releve_point_SpRichness_valid)


# Spatial info
point_geom_experts
names(point_geom_experts)

library(tidyverse)
releve_point_SpRichness_valid

releve_point_SpRichness_valid <- as.data.table(pivot_longer(releve_point_SpRichness_valid, cols = everything()))
releve_point_SpRichness_valid$name <- as.numeric(releve_point_SpRichness_valid$name)
releve_point_SpRichness_valid

releve_point_SpRichness_valid <- merge(releve_point_SpRichness_valid, 
                                       point_geom_experts[, c("POINT_ID", "SURVEY_GRASS_GPS_LAT", "SURVEY_GRASS_GPS_LON")], 
                                       by.x = "name", by.y = "POINT_ID", all.x = TRUE)
releve_point_SpRichness_valid



## Plotting by species richness 

releve_point_SpRichness_valid$value <- as.numeric(releve_point_SpRichness_valid$value)
sort(unique(releve_point_SpRichness_valid$value)) # 3 to 77

releve_point_SpRichness_valid[, sp_richness_class := cut(releve_point_SpRichness_valid$value,
                                                         #breaks = c(0, 1, 8, 17),
                                                         breaks = c(3, 28, 53, 77),
                                                         include.lowest = TRUE,
                                                         right = FALSE)]
releve_point_SpRichness_valid
sort(unique(releve_point_SpRichness_valid$sp_richness_class))

releve_point_SpRichness_valid$value <- as.factor(releve_point_SpRichness_valid$value)


p <- ggplot() +
  geom_sf(data = eur_gisco) +
  geom_point(
    data = releve_point_SpRichness_valid, 
    aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT,
        color = sp_richness_class),
    #size = 0.1
    size = 0.4
  ) +
  
  theme_light() +
  scale_color_viridis(option = "viridis", discrete = TRUE) +
  labs(title = "LUCAS grasslands 2018 (experts): All species richness") + #, x = "TY [°C]", y = "Txxx") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.title = element_text(size = 16)) +
  guides(color = guide_legend("All species richness classes", override.aes = list(size = 2)))

# https://jtr13.github.io/cc21fall2/tutorial-for-scatter-plot-with-marginal-distribution.html
p1 <- ggMarginal(p,
                 aes(colour = sp_richness_class),
                 type = "density", 
                 #type = "histogram", 
                 #type = "densigram", 
                 #color = "darkred",
                 groupColour = TRUE, groupFill = TRUE)
p1

ggsave("point_geom_experts_releve_SpRichness.png", p1, width = 20, height = 20, units = "cm")



### Comparing key sp richness with total richness #### 

releve_point_SpRichness_valid
point_geom_experts#$NUMBER_KEY_SPECIES

releve_point_SpRichness_valid <- merge(releve_point_SpRichness_valid, 
                                       point_geom_experts[, c("POINT_ID", "NUMBER_KEY_SPECIES")], 
                                       by.x = "name", by.y = "POINT_ID", all.x = TRUE)
releve_point_SpRichness_valid

setnames(releve_point_SpRichness_valid, old = c("value", "NUMBER_KEY_SPECIES"), new = c("All_SpRichness", "Key_SpRichness"))
releve_point_SpRichness_valid

releve_point_SpRichness_valid$All_SpRichness <- as.numeric(releve_point_SpRichness_valid$All_SpRichness)
releve_point_SpRichness_valid$Key_SpRichness <- as.numeric(releve_point_SpRichness_valid$Key_SpRichness)

p1 <- ggplot(releve_point_SpRichness_valid,
             aes(x = All_SpRichness, y = Key_SpRichness)) + 
             #aes(x = Key_SpRichness, y = All_SpRichness )) + 
  geom_point()+
  ggpubr::stat_cor(method = "pearson", label.x = 15, label.y = 15) +
  #geom_smooth(method = "lm")
  #coord_fixed(ratio = 1) +
  geom_smooth(method = "lm") +
  theme_light() +
  labs(title = "LUCAS grasslands 2018 (experts): Key Species vs. All Species Richness") +  #, x = "TY [°C]", y = "Txxx") +
  theme(plot.title = element_text(hjust = 0.5))
  
ggsave("point_geom_experts_SpRichness_AllKey.png", p1, width = 20, height = 20, units = "cm")


cor(releve_point_SpRichness_valid[, c(2, 6)], method = c("pearson"))   # 0.77
cor(releve_point_SpRichness_valid[, c(2, 6)], method = c("spearman"))  # 0.78



#





## Buffer to surveyors points ####


point_geom_surveyors

library(sf)

point_surveyors <- st_as_sf(point_geom_surveyors, coords = c("SURVEY_GRASS_GPS_LON", "SURVEY_GRASS_GPS_LAT"), crs = 4326)
point_surveyors

str(point_surveyors)
names(point_surveyors)

sort(unique(point_surveyors$nuts3))
sort(unique(point_surveyors$NUTS0))


# Selecting points from one NUTS3 in FR
#point_surveyors <- point_surveyors[point_surveyors$nuts3 %in% c("FI193"), ]
#point_surveyors

point_surveyors <- point_surveyors[point_surveyors$NUTS0 %in% c("FR"), ]
point_surveyors


# Buffering
sf::st_crs(point_surveyors)

point_surveyors_buf <- st_buffer(point_surveyors, dist = 50)
point_surveyors_buf




## Subsetting plants in the buffers ####

# all plants occurrences
occs_all <- fread(paste0(getwd(), "/../exploring_lucas_data/D5_FFGRCC_gbif_occ/sp_records_20210709.csv"), header = TRUE)

occs_all
names(occs_all)

occs_all_FR <- occs_all[countryCode == "FR", ] 
occs_all_FR

occs_sf <- st_as_sf(as.data.frame(occs_all_FR), coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)#, agr = "constant")
nrow(occs_sf)



point_surveyors_buf

t0 <- Sys.time()
t0
occs_lucas_50 <- st_join(point_surveyors_buf, occs_sf)
Sys.time() - t0


st_write(obj = occs_lucas_50, dsn = "occs_lucas_50.shp", driver = "ESRI Shapefile", append = FALSE)
occs_lucas_50_dt


occs_lucas_50_dt <- as.data.table(occs_lucas_50)
occs_lucas_50_dt

write.csv(occs_lucas_50_dt, "occs_lucas_50_dt.csv", row.names = FALSE)


occs_lucas_50_dt
nrow(occs_lucas_50_dt)  # 322
names(occs_lucas_50_dt)

head(occs_lucas_50_dt$POINT_ID)
length(unique(occs_lucas_50_dt$POINT_ID))  # 320



occs_lucas_50_dt[duplicated(occs_lucas_50_dt$POINT_ID), 1]

occs_lucas_50_dt[occs_lucas_50_dt$POINT_ID == 41002774, ]





## GBIF Occurrences from all species from the releve ####

species_list_releve #<- fread("species_list_releve.csv", sep = ",")

occs_all <- fread(paste0(getwd(), "/../exploring_lucas_data/D5_FFGRCC_gbif_occ/sp_records_20210709.csv"), header = TRUE)
names(occs_all)


sum(unique(occs_all$species) %in% species_list_releve$species_list_releve)  # 1936 (out of 2672)
sum(!species_list_releve$species_list_releve %in% unique(occs_all$species)) # 736 not in GBIF (out of 2672)
sum(!unique(occs_all$species) %in% species_list_releve$species_list_releve) # 12092 (out of 14028)

sps_NotPrenent_GBIF <- species_list_releve$species_list_releve[!species_list_releve$species_list_releve %in% unique(occs_all$species)]
# From the 736 species not present in the GBIF data set, several are taxa at genus or subspecies level.
# Those at species level, they are probably synonyms (e.g. Agrostis tenuis = Agrostis capillaris)
# As we have a good amount of species, at least for now we will be focused on the 1936 species with no issues 
sps_NotPrenent_GBIF[11:20]
unique(occs_all[grepl("Agrostis", occs_all$species), species])



occs_all_releve <- occs_all[species %in% species_list_releve$species_list_releve, ]
occs_all_releve

length(unique(occs_all_releve$species))  # 1936 species
nrow(occs_all_releve)  # 15410591 occurrences






## SDMs for all species ####


### Species from releve ####

releve_valid_points
length(releve_valid_points)

#releve_point_data
#ncol(releve_point_data)

if(sum(is.na(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point[3, ])) == ncol(DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point)) 
  releve_point_sp <- DLV5_lot9_surveydatabase_LUCAS2018_v3_releve_point[-3, ]

releve_point_sp <- releve_point_sp[-c(1:2), ]
ncol(releve_point_sp)
names(releve_point_sp)[1]

is.character(names(releve_point_sp)[1])

releve_point_data_valid <- releve_point_sp[names(releve_point_sp) %in% c(names(releve_point_sp)[1], releve_valid_points)]
ncol(releve_point_data_valid)
releve_point_data_valid


# rarity (more common species -appearing in more points; not more abundant -appearing more often within a point)

releve_point_data_valid$num_points <- apply(releve_point_data_valid[, -1], 1, function(x) sum(!is.na(x)))
releve_point_data_valid
ncol(releve_point_data_valid)
names(releve_point_data_valid)

as.vector(unlist(releve_point_data_valid[1, ]))
as.vector(unlist(releve_point_data_valid[3, ]))


releve_sps_rarity <- as.data.table(releve_point_data_valid[, c("...1", "num_points")])
releve_sps_rarity

setkeyv(releve_sps_rarity, "num_points")
releve_sps_rarity

releve_sps_rarity[num_points == 0]
releve_sps_rarity[num_points == 0]$...1   # 341 species appearing in 0 relev points



releve_sps_rarity <- releve_sps_rarity[order(-rank(num_points))]
releve_sps_rarity[1:100, ]
View(releve_sps_rarity)



### Subset of species from the releve to be modelled ####
# species appearing in 50 or more points of the releve
releve_sps_rarity_50points <- releve_sps_rarity[num_points >= 50]   # 68 species
releve_sps_rarity_50points

releve_sps_rarity_50points_sps <- sort(releve_sps_rarity_50points$...1)
releve_sps_rarity_50points_sps


## gathering points and coordinates

releve_point_sp  # all releve points (also those with coordinates from theoretical lucas point)

releve_point_sp_4modelling <- releve_point_sp %>% 
  filter( ...1 %in% releve_sps_rarity_50points_sps) %>%
  select(-c("type (sub = subspecies, sp = only identified to genus level)", "number of points present"))

releve_point_sp_4modelling


releve_point_sp_4modelling <- releve_point_sp_4modelling %>% 
  pivot_longer(cols = -1) %>%
  #filter(!is.na(value))
  mutate(value = ifelse(is.na(value), 0, 1))

releve_point_sp_4modelling <- as.data.table(releve_point_sp_4modelling)
setnames(releve_point_sp_4modelling, c("...1", "name", "value"), c("species", "POINT_ID", "presence"))
releve_point_sp_4modelling$POINT_ID <- as.numeric(as.character(releve_point_sp_4modelling$POINT_ID))

## merging coordinates
releve_point_sp_4modelling <- merge(releve_point_sp_4modelling,
                                    point_geom_experts[, .SD, .SDcols = c("POINT_ID", "SURVEY_GRASS_GPS_LAT", "SURVEY_GRASS_GPS_LON")],   # coordinates only for harmonized points
                                    all.x = TRUE,
                                    by = "POINT_ID"
                                    )

releve_point_sp_4modelling <- na.omit(releve_point_sp_4modelling)
releve_point_sp_4modelling

# some checks
length(unique(releve_point_sp_4modelling$POINT_ID))
length(unique(releve_point_sp_4modelling$species))

table(releve_point_sp_4modelling[presence == 1, species])
range(table(releve_point_sp_4modelling[presence == 1, species]))  # range of presences: 50 252 points/species
range(table(releve_point_sp_4modelling[presence == 0, species]))  # range of absences: 347 549 points/species
table(releve_point_sp_4modelling[presence == 0, species])


write.csv(releve_point_sp_4modelling, "releve_point_sp_4modelling.csv", quote = FALSE, row.names = FALSE)
releve_point_sp_4modelling <- fread("releve_point_sp_4modelling.csv", header = TRUE)


## Mapping occurrences just for checking

releve_point_sp_4modelling

releve_point_sp_4modelling_pres <- releve_point_sp_4modelling[presence == 1, ]
releve_point_sp_4modelling_pres



# Gisco maps:   https://ropengov.github.io/giscoR/

eur_gisco <- gisco_get_countries(region = "Europe")
eur_gisco
eur_gisco <- st_crop(eur_gisco, xmin = -10.5, xmax = 40, ymin = 33, ymax = 70)

ggplot() +
  geom_sf(data = eur_gisco) +
  geom_point(
    data = releve_point_sp_4modelling_pres, 
    aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT),
    size = 0.1,
    color = "darkgreen"
  ) 




### Corine LC #### 

#### CLC layers for modelling ####
## Producing some relevant land-use variables for modelling

clc_100 <- stack("/eos/jeodpp/data/base/Landcover/EUROPE/CorineLandCover/CLC2018/VER20-b2/Data/GeoTIFF/100m/clc2018_Version_20_b2.tif")
clc_100

clc_100[clc_100$clc2018 %in% c(521, 522, 523)] <- NA   # Removing marine waters


clc_100 <- raster("clc_100_clean.tif")
clc_100

# 231: Pastures, meadows and other permanent grasslands under agricultural use
# 321: Natural grassland


## Aggregating Pastures and to 1km
aggr_fun_1km <- function(x, ...) {     # returns share of grassland at 1km, from 100m grid (0 to 1)
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(231 # Pastures
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_grassl <- aggregate(x = clc_100, 
                            fact = 10,        # from 100m to 1km
                            fun = aggr_fun_1km, 
                            expand = TRUE, 
                            na.rm = TRUE, 
                            filename = "clc_1km_pastures.tif",
                            #filename = "",
                            overwrite = TRUE)
Sys.time() - t0

#clc_1km_pastures <- raster("clc_1km_pastures.tif")
clc_1km_pastures
plot(clc_1km_pastures)



## Aggregating Naturarl grasslands and to 1km
aggr_fun_1km <- function(x, ...) {     # returns share of grassland at 1km, from 100m grid (0 to 1)
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(321 # Natural grasslands
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_NatGrasslands <- aggregate(x = clc_100, 
                            fact = 10,        # from 100m to 1km
                            fun = aggr_fun_1km, 
                            expand = TRUE, 
                            na.rm = TRUE, 
                            filename = "clc_1km_NatGrasslands.tif",
                            #filename = "",
                            overwrite = TRUE)
Sys.time() - t0

#clc_1km_NatGrasslands <- raster("clc_1km_NatGrasslands.tif")
clc_1km_NatGrasslands
plot(clc_1km_NatGrasslands)



## Aggregating arable land to 1km
aggr_fun_1km <- function(x, ...) {     # returns share of arable at 1km, from 100m grid (0 to 1)
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(211,
                             212,
                             213
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_ArableLand <- aggregate(x = clc_100, 
                            fact = 10,        # from 100m to 1km
                            fun = aggr_fun_1km, 
                            expand = TRUE, 
                            na.rm = TRUE, 
                            filename = "clc_1km_ArableLand.tif",
                            #filename = "",
                            overwrite = TRUE)
Sys.time() - t0

#clc_1km_ArableLand <- raster("clc_1km_ArableLand.tif")
clc_1km_ArableLand
plot(clc_1km_ArableLand)




## Aggregating  permanent crops to 1km
aggr_fun_1km <- function(x, ...) {     # returns share of grassland at 1km, from 100m grid (0 to 1)
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(221,
                             222,
                             223
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_PermanentCrops <- aggregate(x = clc_100, 
                            fact = 10,        # from 100m to 1km
                            fun = aggr_fun_1km, 
                            expand = TRUE, 
                            na.rm = TRUE, 
                            filename = "clc_1km_PermanentCrops.tif",
                            #filename = "",
                            overwrite = TRUE)
Sys.time() - t0

#clc_1km_PermanentCrops <- raster("clc_1km_PermanentCrops.tif")
clc_1km_PermanentCrops
plot(clc_1km_PermanentCrops)




## Aggregating  Heterogeneous agricultural areas to 1km
aggr_fun_1km <- function(x, ...) {     # returns share of Heterogeneous agricultural areas at 1km, from 100m grid (0 to 1)
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(241,
                             242,
                             243,
                             244
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_HeterogAgricAreas <- aggregate(x = clc_100, 
                                       fact = 10,        # from 100m to 1km
                                       fun = aggr_fun_1km, 
                                       expand = TRUE, 
                                       na.rm = TRUE, 
                                       filename = "clc_1km_HeterogAgricAreas.tif",
                                       #filename = "",
                                       overwrite = TRUE)
Sys.time() - t0

#clc_1km_HeterogAgricAreas <- raster("clc_1km_HeterogAgricAreas.tif")
clc_1km_HeterogAgricAreas
plot(clc_1km_HeterogAgricAreas)






## Aggregating  Forest areas to 1km
aggr_fun_1km <- function(x, ...) {     
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(311,
                             312,
                             313
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_Forest <- aggregate(x = clc_100, 
                            fact = 10,        # from 100m to 1km
                            fun = aggr_fun_1km, 
                            expand = TRUE, 
                            na.rm = TRUE, 
                            filename = "clc_1km_Forest.tif",
                            #filename = "",
                            overwrite = TRUE)
Sys.time() - t0

#clc_1km_Forest <- raster("clc_1km_Forest.tif")
clc_1km_Forest
plot(clc_1km_Forest)





## Aggregating  Schrubland areas to 1km
aggr_fun_1km <- function(x, ...) {     
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(322,
                             323,
                             324
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_shrublands <- aggregate(x = clc_100, 
                                fact = 10,        # from 100m to 1km
                                fun = aggr_fun_1km, 
                                expand = TRUE, 
                                na.rm = TRUE, 
                                filename = "clc_1km_shrublands.tif",
                                #filename = "",
                                overwrite = TRUE)
Sys.time() - t0

#clc_1km_shrublands <- raster("clc_1km_shrublands.tif")
clc_1km_shrublands
plot(clc_1km_shrublands)





## Aggregating Open spaces with no vegetation to 1km
aggr_fun_1km <- function(x, ...) {     
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(331,
                             332,
                             333,
                             334,
                             335
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_OpenNoVeg <- aggregate(x = clc_100, 
                                fact = 10,        # from 100m to 1km
                                fun = aggr_fun_1km, 
                                expand = TRUE, 
                                na.rm = TRUE, 
                                filename = "clc_1km_OpenNoVeg.tif",
                                #filename = "",
                                overwrite = TRUE)
Sys.time() - t0

#clc_1km_OpenNoVeg <- raster("clc_1km_OpenNoVeg.tif")
clc_1km_OpenNoVeg
plot(clc_1km_OpenNoVeg)







## Aggregating Inland wetlands to 1km
aggr_fun_1km <- function(x, ...) {     
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(411,
                             412
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_InlandWet <- aggregate(x = clc_100, 
                               fact = 10,        # from 100m to 1km
                               fun = aggr_fun_1km, 
                               expand = TRUE, 
                               na.rm = TRUE, 
                               filename = "clc_1km_InlandWet.tif",
                               #filename = "",
                               overwrite = TRUE)
Sys.time() - t0

#clc_1km_InlandWet <- raster("clc_1km_InlandWet.tif")
clc_1km_InlandWet
plot(clc_1km_InlandWet)




## Aggregating Coastal wetlands to 1km
aggr_fun_1km <- function(x, ...) {     
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(421,
                             422,
                             423
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_CoastalWet <- aggregate(x = clc_100, 
                               fact = 10,        # from 100m to 1km
                               fun = aggr_fun_1km, 
                               expand = TRUE, 
                               na.rm = TRUE, 
                               filename = "clc_1km_CoastalWet.tif",
                               #filename = "",
                               overwrite = TRUE)
Sys.time() - t0

#clc_1km_CoastalWet <- raster("clc_1km_CoastalWet.tif")
clc_1km_CoastalWet
plot(clc_1km_CoastalWet)



## Aggregating Inland waters to 1km
aggr_fun_1km <- function(x, ...) {     
  if (all(is.na(x))){
    gr_share <- NA
  }else{
    gr_share <- sum(x %in% c(511,
                             512
                             )#& 
                    #!is.na(x)
                    , na.rm = TRUE) / 100
  }
  return(gr_share)
}



t0 <- Sys.time()
clc_1km_InlandWaters <- aggregate(x = clc_100, 
                               fact = 10,        # from 100m to 1km
                               fun = aggr_fun_1km, 
                               expand = TRUE, 
                               na.rm = TRUE, 
                               filename = "clc_1km_InlandWaters.tif",
                               #filename = "",
                               overwrite = TRUE)
Sys.time() - t0

#clc_1km_InlandWaters <- raster("clc_1km_InlandWaters.tif")
clc_1km_InlandWaters
plot(clc_1km_InlandWaters)





#### Grassland mask ####
## The purpose of this part was to make an European grassland layer (1km grid) to be used as a mask for modelling
## But almost half of the Lucas grassland (2018) releve points were outside this layer because the CLC-100m is not
## fine-grained enough to discriminate small grassland patches within other classes (e.g. forest).
## If trying to include all the lucas releve points in the mask, almost the whole European land was in the mask, 
## making no sense to use it 

run_this <- "no"
if(run_this == "yes"){
  #clc_100 <- stack("/eos/jeodpp/data/base/Landcover/EUROPE/CorineLandCover/CLC2018/VER20-b2/Data/GeoTIFF/100m/clc2018_Version_20_b2.tif")
  #clc_100
  
  clc_100 <- raster("clc_100_clean.tif")
  
  # 231: Pastures, meadows and other permanent grasslands under agricultural use
  # 321: Natural grassland
  
  ## Aggregating grasslands to 1km
  
  aggr_fun_1km <- function(x, ...) {     # returns share of grassland at 1km, from 100m grid (0 to 1)
    if (all(is.na(x))){
      gr_share <- NA
    }else{
      gr_share <- sum(x %in% c(231, # Pastures
                               321 # Natural grassland
      )#& 
      #!is.na(x)
      , na.rm = TRUE) / 100
    }
    return(gr_share)
  }
  
  
  
  t0 <- Sys.time()
  clc_1km_grassl_share <- aggregate(x = clc_100, 
                              fact = 10,        # from 100m to 1km
                              fun = aggr_fun_1km, 
                              expand = TRUE, 
                              na.rm = TRUE, 
                              filename = "clc_1km_grassl_share.tif",
                              #filename = "",
                              overwrite = TRUE)
  Sys.time() - t0
  
  #writeRaster(clc_1km_grassl, "clc_1km_grassl_share.tif", overwrite = TRUE)
  
  clc_1km_grassl <- raster("clc_1km_grassl_share.tif")
  clc_1km_grassl
  
  sort(unique(getValues(clc_1km_grassl)))
  plot(clc_1km_grassl)
  
  
  
  ## Making a mask for the AOI for modelling (grasslands + other 1km pixels with agricultural and/or (semi)natural areas)
  # Some of the lucas grassland points 2018 still might remain out of this mask. These points should be discarded from the study
  
  aggr_fun_1km <- function(x, ...) {     # returns share of grassland at 1km, from 100m grid (0 to 1)
    if (all(is.na(x))){
      gr_share <- NA
    }else{
      gr_share <- sum(x %in% c(231, # Pastures
                               321, # Natural grassland
                               211, 212, 213, 221, 222, 223, 241, 242, 243, 244, # Agricultural areas
                               311, 312, 313, 322, 323, 324, 333#, # Forest and seminatural areas
                               #112, # Discontinuous urban fabric
                               #121, # Industrial or commercial units and public facilities
                               #122, # Road and rail networks and associated land
                               #411, # Inland marshes
                               #511, # Water courses
                               #512  # Water bodies
      )#& 
      #!is.na(x)
      , na.rm = TRUE) / 100
    }
    return(gr_share)
  }
  
  
  t0 <- Sys.time()
  clc_1km_grassl_mask <- aggregate(x = clc_100, 
                                   fact = 10,        # from 100m to 1km
                                   fun = aggr_fun_1km, 
                                   expand = TRUE, 
                                   na.rm = TRUE, 
                                   #filename = "clc_1km_grassl_mask.tif",
                                   filename = "",
                                   overwrite = TRUE)
  Sys.time() - t0
  
  writeRaster(clc_1km_grassl_mask, "clc_1km_grassl_mask.tif", overwrite = TRUE)
  
  clc_1km_grassl_mask <- raster("clc_1km_grassl_mask.tif")
  sort(unique(getValues(clc_1km_grassl_mask$clc_1km_grassl_mask)))
  table(getValues(clc_1km_grassl_mask$clc_1km_grassl_mask))
  
  ## reclassifying from grassland share to grassland/no grassland
  
  # We will use this layer as a mask for the area of interest (AOI) for modelling
  # For this reason, we keep as grasslands all 1km pixels with at least one 100m pixels of relevant classes
  
  clc_1km_grassl_mask[clc_1km_grassl_mask$clc_1km_grassl_mask < 0.01] <- 0
  clc_1km_grassl_mask[clc_1km_grassl_mask$clc_1km_grassl_mask >= 0.01] <- 1
  
  plot(clc_1km_grassl_mask)
  
  
  wkt <- sf::st_crs(4326)[[2]]
  clc_1km_grassl_wgs84 <- projectRaster(clc_1km_grassl,
                                        crs = wkt,
                                        method = "ngb")
  clc_1km_grassl_wgs84
  unique(getValues(clc_1km_grassl))
  unique(getValues(clc_1km_grassl_wgs84))
  
  
  # Plotting with ggplot2
  clc_1km_grassl_pts <- rasterToPoints(clc_1km_grassl_wgs84, spatial = TRUE)
  head(clc_1km_grassl_pts)
  
  clc_1km_grassl_df <- data.frame(clc_1km_grassl_pts)
  head(clc_1km_grassl_df)
  nrow(clc_1km_grassl_df)
  
  clc_1km_grassl_df_1 <- clc_1km_grassl_df %>% mutate(across(c(x, y), round, digits = 4))
  head(clc_1km_grassl_df_1)
  
  
  #jpeg("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/grassland_layer_1km_releves.jpg")
  ggplot() +
    geom_raster(data = clc_1km_grassl_df_1, aes(x, y, fill = clc_1km_grassl,)) +
    scale_fill_continuous(trans = 'reverse') +
    geom_point(
      data = releve_point_sp_4modelling_pres, 
      aes(x = SURVEY_GRASS_GPS_LON, y = SURVEY_GRASS_GPS_LAT),
      size = 0.1,
      color = "red")
  
  #dev.off()
  
  
  
  ## How many releve points fall in CLC grasslands layer??
  
  releve_point_sp_4modelling_pres_sf <- st_as_sf(as.data.frame(releve_point_sp_4modelling_pres), 
                                                 coords = c("SURVEY_GRASS_GPS_LON", "SURVEY_GRASS_GPS_LAT"), 
                                                 crs = 4326)#, agr = "constant")
  releve_point_sp_4modelling_pres_sf
  
  clc_1km_grassl_releve <- as.data.table(extract(clc_1km_grassl_wgs84,
                                                 releve_point_sp_4modelling_pres_sf, 
                                                 sp = TRUE))
  clc_1km_grassl_releve
  
  clc_1km_grassl_releve <- clc_1km_grassl_releve[!duplicated(clc_1km_grassl_releve$POINT_ID), ]
  clc_1km_grassl_releve
  
  # The following was when using only managed grasslands from CLC (class 231)
  #sum(clc_1km_grassl_releve$clc_1km_grassl == 0)  # threshold = 0.3, 454; threshold = 0.2, 412; threshold = 0.01, 331   
  #sum(clc_1km_grassl_releve$clc_1km_grassl == 1)  # threshold = 0.3, 124; threshold = 0.2, 166; threshold = 0.01, 247
  ## even keeping all 1km pixels with only 1 CLC_100m = grassland (class 231), more than half of the 
  ## releve points with occurrences to be modelled are outside CLC grassland
  ## probably there are other CLC classes which also contain grasslands 
  ## (check https://land.copernicus.eu/user-corner/technical-library/corine-land-cover-nomenclature-guidelines/html) 
  
  
  sum(clc_1km_grassl_releve$clc_1km_grassl == 0)  # threshold = 0.01, 275   
  sum(clc_1km_grassl_releve$clc_1km_grassl == 1)  # threshold = 0.01, 303
  
  
  ## Check classes for the points outside grasslands
  clc_1km_grassl_releve_pointOut <- clc_1km_grassl_releve[clc_1km_grassl_releve$clc_1km_grassl == 0, ]
  write.csv(clc_1km_grassl_releve_pointOut, "clc_1km_grassl_releve_pointOut.csv", quote = FALSE, row.names = FALSE)
  
  
  clc_1km_grassl_releve_pointOut_sf <- st_as_sf(as.data.frame(clc_1km_grassl_releve_pointOut), 
                                                coords = c("coords.x1", "coords.x2"), 
                                                crs = 4326)#, agr = "constant")
  
  clc_1km_grassl_releve_pointOut_CLC100 <- as.data.table(extract(clc_100,
                                                                 clc_1km_grassl_releve_pointOut_sf, 
                                                                 sp = TRUE))
  clc_1km_grassl_releve_pointOut_CLC100
  sort(unique(clc_1km_grassl_releve_pointOut_CLC100$clc2018))
  table(clc_1km_grassl_releve_pointOut_CLC100$clc2018)
  
  
} # end of 'run_this'






### Predictors ####

## Bioclimatic variables from WorldClim 
preds_dir <- "/eos/jeodpp/home/users/rotllxa/European_butterflies_SDMs_data/"

worldclim_all <- stack(paste0(preds_dir, "worldclim_all.tif"))
worldclim_all
names(worldclim_all)

worldclim_all_names <- c("wc2.1_30s_bio_1", "wc2.1_30s_bio_10", "wc2.1_30s_bio_11", "wc2.1_30s_bio_12", "wc2.1_30s_bio_13", "wc2.1_30s_bio_14",
                         "wc2.1_30s_bio_15", "wc2.1_30s_bio_16", "wc2.1_30s_bio_17", "wc2.1_30s_bio_18", "wc2.1_30s_bio_19", "wc2.1_30s_bio_2",
                         "wc2.1_30s_bio_3", "wc2.1_30s_bio_4", "wc2.1_30s_bio_5", "wc2.1_30s_bio_6", "wc2.1_30s_bio_7", "wc2.1_30s_bio_8",
                         "wc2.1_30s_bio_9",  "wc2.1_30s_elev") 
worldclim_all_names <- gsub("wc2.1_30s_", "", worldclim_all_names)
names(worldclim_all) <- worldclim_all_names

names(worldclim_all)


## Corine variables
clc_vars <- list.files()[grepl("clc_1km_", list.files())]
clc_vars
clc_vars <- clc_vars[!grepl(".aux", clc_vars)]
clc_vars <- clc_vars[!grepl(".csv", clc_vars)]
clc_vars <- clc_vars[!grepl("mask", clc_vars)]
clc_vars <- clc_vars[!grepl("_grassl", clc_vars)]

clc_vars

clc_vars_all <- stack(clc_vars)
clc_vars_all
names(clc_vars_all)


## Reprojecting Worldclim to LAEA

worldclim_all_laea <- projectRaster(worldclim_all, crs = crs(clc_vars_all),
                                    filename = paste0(preds_dir, "worldclim_all_laea.tif"))

worldclim_all_laea
#writeRaster(worldclim_all_laea, filename = paste0(preds_dir, "worldclim_all_laea.tif"))

worldclim_all_laea <- stack(paste0(preds_dir, "worldclim_all_laea.tif"))



## Cropping rasters
worldclim_all_laea <- crop(worldclim_all_laea, clc_vars_all)
worldclim_all_laea


## Resample
worldclim_all_laea <- resample(worldclim_all_laea, clc_vars_all)

names(worldclim_all_laea) <- worldclim_all_names
worldclim_all_laea
clc_vars_all

#plot(worldclim_all_laea_2[[1]])
#
#writeRaster(worldclim_all_laea_2[[4]], "worldclim_all_laea_2_bio12.tif")
#writeRaster(worldclim_all_laea[[4]], "worldclim_all_laea_bio12.tif")
#writeRaster(worldclim_all[[4]], "worldclim_all_bio12.tif")
#plot(clc_vars_all[[1]])



## WorldClim + Corine

#all_vars <- stack(worldclim_all_laea, clc_vars_all)
#writeRaster(all_vars,"all_vars.grd", format = "raster", overwrite = TRUE)

all_vars <- brick("all_vars.grd")
all_vars

all_vars <- crop(all_vars, extent(2500000, 6000000, 1400000, 5500000))
plot(all_vars[[4]])

names(all_vars)
if("layer.1" %in% names(all_vars))  names(all_vars)[names(all_vars) == "layer.1"] <- "clc_1km_CoastalWet" 
if("layer.2" %in% names(all_vars))  names(all_vars)[names(all_vars) == "layer.2"] <- "clc_1km_InlandWaters" 
if("layer.3" %in% names(all_vars))  names(all_vars)[names(all_vars) == "layer.3"] <- "clc_1km_InlandWet" 
if("layer.4" %in% names(all_vars))  names(all_vars)[names(all_vars) == "layer.4"] <- "clc_1km_OpenNoVeg" 
names(all_vars)

plot(all_vars[[c(1, 28)]])


all_vars <- mask(all_vars, all_vars[[28]])
plot(all_vars[[c(1, 28)]])

all_vars <- mask(all_vars, all_vars[[1]])
plot(all_vars[[c(1, 28)]])

all_vars

writeRaster(all_vars, "all_vars_clean.grd", format = "raster", overwrite = TRUE)

all_vars <- brick("all_vars_clean.grd")


## Removing multicolinarity

library(virtualspecies)

all_vars_NoC <- removeCollinearity(all_vars,
                                   multicollinearity.cutoff = 0.70,
                                   #multicollinearity.cutoff = 0.85,
                                   select.variables = TRUE,  # if TRUE, randomly select one variable of the group. If FALSE, returns a list with the groups
                                   sample.points = TRUE,
                                   nb.points = 10^6,
                                   plot = TRUE)
all_vars_NoC

dev.copy(png,'all_vars_NoC.png')
dev.off()


all_vars_NoC_rstr <- all_vars[[all_vars_NoC]]

writeRaster(all_vars_NoC_rstr, "all_vars_NoC_rstr.grd", format = "raster", overwrite = TRUE)


all_vars_NoC_rstr <- brick("all_vars_NoC_rstr.grd")
names(all_vars_NoC_rstr)



### Modelling with MaxEnt ####


## Predictors
all_vars <- all_vars_NoC_rstr; rm(all_vars_NoC_rstr)

all_vars_data <- as.data.frame(all_vars)
names(all_vars_data)
apply(all_vars_data, 2, function(x) sum(is.na(x)))

all_vars_data <- all_vars_data[complete.cases(all_vars_data), ]
all_vars_data
nrow(all_vars_data)

plot(all_vars[[c(1, 10)]])



## Releve occurrences for (modelling and) validating models
releve_point_sp_4modelling_sf <- st_as_sf(releve_point_sp_4modelling, coords = c("SURVEY_GRASS_GPS_LON", "SURVEY_GRASS_GPS_LAT"), crs = 4326)
releve_point_sp_4modelling_laea <- st_transform(releve_point_sp_4modelling_sf, crs = st_crs(all_vars))
releve_point_sp_4modelling_laea_coords <- data.table(st_coordinates(releve_point_sp_4modelling_laea))
releve_point_sp_4modelling_laea <- cbind(data.table(releve_point_sp_4modelling_laea), releve_point_sp_4modelling_laea_coords)

releve_point_sp_4modelling_laea


releve_point_sp_4modelling_pres_laea <- releve_point_sp_4modelling_laea[presence == 1, ]

length(unique(releve_point_sp_4modelling_pres_laea$POINT_ID))  # 578 releve points
length(unique(releve_point_sp_4modelling_pres_laea$species))  # 68 sp
sort(table(releve_point_sp_4modelling_pres_laea$species))





## GBIF occurrences for modelling

releve_point_sp_4modelling_laea
sort(unique(releve_point_sp_4modelling_laea$species))


sp_releve_sp4model <- gsub("\\.", " ", sort(unique(releve_point_sp_4modelling_laea$species)))

sum(unique(occs_all_releve$species) %in% sp_releve_sp4model)

sp_releve_sp4model[!sp_releve_sp4model %in% unique(occs_all_releve$species)]  # 5 sps from the subset of the releve that have no occurrences in GBIF 


sp_releve_sp4model <- sp_releve_sp4model[sp_releve_sp4model %in% unique(occs_all_releve$species)]  # 63 sps ready to be modelled
sp_releve_sp4model


occs_all_releve_1 <- occs_all_releve[species %in% sp_releve_sp4model]
occs_all_releve_1

sort(unique(occs_all_releve_1$species))
sort(table(occs_all_releve_1$species))    # Sanguisorba minor has only one occurrence in GBIF. Removed
range(table(occs_all_releve_1$species))

occs_all_releve_1 <- occs_all_releve_1[species != "Sanguisorba minor", ]
range(table(occs_all_releve_1$species))   # 1716 - 293287 occurrences


occs_all_releve_1_sf <- st_as_sf(occs_all_releve_1, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
occs_all_releve_1_laea <- st_transform(occs_all_releve_1_sf, crs = st_crs(all_vars))
occs_all_releve_1_laea_coords <- data.table(st_coordinates(occs_all_releve_1_laea))
occs_all_releve_1_laea <- cbind(data.table(occs_all_releve_1_laea), occs_all_releve_1_laea_coords)
occs_all_releve_1_laea


## Background points

bckgr <- dismo::randomPoints(all_vars[[1]], 
                             #n = 10000)
                             n = 15000)
bckgr <- as.data.frame(bckgr)
head(bckgr)
nrow(bckgr)

write.csv(bckgr, "background_points.csv", row.names = FALSE)
bckgr <- read.csv("background_points.csv", header = TRUE)
head(bckgr)



info_models_maxent <- c()
data2save_ReleveValid <- c()

# Threshold to use for converting to presence/absence
# Options: kappa,  spec_sens, no_omission, prevalence, equal_sens_spec, sensitivity
#threshold2use <- "sensitivity"    # deffault 0.9
#threshold2use <- "sensitivity.99"    # sensitivity 0.99
#threshold2use <- "no_omission"    # keeping all presences
#threshold2use <- "prevalence"    # modeled prevalence is closest to observed prevalenc
#threshold2use <- "equal_sens_spec"    # equal sensitivity and specificity
#threshold2use <- "kappa"    #  the threshold at which kappa is highest ("max kappa")
threshold2use <- "spec_sens"    # sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest

Sys.time()

#taxons <- unique(releve_point_sp_4modelling_pres_laea$species)
taxons <- unique(occs_all_releve_1_laea$species)
occs_all <- occs_all_releve_1_laea


for (t in taxons){
  #print(t)
  #t <- taxons[1]
  t0 <- Sys.time()
  print(t0)
  #sps <- spcies[spcies$taxons == t, "sps"]
  sps <- t
  
  print(paste0("running... ", sps))
  
  dir2save_maxent <- paste0("models_maxent_", t, "/")
  if(!dir.exists(paste0("models_maxent_", t))) {
    dir.create(dir2save_maxent)
  }
  
  dir2save_presences <- paste0("pres4modelling", "/")
  if(!dir.exists(paste0("pres4modelling"))) {
    dir.create(dir2save_presences)
  }
  
  occs_i <- occs_all[occs_all$species %in% t, c("X", "Y")]
  occurrences_raw <- nrow(occs_i)
  
  occs_i_shp <- SpatialPointsDataFrame(coords = occs_i[, c("X", "Y")],
                                       data = data.frame(sp = rep(1, nrow(occs_i))),
                                       proj4string = CRS("+init=EPSG:3035"))
  #names(occs_i_shp) <- t
  occs_i_rstr <- rasterize(occs_i_shp, all_vars[[1]], field = "sp", background = 0)
  #names(occs_i_rstr) <- t
  #occs_i_rstr <- occs_i_rstr[[2]]
  occs_i_rstr <- mask(occs_i_rstr, all_vars[[1]])
  #plot(occs_i_rstr)
  
  #assign(paste0(t, "_rstr"), occs_i_rstr)
  #print(sum(getValues(occs_i_rstr) == 1, na.rm = T))
  
  
  ## occurrences for training/testing
  sps_data <- stack(occs_i_rstr, all_vars) 
  sps_data <- as.data.table(as.data.frame(sps_data))
  sps_data[, raster_position := 1:nrow(sps_data)]
  
  # data set for presences
  sps_data_presences <- sps_data[layer == 1, ]
  sps_data_presences <- sps_data_presences[complete.cases(sps_data_presences), ]
  occurrences_1km <- nrow(sps_data_presences)
  rm(sps_data); gc()
  
  # data set for pseudo-absences
  sps_data_absences <- as.data.table(as.data.frame(raster::extract(all_vars, bckgr, cellnumbers = TRUE)))
  sps_data_absences <- sps_data_absences[!sps_data_absences$cells %in% sps_data_presences$raster_position, ]
  names(sps_data_absences)
  
  nrow(sps_data_presences)
  nrow(sps_data_absences)
  
  prop4test <- 0.3
  prop4train <- 1 - prop4test
  
  sps_data_presences_train <- sample_n(sps_data_presences, ceiling(nrow(sps_data_presences) * prop4train))
  sps_data_presences_test <- sps_data_presences[!sps_data_presences$raster_position %in% sps_data_presences_train$raster_position, ]
  
  write.csv(sps_data_presences_train, paste0(dir2save_presences, "/sps_data_presences_train_", t, ".csv"), row.names = FALSE)
  write.csv(sps_data_presences_test, paste0(dir2save_presences, "/sps_data_presences_test_", t, ".csv"), row.names = FALSE)
  write.csv(sps_data_absences, paste0(dir2save_presences, "/sps_data_absences_", t, ".csv"), row.names = FALSE)
  
  #sps_data_presences_train <- fread(paste0(dir2save_presences, "/sps_data_presences_train_", t, ".csv"), header = TRUE)
  #sps_data_presences_test <- fread(paste0(dir2save_presences, "/sps_data_presences_test_", t, ".csv"), header = TRUE)
  #sps_data_absences <- fread(paste0(dir2save_presences, "/sps_data_absences_", t, ".csv"), header = TRUE)  
  
  
  ## Running ENMeval (https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html)
  ## Including a tryCatch to avoid stop process if there's an error because of a bug with "H" transformation or something
  
  library(dismo)
  library(ENMeval)
  
  dir_func <- function(sps_data_presences, all_vars, sps_data_absences, fc){ # to avoid stop modelling if low number of background points or other errors
    res <- tryCatch(
      {
        library(dismo)
        library(ENMeval)
        #modl1 <- ENMevaluate(occs = sps_data_presences[1:3000, .SD, .SDcols = names(all_vars)], 
        modl1 <- ENMevaluate(occs = sps_data_presences_train[, .SD, .SDcols = names(all_vars)], 
                             envs = NULL, 
                             bg = sps_data_absences[, .SD, .SDcols = names(all_vars)], 
                             algorithm = 'maxnet', 
                             #partitions = 'block', 
                             partitions = "testing",
                             #occs.testing = sps_data_presences[3001:3750, .SD, .SDcols = names(all_vars)],  # occurrences for testing; only when partitions = 'testing'
                             occs.testing = sps_data_presences_test[, .SD, .SDcols = names(all_vars)],  # occurrences for testing; only when partitions = 'testing'
                             tune.args = list(
                               fc = fc,
                               rm = c(1, 2, 5)
                               #rm = 1:2
                             ),
                             quiet = TRUE,
                             parallel = TRUE,
                             #parallel = FALSE,
                             numCores = 12
                             #numCores = 4
        )
        
      },
      error = function(con){
        message(con)
        return(NULL)
      }
    )
    if(exists("modl1")){ return(modl1) }else{ return(NULL) }
  } #end of dir_func
  
  
  fc_opts <- list(c("L","LQ","LQH","H"), c("L","LQ","LQH"), c("L","LQ"), "L")
  
  for(fc in fc_opts){
    modl <- dir_func(sps_data_presences, all_vars, sps_data_absences, fc)
    if(!is.null(modl)) break
  }
  
  
  #rm(sps_data_absences); gc()
  modl
  modl@results
  #View(modl@results)
  write.csv(modl@results, file = paste0(dir2save_maxent, "ENMeval_results_", t, ".csv"))
  save(modl, file = paste0(dir2save_maxent, "models_", t, ".RData"))
  #evalplot.stats(e = modl, stats = "or.mtp", color = "fc", x.var = "rm")
  #load(paste0(dir2save_maxent, "models_", t, ".RData"), verbose = TRUE)
  
  occurrences_train <- nrow(modl@occs)
  occurrences_test <- nrow(modl@occs.testing)  # none because cross-validation
  background_points <- nrow(modl@bg)
  
  
  # selecting optimal model
  results <- eval.results(modl)
  results
  #View(results)
  optimal <- results %>% filter(delta.AICc == 0)
  optimal
  if(nrow(optimal) > 1) optimal <- optimal[1, ]
  
  modl_args <- eval.models(modl)[[optimal$tune.args]]
  modl_args$betas
  #str(modl_args)
  
  #dev.off()
  pdf(paste0(dir2save_maxent, "opt_model_RespCurves_", t, ".pdf"))
  plot(modl_args, type = "cloglog")
  # And these are the marginal response curves for the predictor variables wit non-zero 
  # coefficients in our model. We define the y-axis to be the cloglog transformation, which
  # is an approximation of occurrence probability (with assumptions) bounded by 0 and 1
  # (Phillips et al. 2017).
  dev.off()
  
  modl <- modl@models[[optimal$tune.args]]
  gc()
  
  #save(modl, file = paste0(dir2save_maxent, "opt_model_", t, ".RData"))
  
  # making predictions
  #worldclim_all_data <- fread("worldclim_all_data_NoCor.csv", header = TRUE)
  #worldclim_all_data <- fread("worldclim_all_data_NoCor_070.csv", header = TRUE)
  #worldclim_all_data <- fread("worldclim_all_data_NoCor_040.csv", header = TRUE)
  #worldclim_all_data <- fread(paste0(preds_dir, "worldclim_all_data.csv"), header = TRUE)
  #names(worldclim_all_data) <- names(worldclim_all)
  #names(worldclim_all_data) <- gsub("wc2.1_30s_bio_", "worldclim_all.", names(worldclim_all_data))
  #names(worldclim_all_data) <- gsub("wc2.1_30s_elev", "worldclim_all.20", names(worldclim_all_data))
  #names(worldclim_all_data) <- gsub("worldclim_all", "worldclim_all", names(worldclim_all_data))
  #worldclim_all_data <- worldclim_all_data[complete.cases(worldclim_all_data), ]
  
    
  sps_predictions_maxent <- predict(object = modl, 
                                    newdata = all_vars_data, 
                                    clamp = TRUE,
                                    type = c("cloglog")
                                    )
  #rm(worldclim_all_data); gc()
  sps_predictions_maxent <- as.data.table(sps_predictions_maxent)
  head(sps_predictions_maxent)
  range(sps_predictions_maxent)
  nrow(sps_predictions_maxent)
  
  
  all_vars_data0 <- as.data.table(as.data.frame(all_vars[[1]]))
  all_vars_data0$raster_position <- 1:nrow(all_vars_data0)
  
  all_vars_data1 <- all_vars_data0
  all_vars_data1 <- all_vars_data1[complete.cases(all_vars_data1), ]
  
  all_vars_data0 <- all_vars_data0[, .SD, .SDcols = "raster_position"]
  
  all_vars_data1[, predictions := sps_predictions_maxent$V1]
  
  
  all_vars_data0 <- merge(all_vars_data0[, "raster_position", with = FALSE], 
                          all_vars_data1[, .SD, .SDcols = c("raster_position", "predictions")], 
                          by = "raster_position", all.x = TRUE)
  
  #rm(worldclim_all_data1); gc()
  
  sps_preds_rstr <- all_vars[[1]]
  sps_preds_rstr <- setValues(sps_preds_rstr, all_vars_data0$predictions)
  names(sps_preds_rstr) <- "predictions_maxent"
  
  #rm(worldclim_all_data0); gc()
  
  
  #pdf("sps_predictions_maxent_kk.pdf", width = 20, height = 15)
  #par(mfrow = c(1, 2))
  #plot(sps_preds_rstr, zlim = c(0, 1))
  #plot(occs_i_shp, add = TRUE, col = "black")
  #plot(sps_preds_rstr, zlim = c(0, 1))
  #dev.off()
  
  
  #BI_mxnt <- ecospat::ecospat.boyce(fit = sps_preds_rstr,
  #                                  obs = linaria_pres_test_coords, 
  #                                  nclass = 0, 
  #                                  window.w = "default", 
  #                                  res = 100, 
  #                                  PEplot = TRUE)
  
  ## Creating presence/absence map
  # Threshold: minimum presence
  
  #info_models_maxent
  #threshold1 <- min(extract(sps_preds_rstr, occs_i_shp))
  #threshold1 <- quantile(extract(sps_preds_rstr, occs_i_shp), 0.1)#, na.rm = TRUE) # sensitivity = 0.9
  #threshold1
  
  thresholds <- dismo::threshold(dismo::evaluate(raster::extract(sps_preds_rstr, occs_i_shp), 
                                                 raster::extract(sps_preds_rstr, bckgr))) # sensitibity default 0.9
  
  thresholds$sensitivity.99 <- dismo::threshold(dismo::evaluate(raster::extract(sps_preds_rstr, occs_i_shp),
                                                                raster::extract(sps_preds_rstr, bckgr)),
                                                stat = "sensitivity", sensitivity = 0.99) 
  thresholds
  
  #threshold2 <- as.numeric(thresholds$sensitivity)
  #threshold2 <- as.numeric(thresholds$no_omission) # keeping all presences
  threshold2 <- as.numeric(thresholds[names(thresholds) %in% threshold2use])
  threshold_used <- threshold2 #<- 0.09993547
  #threshold2use <- "no_omission"
  #threshold2use <- "sensitivity"
  
  a <- c(0, threshold2, 0)
  b <- c(threshold2, 1, 1)
  thr <- rbind(a, b)
  
  sps_preds_rstr_pres_abs <- reclassify(sps_preds_rstr[["predictions_maxent"]], rcl = thr, filename = '', include.lowest = FALSE, right = TRUE)
  sps_preds_rstr_pres_abs_all <- brick(sps_preds_rstr_pres_abs)
  names(sps_preds_rstr_pres_abs_all) <- c("Pres_Abs_MaxEnt")
  
  #plot(sps_preds_rstr_pres_abs)
  
  #pdf(paste0(dir2save_maxent, "sps_predictions_maxent_", t, ".pdf"), width = 18, height = 15)
  png(paste0(dir2save_maxent, "sps_predictions_maxent_", t, ".png"), width = 20, height = 25, units = "cm", res = 75)
  #par(mar = c(6, 8, 6, 8), oma = c(4,0,8,0))
  par(mar = c(5, 2, 5, 2), oma = c(4,0,8,0))
  par(mfrow = c(2, 2))
  #plot(sps_preds_rstr[["predictions_maxent"]], zlim = c(0, 1), main = "Occurrences (1km)", cex.main = 2, cex.sub = 1.5, legend = FALSE)
  #plot(occs_i_shp, add = TRUE, col = "black")
  plot(occs_i_rstr, col = c("LightGrey", "LightGrey"), main = "Occurrences (1km)", cex.main = 1.2, legend = FALSE)
  occs_i_rstr_centr <- rasterToPoints(occs_i_rstr, fun = function(x){x == 1}, spatial = TRUE)
  plot(occs_i_rstr_centr, add = TRUE, col = "darkred", pch = 3, cex = 0.01)
  plot(sps_preds_rstr[["predictions_maxent"]], zlim = c(0, 1), main = "MaxEnt predictions (cloglog)", cex.main = 1.2, cex.sub = 1.5)
  plot(sps_preds_rstr_pres_abs, main = "Presence-Absence", 
       #sub = paste0("Threshold: '", threshold2use, " (0.99)'"), 
       sub = paste0("Threshold: '", threshold2use, "'"), 
       cex.main = 1.2, cex.sub = 0.9, legend = FALSE)
  title(list(paste0(sps),
             #cex = 4), 
             cex = 2), 
        line = 1, outer = TRUE)
  
  dev.off()
  
  
  
  
  
  #### Validating with releve data
  library(PresenceAbsence)
  
  releve_point_sp_4modelling_laea
  
  releve_point_sp_4modelling_laea_i <- releve_point_sp_4modelling_laea[species %in% gsub(" ", ".", t), ] 
  table(releve_point_sp_4modelling_laea_i$presence)
  releve_point_sp_4modelling_laea_i
  
  
  releve_point_sp_4modelling_laea_i$maxent_preds <- raster::extract(sps_preds_rstr[["predictions_maxent"]], 
                                                                    releve_point_sp_4modelling_laea_i[, .SD, .SDcols = c("X", "Y")]
                                                                    )
  
  releve_point_sp_4modelling_laea_i
  sum(is.na(releve_point_sp_4modelling_laea_i$maxent_preds))  # 20 points with no prediction
  
  
  #obs_pred_conf <- cmx(releve_point_sp_4modelling_laea_i[, .SD, .SDcols = c("POINT_ID", "presence", "maxent_preds")], 
  #                     threshold = threshold_used,
  #                     na.rm = TRUE)
  #obs_pred_conf
  
  data2save_ReleveValid_i <- c()
  
  for(tr in (1:length(thresholds))){
    #print(names(thresholds)[t])
  
    obs_pred_conf_i <- cmx(releve_point_sp_4modelling_laea_i[, .SD, .SDcols = c("POINT_ID", "presence", "maxent_preds")], 
                           threshold = as.numeric(thresholds[tr]),
                           na.rm = TRUE)
    #print(obs_pred_conf_i)
    
    ## Proportion of correctly classified observations
    #print(paste0("Prop Correctly Classified = ", round(pcc(obs_pred_conf_i, st.dev = F), 2)))
    ProportionCorrectlyClassified <- round(pcc(obs_pred_conf_i, st.dev = F), 2)
    #
    ## Sensitivity = true positive rate
    #print(paste0("True Positive Rate = ", sensitivity(obs_pred_conf_i, st.dev = F)))
    TruePositiveRate <- round(sensitivity(obs_pred_conf_i, st.dev = F))
    #
    ## Specificity = true negative rate
    #print(paste0("True Negative Rate = ", round(specificity(obs_pred_conf_i, st.dev = F), 2)))
    TrueNegativeRate <- round(specificity(obs_pred_conf_i, st.dev = F), 2)
    #
    ## Kappa. According to Araujo et al. (2005), Kappa > 0.4 indicate good predictions. 
    #print(paste0("Kappa = ", round(Kappa(obs_pred_conf_i, st.dev = F), 2)))
    Kappa <- round(Kappa(obs_pred_conf_i, st.dev = F), 2)
    #
    ## True skill statistic (sensit+specit-1)
    ## For TSS, we often assume TSS > 0.5 to indicate good predictions
    #print(paste0("TSS = ", round((sensitivity(obs_pred_conf_i, st.dev = F) + 
    #                          specificity(obs_pred_conf_i, st.dev = F) -
    #                          1), 2)))
    TSS <- round((sensitivity(obs_pred_conf_i, st.dev = F) + 
                    specificity(obs_pred_conf_i, st.dev = F) -
                    1), 2)
    #
    #
    #
    #print(" ")
    
    
    data2save_ReleveValid_i_i <- data.frame(species = t, 
                                          threshold = names(thresholds)[tr], 
                                          threshold_value = as.numeric(thresholds[tr]),
                                          ProportionCorrectlyClassified,
                                          TruePositiveRate,
                                          TrueNegativeRate,
                                          Kappa,
                                          TSS)
    
    data2save_ReleveValid_i <- rbind(data2save_ReleveValid_i, data2save_ReleveValid_i_i)
  }
  
  data2save_ReleveValid <- rbind(data2save_ReleveValid, data2save_ReleveValid_i)
  data2save_ReleveValid
  
  write.csv(data2save_ReleveValid, "info_models_maxent_ReleveValid.csv", row.names = FALSE)
  write.csv(data2save_ReleveValid, paste0(dir2save_maxent, "info_models_maxent_ReleveValid.csv"), row.names = FALSE)
  
  
  ####
  
  
  
  running_time <- as.vector(Sys.time() - t0)
  if(exists("data2save")) rm(data2save)
  data2save <- data.frame(species = t, occurrences_raw, occurrences_1km, occurrences_train,
                          occurrences_test, background_points, optimal,
                          thresholds, threshold_used)
  rownames(data2save) <- t
  
  info_models_maxent <- rbind(info_models_maxent, data2save)
  #write.csv(info_models_maxent, "info_models_all_species.csv", row.names = FALSE)
  #write.csv(info_models_maxent, "info_models_all_species_085.csv", row.names = FALSE)
  write.csv(info_models_maxent, "info_models_maxent_all.csv", row.names = FALSE)
  write.csv(info_models_maxent, paste0(dir2save_maxent, "info_models_maxent_all.csv"), row.names = FALSE)
  #info_models_maxent <- fread("info_models_maxent_all.csv")
  #info_models_maxent <- info_models_maxent[-3, ]
  
  print(paste0(t, " run in: ", running_time))
  
  
}  


info_models_maxent_all <- fread("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/info_models_maxent_all.csv", header = TRUE)
info_models_maxent_all


data2save_ReleveValid <- fread("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/info_models_maxent_ReleveValid.csv", header = TRUE)
data2save_ReleveValid


sps_kk <- unique(data2save_ReleveValid$species)#[-1]
sps_kk

sps_kk_i <- sps_kk[3]
sps_kk_i

data2save_ReleveValid[species == sps_kk_i]

data2save_ReleveValid[species == "Rumex obtusifolius"]



max(data2save_ReleveValid$TSS)
max(data2save_ReleveValid$Kappa)

data2save_ReleveValid[TSS >= 0.5]
data2save_ReleveValid[Kappa >= 0.4]




## Comparing with "conventional" validation

data2save_ReleveValid_01 <- read.csv(unz("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/Maxent_01.zip", 
                                      "info_models_maxent_ReleveValid.csv"), 
                                  header = TRUE, sep = ",") %>% data.table

info_models_maxent_all_01 <- read.csv(unz("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/Maxent_01.zip", 
                                       "info_models_maxent_all.csv"), 
                                   header = TRUE, sep = ",") %>% data.table


info_models_maxent_all_01[species == "Bellis perennis", ]
data2save_ReleveValid_01[species == "Bellis perennis", ]

info_models_maxent_all[species == "Galium verum", ]
data2save_ReleveValid[species == "Galium verum", ]

data2save_ReleveValid[species == "Plantago media", ]

