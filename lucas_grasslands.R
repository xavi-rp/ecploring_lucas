

library(data.table)
library(devtools)
#install_github("xavi-rp/PreSPickR", 
#               ref = "v2", 
#               INSTALL_opts = c("--no-multiarch"))  # https://github.com/rstudio/renv/issues/162
library(PreSPickR)

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
point_geom_experts <- paste0(dir_LucasGrassland_Momo, "/estat_e_attr_point_allattr.csv")


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

eur_gisco <- st_crop(eur_gisco, xmin = -10.5, xmax = 40, ymin = 33, ymax = 70)

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


## list of species
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

write.csv(species_list_releve, "species_list_releve.csv", row.names = FALSE, quote = FALSE)
fread("species_list_releve.csv", sep = ",")



apply(releve_point_data, 1, function(x) sum(!is.na(x)))
sum(is.na(apply(releve_point_data, 1, function(x) sum(!is.na(x)))))

all(as.vector(apply(releve_point_data, 1, function(x) sum(!is.na(x)))) == as.numeric(unlist(releve_point[3:nrow(releve_point), 3])))

kk4 <- data.frame(releve_point$...1[-c(1:2)], as.vector(apply(releve_point_data, 1, function(x) sum(!is.na(x)))), as.numeric(unlist(releve_point[3:nrow(releve_point), 3])))
head(kk4)
kk4[kk4$as.vector.apply.releve_point_data..1..function.x..sum..is.na.x.... != kk4$as.numeric.unlist.releve_point.3.nrow.releve_point...3..., ]
kk4[is.na(kk4$as.numeric.unlist.releve_point.3.nrow.releve_point...3...), ]









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

occs_lucas_50


occs_lucas_50_dt <- as.data.table(occs_lucas_50)
occs_lucas_50_dt

nrow(occs_lucas_50_dt)
