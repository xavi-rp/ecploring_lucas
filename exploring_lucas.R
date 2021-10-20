library(rgbif)
library(dplyr)
library(raster)
library(rgdal)
library(rvest)
library(sp)
library(sf)
library(data.table)


library(devtools)
install_github("xavi-rp/PreSPickR", 
               ref = "v2", 
               INSTALL_opts = c("--no-multiarch"))  # https://github.com/rstudio/renv/issues/162
library(PreSPickR)



if(Sys.info()[4] == "D01RI1700308") {
  wd <- "D:/xavi_rp/D5_FFGRCC_gbif_occ/"
}else if(Sys.info()[4] == "S-JRCIPRAP320P") {
  wd <- "D:/rotllxa/D5_FFGRCC_gbif_occ/"
}else if(Sys.info()[4] %in% c("jeodpp-terminal-jd001-03", "jeodpp-terminal-03")) {
  if(!dir.exists("/eos/jeodpp/home/users/rotllxa/exploring_lucas_data/")) 
    dir.create("/eos/jeodpp/home/users/rotllxa/exploring_lucas_data/")
  wd <- "/eos/jeodpp/home/users/rotllxa/exploring_lucas_data/"
  gbif_creds <- "/home/rotllxa/Documents/"
}else{
  wd <- "C:/Users/rotllxa/D5_FFGRCC_occurrences/"
  gbif_creds <- "C:/Users/rotllxa/Documents/"
}

setwd(wd)

# file.edit('~/.Renviron')

## Preliminary checks (GBIF - rgbif) ####

# e.g.
as.data.frame(name_backbone(name='Helianthus annuus', kingdom='plants'))
name_backbone(name='Helianthus annuus', kingdom='plants')$speciesKey


# These are Gymnosperms and Angiosperms
as.data.frame(name_backbone(name='Plantae', kingdom='plants'))
as.data.frame(name_backbone(name='Plantae'))

as.data.frame(name_backbone(name='Tracheophyta'))


as.data.frame(name_backbone(name='Magnoliopsida', kingdom='plants'))
as.data.frame(name_backbone(name='Liliopsida', kingdom='plants'))
as.data.frame(name_backbone(name='Pinopsida', kingdom='plants'))
as.data.frame(name_backbone(name='Gnetopsida', kingdom='plants'))
as.data.frame(name_backbone(name='Cycadopsida', kingdom='plants'))
as.data.frame(name_backbone(name='Ginkgoopsida', kingdom='plants'))



name_backbone(name='Magnoliopsida')$usageKey
name_backbone(name='Chaenorhinum')$usageKey





## Classes to download ####

clss <- c("Magnoliopsida", "Liliopsida", "Pinopsida",  "Gnetopsida", "Cycadopsida", "Ginkgoopsida")
#clss <- c("Liliopsida", "Pinopsida")
#clss <- c("Pinopsida")
clss <- as.data.frame(clss)

write.table(clss, "list_taxons.csv", sep = ", ", row.names = FALSE, col.names = FALSE)





## Downloading with PreSPickR  ####

t0 <- Sys.time()

GetBIF(credentials = paste0(gbif_creds, "/gbif_credentials.RData"),
       taxon_list = clss$clss,
       #taxon_list = "Chaenorhinum",
       #taxon_list = "Ophrys",
       #taxon_list = "Orchis",
       #taxon_list = "list_taxons.csv",
       #taxon_list = "Tracheophyta",
       download_format = "SIMPLE_CSV",
       download_years = c(2000, 2021),
       download_coords = c(-12.69141, 42.71485, 33.4901, 71.9218), #order: xmin, xmax, ymin, ymax
       download_coords_accuracy = c(0, 50),
       rm_dupl = TRUE,
       cols2keep = c("species", "decimalLatitude", "decimalLongitude", #"elevation",
                     "gbifID",
                     "coordinateUncertaintyInMeters",
                     "countryCode", "year", 
                     #"institutionCode",	"collectionCode",
                     #"ownerInstitutionCode",
                     "datasetKey"
                     ),
       out_name = paste0("/D5_FFGRCC_gbif_occ/", "sp_records_", format(Sys.Date(), "%Y%m%d")))

Sys.time() - t0




## Checking results ####

setwd(paste0(getwd(), "/D5_FFGRCC_gbif_occ/"))

## If the download was saved as csv with PreSPickR
#occs_all <- read.csv("sp_records_20210709.csv", header = TRUE)
occs_all <- fread("D5_FFGRCC_gbif_occ/sp_records_20210709.csv", header = TRUE)
#occs_all <- fread("sp_records_20211004_jeodpp_kk.csv", header = TRUE)

## To retrieve the raw data downloaded with PreSPickR::Pres_BIF()
## and create the simpler csv
## This doesn't work in Jeo-Desk... see below for a manual approach

taxon_dir <- "D:/xavi_rp/D5_FFGRCC_gbif_occ/"
taxon_dir <- "/eos/jeodpp/home/users/rotllxa/exploring_lucas_data"
taxons <- clss$clss
data1 <- Prep_BIF(taxon_dir = paste0(taxon_dir, "/"),
                  taxons = taxons,
                  cols2keep = c("species", "decimalLatitude", "decimalLongitude", #"elevation",
                                "gbifID",
                                "coordinateUncertaintyInMeters",
                                "countryCode", "year", 
                                #"institutionCode",	"collectionCode",
                                #"ownerInstitutionCode",
                                "datasetKey"
                  ),
                  #cols2keep = "all",
                  rm_dupl = TRUE)

print(paste0("Saving GBIF data as ", "/sp_records_20210709", ".csv"))
write.csv(data1, file = paste0("sp_records_20210709", ".csv"),
          quote = FALSE, row.names = FALSE)
#data1 <- fread("sp_records_20210709.csv", header = TRUE)
#occs_all <- fread("sp_records_20210709.csv", header = TRUE)



## manual approach to retrieve data

zips <- list.files()[grepl(".zip", list.files())]
data1 <- data.frame()

#for (z in zips[5]){
for (z in zips){
  print(paste0("processing... ", z))
  unzip(z, overwrite = FALSE)
  z1 <- gsub(".zip", "", z)
  data02 <- fread(paste0(z1, ".csv"), header = T)
  print(nrow(data02))
  data1 <- rbind(data1, data02)
}

head(data1)
nrow(data1)
unique(data1$class)
names(data1)


cols2keep = c("species", "decimalLatitude", "decimalLongitude", #"elevation",
              "gbifID",
              "coordinateUncertaintyInMeters",
              "countryCode", "year", 
              #"institutionCode",	"collectionCode",
              #"ownerInstitutionCode",
              "datasetKey")


data1_kk <- data1[, .SD, .SDcols = cols2keep]
    
sum(data1_kk$year == 2018)   
sort(unique(data1_kk$year))
sort(unique(data1_kk$countryCode))
length(unique(data1_kk$countryCode))

length(unique(data1_kk$datasetKey))

paste0(unique(data1_kk$datasetKey), collapse = "', '")


print(paste0("Saving GBIF data as ", getwd(), "/D5_FFGRCC_gbif_occ/sp_records_20211004", ".csv"))
write.csv(data1, file = paste0(getwd(), "/D5_FFGRCC_gbif_occ/sp_records_20211004", ".csv"),
          quote = FALSE, row.names = FALSE)
#



names(data1)
sum(unique(data1$datasetKey) == "7a3679ef-5582-4aaa-81f0-8c2545cafc81")  # Pl@ntNet
#


occs_all <- data1
occs_all <- fread(paste0(getwd(), "/D5_FFGRCC_gbif_occ/sp_records_20210709.csv"), header = TRUE)
nrow(occs_all)
names(occs_all)


sort(unique(occs_all$year))
sum(occs_all$year == 2018)


#cols_order <- c("sp2", "species", "decimalLatitude", "decimalLongitude", #"elevation",
cols_order <- c("species", 
                "decimalLatitude", "decimalLongitude",# "elevation",
                "gbifID",
                #"coordinateUncertaintyInMeters",
                "countryCode", "year"#,
                #"datasetKey"
                )

occs_all <- occs_all[, .SD, .SDcols = cols_order]

View(head(occs_all))
nrow(occs_all)

length(unique(occs_all$species))
nrow(occs_all[occs_all$species == "", ])
View(occs_all[occs_all$species == "", ])

## Removing rows with no species name 
# The lack of species name may be due to the fact that the occurrences are at genus level
occs_all <- occs_all[occs_all$species != "", ]
nrow(occs_all)


range(occs_all$coordinateUncertaintyInMeters)  

range(as.numeric(occs_all$decimalLatitude))
range(as.numeric(occs_all$decimalLongitude))

unique(occs_all$countryCode)
View(occs_all[occs_all$countryCode == "US", ])  # just 15 rows; can be removed
View(occs_all[(occs_all$decimalLatitude - occs_all$decimalLongitude) < 0.001, ])  # 1689 rows; kept!


sort(unique(occs_all$year))
sum(occs_all$year > 2009)
sum(occs_all$year < 2010)

table(occs_all$year)

(sum(occs_all$year >= 2010) * 100) / nrow(occs_all)  # 77% of records are newer than 2010 (included)
(sum(occs_all$year >= 2015) * 100) / nrow(occs_all)  # 52% of records are newer than 2015 (included)



unique(occs_all$institutionCode)
sort(unique(occs_all$datasetName))
unique(occs_all$ownerInstitutionCode)
length(unique(occs_all$datasetKey))


# plantnet "datasetKey":  "7a3679ef-5582-4aaa-81f0-8c2545cafc81"
# inaturalist "datasetKey":  "50c9509d-22c7-4a22-a47d-8c48425ef4a7"

sum(unique(occs_all$datasetKey) == "7a3679ef-5582-4aaa-81f0-8c2545cafc81")  # Pl@ntNet are there

occs_all





## LUCAS points - GPS ####
# LUCAS gps points are those really surveyed; might be different from the theoretical points

library(sf)
library(data.table)

working_year <- 2018

lucas_dir <- "D:/xavi_rp/LUCAS_Momo/LUCAS/output/"

## lucas geometry points (shapefile)
lucas_gps <- st_read(paste0(lucas_dir, "geometry/LUCAS_gps_geom/", "LUCAS_gps_geom.shp"))
lucas_gps

sort(unique(lucas_gps$YEAR))

lucas_gps_dt <- as.data.table(lucas_gps)
str(lucas_gps_dt)

# keeping working year (i.e. 2018)
lucas_gps_18 <- lucas_gps[lucas_gps$YEAR %in% working_year, ]   # 237728 points
lucas_gps_18
rm(lucas_gps)
gc()


# getting coordinates
lucas_gps_18_dt <- as.data.table(lucas_gps_18)
lucas_gps_18$geometry

lucas_gps_18_dt[, ":=" (x = st_coordinates(lucas_gps_18)[, 1], y =  st_coordinates(lucas_gps_18)[, 2])]
lucas_gps_18_dt



## lucas table
lucas_harmo_uf <- fread(paste0(lucas_dir, "table/lucas_harmo_uf.csv"), header = TRUE)
lucas_harmo_uf

View(head(lucas_harmo_uf))
names(lucas_harmo_uf)

lucas_harmo_uf[, 1:4]

sort(unique(lucas_harmo_uf$year))

# question: why the gps shapefile has less points (1144836) than the table 'lucas_harmo_uf' (1351293)


setkeyv(lucas_gps_dt, "POINT_ID")
lucas_gps_dt

lucas_harmo_uf_1_4 <- lucas_harmo_uf[, 1:4]
setkeyv(lucas_harmo_uf_1_4, "point_id")
lucas_harmo_uf_1_4

length(unique(lucas_harmo_uf_1_4$id))   # it seems like 'id' is a unique identifier
nrow(lucas_harmo_uf_1_4)

length(unique(lucas_gps_dt$ID))   # it seems like 'ID' is a unique identifier
nrow(lucas_gps_dt)


## keeping working year (i.e. 2018)
lucas_harmo_uf <- lucas_harmo_uf[year %in% working_year]
unique(lucas_harmo_uf$year)
lucas_harmo_uf[, 1:4]
lucas_harmo_uf
View(as.data.frame(names(lucas_harmo_uf)))


## LC classes

View(head(lucas_harmo_uf[, 28:47]))
kk <- lucas_harmo_uf[, 28:47] 
# 28 <- lc1
# 29 <- lc1_label
pet <- unique(kk[, list(lc1, lc1_label)])
View(pet[order(pet$lc1), ])
unique(lucas_harmo_uf$lc1_spec_label)

# artificial <- A
# forest <- C
# cropland <- B
# cereal <- B1
# root crops <- B2
# oleaginous and industrial crops <- B3
# vegetal and ornamental crops <- B4
# B51 <- Clovers
# B52 <- Lucerne
# B53 <- Other leguminous and mixtures for fodder
# B54 <- Mixed cereals for fodder
# B55 <- Temporary grasslands
# B7 <- Fruit trees
# B81 <- Olive groves
# B82 <- Vineyards
# B83 <- Nurseries
# B84 <- Permanent industrial crops
# BX1 <- Arable land (only pi)
# BX2 <- Permanent crops (only pi)
# C <- Forest
# D <- Shrubland
# E <- Grassland
head(lucas_harmo_uf[, list(point_id, lc1, lc1_label)])
names(lucas_harmo_uf)




## LUCAS points - Occurrences ####

lucas_gps_18_dt
lucas <- lucas_gps_18_dt[, .SD, .SDcols = c("ID", "x", "y")]
lucas
rm(lucas_gps_18_dt)
gc()


occs_all  # for now, keeping occurrences for all years (2000-2021)
#occs <- occs_all[, .SD, .SDcols = c("species", "gbifID", "decimalLongitude", "decimalLatitude")]
#names(occs) <- c("species", "gbifID", "x", "y")
setnames(occs_all, c("decimalLongitude", "decimalLatitude"), c("x", "y"))
occs_all


occs_sf <- st_as_sf(as.data.frame(occs_all), coords = c("x", "y"), crs = 4326)#, agr = "constant")
occs_sf
#rm(occs)
#gc()
#
#
#library(geosphere)
#z <- distHaversine(lucas[sample(1:237728, 10^3), 2:3], occs[sample(1:19384438, 10^4), 2:3])
#head(z)
#length(z)
#
#
#library(rgeos)
#as_Spatial(lucas_gps_18, cast = FALSE)
#
#t0 <- Sys.time()
#dist <- gWithinDistance(spgeom1 = as_Spatial(lucas_gps_18[sample(1:237728, 10^3), ], cast = FALSE), 
#                        spgeom2 = as_Spatial(occs_sf[sample(1:19384438, 10^4), ], cast = FALSE), 
#                        dist = 200, byid = TRUE,
#                        hausdorff = FALSE, densifyFrac = NULL)
#Sys.time() - t0
#
#
#View(head(dist, 100))
#sum(dist)
#
#
#
#
#library(raster)
#t0 <- Sys.time()
#d <- pointDistance(#lucas_gps_18, # cols
#                   p1 = lucas_gps_18[sample(1:237728, 10^3), ], # cols
#                   p2 = occs_sf[sample(1:19384438, 10^5), ],    # rows
#                   lonlat = TRUE)
#Sys.time() - t0
#
#View(head(d))
#nrow(d)
#
#
#d1 <- apply(d < 200, 2, which)
#head(d1)
#str(d1)
#
#
#occ_near <- function(l, o){
#        d <- pointDistance(#lucas_gps_18, # cols
#                l, # cols
#                o,    # rows
#                lonlat = TRUE)
#        d1 <- apply(d < 200, 2, which)
#        return(d1)
#}
#





# 1) buffering over lucas

lucas_sample <- lucas_gps_18[sample(1:237728, 10^3), ]
lucas_sample <- lucas_gps_18
occs_sample <- occs_sf[sample(1:22300974, 10^6), ]
occs_sample <- occs_sf[sample(1:22300974, 22300974/10), ]
occs_sample <- occs_sf

t0 <- Sys.time()
lucas_gps_18_buf <- st_buffer(lucas_sample, dist = 200)
Sys.time() - t0

lucas_gps_18_buf
#lucas_gps_18_buf <- st_read(paste0("lucas_gps_18_buf_200.shp"))


lucas_gps_18[lucas_gps_18$POINT_ID == 36942262, ]$geometry
lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID == 36942262, ]$geometry


plot(as_Spatial(lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID == 36942262, ], cast = FALSE))
plot(as_Spatial(lucas_gps_18[lucas_gps_18$POINT_ID == 36942262, ], cast = FALSE), add = TRUE)



# 2) joining simple features

t0 <- Sys.time()
occs_lucas_200 <- st_join(lucas_gps_18_buf, occs_sample)
Sys.time() - t0

occs_lucas_200
#occs_lucas_200 <- st_read(paste0("occs_lucas_200.shp"))


occs_lucas_200_dt <- as.data.table(occs_lucas_200)
occs_lucas_200_dt
sum(!is.na(occs_lucas_200_dt$species))

occs_lucas_200_dt_NotNA <- occs_lucas_200_dt[!is.na(occs_lucas_200_dt$species)]
length(unique(occs_lucas_200_dt_NotNA$POINT_ID))
# 5631 unique LUCAS with occurrences in their 200m buffer (total of 14470 occurrences)

occs_lucas_200_dt_na <- occs_lucas_200_dt[is.na(occs_lucas_200_dt$species)]
length(unique(occs_lucas_200_dt_na$POINT_ID))


plot(as_Spatial(lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID == 45603620, ], cast = FALSE))
plot(as_Spatial(lucas_gps_18[lucas_gps_18$POINT_ID == 45603620, ], cast = FALSE), add = TRUE)

plot(as_Spatial(occs_sample[occs_sample$species == "Fagus sylvatica", ], cast = FALSE), add = TRUE, col = "red")



## writing out shapefiles
st_write(lucas_gps_18_buf, paste0("lucas_gps_18_buf_200.shp"))
st_write(occs_lucas_200, paste0("occs_lucas_200.shp"))




## some checks
#occs_lucas_200_dt <- fread("occs_lucas_200buffer.csv", header = TRUE)

occs_lucas_200_dt <- as.data.table(occs_lucas_200)
occs_lucas_200_dt
length(unique(occs_lucas_200_dt$POINT_ID))
sum(!is.na(occs_lucas_200_dt$species))

occs_lucas_200_dt_na <- occs_lucas_200_dt[is.na(occs_lucas_200_dt$species)]
length(unique(occs_lucas_200_dt_na$POINT_ID))
# 216138 unique LUCAS (2018) with NO occurrences in their 200m buffer

occs_lucas_200_dt_NotNA <- occs_lucas_200_dt[!is.na(occs_lucas_200_dt$species)]
occs_lucas_200_dt_NotNA[, geometry := NULL]
write.csv(occs_lucas_200_dt_NotNA, "occs_lucas_200buffer.csv", row.names = FALSE)
occs_lucas_200_dt_NotNA <- fread("occs_lucas_200buffer.csv", header = TRUE)

length(unique(occs_lucas_200_dt_NotNA$POINT_ID))
length(unique(occs_lucas_200_dt_NotNA$gbifID))
# 21590 unique LUCAS (2018; over a total of 237728 LUCAS points) with occurrences in their 200m buffer (total of 165446 occurrences)

# How many occurrences do NOT fall into any LUCAS buffer??
22300974 - 165446  # 22135528


#plot(as_Spatial(occs_lucas_200[occs_lucas_200$POINT_ID == 44563818 , ], cast = FALSE))
plot(st_geometry(occs_lucas_200[occs_lucas_200$POINT_ID == 44563818 , ]))

occs2plot <- occs_lucas_200[occs_lucas_200$POINT_ID == 44563818 , ]$gbifID
for(i in occs2plot){
        plot(st_geometry(occs_sample[occs_sample$gbifID %in% i, ]), 
             add = TRUE, 
             col = randomcoloR::randomColor())
}



## LUCAS points with GBIF occurrences in buffer, by country

occs_lucas_200_dt_NotNA
occs_all

# merging with GBIF country code
occs_lucas_200_dt_NotNA[occs_all, on = "gbifID", countryCode := i.countryCode]
occs_lucas_200_dt_NotNA

# checks
occs_lucas_200_dt_NotNA[POINT_ID == 26381954, ]
occs_lucas_200_dt[POINT_ID == 26381954, ]
occs_all[gbifID == "2839575251", ]

# unique lucas per country
occs_lucas_200_dt_NotNA

unique_lcs <- unique(occs_lucas_200_dt_NotNA, by = "POINT_ID")
unique_lcs
unique_lcs_tbl <- as.data.table(table(unique_lcs$countryCode))
library(countrycode)
unique_lcs_tbl[, country := countrycode(as.character(unique_lcs_tbl$V1), "iso2c", "iso.name.en")]
unique_lcs_tbl <- unique_lcs_tbl[, .SD, .SDcols = c("country", "V1", "N")]
setnames(unique_lcs_tbl, c("country", "V1", "N"), c("Country", "countryCode", "N_lucas"))
unique_lcs_tbl <- rbind(unique_lcs_tbl, list("Total", "Total", sum(unique_lcs_tbl$N_lucas)))
unique_lcs_tbl




## plot a map

dev.off()
library(rworldmap)
wrld_map <- getMap()
#worldmap <- rnaturalearth::ne_download(scale = 110,
#                                       type = "countries",
#                                       category = "cultural",
#                                       destdir = tempdir(),
#                                       load = TRUE,
#                                       returnclass = "sf")
#
#library(ggplot2)
#mp2plot <- ggplot() +
#        geom_sf(data = lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID %in% unique_lcs$POINT_ID, ], colour = "red") +
#        geom_sf(data = worldmap, colour = "grey47", fill = NA) +
#        annotation_custom(tableGrob(unique_lcs_tbl[, list(Country, N_lucas)]), xmin = 35, xmax = 50, ymin = 32, ymax = 70)
#
#jpeg(paste0("LucasWithOccurrences.jpg"), width = 23, height = 18, units = "cm", res = 300)
#plot(mp2plot)
#dev.off()



library(gridExtra)
library(grid)

jpeg(paste0("LucasWithOccurrences.jpg"), width = 23, height = 21, units = "cm", res = 300)
par(mar = c(6, 3, 7, 12), bty = "n")
par(xpd = FALSE)
plot(st_geometry(lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID %in% unique_lcs$POINT_ID, ]), col = "red", border = "red")
plot(wrld_map, border = "grey47", add = TRUE)

pushViewport(viewport(x = 0.88, y = 0.45))
thm <- ttheme_default(
        core = list(
                fg_params = list(fontface = c(rep("plain", 50), "bold.italic"),
                                 cex = 0.4),
                bg_params = list(fill = c(rep(c("grey95", "grey90"), length.out = 50)))#,
                #alpha = rep(c(1, 0.5), each = 5))
                ),
        colhead = list(fg_params = list(cex = 0.4)),
        rowhead = list(fg_params = list(cex = 0.4))
)
grid.draw(tableGrob(unique_lcs_tbl[, list(Country, N_lucas)], 
                    theme = thm))

title(main = paste0("LUCAS (2018) With GBIF Occurrences"),
      outer = TRUE,
      line = - 4.3,
      cex.main = 1.5)
mtext("",
      side = 1, line = 3, 
      #at = 5,
      adj = 0,
      #cex = 1)
      cex = 0.8)

dev.off()      
#



## Number of occurrences per LUCAS

occs_lucas_200_dt_NotNA

Noccs_lucas <- occs_lucas_200_dt_NotNA[, list("n_total" = .N), by = .(POINT_ID)]
Noccs_lucas

summary(Noccs_lucas$N)

setkeyv(Noccs_lucas, "n_total")
Noccs_lucas

occs_lucas_200_dt_NotNA[POINT_ID == 42743614]   # DK
length(unique(occs_lucas_200_dt_NotNA[POINT_ID == 42743614]$species)) # 186 unique species
sort(unique(occs_lucas_200_dt_NotNA[POINT_ID == 42743614]$species))

occs_lucas_200_dt_NotNA[POINT_ID == 28822192]   # PT
length(unique(occs_lucas_200_dt_NotNA[POINT_ID == 28822192]$species)) # 666 unique species
sort(unique(occs_lucas_200_dt_NotNA[POINT_ID == 28822192]$species))

occs_lucas_200_dt_NotNA[POINT_ID == 44803630]   # DK
length(unique(occs_lucas_200_dt_NotNA[POINT_ID == 44803630]$species)) # 192 unique species
sort(unique(occs_lucas_200_dt_NotNA[POINT_ID == 44803630]$species))


# Removing duplicated species within the same LUCAS
occs_lucas_200_dt_NotNA
occs_lucas_200_dt_NotNA
occs_lucas_200_dt_NotNA_NotDuplSP <- unique(occs_lucas_200_dt_NotNA, by = c("POINT_ID", "species"))
occs_lucas_200_dt_NotNA_NotDuplSP
length(unique(occs_lucas_200_dt_NotNA_NotDuplSP[POINT_ID == 44803630]$species)) # 192 unique species

Noccs_lucas1 <- occs_lucas_200_dt_NotNA_NotDuplSP[, list("n_unique" = .N), by = .(POINT_ID)]


Noccs_lucas <- Noccs_lucas[Noccs_lucas1, on = "POINT_ID"]
setkeyv(Noccs_lucas, "n_total")
Noccs_lucas

# saving the table
write.csv(Noccs_lucas, "Noccs_lucas.csv", row.names = FALSE)


round(summary(Noccs_lucas$n_unique), 0)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    1       1       1       4       3     666 
table(Noccs_lucas$n_unique)
hist(Noccs_lucas$n_unique, breaks = 20)

round(summary(Noccs_lucas$n_total), 0)
hist(Noccs_lucas$n_total, breaks = 20)

(sum(Noccs_lucas$n_unique == 1) * 100) / nrow(Noccs_lucas) # 51% of LUCAS have 1 occurrence
(sum(Noccs_lucas$n_unique == 2) * 100) / nrow(Noccs_lucas) # 18% of LUCAS have 2 occurrences
(sum(Noccs_lucas$n_unique == 3) * 100) / nrow(Noccs_lucas) # 9% of LUCAS have 3 occurrences
(sum(Noccs_lucas$n_unique == 4) * 100) / nrow(Noccs_lucas) # 5% of LUCAS have 4 occurrences
(sum(Noccs_lucas$n_unique == 5) * 100) / nrow(Noccs_lucas) # 3% of LUCAS have 5 occurrences
(sum(Noccs_lucas$n_unique < 6) * 100) / nrow(Noccs_lucas)  # 86% of LUCAS have 5 or less occurrences
(sum(Noccs_lucas$n_unique < 11) * 100) / nrow(Noccs_lucas) # 93% of LUCAS have 10 or less occurrences
(sum(Noccs_lucas$n_unique < 21) * 100) / nrow(Noccs_lucas) # 96% of LUCAS have 20 or less occurrences




## saving the final csv and a snapshot
occs_lucas_200_dt_NotNA_NotDuplSP
write.csv(occs_lucas_200_dt_NotNA_NotDuplSP, "occs_lucas_200buffer_NotDuplicatedSp.csv", row.names = FALSE)


jpeg(paste0("final_dataset_snapshot.jpg"), width = 20, height = 10, units = "cm", res = 300)
grid.table(head(occs_lucas_200_dt_NotNA_NotDuplSP, 10))
dev.off()



## Some maps

# LUCAS with a certain number (e.g. 10) or more occurrences within the buffer
NoccsInBuff <- 10
clr <- "darkgreen"
NoccsInBuff <- 5
clr <- "blue"


lcs_more10occ <- Noccs_lucas[n_unique >= NoccsInBuff, POINT_ID]
length(lcs_more10occ)

unique_lcs1 <- unique(occs_lucas_200_dt_NotNA, by = "POINT_ID")
unique_lcs1 <- unique_lcs1[POINT_ID %in% lcs_more10occ]
unique_lcs1_tbl <- as.data.table(table(unique_lcs1$countryCode))
#library(countrycode)
unique_lcs1_tbl[, country := countrycode(as.character(unique_lcs1_tbl$V1), "iso2c", "iso.name.en")]
unique_lcs1_tbl <- unique_lcs1_tbl[, .SD, .SDcols = c("country", "V1", "N")]
setnames(unique_lcs1_tbl, c("country", "V1", "N"), c("Country", "countryCode", "N_lucas"))
unique_lcs1_tbl <- rbind(unique_lcs1_tbl, list("Total", "Total", sum(unique_lcs1_tbl$N_lucas)))
unique_lcs1_tbl


jpeg(paste0("Lucas_more", NoccsInBuff, "_Occs.jpg"), width = 23, height = 21, units = "cm", res = 300)
par(mar = c(6, 3, 7, 12), bty = "n")
par(xpd = FALSE)
plot(st_geometry(lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID %in% lcs_more10occ, ]), col = clr, border = clr)
plot(wrld_map, border = "grey47", add = TRUE)

pushViewport(viewport(x = 0.88, y = 0.45))
thm <- ttheme_default(
  core = list(
    fg_params = list(fontface = c(rep("plain", 50), "bold.italic"),
                     cex = 0.4),
    bg_params = list(fill = c(rep(c("grey95", "grey90"), length.out = 50)))#,
    #alpha = rep(c(1, 0.5), each = 5))
  ),
  colhead = list(fg_params = list(cex = 0.4)),
  rowhead = list(fg_params = list(cex = 0.4))
)
grid.draw(tableGrob(unique_lcs1_tbl[, list(Country, N_lucas)], 
                    theme = thm))

title(main = paste0("LUCAS (2018) with ", NoccsInBuff, " or more GBIF Occurrences"),
      outer = TRUE,
      line = - 4.3,
      cex.main = 1.5)
mtext("",
      side = 1, line = 3, 
      #at = 5,
      adj = 0,
      #cex = 1)
      cex = 0.8)

dev.off()      
#



## Adding land cover class info

lucas_lcc <- lucas_harmo_uf[, list(point_id, lc1, lc1_label)]

occs_lucas_200_dt_NotNA_NotDuplSP <- merge(occs_lucas_200_dt_NotNA_NotDuplSP, lucas_lcc, by.x = "POINT_ID", by.y = "point_id", all.x = TRUE)
occs_lucas_200_dt_NotNA_NotDuplSP


## Adding pictures info

cls2keep <- c("point_id",
              "photo_point", "photo_north", "photo_south", "photo_east", "photo_west",           
              "file_path_gisco_point", "file_path_gisco_north", "file_path_gisco_south", "file_path_gisco_east", "file_path_gisco_west")
lucas_photos <- lucas_harmo_uf[, ..cls2keep]
View(head(lucas_photos))

occs_lucas_200_dt_NotNA_NotDuplSP <- merge(occs_lucas_200_dt_NotNA_NotDuplSP, lucas_photos, by.x = "POINT_ID", by.y = "point_id", all.x = TRUE)
occs_lucas_200_dt_NotNA_NotDuplSP


# examples for lucas with 10 or more occurrences
head(occs_lucas_200_dt_NotNA_NotDuplSP[POINT_ID %in% lcs_more10occ & grepl("B", occs_lucas_200_dt_NotNA_NotDuplSP$lc1)])

occs_lucas_200_dt_NotNA_NotDuplSP[POINT_ID == 26601766, c(4, 14)]
occs_lucas_200_dt_NotNA_NotDuplSP[POINT_ID == 26601766, c(4, 14)]$species

#


## saving the final csv and a snapshot
occs_lucas_200_dt_NotNA_NotDuplSP
write.csv(occs_lucas_200_dt_NotNA_NotDuplSP, "occs_lucas_200buffer_NotDuplicatedSp.csv", row.names = FALSE)


jpeg(paste0("final_dataset_snapshot.jpg"), width = 25, height = 10, units = "cm", res = 300)
grid.table(head(occs_lucas_200_dt_NotNA_NotDuplSP[, 1:8], 10))
dev.off()


## plotting one LUCAS with occurrences, as example
jpeg(paste0("example_LucasBuff_species.jpg"), width = 25, height = 15, units = "cm", res = 300)
plot(st_geometry(lucas_gps_18_buf[lucas_gps_18_buf$POINT_ID == 26601766 , ]))
occs2plot <- occs_lucas_200_dt_NotNA_NotDuplSP[occs_lucas_200_dt_NotNA_NotDuplSP$POINT_ID == 26601766 , ]$gbifID
occs2plot_names <- occs_lucas_200_dt_NotNA_NotDuplSP[occs_lucas_200_dt_NotNA_NotDuplSP$POINT_ID == 26601766 , ]$species
cl <- randomcoloR::randomColor(length(occs2plot))
for(i in 1:length(occs2plot)){
  plot(st_geometry(occs_sf[occs_sf$gbifID %in% occs2plot[i], ]), 
       add = TRUE, 
       col = cl[i])
}
legend("right",
       legend = occs2plot_names,
       col = cl,
       fill = cl)
dev.off()
#







## Checks for Catalonia
cat_coords <-  c(xmin = 0, ymin = 40.3, xmax = 3.5, ymax = 43)   # Catalonia


lucas_gps_18_buf_cat <- st_crop(lucas_gps_18_buf, cat_coords) 
occs_lucas_200_cat <- st_crop(occs_lucas_200, cat_coords)



plot(occs_lucas_200_cat)








## Natura200 data ####

library(sf)

Natura2000_end2020_epsg3035 <- st_read("Natura2000_end2020_shp/Natura2000_end2020_epsg3035.shp")
Natura2000_end2020_epsg3035

names(Natura2000_end2020_epsg3035)


natura2000_csv <- list.files("Natura2000_end2020_csv", full.names = TRUE)

for (i in 9:length(natura2000_csv)){
  print(i)
  print(natura2000_csv[i])
  print(names(read.csv(natura2000_csv[i], header = TRUE)))
  print("")
}



## Habitats ("...HABITATCLASS.csv")
i <- 4
natura2000_csv_n4 <- read.csv(natura2000_csv[i], header = TRUE)
head(natura2000_csv_n4)

View(natura2000_csv_n4[natura2000_csv_n4$SITECODE == "CZ0523284", ])
sum(natura2000_csv_n4[natura2000_csv_n4$SITECODE == "CZ0523284", 4])
length(natura2000_csv_n4$SITECODE)

sort(unique(natura2000_csv_n4$HABITATCODE))
# [1] "N01" "N02" "N03" "N04" "N05" "N06" "N07" "N08" "N09" "N10" "N11" "N12" "N13" "N14" "N15"
#[16] "N16" "N17" "N18" "N19" "N20" "N21" "N22" "N23" "N24" "N25" "N26" "N27" "N6"  "N9" 


head(natura2000_csv_n4[natura2000_csv_n4$HABITATCODE == "N9", ], 10)
sum(natura2000_csv_n4$HABITATCODE == "N9")                             # N9 only one case, seems an error
head(natura2000_csv_n4[natura2000_csv_n4$HABITATCODE == "N09", ], 10)
sum(natura2000_csv_n4$HABITATCODE == "N09")
head(natura2000_csv_n4[natura2000_csv_n4$HABITATCODE == "N6", ], 10)    # N6 only one case, seems an error (same SITECODE)



## Habitats ("...HABITATS.csv")
i <- 5
natura2000_csv_n5 <- read.csv(natura2000_csv[i], header = TRUE)
head(natura2000_csv_n5)
names(natura2000_csv_n5)  
nrow(natura2000_csv_n5) # 151305
nrow(natura2000_csv_n4) # 139499

sort(unique(natura2000_csv_n5$HABITATCODE))
length(unique(natura2000_csv_n5$HABITATCODE))  # 237 classes, different from previous csv
                                               # share of each habitat is not provided, only the area


length(natura2000_csv_n5$`?..SITECODE`)
head(natura2000_csv_n5$`?..SITECODE`)
nrow(natura2000_csv_n5[natura2000_csv_n5$`?..SITECODE` == "CZ0523284", ])

sum(!natura2000_csv_n4$SITECODE %in% natura2000_csv_n5$`?..SITECODE`)  # 24000 
sum(!natura2000_csv_n5$`?..SITECODE` %in% natura2000_csv_n4$SITECODE)  #  1208

natura2000_csv_n5[natura2000_csv_n4$SITECODE %in% natura2000_csv_n5$`?..SITECODE`, ][100000, ]
View(head(natura2000_csv_n5[natura2000_csv_n4$SITECODE %in% natura2000_csv_n5$`?..SITECODE`, ], 50))
View((natura2000_csv_n4[natura2000_csv_n4$SITECODE == "ES0000067", ]))
View((natura2000_csv_n5[natura2000_csv_n5$`?..SITECODE` == "ES0000067", ]))
sum(natura2000_csv_n5[natura2000_csv_n5$`?..SITECODE` == "ES0000067", "COVER_HA"])   # 62253.97



## Habitats ("...NATURA2000SITES.csv")
i <- 9
natura2000_csv_n9 <- read.csv(natura2000_csv[i], header = TRUE)
head(natura2000_csv_n9)
names(natura2000_csv_n9)  
nrow(natura2000_csv_n9)  

natura2000_csv_n9[natura2000_csv_n9$SITECODE == "ES0000067", ]
head(natura2000_csv_n9$SITECODE)




## Habitats ("...IMPACT.csv")
i <- 6
natura2000_csv_n6 <- read.csv(natura2000_csv[i], header = TRUE)
head(natura2000_csv_n6)
nrow(natura2000_csv_n6)

unique(natura2000_csv_n6$POLLUTIONCODE)
unique(natura2000_csv_n6$DESCRIPTION)
unique(natura2000_csv_n6$IMPACT_TYPE)
unique(natura2000_csv_n6$IMPACTCODE)





## Habitats ("...SPECIES.csv")
i <- 11
natura2000_csv_n11 <- read.csv(natura2000_csv[i], header = TRUE)
head(natura2000_csv_n11)
names(natura2000_csv_n11)
nrow(natura2000_csv_n11)

length(unique(natura2000_csv_n11$SITECODE))

View(natura2000_csv_n11[natura2000_csv_n11$SITECODE == "ES0000067", ])


unique(natura2000_csv_n11$SPGROUP)
nrow(natura2000_csv_n11[natura2000_csv_n11$SPGROUP == "Plants", ])
View(head(natura2000_csv_n11[natura2000_csv_n11$SPGROUP == "Plants", ]))
View((natura2000_csv_n11[natura2000_csv_n11$SPGROUP == "Plants", ]))

nrow(natura2000_csv_n11[natura2000_csv_n11$SPGROUP == "", ])
View((natura2000_csv_n11[natura2000_csv_n11$SPGROUP == "", ]))



## Habitats ("...OTHERSPECIES.csv")
i <- 10
natura2000_csv_n10 <- read.csv(natura2000_csv[i], header = TRUE)
View(head(natura2000_csv_n10, 50))
names(natura2000_csv_n10)
nrow(natura2000_csv_n10)







## Crop Map ####

list.files("/eos/jeodpp/home/users/rotllxa")
list.files("/eos/jeodpp/home/users/rotllxa/data")
list.files("/storage/rotllxa")
list.files("/storage/rotllxa/Documents")

getwd()
list.files("/home/rotllxa/")
list.files("/home/rotllxa/Documents/")

list.files("/eos/jeodpp/data/base/")
list.files("/eos/jeodpp/data/base/", recursive = TRUE)


list.files("/mnt/cidstorage/cidportal/data/OpenData/EUCROPMAP/")
list.files("/mnt/cidstorage/cidportal/data/OpenData/EUCROPMAP/2018")

list.files("/eos/jeodpp/data/projects/REFOCUS/data/BIODIVERSITY/")
list.files("/eos/jeodpp/data/projects/REFOCUS/data/BIODIVERSITY/DataRestoration")
list.files("/eos/jeodpp/data/projects/REFOCUS/data/BIODIVERSITY/DataRestoration/cropmap_res")


cropmap2018 <- raster("/eos/jeodpp/data/projects/REFOCUS/data/BIODIVERSITY/DataRestoration/cropmap_res/eucropmap_res.tif")  # at 1km
cropmap2018 <- raster("/mnt/cidstorage/cidportal/data/OpenData/EUCROPMAP/2018/EUCROPMAP_2018.tif")  # at 10m
cropmap2018 <- raster("/eos/jeodpp/home/users/rotllxa/exploring_lucas_data/eucropmap_2018/EUCROPMAP_2018.tif")  # at 10m

#writeRaster(cropmap2018, "cropmap2018_10m.tif")

sp::CRS("+init=EPSG:3035")
crs(cropmap2018) <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

cropmap2018
plot(cropmap2018)

cat_coords <-  c(3500000, 4100000, 1850000, 2400000)   # Catalonia (LAEA, m) (xmin, xmax, ymin, ymax)

cropmap2018_cat <- crop(cropmap2018, extent(cat_coords), filename = "cropmap2018_cat.tif")
plot(cropmap2018_cat)

#pixac_v7_byte_masked <- brick("/mnt/cidstorage/cidportal/data/OpenData/EUCROPMAP/2018/pixac_v7_byte_masked.tif")
#pixac_v7_byte_masked


# CropMap classes
#library(readr)
#cropmap_classes_2018 <- read_file("/mnt/cidstorage/cidportal/data/OpenData/EUCROPMAP/2018/EuroCropMap.qml")

crops_categs <- c(100, 211, 212, 213, 214, 215, 216, 217, 218, 219, 221, 222, 223, 230, 231, 232, 233, 240, 250, 290, 300, 500, 600, 800)
crops_names <- c("Artificial", "Common wheat", "Durum wheat", "Barley", "Rye", "Oats", "Maize", "Rice", "Triticale", "Other cereals", "Potatoes", "Sugar beet", "Other root crops", "Other non permanent industrial crops", "Sunflower", "Rape and turnip rape", "Soya", "Dry pulses", "Fodder crops (cereals and leguminous)", "Bare arable land", "Woodland and Shrubland (incl. permanent crops)", "Grasslands", "Bare land", "Wetlands")

cropmap_classes <- data.frame("crop_categ" = crops_categs, "crop_names" = crops_names)
View(cropmap_classes)

# Maiz: 216


## Aggregating maiz to 1km or 10km

#cropmap2018_maiz <- cropmap2018
#cropmap2018_maiz[cropmap2018_maiz$EUCROPMAP_2018 != 216] <- 0

aggr_fun_1km <- function(x, ...) {     # returns share of maize at 1km (0 to 1)
  mz_share <- sum(x == 216, na.rm = TRUE) / 10000
  return(mz_share)
}

cropmap2018_maiz_1km <- aggregate(x = cropmap2018, 
                                  fact = 100,        # 1km
                                  fun = aggr_fun_1km, 
                                  expand = TRUE, 
                                  na.rm = TRUE, 
                                  filename = "cropmap2018_maiz_1km.tif",
                                  overwrite = TRUE)

cropmap2018_maiz_1km



## Maiz weeds ####

weeds_maiz <- read.csv("../weeds/weeds_maize_report_2011.csv", header = TRUE)
head(weeds_maiz)
nrow(weeds_maiz)


occs_all_2018 <- occs_all[occs_all$year == 2018, ]

nrow(occs_all_2018)
sum(occs_all_2018$species == "")
length(unique(occs_all_2018$species))


occs_2018_specie <- unique(occs_all_2018$species)

head(sort(occs_2018_specie), 50)
head(sort(weeds_maiz$Species))

sum(occs_2018_specie %in% weeds_maiz$Species)   # 151 sp
sum(weeds_maiz$Species %in% occs_2018_specie)   # 151 sp
sum(!weeds_maiz$Species %in% occs_2018_specie)   # 53 sp


weeds_maiz_gbib <- sort(weeds_maiz$Species[weeds_maiz$Species %in% occs_2018_specie])
weeds_maiz_not_gbib <- sort(weeds_maiz$Species[!weeds_maiz$Species %in% occs_2018_specie])


occs_all_2018_maiz <- occs_all_2018[occs_all_2018$species %in% weeds_maiz_gbib, ]

nrow(occs_all_2018_maiz)   # 158427 occurrences for 2018
nrow(occs_all_2018)        # over 1591221 in total for 2018




