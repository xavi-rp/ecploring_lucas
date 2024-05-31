
if(Sys.info()[4] == "D01RI1700308") {
  wd <- "D:/xavi_rp/D5_FFGRCC_lucas_grasslands/"
}else if(Sys.info()[4] == "S-JRCIPRAP320P") {
  wd <- "D:/rotllxa/D5_FFGRCC_lucas_grasslands/"
}else if(Sys.info()[4] %in% c("jeodpp-terminal-jd001-03", "jeodpp-terminal-03", "jeodpp-terminal-dev-12", 
                              "jeodpp-terminal-jd002-03", "jeodpp-terminal-jd004-03.cidsn.jrc.it",
                              "jeodpp-terminal-jd001-03.cidsn.jrc.it",
                              "jeodpp-terminal-dev-jd002-12.cidsn.jrc.it")) {
  if(!dir.exists("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/")) 
    dir.create("/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/")
  wd <- "/eos/jeodpp/home/users/rotllxa/lucas_grassland_data/"
  gbif_creds <- "/home/rotllxa/Documents/"
  WhereAmI <- "bdap"
}else if(Sys.info()[4] == "MacBook-MacBook-Pro-de-Xavier.local"){
  wd <- "/Users/xavi_rp/Documents/D5_FFGRCC/lucas_grasslands_data/"
  WhereAmI <- "mac"
}else{
  wd <- "C:/Users/rotllxa/D5_FFGRCC_lucas_grasslands/"
  gbif_creds <- "C:/Users/rotllxa/Documents/"
}

setwd(wd)


library(data.table)
library(sf)
library(tidyverse)




## Final general directory ####

list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/")



## Anonimized directory ####
# The 'anonymised' folder contains EXIF, expert/surveyors surveys and releve files 
# with no geographical information besides country centroid. These are the ones to be 
# published with the data paper.

list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised/")
# "LUCAS_G_2018_EXIF"       "LUCAS_G_2018_KEYSPECIES" "LUCAS_G_2018_RELEVE" 

list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised/LUCAS_G_2018_EXIF")
#  "expertDF_exif.csv"          "LUCAS_2018_Gr_all_EXIF.csv" "surveyorDF_exif.csv"


## 

expertDF_exif <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_EXIF/",
                              "expertDF_exif.csv"))
names(expertDF_exif)
ncol(expertDF_exif)  # 47
nrow(expertDF_exif)  # 2693

surveyorDF_exif <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_EXIF/",
                              "surveyorDF_exif.csv"))
names(surveyorDF_exif)
ncol(surveyorDF_exif)  # 43
nrow(surveyorDF_exif)  # 26513

names(expertDF_exif)[!names(expertDF_exif) %in% names(surveyorDF_exif)] 
names(surveyorDF_exif)[!names(surveyorDF_exif) %in% names(expertDF_exif)] 

unique(surveyorDF_exif$GPSDestDistanceRef)

# removing POINT_ID and adding an ad-hoc unique identifier
head(unique(expertDF_exif$POINT_ID))
head(unique(surveyorDF_exif$POINT_ID))

length(unique(expertDF_exif$POINT_ID))
length(unique(surveyorDF_exif$POINT_ID))

sum(!unique(expertDF_exif$POINT_ID) %in% unique(surveyorDF_exif$POINT_ID)) 
sum(!unique(surveyorDF_exif$POINT_ID) %in% unique(expertDF_exif$POINT_ID))

point_id_all <- c()
point_id_all <- unique(expertDF_exif$POINT_ID)
length(point_id_all)

point_id_all <- c(point_id_all,
                  unique(surveyorDF_exif$POINT_ID)[!unique(surveyorDF_exif$POINT_ID) %in% unique(expertDF_exif$POINT_ID)])
length(point_id_all)
length(unique(point_id_all))
#point_id_all <- data.table("POINT_ID" = point_id_all)
point_id_all


LUCAS_2018_Gr_all_EXIF <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised/LUCAS_G_2018_EXIF/",
                                "LUCAS_2018_Gr_all_EXIF.csv"))
names(LUCAS_2018_Gr_all_EXIF)
ncol(LUCAS_2018_Gr_all_EXIF)  # 52
nrow(LUCAS_2018_Gr_all_EXIF)  # 29205
View(LUCAS_2018_Gr_all_EXIF)  

names(LUCAS_2018_Gr_all_EXIF)[!names(LUCAS_2018_Gr_all_EXIF) %in% names(surveyorDF_exif) &
                                !names(LUCAS_2018_Gr_all_EXIF) %in% names(expertDF_exif) ]   # "survey_type"
unique(LUCAS_2018_Gr_all_EXIF$survey_type)
unique(LUCAS_2018_Gr_all_EXIF$V1)
unique(surveyorDF_exif$V1)
unique(expertDF_exif$V1)


## 

list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs")
#  "estat_e_attr_point_allattr_new_anonym.csv"    "estat_s_attr_point_allattr_new_anonym.csv"

estat_e_attr_point_allattr_new_anonym <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/", 
                                               "estat_e_attr_point_allattr_new_anonym.csv"))

names(estat_e_attr_point_allattr_new_anonym)
head(unique(estat_e_attr_point_allattr_new_anonym$POINT_ID))

sum(!unique(estat_e_attr_point_allattr_new_anonym$POINT_ID) %in% point_id_all)


estat_s_attr_point_allattr_new_anonym <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/", 
                                               "estat_s_attr_point_allattr_new_anonym.csv"))

names(estat_s_attr_point_allattr_new_anonym)
sum(!unique(estat_s_attr_point_allattr_new_anonym$POINT_ID) %in% point_id_all)

point_id_all <- c(point_id_all,
                  unique(estat_s_attr_point_allattr_new_anonym$POINT_ID)[!unique(estat_s_attr_point_allattr_new_anonym$POINT_ID) %in% point_id_all])
length(point_id_all)


##
list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_RELEVE/CSVs/")

releve_data_all <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_RELEVE/CSVs/", 
                                "releve_data_all.csv"), 
                         header = TRUE)
names(releve_data_all)
ncol(releve_data_all) # 730
nrow(releve_data_all) # 2672

names(releve_data_all)[-c(1:2)]
sum(!names(releve_data_all)[-c(1:2)] %in% point_id_all)


# creating anonymised Point_ID
point_id_all <- data.table(POINT_ID = point_id_all, POINT_ID_anonym = paste0("LG18_point_", seq(1:length(point_id_all))))
point_id_all

write.csv(point_id_all, 
          "/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/anonymised_POINT_ID.csv", 
          row.names = FALSE)




## Non-public directory ####
# This folder contains EXIF, expert/surveyors surveys and releve files 
# with geographical information, thus the entire curated data  set

list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/")


##

list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_EXIF/")


expertDF_exif_1 <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/LUCAS_G_2018_EXIF/",
                              "expertDF_exif.csv"))
names(expertDF_exif_1)
ncol(expertDF_exif_1)  # 60
nrow(expertDF_exif_1)  # 2693


surveyorDF_exif_1 <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/LUCAS_G_2018_EXIF/",
                                "surveyorDF_exif.csv"))
names(surveyorDF_exif_1)
ncol(surveyorDF_exif_1)  # 55
nrow(surveyorDF_exif_1)  # 26513


names(expertDF_exif_1)[!names(expertDF_exif_1) %in% names(surveyorDF_exif_1)]
names(surveyorDF_exif_1)[!names(surveyorDF_exif_1) %in% names(expertDF_exif_1)]


##
list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES")
list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/CSVs/")

estat_e_attr_point_allattr_new <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/CSVs/",
                                  "estat_e_attr_point_allattr_new.csv"))

names(estat_e_attr_point_allattr_new)
ncol(estat_e_attr_point_allattr_new)  # 126
nrow(estat_e_attr_point_allattr_new)  # 605


list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/GPKGs/")

estatdb_allattr_e_point <- st_read("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/GPKGs/estatdb_allattr_e_point.gpkg")
estatdb_allattr_e_point
names(estatdb_allattr_e_point)
nrow(estatdb_allattr_e_point) # 605 
ncol(estatdb_allattr_e_point) # 127 (126 + 1, which is the geometries)

estatdb_allattr_s_point <- st_read("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/GPKGs/estatdb_allattr_s_point.gpkg")
estatdb_allattr_s_point
names(estatdb_allattr_s_point)
nrow(estatdb_allattr_s_point) # 2622 
ncol(estatdb_allattr_s_point) # 127 (126 + 1, which is the geometries)


estatdb_allattr_e_transects <- st_read("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/GPKGs/estatdb_allattr_e_transects.gpkg")
estatdb_allattr_e_transects
names(estatdb_allattr_e_transects)
nrow(estatdb_allattr_e_transects) # 78 
ncol(estatdb_allattr_e_transects) # 133 (132 + 1, which is the geometries)


estatdb_allattr_s_transects <- st_read("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_KEYSPECIES/GPKGs/estatdb_allattr_s_transects.gpkg")
estatdb_allattr_s_transects
names(estatdb_allattr_s_transects)
nrow(estatdb_allattr_s_transects) # 2617 
ncol(estatdb_allattr_s_transects) # 133 (132 + 1, which is the geometries)



list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_RELEVE/CSVs/")

releve_data_all <- fread("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/LUCAS_G_2018_RELEVE/CSVs/releve_data_all.csv",
                         header = TRUE)
releve_data_all
names(releve_data_all)
nrow(releve_data_all)
ncol(releve_data_all)


list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/")
l2018_gr_recordDescriptor <- fread("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/l2018_gr_recordDescriptor.csv", 
                                   header = TRUE)
View(l2018_gr_recordDescriptor)
nrow(l2018_gr_recordDescriptor)
ncol(l2018_gr_recordDescriptor)


list.files("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs")

estat_e_attr_point_allattr_new_anonym <- fread("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/estat_e_attr_point_allattr_new_anonym.csv")
names(estat_e_attr_point_allattr_new_anonym)
str(estat_e_attr_point_allattr_new_anonym$X)
str(estat_e_attr_point_allattr_new_anonym)

l2018_gr_recordDescriptor_1 <- fread("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/l2018_gr_recordDescriptor.csv",
                                     header = TRUE)
l2018_gr_recordDescriptor_1
names(l2018_gr_recordDescriptor_1)
nrow(l2018_gr_recordDescriptor_1)
str(l2018_gr_recordDescriptor_1[1, ])

View(l2018_gr_recordDescriptor_1)
View(l2018_gr_recordDescriptor_1[(nrow(l2018_gr_recordDescriptor_1)), ])
View(l2018_gr_recordDescriptor_1[(nrow(l2018_gr_recordDescriptor_1) + 1), ])

l2018_gr_recordDescriptor_1_kk <- data.frame("X", "DECIMAL", "", "", "X coordinates of the NUTS0 centroid", "estat_e_attr_point_allattr_new_anonym.csv,estat_s_attr_point_allattr_new_anonym.csv")
names(l2018_gr_recordDescriptor_1_kk) <- names(l2018_gr_recordDescriptor_1)

l2018_gr_recordDescriptor_1 <- rbind(l2018_gr_recordDescriptor_1, l2018_gr_recordDescriptor_1_kk)

l2018_gr_recordDescriptor_1_kk <- data.frame("Y", "DECIMAL", "", "", "Y coordinates of the NUTS0 centroid", "estat_e_attr_point_allattr_new_anonym.csv,estat_s_attr_point_allattr_new_anonym.csv")
names(l2018_gr_recordDescriptor_1_kk) <- names(l2018_gr_recordDescriptor_1)

l2018_gr_recordDescriptor_1 <- rbind(l2018_gr_recordDescriptor_1, l2018_gr_recordDescriptor_1_kk)


write.csv(l2018_gr_recordDescriptor_1,
          "/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/l2018_gr_recordDescriptor.csv",
          row.names = FALSE)


## Anonymising POINT_ID ####
point_id_all <- fread("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/complete_final/anonymised_POINT_ID.csv",
                      header = TRUE)
point_id_all


expertDF_exif <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_EXIF/",
                              "expertDF_exif.csv"))
names(expertDF_exif)

expertDF_exif <- expertDF_exif %>%
  left_join(point_id_all, by = "POINT_ID") %>% 
  select(!POINT_ID) 

sum(is.na(expertDF_exif$POINT_ID_anonym))
write.csv(expertDF_exif,
          paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_EXIF/",
                 "expertDF_exif.csv"),
          row.names = FALSE)
  


surveyorDF_exif <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_EXIF/",
                                "surveyorDF_exif.csv"))
names(surveyorDF_exif)

surveyorDF_exif <- surveyorDF_exif %>%
  mutate_at("POINT_ID", as.character) %>% 
  left_join(point_id_all, by = "POINT_ID") %>% 
  select(!POINT_ID) 

sum(is.na(surveyorDF_exif$POINT_ID_anonym))

write.csv(surveyorDF_exif,
          paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_EXIF/",
                 "surveyorDF_exif.csv"),
          row.names = FALSE)


#
estat_e_attr_point_allattr_new_anonym <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/", 
                                                      "estat_e_attr_point_allattr_new_anonym.csv"))
names(estat_e_attr_point_allattr_new_anonym)

estat_e_attr_point_allattr_new_anonym <- estat_e_attr_point_allattr_new_anonym %>%
  mutate_at("POINT_ID", as.character) %>% 
  left_join(point_id_all, by = "POINT_ID") %>% 
  select(!POINT_ID) 

sum(is.na(estat_e_attr_point_allattr_new_anonym$POINT_ID_anonym))

write.csv(estat_e_attr_point_allattr_new_anonym,
          paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/", 
                 "estat_e_attr_point_allattr_new_anonym.csv"),
          row.names = FALSE)

#
estat_s_attr_point_allattr_new_anonym <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/", 
                                                      "estat_s_attr_point_allattr_new_anonym.csv"))
names(estat_s_attr_point_allattr_new_anonym)

estat_s_attr_point_allattr_new_anonym <- estat_s_attr_point_allattr_new_anonym %>%
  mutate_at("POINT_ID", as.character) %>% 
  left_join(point_id_all, by = "POINT_ID") %>% 
  select(!POINT_ID) 

sum(is.na(estat_s_attr_point_allattr_new_anonym$POINT_ID_anonym))

write.csv(estat_s_attr_point_allattr_new_anonym,
          paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_KEYSPECIES/CSVs/", 
                 "estat_s_attr_point_allattr_new_anonym.csv"),
          row.names = FALSE)

#
releve_data_all <- fread(paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_RELEVE/CSVs/", 
                                "releve_data_all.csv"), 
                         header = TRUE)
head(names(releve_data_all), 10)
tail(names(releve_data_all), 10)

setnames(releve_data_all, old = point_id_all$POINT_ID, new = point_id_all$POINT_ID_anonym, skip_absent = TRUE)
releve_data_all

write.csv(releve_data_all,
          paste0("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/LUCAS_G_2018_RELEVE/CSVs/", 
                 "releve_data_all.csv"),
          row.names = FALSE)



l2018_gr_recordDescriptor_1 <- fread("/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/l2018_gr_recordDescriptor.csv",
                                   header = TRUE)
View(l2018_gr_recordDescriptor_1)
l2018_gr_recordDescriptor_1[1, ]
unique(l2018_gr_recordDescriptor_1$TYPE)

l2018_gr_recordDescriptor_1_kk <- l2018_gr_recordDescriptor_1[1, ]
l2018_gr_recordDescriptor_1_kk
l2018_gr_recordDescriptor_1_kk$TYPE <- "TEXT"
l2018_gr_recordDescriptor_1_kk$VARIABLE <- "POINT_ID_anonym"
l2018_gr_recordDescriptor_1_kk$DESCRIPTION <- "Anonymised unique identifier"
l2018_gr_recordDescriptor_1_kk$DATASET <- "estatdb_allattr_e_point,estatdb_allattr_e_transects,estatdb_allattr_s_point,estatdb_allattr_s_transects"
l2018_gr_recordDescriptor_1_kk

l2018_gr_recordDescriptor_1 <- rbind(l2018_gr_recordDescriptor_1, l2018_gr_recordDescriptor_1_kk)
l2018_gr_recordDescriptor_1 <- l2018_gr_recordDescriptor_1[-1, ]

write.csv(l2018_gr_recordDescriptor_1,
          "/eos/jeodpp/data/projects/REFOCUS/data/LUCAS2018_Grassland/data/anonymised_final/l2018_gr_recordDescriptor.csv",
          row.names = FALSE)
