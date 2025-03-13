# Loading libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(leaflet)
library(writexl)
library(mapview)

# Place the files in the same folde as the script

# Import data
Occurences_MOTUs <- readxl::read_xlsx(path = "~/Nextcloud/Recherche/2_Plankton/Forams/1-Data_raw/Data_Morard/Biogeography/Table S7.xlsx", sheet = 1, col_types="guess")
Env_data         <- readxl::read_xlsx(path = "~/Nextcloud/Recherche/2_Plankton/Forams/1-Data_raw/Data_Morard/Biogeography/Table S4.xlsx", sheet = 1, col_types="guess")

# Recode Ocean Basin to have only 7 categories

# First fix the biogeography

Biogeography_data <- Env_data %>% select(`Original Sample ID`, Averaged_Latitude, Averaged_Longitude, Province_Spalding_simplified_01) %>% 
  dplyr::rename("Biogeography" = "Province_Spalding_simplified_01", "latitude" = "Averaged_Latitude", "longitude" = "Averaged_Longitude") 

# Create a color palette
pal <- colorFactor(palette = c("#1B9E77", "#7B2E77", "#D95F02", "#7570B3", "#E7298A", "#E6CA00", "#A6AA50", "#1151B1", "#A6A1A1"), domain = Biogeography_data$Biogeography)

# Create a map to check
m <- leaflet(Biogeography_data) %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addCircleMarkers(~longitude, ~latitude, color = ~pal(Biogeography), radius = 1) %>%
  addLegend("bottomright", pal = pal, values = ~Biogeography,
            title = "Biogeography")
m
mapshot(m, file = "~/Nextcloud/Recherche/2_Plankton/Forams/3-Data_processed/Biogeography/Map_biogeography_classification_9regions.pdf")

# Now create an occurrence of MOTUs lvl-2 by region
Occurences_MOTUs_biogeography <- Biogeography_data %>% 
  left_join(Occurences_MOTUs %>% 
              select(`Original Sample ID`, TP_lvl2) %>% 
              filter(TP_lvl2 != "Not_Attributed"), relationship = "many-to-many") %>% unique() %>% drop_na() %>% 
  group_by(Biogeography, TP_lvl2) %>% 
  count() %>% 
  spread(Biogeography,n,fill=0)
              
#Occurences_MOTUs_biogeography[grepl("Neogloboquadrina\\|pachyderma\\|V\\|a", Occurences_MOTUs_biogeography$TP_lvl2),]
# Based on the raw count calculate % of observation
Occurences_MOTUs_biogeography_prop <- Occurences_MOTUs_biogeography %>% 
  column_to_rownames('TP_lvl2') %>% 
  as.matrix() %>% 
  prop.table(1) %>% 
  as.data.frame() %>% 
  rownames_to_column('TP_lvl2') %>% 
  gather(Biogeography, prop_obs, 2:10) 

# Visualize proportion of occurrence by bins of 1%, and add a line to represent the selected threshold
# We start with a threshold at 5%
treshold <- 0.05
ggplot(Occurences_MOTUs_biogeography_prop %>%
         filter(prop_obs>0), 
       aes(x=prop_obs))+
  geom_histogram(binwidth = 0.01, 
                 fill="black", 
                 col="grey")+
  geom_vline(aes(xintercept = treshold), color = "red", linewidth = 1, linetype = "dashed")

#Filter out observation below the treshold and reconstruct the occurence table in presence/absence data
Data_biogeo_filtered <- Occurences_MOTUs_biogeography_prop %>% 
  filter(prop_obs >= treshold) %>% 
  select(-prop_obs) %>% 
  mutate(pres=1) %>% 
  spread(Biogeography,pres,fill=0)

#Export the data
write_xlsx(
  Data_biogeo_filtered %>% as.data.frame(),
  path ="~/Nextcloud/Recherche/2_Plankton/Forams/3-Data_processed/Biogeography/Occurence_Table_9regions.xlsx")

