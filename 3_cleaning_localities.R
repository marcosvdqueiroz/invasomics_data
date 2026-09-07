library(readr)
library(CoordinateCleaner)
library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(countrycode)
library(ggplot2)
library(sf)
library(purrr)
library(rnaturalearth)

#merge all nucleotide summaries into one df (do it just once)

#files = list.files(pattern = "./2_get_info_nucleotide/nucleotide_entrez_summary_clean_.*\\.tsv$",
#                   full.names = TRUE)

#nuc_summary = files %>%
#  lapply(read_tsv, show_col_types = FALSE) %>%
#  bind_rows()

#write_tsv(nuc_summary, './outputs/nucleotide_entrez_summary_clean_merged.tsv')

#From now one, just read the file
nuc_summary = read_tsv('./outputs/2_get_info_nucleotide/nucleotide_entrez_summary_clean_merged.tsv')

#Tier 1 accessions 
  #- Has lat long info
  #- Has gazetter info

tier1 = nuc_summary %>% filter(!is.na(lat_lon), str_detect(country, "^[^:]+:")) %>%
  mutate(tier = 1)
#36,869
write_tsv(tier1, './outputs/3_cleaning_localities/raw_tier1.tsv')

#Tier 2 accessions
  #- Has lat long info
  #- Has only country

tier2 = nuc_summary %>% filter(!is.na(lat_lon), !str_detect(country, "^[^:]+:")) %>%
  mutate(tier = 2)
#9,677
write_tsv(tier2, './outputs/3_cleaning_localities/raw_tier2.tsv')


#Tier 3 accessions
  #-Has lat long info
  #-Has no country info

tier3 = nuc_summary %>% filter(!is.na(lat_lon), is.na(country)) %>%
  mutate(tier = 3)
#1,409
write_tsv(tier3, './outputs/3_cleaning_localities/raw_tier3.tsv')


#Tier 4 accessions
  #-Has gazetter and no lat long info

tier4 = nuc_summary %>% filter(is.na(lat_lon), str_detect(country, "^[^:]+:")) %>%
  mutate(tier = 4)
#71,103
write_tsv(tier4, './outputs/3_cleaning_localities/raw_tier4.tsv')


#Tier 5 accessions
  #-Has only country and no lat long info.

tier5 = nuc_summary %>% filter(is.na(lat_lon), !str_detect(country, "^[^:]+:")) %>%
  mutate(tier = 5)
#85,322
write_tsv(tier5, './outputs/3_cleaning_localities/raw_tier5.tsv')

#To those, let's add the country centroid

country_centroids <- countryref %>%
  select(name, centroid.lon, centroid.lat) %>%
  distinct(name, .keep_all = TRUE)

#Correct some countries (based on crossref df), while excluding historical regions...
tier5_clean <- tier5 %>%
  mutate(
    country = case_when(
      country == "USA" ~ "United States of America",
      country == "Viet Nam" ~ "Vietnam",
      country == "Myanmar" ~ "Myanmar (Burma)",
      country == "Hong Kong" ~ "Hong Kong SAR China",
      country == "Czech Republic" ~ "Czechia",
      country == "Bosnia and Herzegovina" ~ "Bosnia & Herzegovina",
      country == "Sao Tome and Principe" ~ "Sao Tome & Principe",
      country == "North Macedonia" ~ "Macedonia",
      country == "Trinidad and Tobago" ~ "Trinidad & Tobago",
      country == "Democratic Republic of the Congo" ~ "Congo - Kinshasa",
      country == "State of Palestine" ~ "	Palestinian Territories", #FREE PALESTINE!!!
      country == "Republic of the Congo" ~ "Congo - Brazzaville",
      country == "Virgin Islands" ~ "United States of America",
      country == "Micronesia" ~ "Micronesia (Federated States of)",
      country == "Turks and Caicos Islands" ~ "Turks & Caicos Islands",
      country == "Saint Vincent and the Grenadines" ~ "St. Vincent & Grenadines",
      
      TRUE ~ country
    )
  ) %>%
      
      # Remove problematic, non-country, hard to pinpoint and/or historical regions
      filter(!country %in% c(
        "Netherlands Antilles",
        "Kosovo",
        "Yugoslavia",
        "Pacific Ocean",
        "USSR",
        "Indian Ocean",
        "Borneo",
        "Czechoslovakia",
        "Europa Island",
        "not collected",
        "Serbia and Montenegro",
        "Korea"
      )
    )



tier5_clean <- tier5_clean %>%
  left_join(
    country_centroids,
    by = c("country" = "name")
  ) %>%
  mutate(
    lat_dir = if_else(centroid.lat < 0, "S", "N"),
    lon_dir = if_else(centroid.lon < 0, "W", "E"),
    
    lat_lon = if_else(
      !is.na(centroid.lat) & !is.na(centroid.lon),
      sprintf(
        "%.2f %s %.2f %s",
        abs(centroid.lat), lat_dir,
        abs(centroid.lon), lon_dir
      ),
      NA_character_
    )
  ) %>%
  select(-centroid.lat, -centroid.lon, -lat_dir, -lon_dir)

tier5_without_info = tier5_clean %>% 
  filter(is.na(lat_lon))

tier5_clean = tier5_clean %>%
  mutate(
    nums = str_extract_all(lat_lon, "-?\\d+\\.?\\d*"),
    lat = as.numeric(sapply(nums, `[`, 1)),
    lon = as.numeric(sapply(nums, `[`, 2)),
    lat_dir = str_extract(lat_lon, "[NS]"),
    lon_dir = str_extract(lat_lon, "[EW]"),
    decimalLatitude = ifelse(lat_dir == "S", -abs(lat), lat),
    decimalLongitude = ifelse(lon_dir == "W", -abs(lon), lon)
  ) %>%
  select(-nums, -lat, -lon, -lat_dir, -lon_dir) %>%
  filter(!is.na(decimalLatitude))
#74,615
write_tsv(tier5_clean, './outputs/3_cleaning_localities/clean_tier5.tsv')



#Tier 6 accession
  #-Without any georeference... Don't use them. Only calculate to see the percentage of unreferenced accessions

tier6 = nuc_summary %>% filter(is.na(lat_lon), is.na(country)) %>%
  mutate(tier = 6)
#68,608
write_tsv(tier6, './outputs/3_cleaning_localities/raw_tier6.tsv')


#Checking the validity of the lat long values from the original sources, based on the provided localities

tiers_with_latlong = bind_rows(tier1, tier2, tier3)
any(is.na(tiers_with_latlong$lat_lon)) #must be false


#Split the lat long into two columns. Also, standardize them 
tiers_with_latlong <- tiers_with_latlong %>%
  mutate(
    nums = str_extract_all(lat_lon, "-?\\d+\\.?\\d*"),
    lat = as.numeric(sapply(nums, `[`, 1)),
    lon = as.numeric(sapply(nums, `[`, 2)),
    lat_dir = str_extract(lat_lon, "[NS]"),
    lon_dir = str_extract(lat_lon, "[EW]"),
    decimalLatitude = ifelse(lat_dir == "S", -abs(lat), lat),
    decimalLongitude = ifelse(lon_dir == "W", -abs(lon), lon)
  ) %>%
  select(-nums, -lat, -lon, -lat_dir, -lon_dir) %>%
  filter(!is.na(decimalLatitude))

#Let's geolocate the country of tier3 based on the lat lon data

tier3_lat_lon = tiers_with_latlong %>%
  filter(tier == 3) %>%
  select(uid,
         caption,
         decimalLatitude,
         decimalLongitude)

tier3_lat_lon = st_as_sf(
  tier3_lat_lon,
  coords = c('decimalLongitude', 'decimalLatitude'),
  crs = 4326
)
  
world <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

sf::sf_use_s2(FALSE)

points_with_country <- st_join(
  tier3_lat_lon,
  world %>% select(admin)
)

points_with_country <- points_with_country %>%
  st_drop_geometry() %>%
  group_by(uid, caption) %>%
  summarise(
    inferred_country = paste(unique(na.omit(admin)), collapse = "|"),
    .groups = "drop"
  ) %>%
  mutate(
    inferred_country = na_if(inferred_country, "")
  )


#back to the df

tiers_with_latlong <- tiers_with_latlong %>%
  left_join(
    points_with_country,
    by = c("uid", "caption")
  ) %>%
  mutate(
    country = if_else(
      is.na(country),
      inferred_country,
      country
    )
  ) %>%
  select(-inferred_country)
  

#extracting the name of the country

tiers_with_latlong = tiers_with_latlong %>%  
  mutate(
    country_only = str_trim(str_extract(country, "^[^:]+"))
  )

#converting the name of the country to iso codes

tiers_with_latlong$country_code = countrycode(tiers_with_latlong$country_only,
                                 origin = 'country.name',
                                 destination = 'iso3c'
)

#Some accessions came with weird country names. Let's fix it.

#Warning message:
#In countrycode_convert(sourcevar = sourcevar, origin = origin, destination = dest,  :
#  Some values were not matched unambiguously: Kerguelen Archipelago, Kosovo, Micronesia, Virgin Islands

#Kerguelen Archipelago belongs to France
#Micronesia should be Federated States of Micronesia
#Virgin Islands (at least in our datasets) belong to the USA
#Kosovo doesn't have a ISO3-c code

tiers_with_latlong <- tiers_with_latlong %>%
  mutate(
    country = country %>%
      # fix Kerguelen
      str_replace("Kerguelen Archipelago\\s*:", "Kerguelen Archipelago") %>%
      
      # fix Virgin Islands
      str_replace("^Virgin Islands:\\s*", "U.S. Virgin Islands: ") %>%
      
      # fix Micronesia
      str_replace("^Micronesia:?\\s*", "Federated States of Micronesia: ") %>%
      
      # add France prefix when needed
      if_else(
        str_detect(., "Kerguelen Archipelago"),
        paste0("France: ", .),
        .
      )
  )

#Repeat the same code above.

#extracting the name of the country

tiers_with_latlong = tiers_with_latlong %>%  
  mutate(
    country_only = str_trim(str_extract(country, "^[^:]+"))
  )

#converting the name of the country to iso codes

tiers_with_latlong$country_code = countrycode(tiers_with_latlong$country_only,
                                       origin = 'country.name',
                                       destination = 'iso3c',
                                       custom_match = c(Kosovo = 'KSV')
)


#Now, let's check them using coordinatecleaner

library(rnaturalearth)

land = ne_countries(scale = 10, returnclass = 'sf')
land = as(land, 'Spatial')


### CoordinateCleaner ####

flags = clean_coordinates(x = tiers_with_latlong,
                          lon = 'decimalLongitude',
                          lat = 'decimalLatitude',
                          countries = 'country_code',
                          species = 'taxa_accepted',
                          tests = c("capitals", "centroids",
                                    "equal", "zeros", "countries", "institutions", "seas"),
                          seas_ref = land)

summary(flags)
#plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")


#Data cleaning

#1. Flags equal
#Criteria: as longs as there are some gazeteer referencing the locality
flags_equal = flags %>% filter(.equ == FALSE)

#Only one of the flags_equal looks weird - a Canadian accession (uid = 1331390862). The Egyptian accessions looks fine.

flags = flags %>%
  filter(uid != '1331390862')

#2. Flags institutions
#Criteria: if they are cultivated in institutions, they are not in their natural or introduced ranges

flags = flags %>% filter(.inst != FALSE)

#Also, capture those accessions with 'Garden*', 'Collection*' etc in the gazetteer, isolation_source and isolation fields that didn't match for institutions

flags = flags %>% filter(!grepl('Garden*|Collection*|Museum|Arboretum|Arbortum|Botan*|Cultiv*|Nursery*|Experimental|Hort*|Greenhouse', country, ignore.case = TRUE))
flags = flags %>% filter(!grepl('Garden*|Collection*|Museum|Arboretum|Arbortum|Botan*|Cultiv*|Nursery*|Experimental|Hort*|Greenhouse', isolation_source, ignore.case = TRUE))
flags = flags %>% filter(!grepl('Garden*|Collection*|Museum|Arboretum|Arbortum|Botan*|Cultiv*|Nursery*|Experimental|Hort*|Greenhouse', isolate, ignore.case = TRUE))



#3. Flags capitals - Since we are working with invasive species, there's no issue if it's occurring in a capital city, as long as they are not in cultivation

#4. Flags centroids -  Let's keep those, since tier5 (and part of tier4) will also contain country centroids as well.

#5. Flags countries - These can be removed

flags = flags %>% filter(.con != FALSE)

#6. Flags sea - Let's remove these occurrences. Unfortunately, sea-plants might be removed, though...

flags = flags %>% filter(.sea != FALSE)

summary(flags)

#These remained flagged (tier1-3) are suitable for our next steps.

tier1_clean = flags %>% filter(tier == 1) #30,388
tier2_clean = flags %>% filter(tier == 2) #8,764
tier3_clean = flags %>% filter(tier == 3) #1,255


write_tsv(tier1_clean, './outputs/3_cleaning_localities/clean_tier1.tsv')
write_tsv(tier2_clean, './outputs/3_cleaning_localities/clean_tier2.tsv')
write_tsv(tier3_clean, './outputs/3_cleaning_localities/clean_tier3.tsv')

