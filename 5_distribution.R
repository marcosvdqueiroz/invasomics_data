######################
#Capturing the distribution
######################

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


#This script follows the previous ones (3_cleaning_localities.R and 4_georeferencing.R especially)

#First, let's combine the outputs from the 4th script only once. After that, read the merged file

#output_files = Sys.glob('./4_geolocate_metacentrum/output_from_geolocate_*.tsv')
#output_merged = lapply(output_files %>% lapply(read_tsv, show_col_types = F) %>% bind_rows())
#write_tsv(output_merged)

tier4_semiclean = read_tsv('./outputs/4_geolocate_metacentrum/output_from_geolocate_merged.tsv')

################ Filtering the geoloc results ###########################

#Per uid

#1. Remove rows with missing lat/lon, i.e., unable to georeference
#2. Prioritize glcPrecision == 'High'
#3. Within the best precision class, keep highest glcScore
#4. If still tied, keep lowest glcRank
#5. If an uid has only Low precision, add the country centroid.
#6. Mutate the tier column: If geolocalized, tier 4g; If centroid used, tier 4c


country_centroids <- countryref %>%
  select(name, centroid.lon, centroid.lat) %>%
  distinct(name, .keep_all = TRUE)


tier4_clean <- tier4_semiclean %>%
  mutate(
    glcLatitude = as.numeric(glcLatitude),
    glcLongitude = as.numeric(glcLongitude),
    glcScore = as.numeric(glcScore),
    glcRank = as.numeric(glcRank)
  ) %>%
  
  # remove rows without coordinates
  filter(
    !is.na(glcLatitude),
    !is.na(glcLongitude)
  ) %>%
  
  # rank precision
  mutate(
    precision_rank = case_when(
      glcPrecision == "High" ~ 1,
      glcPrecision == "Medium" ~ 2,
      glcPrecision == "Low" ~ 3,
      TRUE ~ 4
    )
  ) %>%
  
  # keep best geolocate result per uid
  group_by(uid) %>%
  arrange(
    precision_rank,
    desc(glcScore),
    glcRank,
    .by_group = TRUE
  ) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  
  # add country centroid
  left_join(
    country_centroids,
    by = c("country_only" = "name")
  ) %>%
  
  # if only Low was available, replace with centroid
  mutate(
    only_low_precision = glcPrecision == "Low",
    
    glcLatitude = if_else(
      only_low_precision & !is.na(centroid.lat),
      centroid.lat,
      glcLatitude
    ),
    
    glcLongitude = if_else(
      only_low_precision & !is.na(centroid.lon),
      centroid.lon,
      glcLongitude
    ),
    
    tier = if_else(
      only_low_precision,
      "4c",
      "4g"
    )
  ) %>%
  select(
    -precision_rank,
    -centroid.lat,
    -centroid.lon,
    -only_low_precision
  )




colnames(tier4_clean)[which(names(tier4_clean) == 'glcLatitude')] = 'decimalLatitude'
colnames(tier4_clean)[which(names(tier4_clean) == 'glcLongitude')] = 'decimalLongitude'

write_tsv(tier4_clean, './outputs/3_cleaning_localities/clean_tier4.tsv')
#56,168

#Geolocalized
nrow(tier4_clean %>% filter(tier == '4g'))
#29,420

#Centroids
nrow(tier4_clean %>% filter(tier == '4c'))
#26,748


###Fetching distribution###

#First, merge our tiers 1, 2 and 3 ("flags", from the 3th script), with the filtered georeferenced tiers 4c and 4g + the centroids from tier5

tier1_clean = read_tsv('./outputs/3_cleaning_localities/clean_tier1.tsv')
tier2_clean = read_tsv('./outputs/3_cleaning_localities/clean_tier2.tsv')
tier3_clean = read_tsv('./outputs/3_cleaning_localities/clean_tier3.tsv')
tier5_clean = read_tsv('./outputs/3_cleaning_localities/clean_tier5.tsv')
tier5_clean$country_only = tier5_clean$country

cols_to_keep <- c(
  "taxa_accepted", "uid", "caption", "title", "genome", "gene",
  "CDS_product", "ref_journal", "ref_authors", "isolate",
  "country", "isolation_source", "altitude", "tier",
  "decimalLatitude", "decimalLongitude", "country_only", "locality")

combined_tiers <- bind_rows(
  tier1_clean %>%
    mutate(tier = as.character(tier)),
  
  tier2_clean %>%
    mutate(tier = as.character(tier)),
  
  tier3_clean %>%
    select(any_of(cols_to_keep)) %>%
    mutate(tier = as.character(tier)),
  
  tier4_clean %>%
    select(any_of(cols_to_keep)) %>%
    mutate(tier = as.character(tier)),
  
  tier5_clean %>%
    select(any_of(cols_to_keep)) %>%
    mutate(tier = as.character(tier))
)
  

write_tsv(combined_tiers, './outputs/5_distribution/combined_tiers.tsv')


#1. POWO
#The powo dataset was download from https://sftp.kew.org/pub/data-repositories/WCVP/ on January 2026

powo_names = read_delim('./powo_df/wcvp_names.csv', delim = '|')
powo_names = powo_names %>% filter(taxon_status == 'Accepted')
powo_dist = read_delim('./powo_df/wcvp_distribution.csv', delim = '|')

#Combine the taxonomic info from powo with our tiers taxa, making sure that we only have unique names. 


taxa_combined_tiers_full_tx_info <- combined_tiers %>%
  left_join(
    powo_names,
    by = c("powo_id_accepted" = "powo_id")
  ) %>%
  distinct(taxa_accepted, .keep_all = TRUE) %>%
  select(
    taxa_accepted,
    family,
    powo_id_accepted,
    plant_name_id
  ) %>%
  rename(powo_id = powo_id_accepted)

no_ids = taxa_combined_tiers_full_tx_info %>% filter(is.na(powo_id))

#Some taxa didn't show any powo id, because they are varieties, subspecies, etc... For most of them we don't know the true taxonomic position, so let's remove them.

taxa_combined_tiers_full_tx_info = taxa_combined_tiers_full_tx_info %>%
  filter(!taxa_accepted %in% no_ids$taxa_accepted)


# -----------------------------
#  Fetching POWO distribution
# -----------------------------




# =========================================================
# 1. Connect species table to POWO names table
# =========================================================
# Here we link:
# taxa_combined_tiers_full_tx_info$powo_id
#            ->
# powo_names$powo_id
#
# This gives us the corresponding plant_name_id,
# which is needed to access the distribution table.

powo_distribution_summary <- taxa_combined_tiers_full_tx_info %>%
  
  
  # =========================================================
# 2. Connect with POWO distribution table
# =========================================================
# plant_name_id is the key linking POWO names
# to geographic distribution information.

left_join(
  powo_dist,
  by = "plant_name_id"
) %>%
  
  
  # =========================================================
# 3. Group records by species
# =========================================================
# One species may have multiple distribution records,
# so we collapse them into a single row per species.

group_by(
  taxa_accepted,
  family,
  powo_id,
  plant_name_id
) %>%
  
  
  # =========================================================
# 4. Summarize native and introduced distributions
# =========================================================
#
# introduced == 0  -> native distribution
# introduced == 1  -> introduced distribution
#
# We:
# - remove NAs
# - keep unique values
# - collapse multiple areas using "|"

summarise(
  
  # -----------------------------
  # Native distribution
  # -----------------------------
  
  native_continent = paste(
    unique(na.omit(continent[introduced == 0])),
    collapse = "|"
  ),
  
  native_region = paste(
    unique(na.omit(region[introduced == 0])),
    collapse = "|"
  ),
  
  native_area = paste(
    unique(na.omit(area[introduced == 0])),
    collapse = "|"
  ),
  
  native_area_code = paste(
    unique(na.omit(area_code_l3[introduced == 0])),
    collapse = "|"
  ),
  
  
  
  # -----------------------------
  # Introduced distribution
  # -----------------------------
  
  introduced_continent = paste(
    unique(na.omit(continent[introduced == 1])),
    collapse = "|"
  ),
  
  introduced_region = paste(
    unique(na.omit(region[introduced == 1])),
    collapse = "|"
  ),
  
  introduced_area = paste(
    unique(na.omit(area[introduced == 1])),
    collapse = "|"
  ),
  
  introduced_area_code = paste(
    unique(na.omit(area_code_l3[introduced == 1])),
    collapse = "|"
  ),
  
  .groups = "drop"
) %>%
  
  
  # =========================================================
# 5. Replace empty strings with NA
# =========================================================
# Some species may have only native OR introduced
# distributions. Empty strings are converted to NA.

mutate(
  across(
    c(
      native_continent,
      native_region,
      native_area,
      native_area_code,
      introduced_continent,
      introduced_region,
      introduced_area,
      introduced_area_code
    ),
    ~ na_if(.x, "")
  )
)

write_tsv(powo_distribution_summary, './outputs/5_distribution/powo_distribution_summary.tsv')



#2. Glonaf Distribution

# Input tables
# taxa_combined_tiers_full_tx_info has:
# taxa_accepted, family_wcvp, powo_id_accepted

gltx <- read_csv("./glonaf_data/glonaf_taxon_wcvp.csv")
glflora     <- read_csv("./glonaf_data/glonaf_flora2.csv")
gldist     <- read_csv("./glonaf_data/glonaf_region.csv")

# Get GloNAF distribution for all taxa in your dataframe

taxa_glonaf_distribution <- taxa_combined_tiers_full_tx_info %>%
  distinct(taxa_accepted, family, powo_id, plant_name_id) %>%
  
  left_join(
    gltx %>%
      select(
        taxon_wcvp_id = id,
        powo_id = powo_id_accepted
      ),
    by = "powo_id"
  ) %>%
  
  left_join(
    glflora %>%
      select(
        glflora_id = id,
        taxon_wcvp_id,
        list_id,
        status,
        region_id
      ),
    by = "taxon_wcvp_id"
  ) %>%
  
  left_join(
    gldist %>%
      select(
        region_id = id,
        region_code = code,
        region_name = name,
        tdwg4_id,
        OBJIDsic
      ),
    by = "region_id"
  ) %>%
  
  select(
    taxa_accepted,
    family,
    powo_id,
    taxon_wcvp_id,
    region_id,
    region_name,
    region_code,
    tdwg4_id,
    OBJIDsic,
    status
  ) %>%
  distinct() %>%
  filter(!is.na(region_id))


# =========================================================
# Identify TDWG codes that are numeric placeholders
# =========================================================
# Example:
# 000-123
# 000-456
#
# These should be replaced by the broadest tdwg possible.

numeric_tdwgs <- grep(
  '^000-*',
  unique(taxa_glonaf_distribution$tdwg4_id),
  value = TRUE
)

taxa_numeric_distribution = taxa_glonaf_distribution %>%
  filter(tdwg4_id %in% numeric_tdwgs)

unique_numeric = unique(taxa_numeric_distribution$tdwg4_id)

tdwg4 = read_delim('./tdwg_df/tblLevel4.txt', delim = '*', locale = readr::locale(encoding = "UTF-8"))

tdwg4_to_tdwg3 <- tdwg4 %>%
  select(
    tdwg4_id = `L4 code`,
    tdwg3_equivalent = `L3 code`,
    iso_code = `L4 ISOcode`
  ) %>%
  distinct()



# ============================================================
# STEP 1. Create manual corrections for problematic regions
# ============================================================
#
# Some GloNAF regions are broader than a single TDWG region,
# absent from TDWG, or correspond to multiple TDWG regions.
#
# The strategy used is:
#
# 1. Use TDWG4 when a direct match exists
# 2. Use TDWG3 when a direct match exists
# 3. Expand large countries into all constituent TDWG3 regions
# 4. Expand archipelagos into all constituent TDWG regions
# 5. Use ISO codes only as a last resort
#
# ============================================================

manual_tdwg_fixes <- tribble(
  ~region_name, ~tdwg_fixed, ~assignment_type,
  
  # ----------------------------------------------------------
  # Large countries split into multiple TDWG3 regions
  # ----------------------------------------------------------
  
  "Brazil", "BZC|BZE|BZL|BZN|BZS", "expanded_country",
  "Australia", "NSW|NTA|QLD|SOA|TAS|VIC|WAU", "expanded_country",
  "Mexico", "MXC|MXE|MXG|MXI|MXN|MXS|MXT", "expanded_country",
  "Argentina", "AGE|AGS|AGW", "expanded_country",
  "Santa Fe", "AGE", "exact_tdwg4",
  "Córdoba", "AGE", "exact_tdwg4",
  "Chile", "CLC|CLN|CLS", "expanded_country",
  
  "China",
  "CHC|CHH|CHI|CHM|CHN|CHQ|CHS|CHT|CHX",
  "expanded_country",
  
  "Canada",
  "ABT|BRC|LAB|MAN|NBR|NSC|NUN|NWT|ONT|PEI|QUE|SAS|YUK|NFL",
  "expanded_country",
  
  "South Africa",
  "CPP|NAT|OFS|TVL",
  "expanded_country",
  
  "New Zealand",
  "NZN|NZS",
  "expanded_country",
  
  "Indonesia",
  "JAW|MOL|SUL|SUM",
  "expanded_country",
  
  "United Kingdom of Great Britain and Northern Ireland (the)",
  "GRB|IRE",
  "expanded_country",
  
  # ----------------------------------------------------------
  # Exact TDWG3 matches
  # ----------------------------------------------------------
  
  "India", "IND", "exact_tdwg3",
  "Japan", "JAP", "exact_tdwg3",
  "Turkey", "TUR", "exact_tdwg3",
  "Yemen", "YEM", "exact_tdwg3",
  "Italy", "ITA", "exact_tdwg3",
  "Portugal", "POR", "exact_tdwg3",
  "Cook Islands", "COO", "exact_tdwg3",
  "Borneo", "BOR", "exact_tdwg3",
  "New Guinea", "NWG", "exact_tdwg3",
  "Korea Peninsula", "KOR", "exact_tdwg3",
  "Netherlands Antilles", "NLA", "exact_tdwg3",
  "Malaysia", "MLY", "exact_tdwg3",
  
  # ----------------------------------------------------------
  # Exact TDWG4 matches
  # ----------------------------------------------------------
  
  "Puducherry", "IND-PO", "exact_tdwg4",
  "Azerbaijan", "TCS-AZ", "exact_tdwg4",
  "Papua New Guinea", "NWG-PN", "exact_tdwg4",
  "Saint-Martin (French part)", "LEE-SM", "exact_tdwg4",
  "St. Martin (island)", "LEE-SM", "exact_tdwg4",
  
  # ----------------------------------------------------------
  # Countries + islands
  # ----------------------------------------------------------
  
  "France (incl. Corsica)",
  "FRA|COR",
  "expanded_country",
  
  "Spain (incl. islands)",
  "SPA|BAL|CNY",
  "expanded_country",
  
  "Ecuador (incl. Galapagos)",
  "ECU|GAL",
  "expanded_country",
  
  "China incl. Taiwan and Hong Kong",
  "CHC|CHH|CHI|CHM|CHN|CHQ|CHS|CHT|CHX|TAI",
  "expanded_country",
  
  # ----------------------------------------------------------
  # Archipelagos
  # ----------------------------------------------------------
  
  "Scattered Islands",
  "MDG",
  "expanded_archipelago",
  
  "French Polynesia",
  "SCI|MRQ|TUA|TUB",
  "expanded_archipelago",
  
  "Kiribati",
  "GIL|LIN|PHX",
  "expanded_archipelago",
  
  "Sao Tome and Principe",
  "GGI-ST|GGI-PR",
  "expanded_archipelago",
  
  "Comoros Archipelago incl. Mayotte",
  "COM-CO|COM-MA",
  "expanded_archipelago",
  
  "Virgin Islands (including British and US Virgin Islands)",
  "LEE-BV|LEE-VI",
  "expanded_archipelago",
  
  "French Subantarctic Islands (Crozet, Kerguelen, St. Paul and Amsterdam",
  "CRZ|KEG|ASP",
  "expanded_archipelago",
  
  # ----------------------------------------------------------
  # Subregions
  # ----------------------------------------------------------
  
  "Baja California Peninsula",
  "MXN-BC|MXN-BS",
  "expanded_region",
  
  "Andaman and Nicobar",
  "AND|NCB",
  "expanded_region",
  
  # ----------------------------------------------------------
  # Russian regions
  # ----------------------------------------------------------
  
  "Ussuri, Russia",
  "PRM",
  "approximated_tdwg3",
  
  "Okhota, Russia",
  "MAG",
  "approximated_tdwg3",
  
  "Lower Zeya, Russia",
  "AMU",
  "approximated_tdwg3",
  
  "Bureya, Russia",
  "AMU",
  "approximated_tdwg3",
  
  "Yamal-Gydan, Russia",
  "WSB",
  "approximated_tdwg3",
  
  "Koryakia, Russia",
  "KAM",
  "approximated_tdwg3",
  
  "Anadyr-Penzhina, Russia",
  "MAG",
  "approximated_tdwg3",
  
  "Russian Far East, temperate part, Russia",
  "AMU|KHA|PRM|SAK|KAM|MAG|KUR",
  "expanded_region",
  
  "Russian Federation (European part)",
  "RUW|RUC|RUS|RUN",
  "approximated_tdwg3",
  
  "Anyui, Russia",
  "KHA|AMU",
  "approximated_tdwg3",
  
  "Krasnoyarsk plus Yamal-Nenetsk, Russia",
  "KRA|WSB",
  "approximated_tdwg3",
  
  "Russian Federation (the)",
  "ALT|AMU|BLR|BLT|BRY|CTA|IRK|KAM|KHA|KRA|KUR|MAG|NCS|PRM|RUC|RUE|RUN|RUS|RUW|SAK|TVA|WSB|YAK",
  "expanded_country",
  
  # ----------------------------------------------------------
  # Broad North America
  # ----------------------------------------------------------
  
  "North Alaska - Yukon Territory",
  "ASK|YUK",
  "expanded_region",
  
  "Hudson Bay - Labrador, Canada",
  "LAB|NUN|QUE",
  "approximated_tdwg3",
  
  "Central Canada",
  "MAN|SAS|ONT",
  "approximated_tdwg3",
  
  # ----------------------------------------------------------
  # USA
  # ----------------------------------------------------------
  
  "United States of America, contiguous (without Alaska and Hawaii)",
  paste(
    c(
      "ALA","ARI","ARK","CAL","CNT","COL","DEL","FLA",
      "GEO","IDA","ILL","INI","IOW","KAN","KTY","LOU",
      "MAI","MRY","MAS","MIC","MIN","MSI","MSO","MNT",
      "NEB","NEV","NWH","NWJ","NWM","NWY","NCA","NDA",
      "OHI","OKL","ORE","PEN","RHO","SCA","SDA","TEN",
      "TEX","UTA","VER","VRG","WAS","WDC","WIS","WVA",
      "WYO"
    ),
    collapse = "|"
  ),
  "expanded_country",
  
  "United States of America (the)",
  paste(
    c(
      "ALA","ARI","ARK","CAL","CNT","COL","DEL","FLA",
      "GEO","IDA","ILL","INI","IOW","KAN","KTY","LOU",
      "MAI","MRY","MAS","MIC","MIN","MSI","MSO","MNT",
      "NEB","NEV","NWH","NWJ","NWM","NWY","NCA","NDA",
      "OHI","OKL","ORE","PEN","RHO","SCA","SDA","TEN",
      "TEX","UTA","VER","VRG","WAS","WDC","WIS","WVA",
      "WYO"
    ),
    collapse = "|"
  ),
  "expanded_country"
)

# ============================================================
# STEP 2. Apply corrections
# ============================================================

taxa_glonaf_distribution_fixed <- taxa_glonaf_distribution %>%
  left_join(
    manual_tdwg_fixes,
    by = "region_name"
  ) %>%
  mutate(
    tdwg_final = case_when(
      !is.na(tdwg_fixed) ~ tdwg_fixed,
      TRUE ~ tdwg4_id
    ),
    
    assignment_type = case_when(
      !is.na(assignment_type) ~ assignment_type,
      tdwg4_id == "000-07" ~ "still_unresolved",
      TRUE ~ "original"
    )
  )

# ============================================================
# STEP 3. Expand rows containing multiple TDWG codes
# ============================================================

taxa_glonaf_distribution_expanded <- taxa_glonaf_distribution_fixed %>%
  separate_rows(
    tdwg_final,
    sep = "\\|"
  ) %>%
  mutate(
    tdwg_final = str_trim(tdwg_final)
  )

# ============================================================
# STEP 4. Build TDWG4 -> TDWG3 lookup table
# ============================================================

tdwg4_lookup <- tdwg4 %>%
  select(
    tdwg4_lookup_id = `L4 code`,
    tdwg3_lookup_id = `L3 code`
  ) %>%
  distinct()

# ============================================================
# STEP 5. Create a universal TDWG3 column
# ============================================================
#
# tdwg3_final will contain:
#
# TDWG4 -> converted to TDWG3
# TDWG3 -> unchanged
# ISO -> unchanged
#
# ============================================================

taxa_glonaf_distribution_expanded <- taxa_glonaf_distribution_expanded %>%
  
  left_join(
    tdwg4_lookup,
    by = c("tdwg_final" = "tdwg4_lookup_id")
  ) %>%
  
  mutate(
    
    tdwg3_final = case_when(
      
      # Original TDWG4
      !is.na(tdwg3_lookup_id) ~ tdwg3_lookup_id,
      
      # Already TDWG3
      str_detect(tdwg_final, "^[A-Z]{3}$") ~ tdwg_final,
      
      # ISO fallback
      TRUE ~ tdwg_final
    ),
    
    tdwg_level_final = case_when(
      str_detect(tdwg_final, "^[A-Z]{3}-") ~ "TDWG4",
      str_detect(tdwg_final, "^[A-Z]{3}$") ~ "TDWG3",
      str_detect(tdwg_final, "^[A-Z]{2}$") ~ "ISO",
      TRUE ~ "OTHER"
    )
  ) %>%
  
  select(-tdwg3_lookup_id)

# ============================================================
# STEP 6. Diagnostics
# ============================================================

# Records still carrying numeric tdwg4 codes
still_unresolved <- taxa_glonaf_distribution_expanded %>%
  filter(tdwg3_final %in% unique_numeric)

# ISO fallback records
iso_fallbacks <- taxa_glonaf_distribution_expanded %>%
  filter(tdwg_level_final == "ISO")

# Count records per TDWG3
tdwg3_summary <- taxa_glonaf_distribution_expanded %>%
  dplyr::count(tdwg3_final, sort = TRUE)

# View diagnostics
still_unresolved
iso_fallbacks
tdwg3_summary







taxa_glonaf_distribution_grouped <- taxa_glonaf_distribution_expanded %>%
  group_by(powo_id) %>%
  summarise(
    
    taxa_accepted = paste(
      sort(unique(na.omit(taxa_accepted))),
      collapse = "|"
    ),
    
    family = paste(
      sort(unique(na.omit(family))),
      collapse = "|"
    ),
    
    region_id = paste(
      sort(unique(na.omit(region_id))),
      collapse = "|"
    ),
    
    region_name = paste(
      sort(unique(na.omit(region_name))),
      collapse = "|"
    ),
    
    region_code = paste(
      sort(unique(na.omit(region_code))),
      collapse = "|"
    ),
    
    tdwg3_final = paste(
      sort(unique(na.omit(tdwg3_final))),
      collapse = "|"
    ),
    
    .groups = "drop"
  )












# Save full long-format distribution table
write_tsv(
  taxa_glonaf_distribution_grouped,
  "./outputs/5_distribution/taxa_glonaf_distribution.csv"
)


#Now, let's compare and merge both distribution data

# =========================================================
# 1. Helper function to split distribution codes
# =========================================================
# This function:
# - handles NA values
# - splits codes separated by "|"
# - removes empty strings
# - removes duplicated codes

split_codes <- function(x) {
  x %>%
    tidyr::replace_na("") %>%
    stringr::str_split("\\|") %>%
    purrr::map(~ unique(.x[.x != ""]))
}


# =========================================================
# 2. Prepare POWO distribution data
# =========================================================
# POWO already contains native and introduced distributions
# at TDWG level 3.

powo_compare <- powo_distribution_summary %>%
  select(
    powo_id,
    taxa_accepted,
    family,
    powo_native_l3 = native_area_code,
    powo_introduced_l3 = introduced_area_code
  ) %>%
  mutate(
    powo_native_l3 = split_codes(powo_native_l3),
    powo_introduced_l3 = split_codes(powo_introduced_l3)
  )


# =========================================================
# 3. Prepare GloNAF distribution data
# =========================================================
# GloNAF distributions are treated as naturalized/invasive
# tdwg3_equivalent should already contain TDWG level 3 codes.

glonaf_compare <- taxa_glonaf_distribution_grouped %>%
  select(
    powo_id,
    glonaf_introduced_l3 = tdwg3_final
  ) %>%
  mutate(
    glonaf_introduced_l3 = split_codes(glonaf_introduced_l3)
  )


# =========================================================
# 4. Join POWO and GloNAF distributions by species
# =========================================================

distribution_comparison <- powo_compare %>%
  left_join(
    glonaf_compare,
    by = "powo_id"
  ) %>%
  
  
  # =========================================================
# 5. Compare GloNAF introduced distribution
#    only against POWO introduced distribution
# =========================================================

mutate(
  
  # Regions where GloNAF overlaps with POWO introduced range
  overlap_glonaf_with_powo_introduced = purrr::map2(
    glonaf_introduced_l3,
    powo_introduced_l3,
    intersect
  ),
  
  # GloNAF introduced regions not reported as introduced in POWO
  glonaf_not_in_powo_introduced = purrr::map2(
    glonaf_introduced_l3,
    powo_introduced_l3,
    setdiff
  ),
  
  # POWO introduced regions missing from GloNAF
  powo_introduced_not_in_glonaf = purrr::map2(
    powo_introduced_l3,
    glonaf_introduced_l3,
    setdiff
  )
) %>%
  
  # =========================================================
# 6. Create combined native and introduced distributions
# =========================================================
# Native distribution:
# - only POWO native records
#
# Introduced distribution:
# - POWO introduced + GloNAF introduced
# - duplicated TDWG3 codes removed

mutate(
  
  combined_native_distribution = powo_native_l3,
  
  combined_introduced_distribution = purrr::map2(
    powo_introduced_l3,
    glonaf_introduced_l3,
    ~ unique(c(.x, .y))
  )
) %>%

  # =========================================================
# 7. Collapse list columns back into readable strings
# =========================================================

mutate(
  across(
    c(
      powo_native_l3,
      powo_introduced_l3,
      glonaf_introduced_l3,
      overlap_glonaf_with_powo_introduced,
      glonaf_not_in_powo_introduced,
      powo_introduced_not_in_glonaf,
      combined_native_distribution,
      combined_introduced_distribution
    ),
    ~ purrr::map_chr(.x, ~ paste(.x, collapse = "|"))
  )
) %>%
  
  
  # =========================================================
# 8. Replace empty strings with NA
# =========================================================

mutate(
  across(
    c(
      powo_native_l3,
      powo_introduced_l3,
      glonaf_introduced_l3,
      overlap_glonaf_with_powo_introduced,
      glonaf_not_in_powo_introduced,
      powo_introduced_not_in_glonaf,
      combined_native_distribution,
      combined_introduced_distribution
    ),
    ~ na_if(.x, "")
  )
)

# =========================================================
# 9. Remove rows lacking native and/or introduced data
# =========================================================
# Keep only species that have BOTH:
# - native distribution
# - introduced distribution

distribution_comparison <- distribution_comparison %>%
  filter(
    !is.na(combined_native_distribution),
    !is.na(combined_introduced_distribution)
  )

# =========================================================
# Check overlap between native and introduced distributions
# =========================================================
# This identifies species where at least one TDWG3 code
# appears in BOTH native and introduced distributions.

distribution_comparison <- distribution_comparison %>%
  
  mutate(
    
    # Split strings back into vectors
    native_codes = str_split(
      combined_native_distribution,
      "\\|"
    ),
    
    introduced_codes = str_split(
      combined_introduced_distribution,
      "\\|"
    ),
    
    
    # Find overlapping TDWG3 codes
    native_introduced_overlap = purrr::map2(
      native_codes,
      introduced_codes,
      intersect
    ),
    
    
    # Count overlaps
    n_native_introduced_overlap =
      purrr::map_int(native_introduced_overlap, length),
    
    
    # TRUE/FALSE if overlap exists
    has_native_introduced_overlap =
      n_native_introduced_overlap > 0
  ) %>%
  
  
  # Collapse overlap list-column into readable string
  mutate(
    native_introduced_overlap =
      purrr::map_chr(
        native_introduced_overlap,
        ~ paste(.x, collapse = "|")
      ),
    
    native_introduced_overlap =
      na_if(native_introduced_overlap, "")
  )

# =========================================================
# Remove introduced regions from native distribution
# =========================================================
# If a TDWG3 code appears in BOTH:
# - combined_native_distribution
# - combined_introduced_distribution
# Then remove the introduced distribution from the native areas

distribution_comparison <- distribution_comparison  %>%
  
  
  # Remove overlapping introduced codes from native range
  mutate(
    corrected_native_distribution = purrr::map2(
      native_codes,
      introduced_codes,
      setdiff
    )
  ) %>%
  
  
  # Collapse back into strings
  mutate(
    corrected_native_distribution = purrr::map_chr(
      corrected_native_distribution,
      ~ paste(unique(.x), collapse = "|")
    ),
    
    corrected_native_distribution =
      na_if(corrected_native_distribution, "")
  ) %>%
  
  # If correction produced empty string but original native
  # existed, keep original value
  mutate(
    corrected_native_distribution = if_else(
      corrected_native_distribution == "" &
        !is.na(combined_native_distribution),
      
      combined_native_distribution,
      
      corrected_native_distribution
    )
  ) %>%
  
  
  # Convert empty strings to NA
  mutate(
    corrected_native_distribution =
      na_if(corrected_native_distribution, "")
  ) %>%
  
  
  # Remove temporary columns
  select(-powo_native_l3,
         -powo_introduced_l3,
         -glonaf_introduced_l3,
         -overlap_glonaf_with_powo_introduced,
         -glonaf_not_in_powo_introduced,
         -powo_introduced_not_in_glonaf,
         -combined_native_distribution,
         -native_introduced_overlap,
         -has_native_introduced_overlap,
         -native_codes,
         -introduced_codes,
         -n_native_introduced_overlap
  )

#Filter those where final combined returns NAs

distribution_comparison = distribution_comparison %>% filter(!is.na(corrected_native_distribution))

write_tsv(distribution_comparison, './outputs/5_distribution/distribution_comparision.tsv')



#####################################################################################################################
#Categorize our georeferenced seqs to its native/introduced areas based on their TDWG location
#####################################################################################################################


library(readr)
library(dplyr)
library(sf)


combined_tiers <- combined_tiers %>%
  mutate(row_id_original = row_number())



# Remove rows without coordinates, but keep them to bind back later
chunk_no_coords <- combined_tiers %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude)) %>%
  mutate(
    ISO_Code = NA_character_,
    Level_4_Na = NA_character_,
    Level4_cod = NA_character_,
    Level4_2 = NA_character_,
    Level3_cod = NA_character_,
    Level2_cod = NA_character_,
    Level1_cod = NA_character_
  )



chunk_coords <- combined_tiers %>%
  filter(!is.na(decimalLatitude),
         !is.na(decimalLongitude))


if (nrow(chunk_coords) > 0) {
  
  sf_use_s2(FALSE)
  
  combined_sf <- chunk_coords %>%
    st_as_sf(
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326,
      remove = FALSE
    )
  
  tdwg <- st_read("./tdwg_df/WGSRPD/level4/level4.shp", quiet = TRUE) %>%
    st_transform(4326) %>%
    st_make_valid()
  
  chunk_joined <- combined_sf %>%
    st_join(
      tdwg %>%
        select(
          ISO_Code,
          Level_4_Na,
          Level4_cod,
          Level4_2,
          Level3_cod,
          Level2_cod,
          Level1_cod
        ),
      left = TRUE
    ) %>%
    st_drop_geometry()
  
} else {
  
  chunk_joined <- chunk_coords %>%
    mutate(
      ISO_Code = NA_character_,
      Level_4_Na = NA_character_,
      Level4_cod = NA_character_,
      Level4_2 = NA_character_,
      Level3_cod = NA_character_,
      Level2_cod = NA_character_,
      Level1_cod = NA_character_
    )
}

chunk_no_coords = chunk_no_coords %>% mutate_at(vars(Level2_cod, Level1_cod), as.numeric)

combined_tiers_tdwg <- bind_rows(chunk_joined, chunk_no_coords) %>%
  arrange(row_id_original)

write_tsv(combined_tiers_tdwg, './outputs/5_distribution/combined_tiers_tdwg.tsv')





# The columns from combined_tiers_tdwg are largely NA at this stage of the pipeline (that's
# exactly why the rescue steps below exist), which is what makes
# read_tsv's type-guessing unreliable for them in the first place --
# force character explicitly above rather than relying on col_guess()
# picking the right type from a mostly-empty sample.


# `row_id_original` is used as a join key in several places below
# (rescue_tdwg_level results, manual TDWG fixes). If the chunked
# upstream processing (5_1_tdwg_spatial_join_chunk.R) ever produces
# overlapping chunks, the same row can end up duplicated with an
# identical row_id_original, which would silently fan out any
# left_join() keyed on it. Check for and drop exact duplicates here,
# and save what was removed for auditing.

dup_row_ids <- combined_tiers_tdwg %>%
  filter(duplicated(row_id_original) | duplicated(row_id_original, fromLast = TRUE))

if (nrow(dup_row_ids) > 0) {
  
  message(
    nrow(dup_row_ids), " row(s) share a duplicated row_id_original ",
    "(", n_distinct(dup_row_ids$row_id_original), " distinct id(s) affected). ",
    "Writing them to duplicated_row_id_original.tsv for review, and keeping ",
    "only the first occurrence of each."
  )
  
  dir.create('./outputs/5_distribution/dup_row_ids', recursive = TRUE, showWarnings = FALSE)
  write_tsv(dup_row_ids, './outputs/5_distribution/dup_row_ids/duplicated_row_id_original.tsv')
  
  combined_tiers_tdwg <- combined_tiers_tdwg %>%
    distinct(row_id_original, .keep_all = TRUE)
}



#check for those without any tdwg


no_tdwg <- combined_tiers_tdwg %>%
  filter(
    is.na(ISO_Code) &
      is.na(Level_4_Na) &
      is.na(Level4_cod) &
      is.na(Level4_2) &
      is.na(Level3_cod) &
      is.na(Level2_cod) &
      is.na(Level1_cod)
  )


#Checking them via QGIS, most of them don't have any tdwg info because their are too close to the coast.
#Let's use a nearest-polygon rescue step to try capture the tdwg info from these data

# -----------------------------
# 3. Convert unmatched rows to sf
# -----------------------------

no_tdwg_sf <- no_tdwg %>%
  filter(
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326,
    remove = FALSE
  )

# -----------------------------
# 5. Function to rescue nearest polygon per TDWG level
# -----------------------------

rescue_tdwg_level <- function(points_sf,
                              tdwg_poly,
                              level,
                              max_dist_m = 10000,
                              target_crs = 3857) {
  
  # Guard against the two conditions that produce the cryptic
  # "zero-length inputs cannot be mixed with those of non-zero
  # length" error out of st_distance(): either input having 0 rows,
  # or tdwg_poly missing a CRS (st_transform then silently produces
  # an sf object sf can't reconcile with points_proj downstream).
  
  if (nrow(points_sf) == 0) {
    return(
      tibble(
        row_id_original = integer(0),
        !!paste0("rescued_tdwg", level, "_ISO_Code") := character(0),
        !!paste0("rescued_tdwg", level, "_Level_4_Na") := character(0),
        !!paste0("rescued_tdwg", level, "_Level4_cod") := character(0),
        !!paste0("rescued_tdwg", level, "_Level4_2") := character(0),
        !!paste0("rescued_tdwg", level, "_Level3_cod") := character(0),
        !!paste0("rescued_tdwg", level, "_Level2_cod") := character(0),
        !!paste0("rescued_tdwg", level, "_Level1_cod") := character(0),
        !!paste0("nearest_tdwg", level, "_dist_m") := numeric(0)
      )
    )
  }
  
  if (nrow(tdwg_poly) == 0) {
    stop(
      "rescue_tdwg_level(): tdwg_poly has 0 features -- check that the ",
      "shapefile it came from was read correctly (all .shp/.shx/.dbf/.prj ",
      "files present, correct path)."
    )
  }
  
  if (is.na(st_crs(tdwg_poly))) {
    stop(
      "rescue_tdwg_level(): tdwg_poly has no CRS. Set st_crs(tdwg_poly) ",
      "explicitly before calling this function."
    )
  }
  
  points_proj <- points_sf %>%
    st_transform(target_crs)
  
  tdwg_proj <- tdwg_poly %>%
    st_transform(target_crs)
  
  nearest_id <- st_nearest_feature(points_proj, tdwg_proj)
  
  if (length(nearest_id) != nrow(points_proj) || anyNA(nearest_id)) {
    stop(
      "rescue_tdwg_level(): st_nearest_feature() returned ",
      length(nearest_id), " match(es) for ", nrow(points_proj),
      " point(s) (NAs present: ", anyNA(nearest_id), "). This points to ",
      "invalid/empty geometries in points_sf or tdwg_poly -- try ",
      "st_make_valid() on both before calling this function."
    )
  }
  
  nearest_dist_m <- st_distance(
    points_proj,
    tdwg_proj[nearest_id, ],
    by_element = TRUE
  )
  
  rescued <- points_proj %>%
    mutate(
      nearest_dist_m_tmp = as.numeric(nearest_dist_m),
      rescued_ISO_Code_tmp   = tdwg_proj$ISO_Code[nearest_id],
      rescued_Level_4_Na_tmp = tdwg_proj$Level_4_Na[nearest_id],
      rescued_Level4_cod_tmp = tdwg_proj$Level4_cod[nearest_id],
      rescued_Level4_2_tmp   = tdwg_proj$Level4_2[nearest_id],
      rescued_Level3_cod_tmp = tdwg_proj$Level3_cod[nearest_id],
      rescued_Level2_cod_tmp = tdwg_proj$Level2_cod[nearest_id],
      rescued_Level1_cod_tmp = tdwg_proj$Level1_cod[nearest_id]
    ) %>%
    st_drop_geometry() %>%
    filter(nearest_dist_m_tmp <= max_dist_m) %>%
    select(
      row_id_original,
      rescued_ISO_Code_tmp,
      rescued_Level_4_Na_tmp,
      rescued_Level4_cod_tmp,
      rescued_Level4_2_tmp,
      rescued_Level3_cod_tmp,
      rescued_Level2_cod_tmp,
      rescued_Level1_cod_tmp,
      nearest_dist_m_tmp
    )
  
  names(rescued) <- c(
    "row_id_original",
    paste0("rescued_tdwg", level, "_ISO_Code"),
    paste0("rescued_tdwg", level, "_Level_4_Na"),
    paste0("rescued_tdwg", level, "_Level4_cod"),
    paste0("rescued_tdwg", level, "_Level4_2"),
    paste0("rescued_tdwg", level, "_Level3_cod"),
    paste0("rescued_tdwg", level, "_Level2_cod"),
    paste0("rescued_tdwg", level, "_Level1_cod"),
    paste0("nearest_tdwg", level, "_dist_m")
  )
  
  rescued
}

# -----------------------------
# 6. Rescue TDWG
# -----------------------------

max_dist_m <- 30000

tdwg = st_read('./tdwg_df/WGSRPD/level4/level4.shp')

# ---- Diagnostics ----
# The "zero-length inputs cannot be mixed with those of non-zero
# length" error from st_distance() inside rescue_tdwg_level() means
# one of the two spatial layers below has 0 rows/features when the
# other doesn't. Print both here so it's obvious which one before
# the function call fails with a much less informative message.

message("no_tdwg_sf: ", nrow(no_tdwg_sf), " point(s), CRS = ", st_crs(no_tdwg_sf)$input)
message("tdwg:       ", nrow(tdwg), " polygon(s), CRS = ", st_crs(tdwg)$input)

if (nrow(tdwg) == 0) {
  stop(
    "tdwg has 0 features after st_read('./tdwg_df/WGSRPD/level4/level4.shp'). ",
    "This usually means the .shp is missing its companion .shx/.dbf/.prj files, ",
    "the path is wrong relative to the current working directory (check getwd()), ",
    "or the download is incomplete/corrupted. Re-download the WGSRPD level4 ",
    "shapefile and confirm all four files (.shp, .shx, .dbf, .prj) are present ",
    "in ./tdwg_df/WGSRPD/level4/ before re-running."
  )
}

if (is.na(st_crs(tdwg))) {
  stop(
    "tdwg has no CRS (likely a missing .prj file). Set it explicitly, e.g. ",
    "st_crs(tdwg) <- 4326 (only if you know the shapefile is really in WGS84), ",
    "before calling rescue_tdwg_level()."
  )
}

if (nrow(no_tdwg_sf) == 0) {
  message(
    "no_tdwg_sf has 0 rows -- nothing to rescue at this step. Skipping the ",
    "rescue_tdwg_level() call and creating an empty rescued_tdwg instead."
  )
  rescued_tdwg <- tibble(
    row_id_original = integer(0),
    rescued_tdwg4_ISO_Code = character(0),
    rescued_tdwg4_Level_4_Na = character(0),
    rescued_tdwg4_Level4_cod = character(0),
    rescued_tdwg4_Level4_2 = character(0),
    rescued_tdwg4_Level3_cod = character(0),
    rescued_tdwg4_Level2_cod = character(0),
    rescued_tdwg4_Level1_cod = character(0),
    nearest_tdwg4_dist_m = numeric(0)
  )
} else {
  rescued_tdwg <- rescue_tdwg_level(no_tdwg_sf, tdwg, level = 4, max_dist_m = max_dist_m)
}

# -----------------------------
# 7. Add rescued matches back to original table
# -----------------------------

# Define column mapping: target columns in combined_tiers_tdwg -> source columns in rescued_tdwg
col_mapping <- c(
  "ISO_Code"   = "rescued_tdwg4_ISO_Code",
  "Level_4_Na" = "rescued_tdwg4_Level_4_Na",
  "Level4_cod" = "rescued_tdwg4_Level4_cod",
  "Level4_2"   = "rescued_tdwg4_Level4_2",
  "Level3_cod" = "rescued_tdwg4_Level3_cod",
  "Level2_cod" = "rescued_tdwg4_Level2_cod",
  "Level1_cod" = "rescued_tdwg4_Level1_cod"
)

combined_tiers_tdwg <- combined_tiers_tdwg %>%
  left_join(
    rescued_tdwg %>% select(row_id_original, all_of(col_mapping)),
    by = "row_id_original",
    suffix = c("", "_new")
  ) %>%
  mutate(across(
    all_of(names(col_mapping)),
    # as.character() on both sides guards against any upstream type
    # mismatch (e.g. a mostly-NA column read_tsv guessed as
    # double/logical instead of character) -- these are TDWG code
    # columns, never meant to be numeric, so this cast is always safe.
    ~ coalesce(
      as.character(get(paste0(cur_column(), "_new"))),
      as.character(.x)
    )
  )) %>%
  select(-ends_with("_new"))



still_no_tdwg <- combined_tiers_tdwg %>%
  filter(
    is.na(ISO_Code) &
      is.na(Level_4_Na) &
      is.na(Level4_cod) &
      is.na(Level4_2) &
      is.na(Level3_cod) &
      is.na(Level2_cod) &
      is.na(Level1_cod)
  )


#These are basically from islands/archipelagos, and the scale of the shapefiles is not enough to capture the data
#Let's add the tdwg manually for those 236 remaining 


unique(still_no_tdwg$country_only)

#French Polynesia
#Society Islands and French Polynesia
sci_oo = tdwg %>% filter(Level_4_Na == 'Society Is.')
#Marquesa Islands
mrq_oo = tdwg %>% filter(Level_4_Na == 'Marquesas')
#Rapa (austral islands, Tubuai)
tub_oo = tdwg %>% filter(Level_4_Na == 'Tubuai Is.')

#São Tomé & Principe
ggi_st = tdwg %>% filter(Level_4_Na == 'São Tomé')

#Comoros
com_co = tdwg %>% filter(Level_4_Na == 'Comoros')

#Bahamas
bah_oo =  tdwg %>% filter(Level_4_Na == 'Bahamas')

#Guinea
gui_oo = tdwg %>% filter(Level_4_Na == 'Guinea')

#Cape Verde
cvi_oo =  tdwg %>% filter(Level_4_Na == 'Cape Verde')

#Kiribati
lin_ki =  tdwg %>% filter(Level_4_Na == 'Kiribati Line Is.')

#Palau
crl_pa = tdwg %>% filter(Level_4_Na == 'Palau')

#Micronesia, Carolina Islands
crl_mf = tdwg %>% filter(Level_4_Na == 'Micronesia Federated States')

#Marshall Islands
mrs_oo = tdwg %>% filter(Level_4_Na == 'Marshall Is.' )

#Australia
#Tasmania
tas_oo = tdwg %>% filter(Level_4_Na == 'Tasmania')
#NSW
nsw_ns = tdwg %>% filter(Level_4_Na == 'New South Wales')

#New Zealand North
nzn_oo = tdwg %>% filter(Level_4_Na == 'New Zealand North')

#Canada
#Nunavut
nun_oo = tdwg %>% filter(Level_4_Na == 'Nunavut')
#Ontario
ont_oo = tdwg %>% filter(Level_4_Na == 'Ontario')

#Chile
jnf_oo = tdwg %>% filter(Level_4_Na == 'Juan Fernández Is.')

#Taiwan/China
tai_oo = tdwg %>% filter(Level_4_Na == 'Taiwan')

#India
and_an = tdwg %>% filter(Level_4_Na == 'Andaman Is.')
ncb_oo = tdwg %>% filter(Level_4_Na == 'Nicobar Is.')

#Japan
jap_hn = tdwg %>% filter(Level_4_Na == 'Honshu')

#Thailand
tha_oo = tdwg %>% filter(Level_4_Na == 'Thailand')

#

tdwg_cols <- c(
  "ISO_Code", "Level_4_Na", "Level4_cod", "Level4_2",
  "Level3_cod", "Level2_cod", "Level1_cod")


tdwg_replacements <- list(
  list(search_col = "country_only", pattern = regex("French Polynesia", ignore_case = TRUE), value = sci_oo),
  list(search_col = "locality",     pattern = regex("Society Islands|Society Is\\.", ignore_case = TRUE), value = sci_oo),
  list(search_col = "locality",     pattern = regex("Marquesas", ignore_case = TRUE), value = mrq_oo),
  list(search_col = "locality",     pattern = regex("Tubuai|Austral Islands|Rapa", ignore_case = TRUE), value = tub_oo),
  list(search_col = "country_only", pattern = regex("Sao Tome & Principe", ignore_case = TRUE), value = ggi_st),
  list(search_col = "country_only", pattern = regex("Comoros", ignore_case = TRUE), value = com_co),
  list(search_col = "country_only", pattern = regex("Bahamas", ignore_case = TRUE), value = bah_oo),
  list(search_col = "country_only", pattern = regex("Guinea", ignore_case = TRUE), value = gui_oo),
  list(search_col = "country_only", pattern = regex("Cape Verde", ignore_case = TRUE), value = cvi_oo),
  list(search_col = "country_only", pattern = regex("Kiribati", ignore_case = TRUE), value = lin_ki),
  list(search_col = "country_only", pattern = regex("Palau", ignore_case = TRUE), value = crl_pa),
  list(search_col = "locality", pattern = regex("Kosrae, Caroline Islands", ignore_case = TRUE), value = crl_mf),
  list(search_col = "country_only", pattern = regex("Marshall Islands", ignore_case = TRUE), value = mrs_oo),
  list(search_col = "locality",     pattern = regex("Tasmania", ignore_case = TRUE), value = tas_oo),
  list(search_col = "locality",     pattern = regex("NSW", ignore_case = TRUE), value = nsw_ns),
  list(search_col = "country",     pattern = regex("New Zealand: North Island", ignore_case = TRUE), value = nzn_oo),
  list(search_col = "locality",     pattern = regex("Nunavut", ignore_case = TRUE), value = nun_oo),
  list(search_col = "locality",     pattern = regex("Ontario", ignore_case = TRUE), value = ont_oo),
  list(search_col = "locality",     pattern = regex("Juan Fernandez", ignore_case = TRUE), value = jnf_oo),
  list(search_col = "locality",     pattern = regex("Taiwan", ignore_case = TRUE), value = tai_oo),
  list(search_col = "locality",     pattern = regex("Andaman", ignore_case = TRUE), value = and_an),
  list(search_col = "locality",     pattern = regex("Nicobar", ignore_case = TRUE), value = ncb_oo),
  list(search_col = "locality",     pattern = regex("Aogashima island", ignore_case = TRUE), value = jap_hn),
  list(search_col = "country_only", pattern = regex("Thailand", ignore_case = TRUE), value = tha_oo)
)

for (r in tdwg_replacements) {
  
  idx <- str_detect(
    still_no_tdwg[[r$search_col]],
    r$pattern
  )
  
  idx[is.na(idx)] <- FALSE
  
  replacement <- st_drop_geometry(r$value[1, tdwg_cols])
  
  for (col in tdwg_cols) {
    replacement[[col]] <- as.character(replacement[[col]])
  }
  
  still_no_tdwg[idx, tdwg_cols] <- replacement
}


combined_tiers_tdwg <- combined_tiers_tdwg %>%
  left_join(
    still_no_tdwg %>%
      select(row_id_original, all_of(tdwg_cols)) %>%
      rename_with(~ paste0("manual_", .x), all_of(tdwg_cols)),
    by = "row_id_original"
  ) %>%
  mutate(
    across(all_of(tdwg_cols), as.character),
    across(starts_with("manual_"), as.character),
    
    ISO_Code   = coalesce(ISO_Code, manual_ISO_Code),
    Level_4_Na = coalesce(Level_4_Na, manual_Level_4_Na),
    Level4_cod = coalesce(Level4_cod, manual_Level4_cod),
    Level4_2   = coalesce(Level4_2, manual_Level4_2),
    Level3_cod = coalesce(Level3_cod, manual_Level3_cod),
    Level2_cod = coalesce(Level2_cod, manual_Level2_cod),
    Level1_cod = coalesce(Level1_cod, manual_Level1_cod)
  ) %>%
  select(-starts_with("manual_"))

#check for NAs 
#Must be empty here

still_no_tdwg = combined_tiers_tdwg %>%
       filter(
         is.na(ISO_Code) &
           is.na(Level_4_Na) &
           is.na(Level4_cod) &
           is.na(Level4_2) &
           is.na(Level3_cod) &
           is.na(Level2_cod) &
           is.na(Level1_cod)
       )



write_tsv(combined_tiers_tdwg, './outputs/5_distribution/combined_tiers_tdwg_clean.tsv')

#####################################################################################################################

#Now, let's see if we can pinpoint if our accessions (combined_tiers_tdwg_clean.tsv) falls between a known native or introduced area, based on
#powo/glonaf datasets

where_they_are <- combined_tiers_tdwg %>%
  left_join(
    distribution_comparison %>%
      select(
        taxa_accepted,
        powo_id,
        corrected_native_distribution,
        combined_introduced_distribution
      ),
    by = "taxa_accepted"
  ) %>%
  mutate(
    native = case_when(
      str_detect(
        corrected_native_distribution,
        paste0("(^|\\|)", Level3_cod, "(\\||$)")
      ) ~ 1,
      TRUE ~ 0
    ),
    
    introduced = case_when(
      str_detect(
        combined_introduced_distribution,
        paste0("(^|\\|)", Level3_cod, "(\\||$)")
      ) ~ 1,
      TRUE ~ 0
    ),
    
    somewhere = case_when(
      native == 0 & introduced == 0 ~ 1,
      TRUE ~ 0
    )
  )

#How many sequences do we have for each distribution?

where_they_are_no_nas = where_they_are %>% filter(!is.na(corrected_native_distribution))
#167,102

where_they_are_nas = where_they_are %>% filter(is.na(corrected_native_distribution))
#4,088 #These are accessions where we couldn't get the native and invasive distribution areas, or they overlapped.

where_summary <- where_they_are_no_nas %>%
  group_by(taxa_accepted) %>%
  summarise(
    n_native = sum(native, na.rm = TRUE),
    n_introduced = sum(introduced, na.rm = TRUE),
    n_somewhere = sum(somewhere, na.rm = TRUE),
    total_records = n(),
    .groups = "drop"
  )

#We need at least 5 accessions total for each species

where_summary_5 = where_summary %>%
  filter(total_records >= 5)

#Filter only those species with 5 or more accessions with known native/introduced distribution

where_summary_5_ni = where_summary %>%
  filter(n_native >= 5,
         n_introduced >= 5)
#747 species

seqs_to_parse = where_they_are %>%
  filter(taxa_accepted %in% where_summary_5_ni$taxa_accepted)
#41,892 accessions

write_tsv(seqs_to_parse, './outputs/seqs_to_parse.tsv')


#####################################################################################################################
# SOMEWHERE DIAGNOSTICS & RESCUE
#####################################################################################################################
#
# `somewhere` (native == 0 & introduced == 0) lumps together several
# different failure modes that call for different fixes:
#
#   (a) precision artifacts   - the point's exact TDWG3 unit isn't
#       listed, but a neighboring TDWG3 unit (within a small radius)
#       IS -- expected mainly for tier 4c/4g (GEOLocate) and tier 5
#       (country-centroid) records, which are the least precise.
#   (b) resolution artifacts  - POWO's L3-level range under-resolves
#       a broader native/introduced region that the accession's
#       own TDWG2/TDWG1 unit is part of.
#   (c) data gaps             - the GloNAF region a taxon is
#       introduced to never resolved to a proper TDWG code (see
#       `still_unresolved` / `iso_fallbacks` above), so it silently
#       dropped out of `combined_introduced_distribution`.
#   (d) taxonomic mismatches  - `taxa_accepted` didn't line up
#       cleanly with POWO's accepted concept, so its native/
#       introduced list is itself incomplete for reasons unrelated
#       to any single accession.
#   (e) genuinely uncertain/undetermined status.
#
# This section flags (a)-(d) explicitly and stratifies by tier and
# by taxon, instead of lumping everything into an unexplained
# "somewhere". Only (a) and (b) are folded back into native/
# introduced automatically (in a NEW `_v2` table, kept alongside the
# original); (c) and (d) need a manual look (edit
# `manual_tdwg_fixes` / check the taxon's POWO linkage) rather than
# a per-accession relabel.
#####################################################################################################################

# ------------------------------------------------------------------
# Small helpers for membership/overlap checks over "|"-collapsed
# code strings, NA-safe.
# ------------------------------------------------------------------

any_overlap <- function(codes_str_a, codes_str_b) {
  if (is.na(codes_str_a) || is.na(codes_str_b) ||
      codes_str_a == "" || codes_str_b == "") {
    return(FALSE)
  }
  any(
    str_split(codes_str_a, "\\|")[[1]] %in%
      str_split(codes_str_b, "\\|")[[1]]
  )
}

code_in_list <- function(code, codes_str) {
  if (is.na(code) || is.na(codes_str) || codes_str == "") {
    return(FALSE)
  }
  code %in% str_split(codes_str, "\\|")[[1]]
}


# ------------------------------------------------------------------
# 1. Tier-stratified diagnostic
# ------------------------------------------------------------------
# Cheap first look: how much of "somewhere" even concentrates in the
# lower-precision tiers (5, 4c) before investing in anything else.

somewhere_by_tier <- where_they_are_no_nas %>%
  group_by(tier) %>%
  summarise(
    n_total = n(),
    n_somewhere = sum(somewhere == 1, na.rm = TRUE),
    pct_somewhere = round(100 * n_somewhere / n_total, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_somewhere))

print(somewhere_by_tier)

write_tsv(somewhere_by_tier, './outputs/5_distribution/somewhere_by_tier.tsv')


# ------------------------------------------------------------------
# 2. Rescue (a): nearby-TDWG3 precision check
# ------------------------------------------------------------------
# For each "somewhere" accession, check whether ANY TDWG4 polygon
# within `rescue_buffer_m` of the point belongs to a TDWG3 unit that
# IS listed in the taxon's native or introduced range. Reuses the
# `tdwg` (level4) polygons already loaded above for the no_tdwg
# rescue step.

rescue_buffer_m <- 50000  # 50 km; tune as needed

somewhere_pts <- where_they_are_no_nas %>%
  filter(
    somewhere == 1,
    !is.na(decimalLatitude),
    !is.na(decimalLongitude)
  )

somewhere_sf <- somewhere_pts %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326,
    remove = FALSE
  )

somewhere_proj <- st_transform(somewhere_sf, 3857)
tdwg_proj <- st_transform(tdwg, 3857)

nearby_idx <- st_is_within_distance(
  somewhere_proj,
  tdwg_proj,
  dist = rescue_buffer_m
)

nearby_l3_codes <- map_chr(nearby_idx, function(idx) {
  if (length(idx) == 0) return(NA_character_)
  codes <- unique(na.omit(tdwg_proj$Level3_cod[idx]))
  if (length(codes) == 0) return(NA_character_)
  paste(codes, collapse = "|")
})

somewhere_pts <- somewhere_pts %>%
  mutate(
    nearby_l3_codes = nearby_l3_codes,
    
    nearby_native_match = map2_lgl(
      nearby_l3_codes, corrected_native_distribution, any_overlap
    ),
    
    nearby_introduced_match = map2_lgl(
      nearby_l3_codes, combined_introduced_distribution, any_overlap
    )
  )

message(
  "Proximity rescue (<=", rescue_buffer_m / 1000, " km): ",
  sum(somewhere_pts$nearby_native_match | somewhere_pts$nearby_introduced_match),
  " / ", nrow(somewhere_pts),
  " 'somewhere' accessions have a nearby TDWG3 unit in the taxon's known range."
)


# ------------------------------------------------------------------
# 3. Rescue (b): coarser TDWG2 / TDWG1 fallback
# ------------------------------------------------------------------
# Aggregate POWO's L3-level native/introduced codes up to L2 and L1
# (using the TDWG hierarchy already present in `tdwg`), then check
# whether the accession's OWN Level2_cod / Level1_cod (from the
# spatial join earlier in this script) falls inside either range.

l3_to_l2l1 <- tdwg %>%
  st_drop_geometry() %>%
  distinct(Level3_cod, Level2_cod, Level1_cod) %>%
  filter(!is.na(Level3_cod))

expand_to_level <- function(l3_codes_str, level_col) {
  if (is.na(l3_codes_str) || l3_codes_str == "") return(NA_character_)
  codes <- str_split(l3_codes_str, "\\|")[[1]]
  higher <- unique(na.omit(
    l3_to_l2l1[[level_col]][l3_to_l2l1$Level3_cod %in% codes]
  ))
  if (length(higher) == 0) return(NA_character_)
  paste(higher, collapse = "|")
}

distribution_comparison_levels <- distribution_comparison %>%
  mutate(
    native_l2 = map_chr(
      corrected_native_distribution, expand_to_level, level_col = "Level2_cod"
    ),
    native_l1 = map_chr(
      corrected_native_distribution, expand_to_level, level_col = "Level1_cod"
    ),
    introduced_l2 = map_chr(
      combined_introduced_distribution, expand_to_level, level_col = "Level2_cod"
    ),
    introduced_l1 = map_chr(
      combined_introduced_distribution, expand_to_level, level_col = "Level1_cod"
    )
  ) %>%
  select(taxa_accepted, powo_id, native_l2, native_l1, introduced_l2, introduced_l1)

somewhere_pts <- somewhere_pts %>%
  left_join(
    distribution_comparison_levels,
    by = c("taxa_accepted", "powo_id")
  ) %>%
  mutate(
    rescued_l2_native      = map2_lgl(Level2_cod, native_l2, code_in_list),
    rescued_l2_introduced  = map2_lgl(Level2_cod, introduced_l2, code_in_list),
    rescued_l1_native      = map2_lgl(Level1_cod, native_l1, code_in_list),
    rescued_l1_introduced  = map2_lgl(Level1_cod, introduced_l1, code_in_list)
  )

message(
  "L2/L1 fallback rescue: ",
  sum(with(
    somewhere_pts,
    rescued_l2_native | rescued_l2_introduced | rescued_l1_native | rescued_l1_introduced
  )),
  " / ", nrow(somewhere_pts),
  " 'somewhere' accessions fall inside a broader native/introduced region."
)


# ------------------------------------------------------------------
# 4. Data-gap flag (c): unresolved / ISO-fallback GloNAF regions
# ------------------------------------------------------------------
# `still_unresolved` and `iso_fallbacks` were already computed above
# while building `taxa_glonaf_distribution_expanded`. Save them for
# manual review (each row is a GloNAF region that never resolved
# cleanly to a TDWG3/4 code), and flag which taxa in `somewhere_pts`
# are affected by at least one such unresolved region.

write_tsv(still_unresolved, './outputs/5_distribution/glonaf_tdwg_still_unresolved.tsv')
write_tsv(iso_fallbacks, './outputs/5_distribution/glonaf_tdwg_iso_fallbacks.tsv')

taxa_with_unresolved_glonaf <- unique(still_unresolved$taxa_accepted)

somewhere_pts <- somewhere_pts %>%
  mutate(
    taxon_has_unresolved_glonaf_region = taxa_accepted %in% taxa_with_unresolved_glonaf
  )


# ------------------------------------------------------------------
# 5. Taxonomic concentration diagnostic (d)
# ------------------------------------------------------------------
# If a small number of species account for a large share of
# "somewhere", that points to a taxon-level (POWO linkage/synonymy)
# problem rather than many independent georeferencing errors.

somewhere_by_taxon <- where_they_are_no_nas %>%
  group_by(taxa_accepted) %>%
  summarise(
    n_total = n(),
    n_somewhere = sum(somewhere == 1, na.rm = TRUE),
    pct_somewhere = round(100 * n_somewhere / n_total, 1),
    .groups = "drop"
  ) %>%
  filter(n_total >= 5) %>%
  arrange(desc(pct_somewhere))

write_tsv(somewhere_by_taxon, './outputs/5_distribution/somewhere_by_taxon.tsv')

message("Top 10 taxa by 'somewhere' rate (>=5 accessions):")
print(head(somewhere_by_taxon, 10))


# ------------------------------------------------------------------
# 6. Combine into a single somewhere_reason flag
# ------------------------------------------------------------------
# Precedence: an accession can in principle satisfy more than one
# rescue at once; proximity (a) is checked first since it directly
# addresses the accession's own coordinates, then the coarser
# resolution fallback (b), then the taxon-level data-gap flag (c).
# Anything matching BOTH native and introduced under the same rescue
# is left "ambiguous" rather than guessed.

somewhere_pts <- somewhere_pts %>%
  mutate(
    proximity_ambiguous = nearby_native_match & nearby_introduced_match,
    level_ambiguous = (rescued_l2_native | rescued_l1_native) &
      (rescued_l2_introduced | rescued_l1_introduced),
    
    somewhere_reason = case_when(
      proximity_ambiguous ~ "ambiguous_rescue_nearby_tdwg",
      nearby_native_match | nearby_introduced_match ~ "precision_artifact_nearby_tdwg",
      level_ambiguous ~ "ambiguous_rescue_l2_l1",
      rescued_l2_native | rescued_l1_native |
        rescued_l2_introduced | rescued_l1_introduced ~ "resolution_artifact_l2_l1",
      taxon_has_unresolved_glonaf_region ~ "data_gap_unresolved_glonaf_region",
      TRUE ~ "unexplained"
    ),
    
    # direction only meaningful for the two non-ambiguous rescue reasons
    rescued_direction = case_when(
      somewhere_reason == "precision_artifact_nearby_tdwg" & nearby_native_match ~ "native",
      somewhere_reason == "precision_artifact_nearby_tdwg" & nearby_introduced_match ~ "introduced",
      somewhere_reason == "resolution_artifact_l2_l1" &
        (rescued_l2_native | rescued_l1_native) ~ "native",
      somewhere_reason == "resolution_artifact_l2_l1" &
        (rescued_l2_introduced | rescued_l1_introduced) ~ "introduced",
      TRUE ~ NA_character_
    )
  )

write_tsv(
  somewhere_pts %>%
    select(
      taxa_accepted, uid, tier,
      decimalLatitude, decimalLongitude,
      Level3_cod, Level2_cod, Level1_cod,
      nearby_native_match, nearby_introduced_match,
      rescued_l2_native, rescued_l2_introduced,
      rescued_l1_native, rescued_l1_introduced,
      taxon_has_unresolved_glonaf_region,
      somewhere_reason, rescued_direction
    ),
  './outputs/5_distribution/somewhere_diagnostics.tsv'
)

message("somewhere_reason breakdown:")
print(count(somewhere_pts, somewhere_reason, sort = TRUE))


# ------------------------------------------------------------------
# 7. Fold non-ambiguous rescues back into native/introduced
# ------------------------------------------------------------------
# Only "precision_artifact_nearby_tdwg" and "resolution_artifact_l2_l1"
# are folded back in, each into the specific direction
# (`rescued_direction`) it actually matched. Ambiguous rescues, data
# gaps, and unexplained cases are left as "somewhere" -- those need a
# closer look (better georeferencing, editing
# `manual_tdwg_fixes`/GloNAF linkage, or checking the taxon's POWO
# concept), not an automatic relabel.
#
# IMPORTANT: match on the COMBINATION of `row_id_original` AND
# `powo_id`, not on either alone. `where_they_are` is built via
# `left_join(distribution_comparison %>% select(taxa_accepted,
# powo_id, ...), by = "taxa_accepted")`, which fans a single
# accession out into multiple rows whenever a taxa_accepted string
# matches more than one powo_id -- and BOTH `uid` and
# `row_id_original` get carried through unchanged into every fanned-
# out copy, so neither is unique on its own once that happens. What
# IS unique per fanned-out row is (row_id_original, powo_id): each
# copy has a different powo_id (that's literally what's being fanned
# out on) paired with the same row_id_original. A left_join on both
# columns together -- instead of a single-column %in% check -- is
# the correct way to target only the specific row that was rescued.

diagnosed_fanout <- where_they_are_no_nas %>%
  count(row_id_original, name = "n_rows_for_this_accession") %>%
  filter(n_rows_for_this_accession > 1)

if (nrow(diagnosed_fanout) > 0) {
  message(
    nrow(diagnosed_fanout), " row_id_original value(s) appear more than ",
    "once in where_they_are_no_nas (", sum(diagnosed_fanout$n_rows_for_this_accession),
    " rows total) -- i.e. that many accessions matched more than one ",
    "powo_id for their taxa_accepted string in distribution_comparison. ",
    "This is the taxa_accepted -> powo_id fan-out; see ",
    "accession_powo_fanout.tsv. Worth checking WHY those names match ",
    "multiple powo_id (e.g. distribution_comparison / powo_names not ",
    "deduplicated per accepted name)."
  )
  write_tsv(
    where_they_are_no_nas %>%
      semi_join(diagnosed_fanout, by = "row_id_original") %>%
      select(taxa_accepted, uid, row_id_original, powo_id,
             corrected_native_distribution, combined_introduced_distribution) %>%
      arrange(row_id_original),
    './outputs/5_distribution/accession_powo_fanout.tsv'
  )
}

rescued_flags <- somewhere_pts %>%
  filter(!is.na(rescued_direction)) %>%
  distinct(row_id_original, powo_id, rescued_direction)

where_they_are_v2 <- where_they_are_no_nas %>%
  left_join(rescued_flags, by = c("row_id_original", "powo_id")) %>%
  mutate(
    native = if_else(rescued_direction == "native", 1L, native),
    introduced = if_else(rescued_direction == "introduced", 1L, introduced),
    somewhere = if_else(!is.na(rescued_direction), 0L, somewhere)
  ) %>%
  select(-rescued_direction)

# Sanity check: this must never fire. If it does, something upstream
# (not this rescue step) is producing rows where native/introduced
# were already both 1 before the fold-back -- worth investigating
# `where_they_are` itself rather than this block.
impossible_v2 <- where_they_are_v2 %>%
  group_by(taxa_accepted) %>%
  summarise(
    n_native = sum(native, na.rm = TRUE),
    n_introduced = sum(introduced, na.rm = TRUE),
    n_somewhere = sum(somewhere, na.rm = TRUE),
    total_records = n(),
    .groups = "drop"
  ) %>%
  filter((n_native + n_introduced + n_somewhere) > total_records)

if (nrow(impossible_v2) > 0) {
  warning(
    nrow(impossible_v2), " taxa still show n_native + n_introduced + ",
    "n_somewhere > total_records even after keying the fold-back on ",
    "(row_id_original, powo_id). Writing them to ",
    "impossible_distribution_counts.tsv for review -- this points to a ",
    "PRE-EXISTING issue upstream of this rescue step (e.g. a single row ",
    "of where_they_are_no_nas already had native == 1 AND introduced == 1 ",
    "before this script ever ran the rescue logic)."
  )
  write_tsv(impossible_v2, './outputs/5_distribution/impossible_distribution_counts.tsv')
} else {
  message("No impossible native+introduced+somewhere > total_records counts in where_they_are_v2.")
}

where_summary_v2 <- where_they_are_v2 %>%
  group_by(taxa_accepted) %>%
  summarise(
    n_native = sum(native, na.rm = TRUE),
    n_introduced = sum(introduced, na.rm = TRUE),
    n_somewhere = sum(somewhere, na.rm = TRUE),
    total_records = n(),
    .groups = "drop"
  )

where_summary_5_ni_v2 <- where_summary_v2 %>%
  filter(n_native >= 5, n_introduced >= 5)

seqs_to_parse_v2 <- where_they_are_v2 %>%
  filter(taxa_accepted %in% where_summary_5_ni_v2$taxa_accepted)

#write_tsv(seqs_to_parse_v2, './outputs/5_distribution/seqs_to_parse_v2.tsv')

message(
  "seqs_to_parse:     ", nrow(seqs_to_parse), " rows / ",
  n_distinct(seqs_to_parse$taxa_accepted), " taxa\n",
  "seqs_to_parse_v2:  ", nrow(seqs_to_parse_v2), " rows / ",
  n_distinct(seqs_to_parse_v2$taxa_accepted), " taxa",
  " (after folding back precision/resolution 'somewhere' rescues)"
)

#get all unique values from v2 that it's not in seqs_to_parse

unique_v2 = anti_join(seqs_to_parse_v2, seqs_to_parse, by = 'uid')


final_seqs_to_parse = rbind(seqs_to_parse, unique_v2)
#42,375 rows

write_tsv(final_seqs_to_parse, './outputs/5_distribution/final_seqs_to_parse.tsv')
