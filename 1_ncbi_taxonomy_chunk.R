#1_ncbi_taxonomy_chunk.R

##PART A - MATCHING THE TAXONOMIC BACKBONE OF GLONAF WITH NCBI'S

library(dplyr)
library(purrr)
library(tibble)
library(readr)
library(rentrez)

# -----------------------------
# 1. Arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

start_i <- as.integer(args[1])
end_i   <- as.integer(args[2])

if (is.na(start_i) || is.na(end_i)) {
  stop("Usage: Rscript ncbi_taxonomy_chunk.R <start_i> <end_i>")
}

# -----------------------------
# 2. NCBI settings
# -----------------------------
ncbi_sleep <- 0.34

# Optional, but recommended
Sys.setenv(ENTREZ_KEY = "9263da97822c1746429a7dce2aa901ee1708")
entrez_email <- "marvin.danque@gmail.com"

# -----------------------------
# 3. Input data
# -----------------------------
glo_wcvp <- read_csv("glonaf_taxon_wcvp.csv", col_names = TRUE)
glo_flora <- read_csv("glonaf_flora2.csv", col_names = TRUE)

glonaf_f <- glo_flora %>%
  left_join(glo_wcvp,join_by(taxon_wcvp_id==id)) %>%
  distinct(accepted_plant_name_id,.keep_all = TRUE) %>%
  select(accepted_plant_name_id,taxa_accepted, powo_id_accepted) %>%
  filter(!is.na(accepted_plant_name_id))
#16,429 unique taxa

all_species <- as.vector(glonaf_f$taxa_accepted)

# Prevent end_i from exceeding number of species
end_i <- min(end_i, length(all_species))

species_chunk <- all_species[start_i:end_i]

message("Running species ", start_i, " to ", end_i)
message("Total species in this chunk: ", length(species_chunk))

# -----------------------------
# 4. Helper functions
# -----------------------------
safe_entrez_search <- function(db, term, retmax = 5000) {
  Sys.sleep(ncbi_sleep)
  tryCatch(
    entrez_search(db = db, term = term, retmax = retmax, use_history = TRUE),
    error = function(e) {
      message("Search error for term: ", term)
      message("Error: ", e$message)
      NULL
    }
  )
}

safe_entrez_summary <- function(db, id) {
  Sys.sleep(ncbi_sleep)
  tryCatch(
    entrez_summary(db = db, id = id),
    error = function(e) {
      message("Summary error for id: ", id)
      message("Error: ", e$message)
      NULL
    }
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# -----------------------------
# 5. Find NCBI Taxonomy match
# -----------------------------
get_taxonomy_match <- function(sp_name) {
  
  message("Searching: ", sp_name)
  
  term <- paste0('"', sp_name, '"[Scientific Name]')
  tx <- safe_entrez_search(db = "taxonomy", term = term, retmax = 10)
  
  if (is.null(tx) || is.null(tx$count) || tx$count == 0) {
    return(tibble(
      taxa_accepted = sp_name,
      ncbi_taxid = NA_character_,
      ncbi_scientific_name = NA_character_,
      rank = NA_character_,
      taxonomy_found = FALSE
    ))
  }
  
  tx_sum <- safe_entrez_summary(db = "taxonomy", id = tx$ids[1])
  
  if (is.null(tx_sum)) {
    return(tibble(
      taxa_accepted = sp_name,
      ncbi_taxid = as.character(tx$ids[1]),
      ncbi_scientific_name = NA_character_,
      rank = NA_character_,
      taxonomy_found = FALSE
    ))
  }
  
  tibble(
    taxa_accepted = sp_name,
    ncbi_taxid = as.character(tx$ids[1]),
    ncbi_scientific_name = tx_sum$scientificname %||% NA_character_,
    rank = tx_sum$rank %||% NA_character_,
    taxonomy_found = TRUE
  )
}

### Search for taxonomic matches along the NCBI db


taxonomy_tbl <- map_dfr(species_chunk, get_taxonomy_match)

matched_taxonomy_tbl <- taxonomy_tbl %>%
  filter(taxonomy_found, !is.na(ncbi_taxid), !is.na(ncbi_scientific_name))

unmatched_taxonomy_tbl <- taxonomy_tbl %>%
  filter(!taxonomy_found)

glonaf_matched <- glonaf_f %>%
  inner_join(
    matched_taxonomy_tbl %>%
      dplyr::select(taxa_accepted, ncbi_taxid),
    by = "taxa_accepted") %>%
  select(!taxa_accepted)

# -----------------------------
# 6. Save chunk outputs
# -----------------------------
out_prefix <- paste0(start_i, "_", end_i)

write_tsv(unmatched_taxonomy_tbl, paste0("unmatched_taxonomy_", out_prefix, ".tsv"))
write_tsv(glonaf_matched, paste0("matched_taxonomy_", out_prefix, ".tsv"))

message("Finished chunk: ", start_i, "-", end_i)

### PART B - CHECKING THE UNMATCHED TAXONOMIC BACKBONE.

### Only run this part after the part A is complete.

## First, merge all the matched/unmatched chunks into single files.

#unmatched_files = Sys.glob('unmatched_taxonomy_*.tsv')
#matched_files = Sys.glob('matched_taxonomy_*.tsv')

#unmatched_merged = unmatched_files %>%
#  lapply(read_tsv, show_col_types = FALSE) %>%
#  bind_rows()

#matched_merged = matched_files %>%
#  lapply(read_tsv, show_col_types = FALSE) %>%
#  bind_rows()

#Just to be sure that the names from the unmatched table don't have any misspelling, let's compare them with
#those from the powo backbone
powo_names = read_delim('wcvp_names.csv', delim =  '|')

unmatched_merged_pw = unmatched_merged %>% 
  left_join(glonaf_f, by = 'taxa_accepted') %>%
  select(taxa_accepted, accepted_plant_name_id, powo_id_accepted)
  
unmatched_merged_different = unmatched_merged_pw %>% left_join(powo_names, join_by(powo_id_accepted==powo_id)) %>%
  select(taxa_accepted, taxon_name, accepted_plant_name_id.x, powo_id_accepted) %>%
  mutate(comparision = case_when(
    is.na(taxon_name) ~ 'no_powo_match',
    taxa_accepted == taxon_name ~ 'same',
    taxa_accepted != taxon_name ~ 'different'
  )) %>%
  filter(comparision == 'different')

#These couldn't be matched to ncbi's taxonomy list
unmatched_merged_same <- unmatched_merged_pw %>%
  left_join(powo_names, join_by(powo_id_accepted == powo_id)) %>%
  select(taxa_accepted, taxon_name, accepted_plant_name_id.x, powo_id_accepted) %>%
  mutate(comparision = case_when(
    is.na(taxon_name) ~ 'no_powo_match',
    taxa_accepted == taxon_name ~ 'same',
    taxa_accepted != taxon_name ~ 'different'
  )) %>%
  filter(comparision == 'same' | comparision == 'no_powo_match')


#Check again with the powo name

species_chunk = as.vector(unmatched_merged_different$taxon_name)

taxonomy_tbl2 <- map_dfr(species_chunk, get_taxonomy_match)

matched_taxonomy_tbl2 <- taxonomy_tbl2 %>%
  filter(taxonomy_found, !is.na(ncbi_taxid), !is.na(ncbi_scientific_name))

#These also couldn't be matched to ncbi's taxonomy list
unmatched_taxonomy_tbl2 <- taxonomy_tbl2 %>%
  filter(!taxonomy_found)

#Let's merge what couldn't truly be matched

truly_unmatched <- bind_rows(
  unmatched_merged_same %>% select(name = taxa_accepted),
  unmatched_taxonomy_tbl2 %>% select(name = taxa_accepted)
) %>%
  distinct(.keep_all = TRUE)

#write_tsv(truly_unmatched, './1_ncbi_taxonomy_chunk/truly_unmatched.tsv')
#2,432 unique taxa

#These are the saved recovered matches
glonaf_matched2 <- unmatched_merged_different %>%
  inner_join(
    matched_taxonomy_tbl2 %>%
      dplyr::select(taxa_accepted, ncbi_taxid),
    join_by(taxon_name==taxa_accepted)) %>%
  select(taxon_name, accepted_plant_name_id.x, ncbi_taxid, powo_id_accepted)

#Let's add them back to the matched names

original_glonaf_names = unmatched_merged_different %>%
  left_join(matched_taxonomy_tbl2, join_by(taxon_name==taxa_accepted)) %>%
  filter(taxonomy_found == TRUE) %>%
  select(accepted_plant_name_id.x, ncbi_taxid)

truly_matched <- bind_rows(
  matched_merged %>%
    select(accepted_plant_name_id, ncbi_taxid) %>%
    mutate(ncbi_taxid = as.character(ncbi_taxid)),
  original_glonaf_names %>%
    select(accepted_plant_name_id = accepted_plant_name_id.x, ncbi_taxid) %>%
    mutate(ncbi_taxid = as.character(ncbi_taxid))
)

#write_tsv(truly_matched, './1_ncbi_taxonomy_chunk/truly_matched.tsv')
#13,997 unique taxa
