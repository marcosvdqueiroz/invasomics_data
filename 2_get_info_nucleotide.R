# ============================================================
# Search NCBI nucleotide records and extract locality metadata
# Faster version using batch efetch
# ============================================================

library(rentrez)   # Access NCBI Entrez
library(dplyr)     # Data manipulation
library(purrr)     # map/map_dfr functions
library(tibble)    # tibble data frames
library(readr)     # read_tsv/write_tsv
library(tidyr)     # unnest/pivot_wider
library(stringr)   # string manipulation
library(xml2)      # XML parsing

# -----------------------------
# 1. Prepare taxa list
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript get_info_nucleotide.R START_ROW END_ROW CHUNK_ID")
}

start_row <- as.integer(args[1])
end_row   <- as.integer(args[2])
chunk_id  <- as.integer(args[3])

message("Running chunk: ", chunk_id)
message("Rows: ", start_row, " to ", end_row)


glonaf_matched = read_tsv("./1_ncbi_taxonomy_chunk/truly_matched.tsv", show_col_types = FALSE)

glo_wcvp <- read_csv("glonaf_taxon_wcvp.csv", col_names = TRUE)
glo_flora <- read_csv("glonaf_flora2.csv", col_names = TRUE)

glonaf_f <- glo_flora %>%
  left_join(glo_wcvp,join_by(taxon_wcvp_id==id)) %>%
  distinct(accepted_plant_name_id,.keep_all = TRUE) %>%
  select(accepted_plant_name_id,taxa_accepted, powo_id_accepted) %>%
  filter(!is.na(accepted_plant_name_id))

glonaf_matched_f = glonaf_matched %>%
  left_join(glonaf_f, by='accepted_plant_name_id') %>%
  select(accepted_plant_name_id, powo_id_accepted, ncbi_taxid, taxa_accepted)
  


taxa <- glonaf_matched_f$taxa_accepted %>%
  unique() %>%
  na.omit() %>%
  as.character()

end_row <- min(end_row, length(taxa))

taxa <- taxa[start_row:end_row]

message("Number of taxa in this job: ", length(taxa))

# -----------------------------
# 2. Helper functions
# -----------------------------

# Safely extract a field from an Entrez summary object
safe_extract <- function(x, field) {
  
  # Return NA if object is missing
  if (is.null(x) || length(x) == 0) {
    return(NA_character_)
  }
  
  # Return NA if requested field is absent
  if (!field %in% names(x)) {
    return(NA_character_)
  }
  
  # Extract value
  value <- x[[field]]
  
  # Return NA if value is empty
  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }
  
  # Collapse multiple values into one string
  paste(value, collapse = "|")
}

# Extract gene/CDS product information from GenBank XML
extract_feature_products <- function(record, feature_key, qualifier_name) {
  
  # Find all GenBank features with the requested feature key
  features <- xml_find_all(
    record,
    paste0(".//GBFeature[GBFeature_key='", feature_key, "']")
  )
  
  # Return NA if no such feature exists
  if (length(features) == 0) {
    return(NA_character_)
  }
  
  # Extract qualifier values from each feature
  values <- map_chr(features, function(feature) {
    
    # Find qualifier value inside the feature
    val <- xml_find_all(
      feature,
      paste0(
        ".//GBQualifier[GBQualifier_name='",
        qualifier_name,
        "']/GBQualifier_value"
      )
    ) %>%
      xml_text()
    
    # Return NA if qualifier is absent
    if (length(val) == 0) {
      NA_character_
    } else {
      paste(val, collapse = "|")
    }
  })
  
  # Remove missing values
  values <- values[!is.na(values)]
  
  # Collapse unique values
  if (length(values) == 0) {
    NA_character_
  } else {
    paste(unique(values), collapse = "|")
  }
}

# Extract journal or authors from the first GenBank reference
extract_reference_info <- function(record, field) {
  
  # Find first reference block
  ref <- xml_find_first(record, ".//GBReference")
  
  # Return NA if no reference exists
  if (length(ref) == 0 || inherits(ref, "xml_missing")) {
    return(NA_character_)
  }
  
  # Authors are stored as repeated GBAuthor nodes
  if (field == "authors") {
    vals <- xml_find_all(ref, ".//GBAuthor") %>%
      xml_text()
  } else {
    vals <- xml_find_all(ref, paste0(".//", field)) %>%
      xml_text()
  }
  
  # Collapse reference values
  if (length(vals) == 0) {
    NA_character_
  } else {
    paste(vals, collapse = "|")
  }
}

# -----------------------------
# 3. Batch XML metadata extraction
# -----------------------------

# Fetch and parse XML metadata for many UIDs at once
get_nucleotide_xml_metadata_batch <- function(uids, batch_size = 200) {
  
  # Split UID vector into batches to avoid very large NCBI requests
  uid_batches <- split(uids, ceiling(seq_along(uids) / batch_size))
  
  # Process each batch and combine results
  map_dfr(uid_batches, function(uid_batch) {
    
    # Fetch GenBank XML for all UIDs in this batch
    xml_txt <- tryCatch(
      entrez_fetch(
        db = "nucleotide",
        id = uid_batch,
        rettype = "gb",
        retmode = "xml"
      ),
      error = function(e) NA_character_
    )
    
    # Pause to respect NCBI rate limits
    Sys.sleep(0.34)
    
    # If fetch failed, return empty metadata for all UIDs in batch
    if (is.na(xml_txt)) {
      return(tibble(
        uid = uid_batch,
        gene = NA_character_,
        CDS_product = NA_character_,
        ref_journal = NA_character_,
        ref_authors = NA_character_
      ))
    }
    
    # Parse XML using HUGE option to avoid long text-node errors
    xml_doc <- tryCatch(
      read_xml(xml_txt, options = "HUGE"),
      error = function(e) NULL
    )
    
    # If XML parsing failed, return empty metadata
    if (is.null(xml_doc)) {
      return(tibble(
        uid = uid_batch,
        gene = NA_character_,
        CDS_product = NA_character_,
        ref_journal = NA_character_,
        ref_authors = NA_character_
      ))
    }
    
    # Find each individual GenBank sequence record
    records <- xml_find_all(xml_doc, ".//GBSeq")
    
    # If no records were parsed, return empty metadata
    if (length(records) == 0) {
      return(tibble(
        uid = uid_batch,
        gene = NA_character_,
        CDS_product = NA_character_,
        ref_journal = NA_character_,
        ref_authors = NA_character_
      ))
    }
    
    # Parse each GenBank record
    parsed_records <- map_dfr(records, function(rec) {
      
      # Extract accession from XML
      accession <- xml_find_first(rec, ".//GBSeq_primary-accession") %>%
        xml_text()
      
      # Extract gene names, CDS products, and reference metadata
      tibble(
        accession = accession,
        gene = extract_feature_products(rec, "gene", "gene"),
        CDS_product = extract_feature_products(rec, "CDS", "product"),
        ref_journal = extract_reference_info(rec, "GBReference_journal"),
        ref_authors = extract_reference_info(rec, "authors")
      )
    })
    
    # Use entrez_summary caption/accession matching later, so keep accession here
    parsed_records
  })
}

# -----------------------------
# 4. Query NCBI nucleotide
# -----------------------------

# Search NCBI nucleotide for one taxon and return summary + XML metadata
get_nucleotide_summary <- function(taxon, retmax = 500) {
  
  # Print progress message
  message("Searching: ", taxon)
  
  # Search nucleotide database for records assigned to the taxon
  my_ids <- tryCatch(
    entrez_search(
      db = "nucleotide",
      term = paste0("\"", taxon, "\"[Organism]"),
      retmax = retmax
    ),
    error = function(e) NULL
  )
  
  # Pause to respect NCBI rate limits
  Sys.sleep(0.34)
  
  # If search failed or returned no IDs, return one empty row
  if (is.null(my_ids) || length(my_ids$ids) == 0) {
    return(tibble(
      taxa_accepted = taxon,
      uid = NA_character_,
      caption = NA_character_,
      title = NA_character_,
      subtype = NA_character_,
      subname = NA_character_,
      genome = NA_character_,
      gene = NA_character_,
      CDS_product = NA_character_,
      ref_journal = NA_character_,
      ref_authors = NA_character_
    ))
  }
  
  # Retrieve Entrez summary metadata for all IDs
  my_hits <- tryCatch(
    entrez_summary(
      db = "nucleotide",
      id = my_ids$ids
    ),
    error = function(e) NULL
  )
  
  # Pause to respect NCBI rate limits
  Sys.sleep(0.34)
  
  # If summary retrieval failed, return one empty row
  if (is.null(my_hits)) {
    return(tibble(
      taxa_accepted = taxon,
      uid = NA_character_,
      caption = NA_character_,
      title = NA_character_,
      subtype = NA_character_,
      subname = NA_character_,
      genome = NA_character_,
      gene = NA_character_,
      CDS_product = NA_character_,
      ref_journal = NA_character_,
      ref_authors = NA_character_
    ))
  }
  
  # Convert Entrez summary records into a table
  summary_tbl <- map_dfr(seq_along(my_hits), function(j) {
    
    # Extract one summary hit
    hit <- my_hits[[j]]
    
    # Store selected summary fields
    tibble(
      taxa_accepted = taxon,
      uid = safe_extract(hit, "uid"),
      caption = safe_extract(hit, "caption"),
      title = safe_extract(hit, "title"),
      subtype = safe_extract(hit, "subtype"),
      subname = safe_extract(hit, "subname"),
      genome = safe_extract(hit, "genome")
    )
  })
  
  # Extract valid UIDs for XML batch fetch
  uids_to_fetch <- summary_tbl %>%
    filter(!is.na(uid), uid != "") %>%
    pull(uid) %>%
    unique()
  
  # If there are no UIDs, add empty XML metadata columns
  if (length(uids_to_fetch) == 0) {
    return(
      summary_tbl %>%
        mutate(
          gene = NA_character_,
          CDS_product = NA_character_,
          ref_journal = NA_character_,
          ref_authors = NA_character_
        )
    )
  }
  
  # Fetch XML metadata in batches instead of one UID at a time
  xml_metadata <- get_nucleotide_xml_metadata_batch(
    uids = uids_to_fetch,
    batch_size = 200
  )
  
  # If XML metadata failed completely, create empty accession-based table
  if (nrow(xml_metadata) == 0 || !"accession" %in% names(xml_metadata)) {
    xml_metadata <- tibble(
      accession = summary_tbl$caption,
      gene = NA_character_,
      CDS_product = NA_character_,
      ref_journal = NA_character_,
      ref_authors = NA_character_
    )
  }
  
  # Join XML metadata to summary metadata using accession/caption
  summary_tbl %>%
    left_join(xml_metadata, by = c("caption" = "accession"))
}

# -----------------------------
# 5. Run query
# -----------------------------

# Run search for all taxa and combine all results
nucleotide_summary <- map_dfr(
  taxa,
  ~ get_nucleotide_summary(.x, retmax = 500)
)

# Save raw summary table (if needed)
write_tsv(
  nucleotide_summary,
  paste0("./2_get_info_nucleotide/nucleotide_entrez_summary_raw_", chunk_id, ".tsv")
)

# -----------------------------
# 6. Parse subtype/subname metadata
# -----------------------------

# Convert subtype/subname fields into long-format metadata
metadata_long <- nucleotide_summary %>%
  
  # Add row ID so each sequence remains uniquely identifiable
  mutate(row_id = row_number()) %>%
  
  # Split subtype and subname fields into lists using "|"
  mutate(
    subtype_list = strsplit(as.character(subtype), "\\|"),
    subname_list = strsplit(as.character(subname), "\\|")
  ) %>%
  
  # Pair subtype keys with subname values
  mutate(
    metadata = map2(subtype_list, subname_list, function(keys, values) {
      
      # Return empty metadata if keys or values are missing
      if (length(keys) == 0 || length(values) == 0 ||
          all(is.na(keys)) || all(is.na(values))) {
        return(tibble(
          metadata_field = NA_character_,
          metadata_value = NA_character_
        ))
      }
      
      # Use the maximum length in case key/value vectors differ
      n <- max(length(keys), length(values))
      
      # Create key-value metadata table
      tibble(
        metadata_field = keys[seq_len(n)],
        metadata_value = values[seq_len(n)]
      )
    })
  ) %>%
  
  # Remove temporary list columns
  select(-subtype_list, -subname_list) %>%
  
  # Expand nested metadata tibbles
  unnest(metadata) %>%
  
  # Remove empty metadata fields
  filter(!is.na(metadata_field), metadata_field != "") %>%
  
  # Standardize metadata field names
  mutate(
    metadata_field = str_to_lower(metadata_field),
    metadata_field = str_replace_all(metadata_field, "[^A-Za-z0-9]+", "_"),
    metadata_field = str_replace_all(metadata_field, "^_|_$", "")
  )

# -----------------------------
# 7. Keep locality-related metadata
# -----------------------------

# Define metadata fields related to geography/locality
locality_terms <- c(
  "country",
  "lat_lon",
  "lat_long",
  "latitude_longitude",
  "latitude_and_longitude",
  "geo_loc_name",
  "geographic_location",
  "locality",
  "location",
  "collection_site",
  "collection_location",
  "isolation_source",
  "isolate",
  "altitude",
  "elevation",
  "depth"
)

# Keep only locality metadata
metadata_locality <- metadata_long %>%
  filter(metadata_field %in% locality_terms)

# -----------------------------
# 8. Convert locality metadata to wide format
# -----------------------------

# Convert locality metadata from long to wide format
nucleotide_summary_clean <- metadata_locality %>%
  
  # Keep identifiers, XML fields, and locality key-value columns
  select(
    row_id,
    taxa_accepted,
    uid,
    caption,
    title,
    genome,
    gene,
    CDS_product,
    ref_journal,
    ref_authors,
    metadata_field,
    metadata_value
  ) %>%
  
  # Spread metadata fields into columns
  pivot_wider(
    id_cols = c(
      row_id,
      taxa_accepted,
      uid,
      caption,
      title,
      genome,
      gene,
      CDS_product,
      ref_journal,
      ref_authors
    ),
    names_from = metadata_field,
    values_from = metadata_value,
    values_fn = ~ paste(unique(na.omit(.x)), collapse = "; ")
  ) %>%
  
  # Remove temporary row ID
  select(-row_id)


#Join important columns from taxonomy

nucleotide_summary_clean = nucleotide_summary_clean %>%
  left_join(
    glonaf_matched_f %>%
      select(taxa_accepted, accepted_plant_name_id, powo_id_accepted) %>%
      distinct(),
    by = 'taxa_accepted'
  )

# -----------------------------
# 9. Save final table
# -----------------------------

# Save cleaned locality table
write_tsv(
  nucleotide_summary_clean,
  paste0("./2_get_info_nucleotide/nucleotide_entrez_summary_clean_", chunk_id, ".tsv")
)