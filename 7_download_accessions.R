# ============================================================
# 7. EXTRACT MARKER-SPECIFIC FASTA SEQUENCES FROM NCBI/GENBANK
# ============================================================
#
# Input: outputs/6_parsing_accessions/seqs_to_download.tsv
# produced by the structured-annotation version of script 6.
#
# Core rationale:
#   * group by the canonical `marker`
#   * use GenBank FEATURES coordinates whenever available;
#   * respect join()/order()/complement() locations;
#   * fetch each accession only once (and reuse script 6's XML cache);
#   * use whole-record sequence only as a conservative fallback when the
#     accession itself is clearly a marker/amplicon record;
#   * never silently put a whole plastome/mitogenome/chromosome into a
#     marker-specific FASTA;
#   * record extraction provenance and failures for auditing.
#
# Output tree:
#   outputs/7_download_accessions/fasta_by_taxon/
#     Species_name/
#       Native|Introduced|Somewhere/
#         Nuclear|Plastidial|Mitochondrial/
#           marker.fasta
#
# Required packages:
#   dplyr, tidyr, stringr, stringi, purrr, readr, rentrez,
#   xml2, tibble, fs, Biostrings
# ============================================================

library(dplyr)
library(tidyr)
library(stringr)
library(stringi)
library(purrr)
library(readr)
library(rentrez)
library(xml2)
library(tibble)
library(fs)
library(Biostrings)

# ============================================================
# 1. SETTINGS
# ============================================================

input_file <- './outputs/6_parsing_accessions/seqs_to_download.tsv'
output_root <- './outputs/7_download_accessions'
fasta_root <- file.path(output_root, 'fasta_by_taxon')

# Reuse the XML cache created by script 6. This is important because GenBank
# XML contains both the nucleotide sequence and the feature annotation.
xml_cache_dir <- './outputs/6_parsing_accessions/genbank_xml_cache'

# If a cache entry is absent/corrupt, script 7 can refetch it.
max_fetch_attempts <- 5
request_delay_no_key <- 0.36
request_delay_with_key <- 0.12

# Whole-record fallback is intentionally conservative. A record longer than
# this is not treated as a marker amplicon unless its marker source is explicit
# and the title clearly describes the requested region.
max_whole_record_fallback_length <- 12000L

max_species_name_length <- 80
max_distribution_name_length <- 20
max_genome_name_length <- 30
max_marker_name_length <- 100

# Do NOT hard-code an NCBI API key in this script.
# In the shell, use for example:
#   export ENTREZ_KEY='...'
#   export ENTREZ_EMAIL='you@example.org'
if (nzchar(Sys.getenv('ENTREZ_KEY'))) {
  rentrez::set_entrez_key(Sys.getenv('ENTREZ_KEY'))
}

dir_create(output_root, recurse = TRUE)
dir_create(fasta_root, recurse = TRUE)
dir_create(xml_cache_dir, recurse = TRUE)

# ============================================================
# 2. PATH / FASTA HELPERS
# ============================================================

sanitize_path_component <- function(x, max_length = 100) {
  x <- as.character(x)
  x <- str_squish(x)
  
  x <- str_replace_all(
    x,
    c(
      'α'='alpha', 'β'='beta', 'γ'='gamma', 'δ'='delta', 'ε'='epsilon',
      'κ'='kappa', 'λ'='lambda', 'μ'='mu', 'ω'='omega',
      'Α'='Alpha', 'Β'='Beta', 'Γ'='Gamma', 'Δ'='Delta', 'Ε'='Epsilon',
      'Κ'='Kappa', 'Λ'='Lambda', 'Μ'='Mu', 'Ω'='Omega'
    )
  )
  
  x <- stringi::stri_trans_general(x, 'Latin-ASCII')
  x <- str_replace_all(x, '[/\\\\:*?"<>|]', '_')
  x <- str_replace_all(x, '[[:space:]]+', '_')
  x <- str_replace_all(x, '[^A-Za-z0-9._-]', '_')
  x <- str_replace_all(x, '_+', '_')
  x <- str_remove_all(x, '^[_\\.-]+|[_\\.-]+$')
  x[is.na(x) | x == ''] <- 'unknown'
  str_sub(x, 1, max_length)
}

make_marker_filename <- function(marker, max_length = 100) {
  paste0(sanitize_path_component(marker, max_length), '.fasta')
}

wrap_sequence <- function(x, width = 80L) {
  if (is.na(x) || !nzchar(x)) return(character(0))
  starts <- seq.int(1L, nchar(x), by = width)
  vapply(starts, function(i) substr(x, i, min(i + width - 1L, nchar(x))), character(1))
}

make_fasta_record <- function(header, sequence) {
  paste(c(paste0('>', header), wrap_sequence(sequence)), collapse = '\n')
}

# ============================================================
# 3. INPUT
# ============================================================

seqs_to_download <- read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(uid = as.character(uid))

required_columns <- c(
  'taxa_accepted', 'native', 'introduced', 'somewhere',
  'genome', 'marker', 'marker_source', 'region_type', 'uid',
  'feature_location', 'feature_start', 'feature_end', 'strand', 'title'
)

missing_columns <- setdiff(required_columns, names(seqs_to_download))
if (length(missing_columns) > 0) {
  stop(
    'seqs_to_download.tsv is missing columns required by the structured ',
    'script 7: ', paste(missing_columns, collapse = ', '),
    '\nMake sure script 6_parsing_accessions_structured.R was run first.'
  )
}

# Exactly one distribution category is required for directory assignment.
download_table_all <- seqs_to_download %>%
  mutate(
    taxa_accepted = str_squish(as.character(taxa_accepted)),
    genome = str_to_lower(str_squish(as.character(genome))),
    marker = str_squish(as.character(marker)),
    marker_source = as.character(marker_source),
    region_type = as.character(region_type),
    feature_location = na_if(str_squish(as.character(feature_location)), ''),
    title = as.character(title),
    native = coalesce(as.integer(native), 0L),
    introduced = coalesce(as.integer(introduced), 0L),
    somewhere = coalesce(as.integer(somewhere), 0L),
    n_positive_distribution_flags = native + introduced + somewhere,
    distribution_status = case_when(
      native == 1L & introduced == 0L & somewhere == 0L ~ 'Native',
      native == 0L & introduced == 1L & somewhere == 0L ~ 'Introduced',
      native == 0L & introduced == 0L & somewhere == 1L ~ 'Somewhere',
      TRUE ~ NA_character_
    )
  )

invalid_distribution_rows <- download_table_all %>%
  filter(is.na(distribution_status))

write_tsv(
  invalid_distribution_rows,
  file.path(output_root, 'invalid_distribution_assignments.tsv')
)

if (nrow(invalid_distribution_rows) > 0) {
  warning(
    nrow(invalid_distribution_rows),
    ' row(s) did not have exactly one distribution assignment and were excluded.'
  )
}

download_table <- download_table_all %>%
  filter(
    !is.na(distribution_status),
    !is.na(taxa_accepted), taxa_accepted != '',
    genome %in% c('nuclear', 'plastidial', 'mitochondrial'),
    !is.na(marker), marker != '',
    !is.na(uid), uid != ''
  ) %>%
  mutate(
    genome = recode(
      genome,
      nuclear='Nuclear',
      plastidial='Plastidial',
      mitochondrial='Mitochondrial'
    )
  ) %>%
  distinct(
    taxa_accepted, distribution_status, genome, marker, uid,
    .keep_all = TRUE
  )

write_tsv(download_table, file.path(output_root, 'extraction_manifest.tsv'))

# ============================================================
# 4. FETCH / READ ONE GENBANK XML RECORD
# ============================================================

xml_cache_file <- function(uid) file.path(xml_cache_dir, paste0(uid, '.xml'))

xml_has_record <- function(x) {
  !is.na(x) && nchar(x) > 50 && str_detect(x, '<GBSeq')
}

get_genbank_xml <- function(uid) {
  cache_file <- xml_cache_file(uid)
  
  if (file_exists(cache_file) && file_size(cache_file) > 50) {
    cached <- tryCatch(read_file(cache_file), error = function(e) NA_character_)
    if (xml_has_record(cached)) {
      return(list(xml = cached, source = 'script6_xml_cache', error = NA_character_))
    }
  }
  
  last_error <- NA_character_
  
  for (attempt in seq_len(max_fetch_attempts)) {
    ans <- tryCatch(
      rentrez::entrez_fetch(
        db='nuccore', id=uid, rettype='gb', retmode='xml'
      ),
      error = function(e) {
        last_error <<- conditionMessage(e)
        NA_character_
      }
    )
    
    if (xml_has_record(ans)) {
      try(write_file(ans, cache_file), silent = TRUE)
      Sys.sleep(if (nzchar(Sys.getenv('ENTREZ_KEY'))) {
        request_delay_with_key
      } else {
        request_delay_no_key
      })
      
      return(list(xml = ans, source = 'refetched_genbank_xml', error = NA_character_))
    }
    
    Sys.sleep(min(2^attempt, 30))
  }
  
  list(xml = NA_character_, source = 'fetch_failed', error = last_error)
}

# Parse the record-level sequence and useful identifiers from GenBank XML.
parse_genbank_record <- function(xml_text, uid) {
  if (!xml_has_record(xml_text)) {
    return(tibble(
      uid=uid, gb_accession=NA_character_, gb_accession_version=NA_character_,
      gb_definition=NA_character_, sequence=NA_character_, sequence_length=NA_integer_
    ))
  }
  
  doc <- read_xml(xml_text)
  gb <- xml_find_first(doc, './/GBSeq')
  
  get_text <- function(path) {
    node <- xml_find_first(gb, path)
    if (inherits(node, 'xml_missing')) NA_character_ else xml2::xml_text(node)
  }
  
  sequence <- str_to_upper(get_text('./GBSeq_sequence'))
  sequence <- str_replace_all(coalesce(sequence, ''), '[^ACGTNRYKMSWBDHV-]', '')
  if (sequence == '') sequence <- NA_character_
  
  tibble(
    uid = uid,
    gb_accession = get_text('./GBSeq_primary-accession'),
    gb_accession_version = get_text('./GBSeq_accession-version'),
    gb_definition = get_text('./GBSeq_definition'),
    sequence = sequence,
    sequence_length = ifelse(is.na(sequence), NA_integer_, nchar(sequence))
  )
}

# ============================================================
# 5. GENBANK LOCATION PARSER
# ============================================================
#
# Supports the common local GenBank location forms needed here:
#   123..456
#   complement(123..456)
#   join(123..200,300..456)
#   complement(join(123..200,300..456))
#   order(...)
#   <123..>456
#   123
#
# Remote locations (e.g. OTHERACC:1..100) are rejected because those bases
# are not necessarily contained in the current accession sequence.
# ============================================================

strip_outer_call <- function(x, fun) {
  prefix <- paste0(fun, '(')
  if (!str_starts(str_to_lower(x), fixed(prefix)) || !str_ends(x, fixed(')'))) return(NULL)
  substr(x, nchar(prefix) + 1L, nchar(x) - 1L)
}

split_top_level_commas <- function(x) {
  chars <- strsplit(x, '', fixed = TRUE)[[1]]
  depth <- 0L
  starts <- 1L
  pieces <- character(0)
  
  for (i in seq_along(chars)) {
    if (chars[i] == '(') depth <- depth + 1L
    if (chars[i] == ')') depth <- depth - 1L
    if (chars[i] == ',' && depth == 0L) {
      pieces <- c(pieces, substr(x, starts, i - 1L))
      starts <- i + 1L
    }
  }
  
  c(pieces, substr(x, starts, nchar(x))) %>% str_squish()
}

parse_location_node <- function(location) {
  x <- str_replace_all(str_squish(location), '\\s+', '')
  
  # Reject remote accession coordinates.
  if (str_detect(x, '(?:^|[,(])[A-Za-z][A-Za-z0-9_.]*:')) {
    return(list(ok=FALSE, type='remote', message='remote accession location'))
  }
  
  inner <- strip_outer_call(x, 'complement')
  if (!is.null(inner)) {
    child <- parse_location_node(inner)
    if (!isTRUE(child$ok)) return(child)
    return(list(ok=TRUE, type='complement', child=child))
  }
  
  for (fun in c('join', 'order')) {
    inner <- strip_outer_call(x, fun)
    if (!is.null(inner)) {
      parts <- split_top_level_commas(inner)
      children <- lapply(parts, parse_location_node)
      if (any(!vapply(children, function(z) isTRUE(z$ok), logical(1)))) {
        bad <- children[[which(!vapply(children, function(z) isTRUE(z$ok), logical(1)))[1]]]
        return(bad)
      }
      return(list(ok=TRUE, type=fun, children=children))
    }
  }
  
  clean <- str_replace_all(x, '[<>]', '')
  
  if (str_detect(clean, '^\\d+\\.\\.\\d+$')) {
    nums <- as.integer(str_split(clean, '\\.\\.', simplify=TRUE))
    return(list(ok=TRUE, type='interval', start=nums[1], end=nums[2]))
  }
  
  # Some INSDC locations use 123^124 for a between-base site. This does not
  # represent an extractable nucleotide interval, so fail explicitly.
  if (str_detect(clean, '^\\d+\\^\\d+$')) {
    return(list(ok=FALSE, type='between', message='between-base location'))
  }
  
  if (str_detect(clean, '^\\d+$')) {
    pos <- as.integer(clean)
    return(list(ok=TRUE, type='interval', start=pos, end=pos))
  }
  
  list(ok=FALSE, type='unsupported', message=paste('unsupported location:', location))
}

extract_location_node <- function(node, full_sequence) {
  if (!isTRUE(node$ok)) stop(node$message)
  
  if (node$type == 'interval') {
    if (is.na(node$start) || is.na(node$end) ||
        node$start < 1L || node$end < node$start ||
        node$end > nchar(full_sequence)) {
      stop('feature coordinates outside accession sequence')
    }
    return(substr(full_sequence, node$start, node$end))
  }
  
  if (node$type %in% c('join', 'order')) {
    return(paste0(vapply(
      node$children,
      extract_location_node,
      full_sequence=full_sequence,
      FUN.VALUE=character(1)
    ), collapse=''))
  }
  
  if (node$type == 'complement') {
    child_seq <- extract_location_node(node$child, full_sequence)
    return(as.character(reverseComplement(DNAString(child_seq))))
  }
  
  stop('unsupported parsed location node')
}

extract_by_feature_location <- function(full_sequence, feature_location) {
  if (is.na(full_sequence) || !nzchar(full_sequence)) {
    return(list(sequence=NA_character_, ok=FALSE, error='record sequence missing'))
  }
  if (is.na(feature_location) || !nzchar(feature_location)) {
    return(list(sequence=NA_character_, ok=FALSE, error='feature location missing'))
  }
  
  parsed <- parse_location_node(feature_location)
  if (!isTRUE(parsed$ok)) {
    return(list(sequence=NA_character_, ok=FALSE, error=parsed$message))
  }
  
  ans <- tryCatch(
    extract_location_node(parsed, full_sequence),
    error=function(e) structure(NA_character_, extraction_error=conditionMessage(e))
  )
  
  if (is.na(ans) || !nzchar(ans)) {
    return(list(
      sequence=NA_character_, ok=FALSE,
      error=coalesce(attr(ans, 'extraction_error'), 'empty extracted sequence')
    ))
  }
  
  list(sequence=ans, ok=TRUE, error=NA_character_)
}

# ============================================================
# 6. CONSERVATIVE WHOLE-RECORD FALLBACK
# ============================================================
#
# This fallback is for marker-specific amplicon/sequence records that have no
# usable feature coordinates. It is NOT used for generic genomic/organellar
# genome records.
# ============================================================

escape_regex <- function(x) {
  str_replace_all(x, '([.\\^$|()\\[\\]{}*+?\\\\-])', '\\\\\\1')
}

marker_is_explicit_in_title <- function(marker, title) {
  if (is.na(marker) || is.na(title)) return(FALSE)
  
  # Match the canonical marker literally after normalizing dash characters.
  m <- str_replace_all(marker, '[–—]', '-')
  t <- str_replace_all(title, '[–—]', '-')
  
  str_detect(t, fixed(m, ignore_case=TRUE))
}

whole_record_fallback_allowed <- function(row, record_length) {
  if (is.na(record_length) || record_length < 1L ||
      record_length > max_whole_record_fallback_length) return(FALSE)
  
  title <- coalesce(row$title, '')
  source <- coalesce(row$marker_source, '')
  
  # Never treat obvious whole-genome/assembly records as marker amplicons.
  if (str_detect(
    title,
    regex(
      'complete genome|complete sequence|whole genome|genome assembly|chromosome|plastome|mitogenome',
      ignore_case=TRUE
    )
  )) return(FALSE)
  
  # Title-derived markers are acceptable only when the marker is explicit in
  # the title. Input-field markers are also allowed if the title explicitly
  # contains the same marker.
  if (str_detect(source, '^title_') && marker_is_explicit_in_title(row$marker, title)) {
    return(TRUE)
  }
  
  if (source == 'input_gene_field' && marker_is_explicit_in_title(row$marker, title)) {
    return(TRUE)
  }
  
  FALSE
}

# ============================================================
# 7. EXTRACT ONE ACCESSION x MARKER ROW
# ============================================================

extract_one_row <- function(row, record) {
  has_feature_location <- !is.na(row$feature_location) && nzchar(row$feature_location)
  
  if (has_feature_location) {
    extracted <- extract_by_feature_location(record$sequence, row$feature_location)
    
    if (isTRUE(extracted$ok)) {
      return(tibble(
        extracted_sequence=extracted$sequence,
        extraction_method='feature_coordinates',
        extraction_status='success',
        extraction_error=NA_character_
      ))
    }
    
    # Important: feature-annotated rows do not silently fall back to the whole
    # record if their coordinates fail. That would contaminate marker FASTAs.
    return(tibble(
      extracted_sequence=NA_character_,
      extraction_method='feature_coordinates',
      extraction_status='failed',
      extraction_error=extracted$error
    ))
  }
  
  if (whole_record_fallback_allowed(row, record$sequence_length)) {
    return(tibble(
      extracted_sequence=record$sequence,
      extraction_method='whole_record_fallback',
      extraction_status='success',
      extraction_error=NA_character_
    ))
  }
  
  tibble(
    extracted_sequence=NA_character_,
    extraction_method='none',
    extraction_status='failed',
    extraction_error='no usable feature coordinates and whole-record fallback not justified'
  )
}

# ============================================================
# 8. FETCH EACH UNIQUE ACCESSION ONCE
# ============================================================

unique_uids <- sort(unique(download_table$uid))
message('Unique accessions to process: ', length(unique_uids))

record_list <- vector('list', length(unique_uids))
fetch_log_list <- vector('list', length(unique_uids))

for (i in seq_along(unique_uids)) {
  uid <- unique_uids[[i]]
  
  if (i == 1L || i %% 100L == 0L || i == length(unique_uids)) {
    message('Reading/fetching accession ', i, '/', length(unique_uids), ': ', uid)
  }
  
  fetched <- get_genbank_xml(uid)
  rec <- parse_genbank_record(fetched$xml, uid)
  
  record_list[[i]] <- rec
  fetch_log_list[[i]] <- tibble(
    uid=uid,
    record_source=fetched$source,
    fetch_error=fetched$error,
    record_sequence_length=rec$sequence_length
  )
}

records <- bind_rows(record_list)
fetch_log <- bind_rows(fetch_log_list)

write_tsv(fetch_log, file.path(output_root, 'accession_fetch_log.tsv'))

# ============================================================
# 9. EXTRACT ALL REQUESTED MARKERS
# ============================================================

rows_with_record <- download_table %>%
  left_join(records, by='uid')

extraction_results <- map_dfr(seq_len(nrow(rows_with_record)), function(i) {
  row <- rows_with_record[i, , drop=FALSE]
  
  if (is.na(row$sequence[[1]]) || !nzchar(row$sequence[[1]])) {
    result <- tibble(
      extracted_sequence=NA_character_,
      extraction_method='none',
      extraction_status='failed',
      extraction_error='GenBank record sequence unavailable'
    )
  } else {
    result <- extract_one_row(row, row)
  }
  
  bind_cols(row, result)
}) %>%
  mutate(
    extracted_length = ifelse(
      is.na(extracted_sequence),
      NA_integer_,
      nchar(extracted_sequence)
    ),
    # A transparent identifier: UID + marker + extraction mode.
    fasta_id = paste(uid, sanitize_path_component(marker, 50), sep='|'),
    fasta_header = paste0(
      fasta_id,
      '|taxon=', sanitize_path_component(taxa_accepted, 100),
      '|status=', distribution_status,
      '|genome=', genome,
      '|marker_source=', sanitize_path_component(marker_source, 60),
      '|extract=', extraction_method,
      ifelse(
        !is.na(coalesce(accession_version, gb_accession_version)),
        paste0('|acc=', coalesce(accession_version, gb_accession_version)),
        ''
      )
    )
  )

# Never keep duplicate accession x marker x distribution records.
extraction_results <- extraction_results %>%
  arrange(taxa_accepted, distribution_status, genome, marker, uid) %>%
  distinct(taxa_accepted, distribution_status, genome, marker, uid, .keep_all=TRUE)

# Full provenance table, minus the potentially huge sequence strings.
extraction_audit <- extraction_results %>%
  select(-sequence, -extracted_sequence)

write_tsv(extraction_audit, file.path(output_root, 'extraction_audit.tsv'))

failed_extractions <- extraction_audit %>%
  filter(extraction_status != 'success')
write_tsv(failed_extractions, file.path(output_root, 'failed_extractions.tsv'))

fallback_extractions <- extraction_audit %>%
  filter(extraction_method == 'whole_record_fallback')
write_tsv(fallback_extractions, file.path(output_root, 'whole_record_fallbacks.tsv'))

# ============================================================
# 10. WRITE MARKER-SPECIFIC FASTA FILES
# ============================================================

successful <- extraction_results %>%
  filter(
    extraction_status == 'success',
    !is.na(extracted_sequence),
    extracted_sequence != ''
  )

fasta_groups <- successful %>%
  group_by(taxa_accepted, distribution_status, genome, marker) %>%
  group_split()

fasta_log <- map_dfr(fasta_groups, function(g) {
  species <- g$taxa_accepted[[1]]
  distribution <- g$distribution_status[[1]]
  genome <- g$genome[[1]]
  marker <- g$marker[[1]]
  
  target_dir <- file.path(
    fasta_root,
    sanitize_path_component(species, max_species_name_length),
    sanitize_path_component(distribution, max_distribution_name_length),
    sanitize_path_component(genome, max_genome_name_length)
  )
  dir_create(target_dir, recurse=TRUE)
  
  target_file <- file.path(
    target_dir,
    make_marker_filename(marker, max_marker_name_length)
  )
  
  # Deterministic order makes reruns reproducible.
  g <- g %>% arrange(uid)
  
  fasta_text <- map2_chr(
    g$fasta_header,
    g$extracted_sequence,
    make_fasta_record
  )
  
  writeLines(paste(fasta_text, collapse='\n'), target_file, useBytes=TRUE)
  
  tibble(
    taxa_accepted=species,
    distribution_status=distribution,
    genome=genome,
    marker=marker,
    n_sequences=nrow(g),
    n_feature_coordinates=sum(g$extraction_method == 'feature_coordinates'),
    n_whole_record_fallback=sum(g$extraction_method == 'whole_record_fallback'),
    min_sequence_length=min(g$extracted_length, na.rm=TRUE),
    median_sequence_length=median(g$extracted_length, na.rm=TRUE),
    max_sequence_length=max(g$extracted_length, na.rm=TRUE),
    output_file=target_file
  )
})

write_tsv(fasta_log, file.path(output_root, 'fasta_manifest.tsv'))

# ============================================================
# 11. POST-EXTRACTION 5/5 AUDIT
# ============================================================
#
# Script 6 applied >=5 native and >=5 introduced BEFORE sequence extraction.
# Some records can fail extraction, so check the threshold again on sequences
# that were actually written. This does not delete files automatically; it
# identifies groups that fell below the desired threshold after QC/extraction.
# ============================================================

post_extraction_counts <- successful %>%
  group_by(taxa_accepted, genome, marker) %>%
  summarise(
    n_native=n_distinct(uid[distribution_status == 'Native']),
    n_introduced=n_distinct(uid[distribution_status == 'Introduced']),
    n_somewhere=n_distinct(uid[distribution_status == 'Somewhere']),
    .groups='drop'
  ) %>%
  mutate(
    passes_5_5_after_extraction = n_native >= 5L & n_introduced >= 5L
  )

write_tsv(
  post_extraction_counts,
  file.path(output_root, 'post_extraction_marker_counts.tsv')
)

write_tsv(
  post_extraction_counts %>% filter(!passes_5_5_after_extraction),
  file.path(output_root, 'markers_failing_5_5_after_extraction.tsv')
)

# ============================================================
# 12. SUMMARY
# ============================================================

summary_table <- extraction_results %>%
  count(extraction_status, extraction_method, name='n_accession_marker_rows') %>%
  arrange(extraction_status, extraction_method)

print(summary_table, n=Inf)

message('')
message('Marker-specific FASTA extraction finished.')
message('Requested accession-marker rows: ', nrow(extraction_results))
message('Successful extractions: ', sum(extraction_results$extraction_status == 'success'))
message('  from feature coordinates: ', sum(extraction_results$extraction_method == 'feature_coordinates'))
message('  whole-record fallbacks: ', sum(extraction_results$extraction_method == 'whole_record_fallback'))
message('Failed extractions: ', sum(extraction_results$extraction_status != 'success'))
message('FASTA files written: ', nrow(fasta_log))
message('Groups still passing >=5 native / >=5 introduced after extraction: ',
        sum(post_extraction_counts$passes_5_5_after_extraction, na.rm=TRUE),
        '/', nrow(post_extraction_counts))
message('Output directory: ', normalizePath(output_root, mustWork=FALSE))
