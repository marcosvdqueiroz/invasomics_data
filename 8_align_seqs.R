# ============================================================
# SCRIPT 8
# Sequence QC + joint Native/Introduced alignment
# ============================================================
#
# This script follows:
#   6_parsing_accessions.R
#   7_download_accessions.R
#
# Rationale:
#   * use the structured marker information produced by scripts 6/7;
#   * use script 7's fasta_manifest.tsv as the metadata source of truth;
#   * QC the actual extracted marker sequences before alignment;
#   * align Native + Introduced sequences TOGETHER for each
#       taxa_accepted x genome x marker combination;
#   * derive Native and Introduced subsets FROM THE SAME ALIGNMENT;
#   * retain a raw joint alignment and a site-coverage-filtered alignment;
#   * write explicit sequence-QC, alignment-QC, eligibility, and failure logs.
#
# IMPORTANT:
# Identical nucleotide sequences are NOT deduplicated. Different accessions
# carrying the same haplotype must remain separate for downstream diversity
# and haplotype-frequency analyses.
#
# Required packages:
#   mamba install -c conda-forge -c bioconda \
#       r-readr r-dplyr r-tidyr r-stringr r-purrr r-fs \
#       bioconductor-biostrings bioconductor-decipher mafft
#
# MAFFT is the default alignment backend. DECIPHER can be selected below.
# ============================================================


# ============================================================
# 1. PACKAGES
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(fs)
library(Biostrings)
library(DECIPHER)


# ============================================================
# 2. SETTINGS
# ============================================================

script7_root <- "./outputs/7_download_accessions"
fasta_manifest_file <- file.path(script7_root, "fasta_manifest.tsv")

output_root <- "./outputs/8_align_seqs"
aligned_root <- file.path(output_root, "aligned_by_taxon")

# Minimum number of accessions required in EACH distribution category.
min_seqs <- 5L

# Alignment backend:
#   "mafft"    -> recommended/default
#   "decipher" -> pure R fallback
alignment_method <- "mafft"

# MAFFT executable. Change this if MAFFT is not on PATH.
mafft_executable <- "mafft"

# ------------------------------------------------------------
# Sequence-QC settings
# ------------------------------------------------------------
#
# "Hard" failures are always excluded:
#   * unreadable FASTA records
#   * zero-length sequences
#
# The following criteria are QC FLAGS by default. They are written to the
# audit table but are NOT automatically excluded unless
# exclude_flagged_sequences <- TRUE.
#
# This is deliberate: GenBank markers often differ in amplicon boundaries,
# so a conservative audit-first approach avoids silently discarding data.

minimum_absolute_length <- 50L
minimum_relative_length <- 0.50
maximum_ambiguous_fraction <- 0.10

exclude_flagged_sequences <- FALSE

# ------------------------------------------------------------
# Alignment-site filtering
# ------------------------------------------------------------
#
# Keep a raw joint alignment AND a coverage-filtered alignment.
# A site is retained in the filtered alignment when at least this fraction
# of sequences contains an unambiguous A/C/G/T nucleotide at that site.

minimum_site_coverage <- 0.80

# Path component limits.
max_species_name_length <- 80
max_genome_name_length <- 30
max_marker_name_length <- 100


dir_create(output_root, recurse = TRUE)
dir_create(aligned_root, recurse = TRUE)


# ============================================================
# 3. HELPERS
# ============================================================

sanitize_path_component <- function(x, max_length = 100) {

  x <- as.character(x)
  x <- str_squish(x)

  x <- stringi::stri_trans_general(x, "Latin-ASCII")

  x <- str_replace_all(
    x,
    "[/\\\\:*?\"<>|]",
    "_"
  )

  x <- str_replace_all(x, "[[:space:]]+", "_")
  x <- str_replace_all(x, "[^A-Za-z0-9._-]", "_")
  x <- str_replace_all(x, "_+", "_")
  x <- str_remove_all(x, "^[_\\.-]+|[_\\.-]+$")

  x[is.na(x) | x == ""] <- "unknown"

  str_sub(x, 1, max_length)
}


resolve_input_path <- function(path, base_dir = ".") {

  path <- as.character(path)

  if (is.na(path) || path == "") {
    return(NA_character_)
  }

  # Already-valid path.
  if (file_exists(path)) {
    return(path_abs(path))
  }

  # Try relative to supplied base directory.
  candidate <- file.path(base_dir, path)

  if (file_exists(candidate)) {
    return(path_abs(candidate))
  }

  # Return original path for transparent downstream failure logging.
  path
}


extract_uid_from_header <- function(x) {

  # Script 7 writes:
  # >UID|marker|taxon=...|status=...
  str_extract(x, "^[^|]+")
}


read_fasta_safely <- function(path) {

  tryCatch(
    readDNAStringSet(path),
    error = function(e) {
      message(
        "Could not read FASTA: ",
        path,
        " -- ",
        conditionMessage(e)
      )
      NULL
    }
  )
}


# ============================================================
# 4. READ SCRIPT-7 FASTA MANIFEST
# ============================================================

if (!file_exists(fasta_manifest_file)) {
  stop(
    "Could not find script-7 FASTA manifest:\n",
    fasta_manifest_file,
    "\nRun the structured script 7 first."
  )
}

fasta_manifest <- read_tsv(
  fasta_manifest_file,
  show_col_types = FALSE
)

required_manifest_columns <- c(
  "taxa_accepted",
  "distribution_status",
  "genome",
  "marker",
  "n_sequences",
  "output_file"
)

missing_manifest_columns <- setdiff(
  required_manifest_columns,
  names(fasta_manifest)
)

if (length(missing_manifest_columns) > 0) {
  stop(
    "fasta_manifest.tsv is missing required columns: ",
    paste(missing_manifest_columns, collapse = ", ")
  )
}


fasta_manifest <- fasta_manifest %>%

  mutate(
    taxa_accepted = str_squish(as.character(taxa_accepted)),
    distribution_status = as.character(distribution_status),
    genome = as.character(genome),
    marker = str_squish(as.character(marker)),

    fasta_path = map_chr(
      output_file,
      resolve_input_path,
      base_dir = script7_root
    )
  ) %>%

  filter(
    distribution_status %in% c("Native", "Introduced"),
    !is.na(taxa_accepted),
    taxa_accepted != "",
    !is.na(genome),
    genome != "",
    !is.na(marker),
    marker != ""
  )


# ============================================================
# 5. CHECK MANIFEST FILES
# ============================================================

missing_fasta_files <- fasta_manifest %>%
  filter(
    is.na(fasta_path) |
      !file_exists(fasta_path)
  )

write_tsv(
  missing_fasta_files,
  file.path(output_root, "missing_fasta_files.tsv")
)

if (nrow(missing_fasta_files) > 0) {
  warning(
    nrow(missing_fasta_files),
    " FASTA file(s) listed by script 7 could not be found. ",
    "They are listed in missing_fasta_files.tsv."
  )
}

fasta_manifest_valid <- fasta_manifest %>%
  filter(
    !is.na(fasta_path),
    file_exists(fasta_path)
  )


# ============================================================
# 6. READ EVERY FASTA AND BUILD ONE SEQUENCE-LEVEL TABLE
# ============================================================
#
# We intentionally use the script-7 manifest for biological metadata rather
# than reconstructing species/genome/marker from directory names.
# ============================================================

read_manifest_fasta <- function(row) {

  seqs <- read_fasta_safely(row$fasta_path)

  if (is.null(seqs)) {

    return(
      tibble(
        taxa_accepted = row$taxa_accepted,
        distribution_status = row$distribution_status,
        genome = row$genome,
        marker = row$marker,
        fasta_path = row$fasta_path,
        fasta_header = NA_character_,
        uid = NA_character_,
        sequence = NA_character_,
        readable = FALSE
      )
    )
  }

  tibble(
    taxa_accepted = row$taxa_accepted,
    distribution_status = row$distribution_status,
    genome = row$genome,
    marker = row$marker,
    fasta_path = row$fasta_path,
    fasta_header = names(seqs),
    uid = extract_uid_from_header(names(seqs)),
    sequence = as.character(seqs),
    readable = TRUE
  )
}


sequence_table <- map_dfr(
  seq_len(nrow(fasta_manifest_valid)),
  function(i) {
    read_manifest_fasta(fasta_manifest_valid[i, , drop = FALSE])
  }
)


# ============================================================
# 7. BASIC SEQUENCE QC
# ============================================================

sequence_table <- sequence_table %>%

  mutate(
    sequence = str_to_upper(sequence),

    sequence_length = if_else(
      readable & !is.na(sequence),
      nchar(sequence),
      NA_integer_
    ),

    n_acgt = if_else(
      readable & !is.na(sequence),
      str_count(sequence, "[ACGT]"),
      NA_integer_
    ),

    ambiguous_fraction = case_when(
      !readable ~ NA_real_,
      is.na(sequence_length) ~ NA_real_,
      sequence_length == 0L ~ NA_real_,
      TRUE ~ 1 - (n_acgt / sequence_length)
    ),

    hard_qc_pass =
      readable &
      !is.na(sequence) &
      !is.na(sequence_length) &
      sequence_length > 0L
  )


# Median length is calculated across Native + Introduced together for the
# same species x genome x marker.

sequence_table <- sequence_table %>%

  group_by(
    taxa_accepted,
    genome,
    marker
  ) %>%

  mutate(
    marker_median_length = median(
      sequence_length[hard_qc_pass],
      na.rm = TRUE
    ),

    relative_length = case_when(
      !hard_qc_pass ~ NA_real_,
      is.na(marker_median_length) ~ NA_real_,
      marker_median_length <= 0 ~ NA_real_,
      TRUE ~ sequence_length / marker_median_length
    )
  ) %>%

  ungroup()


# Protect against median(numeric(0), na.rm=TRUE) -> NA/warning-like edge cases.
sequence_table <- sequence_table %>%
  mutate(
    marker_median_length = if_else(
      is.finite(marker_median_length),
      marker_median_length,
      NA_real_
    )
  )


sequence_table <- sequence_table %>%

  mutate(
    flag_short_absolute =
      hard_qc_pass &
      sequence_length < minimum_absolute_length,

    flag_short_relative =
      hard_qc_pass &
      !is.na(relative_length) &
      relative_length < minimum_relative_length,

    flag_high_ambiguity =
      hard_qc_pass &
      !is.na(ambiguous_fraction) &
      ambiguous_fraction > maximum_ambiguous_fraction,

    n_qc_flags =
      coalesce(as.integer(flag_short_absolute), 0L) +
      coalesce(as.integer(flag_short_relative), 0L) +
      coalesce(as.integer(flag_high_ambiguity), 0L),

    qc_flagged = n_qc_flags > 0L,

    qc_pass_for_alignment = case_when(
      !hard_qc_pass ~ FALSE,
      exclude_flagged_sequences ~ !qc_flagged,
      TRUE ~ TRUE
    )
  )


# ============================================================
# 8. DUPLICATE-UID AUDIT
# ============================================================
#
# Duplicate nucleotide sequences are allowed and retained.
# Duplicate accession UIDs within the SAME biological group are not expected.
# ============================================================

duplicate_uid_audit <- sequence_table %>%

  filter(
    hard_qc_pass,
    !is.na(uid),
    uid != ""
  ) %>%

  count(
    taxa_accepted,
    distribution_status,
    genome,
    marker,
    uid,
    name = "n_rows"
  ) %>%

  filter(n_rows > 1L)


write_tsv(
  duplicate_uid_audit,
  file.path(output_root, "duplicate_uid_audit.tsv")
)


# Keep one copy of an accidentally repeated UID inside a group.
# This removes duplicated records, NOT repeated haplotypes from different UIDs.

sequence_table <- sequence_table %>%

  arrange(
    taxa_accepted,
    distribution_status,
    genome,
    marker,
    uid
  ) %>%

  distinct(
    taxa_accepted,
    distribution_status,
    genome,
    marker,
    uid,
    .keep_all = TRUE
  )


write_tsv(
  sequence_table %>%
    select(-sequence),
  file.path(output_root, "sequence_qc.tsv")
)


# ============================================================
# 9. FINAL 5/5 ELIGIBILITY AFTER SEQUENCE QC
# ============================================================

qc_counts <- sequence_table %>%

  filter(qc_pass_for_alignment) %>%

  group_by(
    taxa_accepted,
    genome,
    marker
  ) %>%

  summarise(
    n_native = n_distinct(
      uid[distribution_status == "Native"]
    ),

    n_introduced = n_distinct(
      uid[distribution_status == "Introduced"]
    ),

    .groups = "drop"
  ) %>%

  mutate(
    passes_final_5_5 =
      n_native >= min_seqs &
      n_introduced >= min_seqs
  )


write_tsv(
  qc_counts,
  file.path(output_root, "post_qc_marker_counts.tsv")
)

write_tsv(
  qc_counts %>%
    filter(!passes_final_5_5),
  file.path(output_root, "markers_failing_5_5_after_qc.tsv")
)


eligible_markers <- qc_counts %>%
  filter(passes_final_5_5)

write_tsv(
  eligible_markers,
  file.path(output_root, "alignment_eligible_markers.tsv")
)


message(
  nrow(eligible_markers),
  " species × genome × marker combinations pass the final ",
  min_seqs,
  "/",
  min_seqs,
  " criterion."
)


# ============================================================
# 10. ALIGNMENT HELPERS
# ============================================================

make_joint_sequences <- function(data) {

  data <- data %>%
    filter(qc_pass_for_alignment) %>%
    arrange(distribution_status, uid)

  seqs <- DNAStringSet(data$sequence)

  # Names from script 7 already contain provenance and distribution status.
  # Make them unique defensively.
  names(seqs) <- make.unique(data$fasta_header, sep = "__dup")

  list(
    seqs = seqs,
    metadata = data %>%
      mutate(alignment_name = names(seqs))
  )
}


align_with_decipher <- function(seqs) {

  tryCatch(
    DECIPHER::AlignSeqs(
      seqs,
      verbose = FALSE
    ),
    error = function(e) {
      attr(e, "alignment_backend") <- "decipher"
      stop(e)
    }
  )
}


align_with_mafft <- function(seqs) {

  mafft_path <- Sys.which(mafft_executable)

  if (mafft_path == "") {
    stop(
      "MAFFT executable was not found on PATH: ",
      mafft_executable
    )
  }

  temp_in <- tempfile(fileext = ".fasta")
  temp_out <- tempfile(fileext = ".fasta")
  temp_err <- tempfile(fileext = ".log")

  on.exit(
    unlink(
      c(temp_in, temp_out, temp_err),
      force = TRUE
    ),
    add = TRUE
  )

  writeXStringSet(
    seqs,
    filepath = temp_in,
    format = "fasta"
  )

  status <- system2(
    command = mafft_path,
    args = c(
      "--auto",
      "--quiet",
      shQuote(temp_in)
    ),
    stdout = temp_out,
    stderr = temp_err
  )

  if (!identical(status, 0L)) {

    error_text <- if (file_exists(temp_err)) {
      paste(readLines(temp_err, warn = FALSE), collapse = "\n")
    } else {
      ""
    }

    stop(
      "MAFFT failed with exit status ",
      status,
      if (nzchar(error_text)) paste0("\n", error_text) else ""
    )
  }

  aligned <- readDNAStringSet(temp_out)

  # Preserve input order rather than MAFFT output ordering.
  aligned <- aligned[names(seqs)]

  aligned
}


align_joint_sequences <- function(seqs) {

  if (length(seqs) < 2L) {
    stop("At least two sequences are required for alignment.")
  }

  if (alignment_method == "mafft") {
    return(align_with_mafft(seqs))
  }

  if (alignment_method == "decipher") {
    return(align_with_decipher(seqs))
  }

  stop(
    "Unknown alignment_method: ",
    alignment_method,
    ". Use 'mafft' or 'decipher'."
  )
}


# ============================================================
# 11. ALIGNMENT QC HELPERS
# ============================================================

alignment_matrix <- function(aligned) {

  if (length(aligned) == 0L) {
    return(matrix(character(0), nrow = 0, ncol = 0))
  }

  do.call(
    rbind,
    strsplit(
      as.character(aligned),
      split = "",
      fixed = TRUE
    )
  )
}


calculate_alignment_qc <- function(aligned) {

  mat <- alignment_matrix(aligned)

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    return(
      tibble(
        n_sequences = length(aligned),
        alignment_length = 0L,
        gap_fraction = NA_real_,
        ambiguous_fraction = NA_real_,
        mean_site_coverage = NA_real_,
        n_variable_sites = 0L,
        n_parsimony_informative_sites = 0L
      )
    )
  }

  mat <- toupper(mat)

  canonical <- mat %in% c("A", "C", "G", "T")
  gaps <- mat == "-"
  ambiguous <- !(canonical | gaps)

  site_coverage <- colSums(
    matrix(
      canonical,
      nrow = nrow(mat),
      ncol = ncol(mat)
    )
  ) / nrow(mat)

  site_stats <- map_dfr(
    seq_len(ncol(mat)),
    function(j) {

      x <- mat[, j]
      x <- x[x %in% c("A", "C", "G", "T")]

      counts <- table(x)

      tibble(
        variable = length(counts) >= 2L,
        parsimony_informative =
          sum(counts >= 2L) >= 2L
      )
    }
  )

  tibble(
    n_sequences = nrow(mat),
    alignment_length = ncol(mat),

    gap_fraction = mean(gaps),

    ambiguous_fraction = mean(ambiguous),

    mean_site_coverage = mean(site_coverage),

    n_variable_sites = sum(site_stats$variable),

    n_parsimony_informative_sites =
      sum(site_stats$parsimony_informative)
  )
}


coverage_filter_alignment <- function(
  aligned,
  minimum_coverage = 0.80
) {

  mat <- alignment_matrix(aligned)

  if (nrow(mat) == 0L || ncol(mat) == 0L) {
    return(NULL)
  }

  mat <- toupper(mat)

  canonical <- mat %in% c("A", "C", "G", "T")

  canonical <- matrix(
    canonical,
    nrow = nrow(mat),
    ncol = ncol(mat)
  )

  site_coverage <- colSums(canonical) / nrow(mat)

  keep_sites <- site_coverage >= minimum_coverage

  if (!any(keep_sites)) {
    return(NULL)
  }

  filtered_strings <- apply(
    mat[, keep_sites, drop = FALSE],
    1,
    paste0,
    collapse = ""
  )

  filtered <- DNAStringSet(filtered_strings)
  names(filtered) <- rownames(mat)

  filtered
}


# ============================================================
# 12. OUTPUT HELPERS
# ============================================================

write_alignment_set <- function(
  aligned,
  metadata,
  marker_dir,
  marker_file_stem,
  suffix
) {

  all_file <- file.path(
    marker_dir,
    paste0(
      marker_file_stem,
      "_",
      suffix,
      "_all.fasta"
    )
  )

  native_file <- file.path(
    marker_dir,
    paste0(
      marker_file_stem,
      "_",
      suffix,
      "_native.fasta"
    )
  )

  introduced_file <- file.path(
    marker_dir,
    paste0(
      marker_file_stem,
      "_",
      suffix,
      "_introduced.fasta"
    )
  )

  writeXStringSet(
    aligned,
    filepath = all_file,
    format = "fasta"
  )

  native_names <- metadata %>%
    filter(distribution_status == "Native") %>%
    pull(alignment_name)

  introduced_names <- metadata %>%
    filter(distribution_status == "Introduced") %>%
    pull(alignment_name)

  native_aligned <- aligned[
    names(aligned) %in% native_names
  ]

  introduced_aligned <- aligned[
    names(aligned) %in% introduced_names
  ]

  writeXStringSet(
    native_aligned,
    filepath = native_file,
    format = "fasta"
  )

  writeXStringSet(
    introduced_aligned,
    filepath = introduced_file,
    format = "fasta"
  )

  list(
    all = all_file,
    native = native_file,
    introduced = introduced_file
  )
}


# ============================================================
# 13. ALIGN ONE SPECIES × GENOME × MARKER
# ============================================================

align_one_marker <- function(
  species,
  genome,
  marker
) {

  marker_data <- sequence_table %>%

    semi_join(
      eligible_markers %>%
        filter(
          .data$taxa_accepted == .env$species,
          .data$genome == .env$genome,
          .data$marker == .env$marker
        ),
      by = c(
        "taxa_accepted",
        "genome",
        "marker"
      )
    ) %>%

    filter(
      taxa_accepted == species,
      .data$genome == genome,
      .data$marker == marker,
      distribution_status %in% c(
        "Native",
        "Introduced"
      ),
      qc_pass_for_alignment
    )

  n_native <- n_distinct(
    marker_data$uid[
      marker_data$distribution_status == "Native"
    ]
  )

  n_introduced <- n_distinct(
    marker_data$uid[
      marker_data$distribution_status == "Introduced"
    ]
  )

  species_dir <- sanitize_path_component(
    species,
    max_species_name_length
  )

  genome_dir <- sanitize_path_component(
    genome,
    max_genome_name_length
  )

  marker_dir_name <- sanitize_path_component(
    marker,
    max_marker_name_length
  )

  marker_dir <- file.path(
    aligned_root,
    species_dir,
    genome_dir,
    marker_dir_name
  )

  dir_create(marker_dir, recurse = TRUE)

  marker_stem <- marker_dir_name

  joint <- make_joint_sequences(marker_data)

  alignment_error <- NA_character_

  aligned_raw <- tryCatch(
    align_joint_sequences(joint$seqs),
    error = function(e) {
      alignment_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(aligned_raw)) {

    return(
      tibble(
        taxa_accepted = species,
        genome = genome,
        marker = marker,
        n_native = n_native,
        n_introduced = n_introduced,
        alignment_method = alignment_method,
        alignment_status = "failed",
        alignment_error = alignment_error,
        raw_all_file = NA_character_,
        raw_native_file = NA_character_,
        raw_introduced_file = NA_character_,
        filtered_all_file = NA_character_,
        filtered_native_file = NA_character_,
        filtered_introduced_file = NA_character_,
        raw_alignment_length = NA_integer_,
        filtered_alignment_length = NA_integer_,
        n_sites_retained = NA_integer_,
        fraction_sites_retained = NA_real_
      )
    )
  }

  # ----------------------------------------------------------
  # Raw joint alignment
  # ----------------------------------------------------------

  raw_files <- write_alignment_set(
    aligned = aligned_raw,
    metadata = joint$metadata,
    marker_dir = marker_dir,
    marker_file_stem = marker_stem,
    suffix = "raw_aligned"
  )

  raw_qc <- calculate_alignment_qc(
    aligned_raw
  )


  # ----------------------------------------------------------
  # Coverage-filtered alignment
  # ----------------------------------------------------------

  aligned_filtered <- coverage_filter_alignment(
    aligned_raw,
    minimum_coverage = minimum_site_coverage
  )

  if (is.null(aligned_filtered)) {

    filtered_files <- list(
      all = NA_character_,
      native = NA_character_,
      introduced = NA_character_
    )

    filtered_qc <- tibble(
      n_sequences = length(aligned_raw),
      alignment_length = 0L,
      gap_fraction = NA_real_,
      ambiguous_fraction = NA_real_,
      mean_site_coverage = NA_real_,
      n_variable_sites = 0L,
      n_parsimony_informative_sites = 0L
    )

    filtering_status <- "no sites passed coverage filter"

  } else {

    filtered_files <- write_alignment_set(
      aligned = aligned_filtered,
      metadata = joint$metadata,
      marker_dir = marker_dir,
      marker_file_stem = marker_stem,
      suffix = paste0(
        "coverage",
        round(minimum_site_coverage * 100),
        "_aligned"
      )
    )

    filtered_qc <- calculate_alignment_qc(
      aligned_filtered
    )

    filtering_status <- "complete"
  }


  # ----------------------------------------------------------
  # Save accession metadata used in this alignment
  # ----------------------------------------------------------

  marker_metadata_file <- file.path(
    marker_dir,
    paste0(
      marker_stem,
      "_sequence_metadata.tsv"
    )
  )

  joint$metadata %>%
    select(
      taxa_accepted,
      distribution_status,
      genome,
      marker,
      uid,
      fasta_header,
      alignment_name,
      sequence_length,
      marker_median_length,
      relative_length,
      ambiguous_fraction,
      flag_short_absolute,
      flag_short_relative,
      flag_high_ambiguity,
      qc_flagged,
      qc_pass_for_alignment
    ) %>%
    write_tsv(marker_metadata_file)


  # ----------------------------------------------------------
  # Alignment-level QC
  # ----------------------------------------------------------

  alignment_qc_file <- file.path(
    marker_dir,
    paste0(
      marker_stem,
      "_alignment_qc.tsv"
    )
  )

  bind_rows(
    raw_qc %>%
      mutate(
        alignment_version = "raw"
      ),

    filtered_qc %>%
      mutate(
        alignment_version = paste0(
          "coverage_",
          minimum_site_coverage
        )
      )
  ) %>%

    mutate(
      taxa_accepted = species,
      genome = genome,
      marker = marker,
      .before = 1
    ) %>%

    write_tsv(alignment_qc_file)


  raw_length <- raw_qc$alignment_length[[1]]
  filtered_length <- filtered_qc$alignment_length[[1]]

  tibble(
    taxa_accepted = species,
    genome = genome,
    marker = marker,

    n_native = n_native,
    n_introduced = n_introduced,

    alignment_method = alignment_method,
    alignment_status = "complete",
    alignment_error = NA_character_,

    filtering_status = filtering_status,

    raw_all_file = raw_files$all,
    raw_native_file = raw_files$native,
    raw_introduced_file = raw_files$introduced,

    filtered_all_file = filtered_files$all,
    filtered_native_file = filtered_files$native,
    filtered_introduced_file = filtered_files$introduced,

    sequence_metadata_file = marker_metadata_file,
    alignment_qc_file = alignment_qc_file,

    raw_alignment_length = raw_length,
    filtered_alignment_length = filtered_length,

    n_sites_retained = filtered_length,

    fraction_sites_retained = if_else(
      !is.na(raw_length) &
        raw_length > 0L &
        !is.na(filtered_length),
      filtered_length / raw_length,
      NA_real_
    )
  )
}


# ============================================================
# 14. RUN ALL ELIGIBLE ALIGNMENTS
# ============================================================

if (nrow(eligible_markers) == 0L) {

  warning(
    "No species × genome × marker combinations pass the final ",
    min_seqs,
    "/",
    min_seqs,
    " criterion."
  )

  alignment_manifest <- tibble(
    taxa_accepted = character(),
    genome = character(),
    marker = character(),
    n_native = integer(),
    n_introduced = integer(),
    alignment_method = character(),
    alignment_status = character(),
    alignment_error = character(),
    filtering_status = character(),
    raw_all_file = character(),
    raw_native_file = character(),
    raw_introduced_file = character(),
    filtered_all_file = character(),
    filtered_native_file = character(),
    filtered_introduced_file = character(),
    sequence_metadata_file = character(),
    alignment_qc_file = character(),
    raw_alignment_length = integer(),
    filtered_alignment_length = integer(),
    n_sites_retained = integer(),
    fraction_sites_retained = double()
  )

} else {

  alignment_manifest <- pmap_dfr(
    eligible_markers %>%
      select(
        taxa_accepted,
        genome,
        marker
      ),
    function(taxa_accepted, genome, marker) {
      align_one_marker(
        species = taxa_accepted,
        genome = genome,
        marker = marker
      )
    }
  )
}


# ============================================================
# 15. SAVE MASTER ALIGNMENT MANIFEST
# ============================================================

write_tsv(
  alignment_manifest,
  file.path(
    output_root,
    "alignment_manifest.tsv"
  )
)


alignment_failures <- alignment_manifest %>%
  filter(
    alignment_status != "complete"
  )

write_tsv(
  alignment_failures,
  file.path(
    output_root,
    "alignment_failures.tsv"
  )
)


# ============================================================
# 16. MASTER ALIGNMENT-QC TABLE
# ============================================================

alignment_qc_master <- map_dfr(
  alignment_manifest$alignment_qc_file[
    !is.na(alignment_manifest$alignment_qc_file) &
      file_exists(alignment_manifest$alignment_qc_file)
  ],
  ~ read_tsv(
    .x,
    show_col_types = FALSE
  )
)

write_tsv(
  alignment_qc_master,
  file.path(
    output_root,
    "alignment_qc_summary.tsv"
  )
)


# ============================================================
# 17. SUMMARY
# ============================================================

message("")
message("Alignment stage finished.")
message("")

message(
  "Sequence rows read: ",
  nrow(sequence_table)
)

message(
  "QC-flagged sequence rows: ",
  sum(
    sequence_table$qc_flagged,
    na.rm = TRUE
  )
)

message(
  "QC flags excluded from alignment: ",
  ifelse(
    exclude_flagged_sequences,
    "YES",
    "NO"
  )
)

message(
  "Eligible species × genome × marker groups: ",
  nrow(eligible_markers)
)

if (nrow(alignment_manifest) > 0L) {

  message(
    "Successful joint alignments: ",
    sum(
      alignment_manifest$alignment_status == "complete",
      na.rm = TRUE
    )
  )

  message(
    "Failed joint alignments: ",
    sum(
      alignment_manifest$alignment_status != "complete",
      na.rm = TRUE
    )
  )
}

message(
  "Alignment backend: ",
  alignment_method
)

message(
  "Coverage-filter threshold: ",
  minimum_site_coverage
)

message(
  "Output directory: ",
  path_abs(output_root)
)
