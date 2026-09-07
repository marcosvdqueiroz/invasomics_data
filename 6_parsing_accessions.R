library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(rentrez)
library(xml2)
library(tibble)
library(fs)

# ============================================================
# 6. PARSE GENBANK ACCESSIONS USING STRUCTURED ANNOTATION FIRST
# ============================================================
#
# Evidence hierarchy used here:
#   1. GenBank FEATURES + /organelle qualifier
#   2. Existing structured Entrez genome/gene fields in the input table
#   3. Explicit information in the GenBank title
#   4. Conservative marker-based genome inference
#   5. unresolved
#
# Important design changes:
#   * one accession can yield MANY marker rows;
#   * 'marker' replaces the old overloaded 'gene' grouping concept;
#   * marker provenance and genome provenance are retained;
#   * feature coordinates/location are retained for later subsequence extraction;
#   * title parsing is now a rescue layer, not the primary classifier;
#   * 5/5 thresholds count DISTINCT accessions per marker.
#
# This script deliberately does NOT force a flanking gene into an intergenic
# spacer family. For example, psbA remains psbA, while trnH-psbA is retained
# as an intergenic_spacer only when it is explicitly annotated/described as one.
# ============================================================

input_file <- './outputs/5_distribution/final_seqs_to_parse.tsv'
out_dir    <- './outputs/6_parsing_accessions'
cache_dir  <- file.path(out_dir, 'genbank_xml_cache')

dir_create(out_dir, recurse = TRUE)
dir_create(cache_dir, recurse = TRUE)

# NCBI recommends identifying requests when possible.
# Set these in your shell if desired:
#   export ENTREZ_EMAIL='xx'
#   export ENTREZ_KEY='ENTREZ_KEY="yyy"'
# rentrez automatically uses ENTREZ_KEY when available.

final_seqs_to_parse <- read_tsv(input_file, show_col_types = FALSE) %>%
  mutate(
    uid = as.character(uid),
    title = as.character(title),
    caption = as.character(caption),
    input_genome = as.character(genome),
    input_gene = as.character(gene)
  )

# ------------------------------------------------------------
# Remove whole-genome / assembly records that are treated elsewhere
# ------------------------------------------------------------
final_seqs_to_parse <- final_seqs_to_parse %>%
  filter(
    is.na(input_genome) | input_genome != 'chromosome',
    !str_detect(coalesce(title, ''), regex('\\bchromosome\\b', ignore_case = TRUE)),
    !str_detect(coalesce(title, ''), regex('complete genome', ignore_case = TRUE)),
    #!str_detect(coalesce(title, ''), regex('chloroplast, complete sequence', ignore_case = TRUE)),
    #!str_detect(coalesce(title, ''), regex('chloroplast, partial genome', ignore_case = TRUE)),
    #!str_detect(coalesce(title, ''), regex('plastid, partial genome', ignore_case = TRUE)),
    !(is.na(input_gene) & str_detect(coalesce(title, ''), regex('genome assembly', ignore_case = TRUE))),
    !str_detect(coalesce(title, ''), regex('hypothetical', ignore_case = TRUE)),
    !str_detect(coalesce(title, ''), regex('microsatellite', ignore_case = TRUE))
  )

# ============================================================
# 1. NORMALIZATION HELPERS
# ============================================================

normalize_genome_value <- function(x) {
  x <- str_to_lower(str_squish(as.character(x)))
  x[x == ''] <- NA_character_

  case_when(
    is.na(x) ~ NA_character_,
    str_detect(x, 'chloroplast|plastid|plastome') ~ 'plastidial',
    str_detect(x, 'mitochond') ~ 'mitochondrial',
    str_detect(x, 'nucleomorph|nuclear|genomic') ~ 'nuclear',
    TRUE ~ NA_character_
  )
}

# Canonical marker names. This function standardizes strings that have already
# been identified as marker candidates; it does NOT decide whether they are
# biologically valid marker candidates.
standardize_marker_name <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x == ''] <- NA_character_
  x <- str_replace_all(x, '[–—]', '-')
  x <- str_replace_all(x, '\\s*-\\s*', '-')

  # Written tRNA products -> symbols.
  trna_names <- c(
    'trna-ala'='trnA', 'trna-arg'='trnR', 'trna-asn'='trnN',
    'trna-asp'='trnD', 'trna-cys'='trnC', 'trna-gln'='trnQ',
    'trna-glu'='trnE', 'trna-gly'='trnG', 'trna-his'='trnH',
    'trna-ile'='trnI', 'trna-leu'='trnL', 'trna-lys'='trnK',
    'trna-met'='trnM', 'trna-phe'='trnF', 'trna-pro'='trnP',
    'trna-ser'='trnS', 'trna-thr'='trnT', 'trna-trp'='trnW',
    'trna-tyr'='trnY', 'trna-val'='trnV', 'trna-fmet'='trnfM'
  )

  for (a in names(trna_names)) {
    x <- str_replace_all(x, regex(paste0('\\b', a, '\\b'), ignore_case = TRUE), trna_names[[a]])
  }

  # Common aliases and spelling variants.
  aliases <- c(
    'rbcl'='rbcL', 'matk'='matK', 'nadhf'='ndhF',
    'coi'='cox1', 'coxi'='cox1', 'cytb'='cob',
    'rpoc1'='rpoC1', 'rpoc2'='rpoC2',
    'matr'='matR', 'mttb'='mttB',
    'ccmb'='ccmB', 'ccmc'='ccmC', 'ccmfc'='ccmFc', 'ccmfn'='ccmFn',
    'accd'='accD', 'ccsa'='ccsA', 'cema'='cemA', 'infa'='infA', 'clpp'='clpP',
    'its'='ITS', 'its1'='ITS1', 'its2'='ITS2', 'ets'='ETS',
    'rpb2'='RPB2', 'leafy'='LFY', 'fae1'='FAE1', 'g3pdh'='G3PDH'
  )

  low <- str_to_lower(x)
  exact_hit <- aliases[low]
  use_alias <- !is.na(exact_hit)
  x[use_alias] <- unname(exact_hit[use_alias])

  # Plastid families.
  x <- str_replace_all(x, regex('^ndh([a-k])$', ignore_case = TRUE), function(z) paste0('ndh', str_to_upper(str_sub(z,-1))))
  x <- str_replace_all(x, regex('^psa([a-z])$', ignore_case = TRUE), function(z) paste0('psa', str_to_upper(str_sub(z,-1))))
  x <- str_replace_all(x, regex('^(?:psb|pbs)([a-z])$', ignore_case = TRUE), function(z) paste0('psb', str_to_upper(str_sub(z,-1))))
  x <- str_replace_all(x, regex('^pet([a-z])$', ignore_case = TRUE), function(z) paste0('pet', str_to_upper(str_sub(z,-1))))
  x <- str_replace_all(x, regex('^atp([a-i])$', ignore_case = TRUE), function(z) paste0('atp', str_to_upper(str_sub(z,-1))))
  x <- str_replace_all(x, regex('^rpo([a-c])([12]?)$', ignore_case = TRUE), function(z) {
    m <- str_match(z, regex('^rpo([a-c])([12]?)$', ignore_case = TRUE))
    paste0('rpo', str_to_upper(m[,2]), m[,3])
  })
  x <- str_replace_all(x, regex('^rpl(\\d+)$', ignore_case = TRUE), function(z) paste0('rpl', str_extract(z,'\\d+')))
  x <- str_replace_all(x, regex('^rps(\\d+)$', ignore_case = TRUE), function(z) paste0('rps', str_extract(z,'\\d+')))
  x <- str_replace_all(x, regex('^ycf(\\d+)([a-z]?)$', ignore_case = TRUE), function(z) {
    m <- str_match(z, regex('^ycf(\\d+)([a-z]?)$', ignore_case = TRUE))
    paste0('ycf', m[,2], str_to_lower(m[,3]))
  })

  # Mitochondrial families.
  x <- str_replace_all(x, regex('^nad([1-9])l$', ignore_case = TRUE), function(z) paste0('nad', str_extract(z,'\\d+'), 'L'))
  x <- str_replace_all(x, regex('^nad([1-9])$', ignore_case = TRUE), function(z) paste0('nad', str_extract(z,'\\d+')))
  x <- str_replace_all(x, regex('^cox([1-3])$', ignore_case = TRUE), function(z) paste0('cox', str_extract(z,'\\d+')))

  # tRNA symbols and anticodons.
  x <- str_replace_all(x, regex('^trn([a-z])$', ignore_case = TRUE), function(z) paste0('trn', str_to_upper(str_sub(z,-1))))
  x <- str_replace_all(x, regex('^trnfm$', ignore_case = TRUE), 'trnfM')
  x <- str_replace_all(x, regex('^(trn(?:fM|[A-Z]))\\s*\\(([ACGTU]{3})\\)$', ignore_case = TRUE), function(z) {
    m <- str_match(z, regex('^(trn(?:fM|[A-Z]))\\s*\\(([ACGTU]{3})\\)$', ignore_case = TRUE))
    trn <- ifelse(str_to_lower(m[,2]) == 'trnfm', 'trnfM', paste0('trn', str_to_upper(str_sub(m[,2],-1))))
    paste0(trn, '-', str_to_upper(m[,3]))
  })

  # Ribosomal labels.
  x <- case_when(
    str_detect(x, regex('^5\\.8S(?: ribosomal RNA| rRNA)?$', ignore_case=TRUE)) ~ '5.8S rRNA',
    str_detect(x, regex('^18S(?: ribosomal RNA| rRNA)?$', ignore_case=TRUE)) ~ '18S rRNA',
    str_detect(x, regex('^26S(?: ribosomal RNA| rRNA)?$', ignore_case=TRUE)) ~ '26S rRNA',
    str_detect(x, regex('^28S(?: ribosomal RNA| rRNA)?$', ignore_case=TRUE)) ~ '28S rRNA',
    str_detect(x, regex('^large subunit ribosomal RNA$|^LSU(?: rRNA)?$', ignore_case=TRUE)) ~ 'LSU rRNA',
    TRUE ~ x
  )

  x
}

marker_symbol_pattern <- paste0(
  '(?:',
  'rbcL|matK|',
  'ndh[A-K]|nad[1-9](?:L)?|',
  'ps[ab][A-Z]|pet[A-Z]|atp[A-I]|',
  'rpl\\d+|rps\\d+|rpo[A-C][12]?|',
  'ycf\\d+[a-z]?|clpP\\d?|accD|cemA|ccsA|infA|',
  'trn(?:fM|[A-Z])(?:\\([ACGU]{3}\\)|-[ACGU]{3})?|',
  'cox[1-3]|cob|matR|mttB|ccmB|ccmC|ccmFc|ccmFn|',
  'ITS1|ITS2|ITS|ETS|RPB2|LFY',
  ')'
)

# ============================================================
# 2. GENBANK XML FETCHING WITH CACHE
# ============================================================

safe_cache_name <- function(uid) file.path(cache_dir, paste0(uid, '.xml'))

fetch_one_genbank_xml <- function(uid, max_attempts = 5) {
  cache_file <- safe_cache_name(uid)
  if (file_exists(cache_file) && file_size(cache_file) > 50) {
    return(read_file(cache_file))
  }

  for (attempt in seq_len(max_attempts)) {
    ans <- tryCatch(
      entrez_fetch(db='nuccore', id=uid, rettype='gb', retmode='xml'),
      error = function(e) NULL
    )

    if (!is.null(ans) && nchar(ans) > 50 && str_detect(ans, '<GBSeq')) {
      write_file(ans, cache_file)
      Sys.sleep(ifelse(nzchar(Sys.getenv('ENTREZ_KEY')), 0.12, 0.36))
      return(ans)
    }

    Sys.sleep(min(2^attempt, 30))
  }

  NA_character_
}

# ============================================================
# 3. XML FEATURE PARSING
# ============================================================

qualifier_table <- function(feature_node) {
  qs <- xml_find_all(feature_node, './/GBQualifier')
  if (length(qs) == 0) return(tibble(name=character(), value=character()))

  tibble(
    name = map_chr(qs, ~ xml_text(xml_find_first(.x, './GBQualifier_name'))),
    value = map_chr(qs, ~ xml_text(xml_find_first(.x, './GBQualifier_value')))
  )
}

get_qual <- function(qtab, name) {
  vals <- qtab$value[str_to_lower(qtab$name) == str_to_lower(name)]
  vals <- vals[!is.na(vals) & vals != '']
  if (length(vals) == 0) NA_character_ else paste(unique(vals), collapse=' | ')
}

parse_feature_coordinates <- function(feature_node) {
  location <- xml_text(xml_find_first(feature_node, './GBFeature_location'))
  intervals <- xml_find_all(feature_node, './/GBInterval')

  nums <- integer(0)
  if (length(intervals) > 0) {
    for (iv in intervals) {
      vals <- c(
        xml_text(xml_find_first(iv, './GBInterval_from')),
        xml_text(xml_find_first(iv, './GBInterval_to')),
        xml_text(xml_find_first(iv, './GBInterval_point'))
      )
      vals <- suppressWarnings(as.integer(vals[vals != '']))
      nums <- c(nums, vals[!is.na(vals)])
    }
  }

  tibble(
    feature_location = ifelse(is.na(location) || location == '', NA_character_, location),
    feature_start = ifelse(length(nums) == 0, NA_integer_, min(nums)),
    feature_end = ifelse(length(nums) == 0, NA_integer_, max(nums)),
    strand = ifelse(!is.na(location) && str_detect(location, '^complement\\('), '-', '+')
  )
}

organelle_to_genome <- function(organelle) {
  x <- str_to_lower(coalesce(organelle, ''))
  case_when(
    str_detect(x, 'plastid|chloroplast') ~ 'plastidial',
    str_detect(x, 'mitochond') ~ 'mitochondrial',
    TRUE ~ NA_character_
  )
}

infer_region_type <- function(feature_key, marker, text='') {
  fk <- str_to_lower(coalesce(feature_key, ''))
  z <- str_to_lower(paste(coalesce(marker,''), coalesce(text,'')))

  case_when(
    str_detect(z, 'intergenic spacer|spacer region') ~ 'intergenic_spacer',
    str_detect(z, '\\bintron\\b') ~ 'intron',
    fk == 'trna' | str_detect(coalesce(marker,''), '^trn') ~ 'tRNA',
    fk == 'rrna' ~ 'rRNA',
    str_detect(z, 'its1|its2|internal transcribed spacer|external transcribed spacer|\\bets\\b') ~ 'ribosomal_region',
    fk %in% c('cds','gene') ~ 'coding_gene',
    fk %in% c('misc_feature','misc_rna','regulatory') ~ 'other_region',
    TRUE ~ 'unknown'
  )
}

# Detect a marker from structured feature qualifiers. Preference is gene ->
# standard_name -> product/note parsing. The raw values are retained separately.
marker_from_feature <- function(feature_key, gene_q, product_q, note_q, standard_q) {
  # 1) /gene is usually the cleanest source.
  if (!is.na(gene_q) && gene_q != '') {
    # Some records place several gene names in one qualifier.
    parts <- str_split(gene_q, '\\s*[|;,]\\s*')[[1]]
    parts <- parts[parts != '']
    return(tibble(
      marker_raw = parts,
      marker = standardize_marker_name(parts),
      marker_source = 'feature_gene'
    ))
  }

  # 2) /standard_name.
  if (!is.na(standard_q) && standard_q != '') {
    return(tibble(
      marker_raw = standard_q,
      marker = standardize_marker_name(standard_q),
      marker_source = 'feature_standard_name'
    ))
  }

  text <- str_squish(paste(coalesce(product_q,''), coalesce(note_q,'')))
  if (text == '') return(tibble(marker_raw=character(), marker=character(), marker_source=character()))

  # 3) Explicit intergenic spacer in product/note.
  spacer <- str_match(
    text,
    regex(
      paste0('\\b(', marker_symbol_pattern, ')\\s*[-–—]\\s*(', marker_symbol_pattern,
             ')\\s+(?:intergenic\\s+spacer(?:\\s+region)?|spacer(?:\\s+region)?)\\b'),
      ignore_case=TRUE
    )
  )
  if (!is.na(spacer[1,2]) && !is.na(spacer[1,3])) {
    m <- paste0(standardize_marker_name(spacer[1,2]), '-', standardize_marker_name(spacer[1,3]))
    return(tibble(marker_raw=spacer[1,1], marker=m, marker_source='feature_spacer_text'))
  }

  # 4) Ribosomal regions described by words.
  ribo <- case_when(
    str_detect(text, regex('internal transcribed spacer 1.*internal transcribed spacer 2|ITS1.*ITS2', ignore_case=TRUE)) ~ 'ITS1-5.8S-ITS2',
    str_detect(text, regex('internal transcribed spacer 1|\\bITS1\\b', ignore_case=TRUE)) ~ 'ITS1',
    str_detect(text, regex('internal transcribed spacer 2|\\bITS2\\b', ignore_case=TRUE)) ~ 'ITS2',
    str_detect(text, regex('external transcribed spacer|\\bETS\\b', ignore_case=TRUE)) ~ 'ETS',
    str_detect(text, regex('5\\.8S ribosomal RNA', ignore_case=TRUE)) ~ '5.8S rRNA',
    str_detect(text, regex('18S ribosomal RNA', ignore_case=TRUE)) ~ '18S rRNA',
    str_detect(text, regex('26S ribosomal RNA', ignore_case=TRUE)) ~ '26S rRNA',
    str_detect(text, regex('28S ribosomal RNA|large subunit ribosomal RNA', ignore_case=TRUE)) ~ '28S rRNA',
    TRUE ~ NA_character_
  )
  if (!is.na(ribo)) {
    return(tibble(marker_raw=text, marker=ribo, marker_source='feature_product'))
  }

  # 5) Written tRNA product.
  trna <- str_match(text, regex('\\btRNA-([A-Za-z]{3})\\b', ignore_case=TRUE))
  if (!is.na(trna[1,2])) {
    tmap <- c(Ala='trnA',Arg='trnR',Asn='trnN',Asp='trnD',Cys='trnC',Gln='trnQ',
              Glu='trnE',Gly='trnG',His='trnH',Ile='trnI',Leu='trnL',Lys='trnK',
              Met='trnM',Phe='trnF',Pro='trnP',Ser='trnS',Thr='trnT',Trp='trnW',
              Tyr='trnY',Val='trnV')
    aa <- str_to_title(trna[1,2])
    if (aa %in% names(tmap)) {
      return(tibble(marker_raw=trna[1,1], marker=unname(tmap[aa]), marker_source='feature_product'))
    }
  }

  # 6) Recognizable marker symbol anywhere in product/note.
  syms <- str_extract_all(text, regex(marker_symbol_pattern, ignore_case=TRUE))[[1]]
  syms <- unique(standardize_marker_name(syms))
  syms <- syms[!is.na(syms) & syms != '']
  if (length(syms) > 0) {
    return(tibble(marker_raw=syms, marker=syms, marker_source='feature_text_symbol'))
  }

  tibble(marker_raw=character(), marker=character(), marker_source=character())
}

parse_genbank_xml <- function(xml_text_string, requested_uid) {
  empty_record <- tibble(
    uid=character(), accession=character(), accession_version=character(),
    record_definition=character(), organelle=character(), structured_genome=character()
  )
  empty_features <- tibble(
    uid=character(), feature_key=character(), marker_raw=character(), marker=character(),
    marker_source=character(), region_type=character(), feature_location=character(),
    feature_start=integer(), feature_end=integer(), strand=character(),
    gene_qualifier=character(), product_qualifier=character(), note_qualifier=character()
  )

  if (is.na(xml_text_string) || xml_text_string == '') {
    return(list(record=empty_record, features=empty_features))
  }

  doc <- tryCatch(read_xml(xml_text_string), error=function(e) NULL)
  if (is.null(doc)) return(list(record=empty_record, features=empty_features))

  gb <- xml_find_first(doc, './/GBSeq')
  if (inherits(gb, 'xml_missing')) return(list(record=empty_record, features=empty_features))

  accession <- xml_text(xml_find_first(gb, './GBSeq_primary-accession'))
  accession_version <- xml_text(xml_find_first(gb, './GBSeq_accession-version'))
  definition <- xml_text(xml_find_first(gb, './GBSeq_definition'))

  features <- xml_find_all(gb, './/GBFeature')

  # Source qualifiers provide the strongest genome evidence.
  organelle <- NA_character_
  if (length(features) > 0) {
    source_idx <- which(map_chr(features, ~ xml_text(xml_find_first(.x, './GBFeature_key'))) == 'source')
    if (length(source_idx) > 0) {
      qsource <- qualifier_table(features[[source_idx[1]]])
      organelle <- get_qual(qsource, 'organelle')
    }
  }

  record_tbl <- tibble(
    uid = as.character(requested_uid),
    accession = accession,
    accession_version = accession_version,
    record_definition = definition,
    organelle = organelle,
    structured_genome = organelle_to_genome(organelle)
  )

  if (length(features) == 0) return(list(record=record_tbl, features=empty_features))

  feat_tbl <- map_dfr(features, function(ft) {
    key <- xml_text(xml_find_first(ft, './GBFeature_key'))
    if (is.na(key) || key == 'source') return(tibble())

    q <- qualifier_table(ft)
    gene_q <- get_qual(q, 'gene')
    product_q <- get_qual(q, 'product')
    note_q <- get_qual(q, 'note')
    standard_q <- get_qual(q, 'standard_name')

    markers <- marker_from_feature(key, gene_q, product_q, note_q, standard_q)
    if (nrow(markers) == 0) return(tibble())

    coords <- parse_feature_coordinates(ft)

    markers %>%
      mutate(
        uid = as.character(requested_uid),
        feature_key = key,
        region_type = map2_chr(marker, marker_raw, ~ infer_region_type(key, .x, paste(.y, note_q, product_q))),
        feature_location = coords$feature_location,
        feature_start = coords$feature_start,
        feature_end = coords$feature_end,
        strand = coords$strand,
        gene_qualifier = gene_q,
        product_qualifier = product_q,
        note_qualifier = note_q
      ) %>%
      select(uid, feature_key, marker_raw, marker, marker_source, region_type,
             feature_location, feature_start, feature_end, strand,
             gene_qualifier, product_qualifier, note_qualifier)
  })

  list(record=record_tbl, features=feat_tbl)
}

# ============================================================
# 4. FETCH + PARSE ALL RECORDS
# ============================================================

message('Fetching/parsing GenBank XML for ', nrow(final_seqs_to_parse), ' accessions...')

parsed <- map(final_seqs_to_parse$uid, function(id) {
  xml <- fetch_one_genbank_xml(id)
  parse_genbank_xml(xml, id)
})

record_metadata <- map_dfr(parsed, 'record')
accession_features <- map_dfr(parsed, 'features')

write_tsv(record_metadata, file.path(out_dir, 'genbank_record_metadata.tsv'))
write_tsv(accession_features, file.path(out_dir, 'genbank_accession_features.tsv'))

empty_marker_table <- function() {
  tibble(
    marker_raw=character(),
    marker=character(),
    marker_source=character(),
    region_type=character()
  )
}

# ============================================================
# 5. TITLE FALLBACK PARSER
# ============================================================
# This parser is intentionally much smaller than the old giant case_when().
# It is only used to rescue markers absent from structured feature parsing.

parse_title_markers <- function(title) {
  x <- str_squish(as.character(title))
  if (is.na(x) || x == '') return(empty_marker_table())

  # A) Explicit intergenic spacer. Highest priority: don't reinterpret its
  # flanking gene symbols as independent genes from title text alone.
  spacer <- str_match(
    x,
    regex(
      paste0('\\b(', marker_symbol_pattern, ')\\s*[-–—]\\s*(', marker_symbol_pattern,
             ')\\s+(?:intergenic\\s+spacer(?:\\s+region)?|spacer(?:\\s+region)?)\\b'),
      ignore_case=TRUE
    )
  )
  if (!is.na(spacer[1,2]) && !is.na(spacer[1,3])) {
    marker <- paste0(standardize_marker_name(spacer[1,2]), '-', standardize_marker_name(spacer[1,3]))
    return(tibble(
      marker_raw = spacer[1,1], marker=marker,
      marker_source='title_explicit_spacer', region_type='intergenic_spacer'
    ))
  }

  # B) ITS/ribosomal composite regions.
  has_its1 <- str_detect(x, regex('internal transcribed spacer 1|\\bITS1\\b', ignore_case=TRUE))
  has_its2 <- str_detect(x, regex('internal transcribed spacer 2|\\bITS2\\b', ignore_case=TRUE))
  has_58s  <- str_detect(x, regex('\\b5\\.8S\\b', ignore_case=TRUE))

  if (has_its1 && has_its2) {
    return(tibble(marker_raw='ITS1 + 5.8S + ITS2', marker='ITS1-5.8S-ITS2',
                  marker_source='title_ribosomal_region', region_type='ribosomal_region'))
  }
  if (has_its1 && has_58s) {
    return(tibble(marker_raw='ITS1 + 5.8S', marker='ITS1-5.8S',
                  marker_source='title_ribosomal_region', region_type='ribosomal_region'))
  }
  if (has_58s && has_its2) {
    return(tibble(marker_raw='5.8S + ITS2', marker='5.8S-ITS2',
                  marker_source='title_ribosomal_region', region_type='ribosomal_region'))
  }

  ribo <- tibble(
    pattern=c(
      'internal transcribed spacer 1|\\bITS1\\b',
      'internal transcribed spacer 2|\\bITS2\\b',
      'internal transcribed spacer|\\bITS\\b',
      'external transcribed spacer|\\bETS\\b',
      '\\b5\\.8S\\b', '\\b18S\\b', '\\b26S\\b', '\\b28S\\b|large subunit ribosomal RNA|\\bLSU\\b'
    ),
    marker=c('ITS1','ITS2','ITS','ETS','5.8S rRNA','18S rRNA','26S rRNA','28S rRNA')
  )
  rhit <- ribo %>% filter(map_lgl(pattern, ~ str_detect(x, regex(.x, ignore_case=TRUE))))
  if (nrow(rhit) > 0) {
    return(rhit %>% transmute(marker_raw=marker, marker=marker,
                              marker_source='title_ribosomal_region', region_type='ribosomal_region'))
  }

  # C) Generic explicit marker symbols. Return ALL unique symbols, not first hit.
  symbols <- str_extract_all(x, regex(marker_symbol_pattern, ignore_case=TRUE))[[1]]
  symbols <- unique(standardize_marker_name(symbols))
  symbols <- symbols[!is.na(symbols) & symbols != '']

  if (length(symbols) > 0) {
    return(tibble(
      marker_raw=symbols,
      marker=symbols,
      marker_source='title_symbol',
      region_type=case_when(
        str_detect(symbols, '^trn') ~ 'tRNA',
        TRUE ~ 'coding_gene'
      )
    ))
  }

  # D) A small set of common full-name descriptions that do not expose the symbol.
  descriptions <- tribble(
    ~pattern, ~marker,
    'ribulose[ -]1,5-bi?sphosphate carboxylase(?:/oxygenase)? large subunit|ribulose-1,5-bisphosphate carboxylase oxygenase', 'rbcL',
    'maturase K', 'matK',
    'RNA polymerase beta[′\' ]? subunit C1', 'rpoC1',
    'RNA polymerase beta[′\' ]{0,2} subunit C2', 'rpoC2',
    'RNA polymerase beta(?: chain| subunit)', 'rpoB',
    'ribosomal protein L16', 'rpl16',
    'ribosomal protein S16', 'rps16',
    'acetyl-coenzyme A carboxylase', 'accD'
  )

  dhit <- descriptions %>%
    filter(map_lgl(pattern, ~ str_detect(x, regex(.x, ignore_case=TRUE))))

  if (nrow(dhit) > 0) {
    return(dhit %>%
      transmute(marker_raw=marker, marker=marker,
                marker_source='title_description', region_type='coding_gene') %>%
      distinct(marker, .keep_all=TRUE))
  }

  empty_marker_table()
}

# ============================================================
# 6. INPUT GENE FIELD AS SECONDARY STRUCTURED EVIDENCE
# ============================================================

parse_input_gene <- function(gene) {
  g <- str_squish(as.character(gene))
  if (is.na(g) || g == '') return(empty_marker_table())

  # Pipe/semicolon/comma are treated as multiple explicit annotations.
  # Hyphens are NOT split here because they may represent a real spacer/region.
  parts <- str_split(g, '\\s*[|;,]\\s*')[[1]]
  parts <- parts[parts != '']
  parts <- unique(standardize_marker_name(parts))

  tibble(
    marker_raw=parts,
    marker=parts,
    marker_source='input_gene_field',
    region_type=case_when(
      str_detect(parts, regex('ITS|ETS|5\\.8S|18S|26S|28S|LSU', ignore_case=TRUE)) ~ 'ribosomal_region',
      str_detect(parts, '^trn') ~ 'tRNA',
      TRUE ~ 'coding_gene'
    )
  )
}

# ============================================================
# 7. BUILD MARKER CANDIDATES FROM ALL EVIDENCE SOURCES
# ============================================================

feature_candidates <- accession_features %>%
  select(uid, marker_raw, marker, marker_source, region_type,
         feature_key, feature_location, feature_start, feature_end, strand,
         gene_qualifier, product_qualifier, note_qualifier) %>%
  mutate(evidence_rank=1L)

input_candidates <- final_seqs_to_parse %>%
  transmute(uid, marker_info=map(input_gene, parse_input_gene)) %>%
  unnest(marker_info) %>%
  mutate(
    feature_key=NA_character_, feature_location=NA_character_,
    feature_start=NA_integer_, feature_end=NA_integer_, strand=NA_character_,
    gene_qualifier=NA_character_, product_qualifier=NA_character_, note_qualifier=NA_character_,
    evidence_rank=2L
  )

title_candidates <- final_seqs_to_parse %>%
  transmute(uid, marker_info=map(title, parse_title_markers)) %>%
  unnest(marker_info) %>%
  mutate(
    feature_key=NA_character_, feature_location=NA_character_,
    feature_start=NA_integer_, feature_end=NA_integer_, strand=NA_character_,
    gene_qualifier=NA_character_, product_qualifier=NA_character_, note_qualifier=NA_character_,
    evidence_rank=3L
  )

all_marker_candidates <- bind_rows(feature_candidates, input_candidates, title_candidates) %>%
  filter(!is.na(marker), marker != '') %>%
  mutate(marker=standardize_marker_name(marker))

# Select one best evidence row per accession x marker, prioritizing FEATURES.
# The complete candidate table is also written for auditability.
accession_markers_core <- all_marker_candidates %>%
  arrange(uid, marker, evidence_rank) %>%
  group_by(uid, marker) %>%
  slice_head(n=1) %>%
  ungroup()

# ============================================================
# 8. GENOME CLASSIFICATION WITH PROVENANCE
# ============================================================

infer_genome_from_title <- function(title) {
  x <- coalesce(as.character(title), '')
  case_when(
    str_detect(x, regex('chloroplast|plastid|plastome', ignore_case=TRUE)) ~ 'plastidial',
    str_detect(x, regex('mitochondri(?:on|al)', ignore_case=TRUE)) ~ 'mitochondrial',
    str_detect(x, regex('nuclear|genomic DNA', ignore_case=TRUE)) ~ 'nuclear',
    TRUE ~ NA_character_
  )
}

# Conservative marker inference: only marker groups that are strongly associated
# with one compartment are used. Ambiguous ribosomal/protein names remain NA.
infer_genome_from_marker <- function(marker) {
  m <- coalesce(marker, '')
  case_when(
    str_detect(m, regex('^(rbcL|matK|ndh[A-K]|psa[A-Z]|psb[A-Z]|pet[A-Z]|rpo[A-C][12]?|accD|cemA|ccsA|infA|ycf\\d+)$', ignore_case=TRUE)) ~ 'plastidial',
    str_detect(m, regex('^(cox[1-3]|cob|nad[1-9]L?|matR|mttB|ccmB|ccmC|ccmFc|ccmFn)$', ignore_case=TRUE)) ~ 'mitochondrial',
    str_detect(m, regex('^(ITS|ITS1|ITS2|ITS1-5\\.8S|5\\.8S-ITS2|ITS1-5\\.8S-ITS2|ETS|RPB2|LFY)$', ignore_case=TRUE)) ~ 'nuclear',
    TRUE ~ NA_character_
  )
}

base_metadata <- final_seqs_to_parse %>%
  left_join(record_metadata, by='uid') %>%
  mutate(
    input_genome_normalized = normalize_genome_value(input_genome),
    title_genome = map_chr(title, infer_genome_from_title)
  )

accession_markers <- accession_markers_core %>%
  left_join(base_metadata, by='uid') %>%
  mutate(
    inferred_marker_genome = map_chr(marker, infer_genome_from_marker),
    genome = coalesce(
      structured_genome,
      input_genome_normalized,
      title_genome,
      inferred_marker_genome
    ),
    genome_source = case_when(
      !is.na(structured_genome) ~ 'feature_organelle',
      !is.na(input_genome_normalized) ~ 'entrez_input_field',
      !is.na(title_genome) ~ 'title_explicit',
      !is.na(inferred_marker_genome) ~ 'marker_inference',
      TRUE ~ 'unresolved'
    )
  ) %>%
  select(
    taxa_accepted, uid, caption, accession, accession_version, title,
    genome, genome_source, structured_genome, organelle, input_genome,
    marker, marker_raw, marker_source, region_type,
    feature_key, feature_location, feature_start, feature_end, strand,
    gene_qualifier, product_qualifier, note_qualifier,
    native, introduced, somewhere,
    everything()
  ) %>%
  distinct(uid, marker, .keep_all=TRUE)

# ============================================================
# 9. AUDITS
# ============================================================

# Records for which no marker could be recovered from FEATURES, input gene,
# or title fallback.
unresolved_markers <- final_seqs_to_parse %>%
  anti_join(accession_markers %>% distinct(uid), by='uid') %>%
  select(taxa_accepted, uid, caption, title, input_genome, input_gene, everything())

# Marker records whose genome remains unresolved.
unresolved_genomes <- accession_markers %>%
  filter(is.na(genome)) %>%
  select(taxa_accepted, uid, caption, title, marker, marker_source,
         structured_genome, input_genome, organelle, genome_source)

# Conflicts are not silently reconciled. They are reported whenever multiple
# non-NA evidence sources disagree about genome identity.
genome_conflicts <- accession_markers %>%
  rowwise() %>%
  mutate(
    genome_evidence_values = paste(unique(na.omit(c(
      structured_genome, input_genome_normalized, title_genome, inferred_marker_genome
    ))), collapse=' | '),
    genome_n_distinct = n_distinct(na.omit(c(
      structured_genome, input_genome_normalized, title_genome, inferred_marker_genome
    )))
  ) %>%
  ungroup() %>%
  filter(genome_n_distinct > 1) %>%
  select(taxa_accepted, uid, caption, title, marker,
         structured_genome, input_genome_normalized, title_genome,
         inferred_marker_genome, genome, genome_source, genome_evidence_values)

# Feature-vs-title/input disagreements about marker identity are retained here.
marker_candidate_audit <- all_marker_candidates %>%
  left_join(final_seqs_to_parse %>% select(uid, taxa_accepted, caption, title), by='uid') %>%
  arrange(uid, evidence_rank, marker)

# ============================================================
# 10. COUNT DISTINCT ACCESSIONS AND APPLY 5/5 FILTER
# ============================================================

marker_counts <- accession_markers %>%
  filter(!is.na(genome), !is.na(marker)) %>%
  group_by(taxa_accepted, genome, marker) %>%
  summarise(
    n_native = n_distinct(uid[native == 1]),
    n_introduced = n_distinct(uid[introduced == 1]),
    n_somewhere = n_distinct(uid[somewhere == 1]),
    n_accessions = n_distinct(uid),
    n_feature_annotated = n_distinct(uid[str_detect(marker_source, '^feature_')]),
    .groups='drop'
  )

eligible_markers <- marker_counts %>%
  filter(n_native >= 5, n_introduced >= 5)

seqs_to_download <- accession_markers %>%
  semi_join(eligible_markers, by=c('taxa_accepted','genome','marker')) %>%
  arrange(taxa_accepted, genome, marker, uid)

# ============================================================
# 11. WRITE OUTPUTS
# ============================================================

write_tsv(accession_markers,       file.path(out_dir, 'accession_markers.tsv'))
write_tsv(marker_counts,           file.path(out_dir, 'marker_counts.tsv'))
write_tsv(eligible_markers,        file.path(out_dir, 'eligible_markers.tsv'))
write_tsv(seqs_to_download,        file.path(out_dir, 'seqs_to_download.tsv'))
write_tsv(unresolved_markers,      file.path(out_dir, 'unresolved_markers.tsv'))
write_tsv(unresolved_genomes,      file.path(out_dir, 'unresolved_genomes.tsv'))
write_tsv(genome_conflicts,        file.path(out_dir, 'genome_conflicts.tsv'))
write_tsv(marker_candidate_audit,  file.path(out_dir, 'marker_candidate_audit.tsv'))
write_tsv(all_marker_candidates,   file.path(out_dir, 'all_marker_candidates.tsv'))

# A compact accession-level audit table preserving the original fields plus
# structured GenBank metadata.
seqs_to_parse_good <- base_metadata %>%
  mutate(
    genome = coalesce(structured_genome, input_genome_normalized, title_genome),
    genome_source = case_when(
      !is.na(structured_genome) ~ 'feature_organelle',
      !is.na(input_genome_normalized) ~ 'entrez_input_field',
      !is.na(title_genome) ~ 'title_explicit',
      TRUE ~ 'unresolved'
    )
  )

write_tsv(seqs_to_parse_good, file.path(out_dir, 'seqs_to_parse_good.tsv'))

# ============================================================
# 12. SUMMARY TO CONSOLE
# ============================================================

message('Done.')
message('Accessions in input: ', n_distinct(final_seqs_to_parse$uid))
message('Accessions with >=1 marker: ', n_distinct(accession_markers$uid))
message('Unique accession-marker combinations: ', nrow(accession_markers))
message('Feature-derived accession-marker combinations: ',
        sum(str_detect(accession_markers$marker_source, '^feature_'), na.rm=TRUE))
message('Unresolved marker accessions: ', n_distinct(unresolved_markers$uid))
message('Genome conflicts requiring audit: ', nrow(genome_conflicts))
message('Eligible species/genome/marker groups (>=5 native and >=5 introduced): ', nrow(eligible_markers))
