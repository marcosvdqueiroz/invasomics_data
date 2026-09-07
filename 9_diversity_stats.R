# ============================================================
# SCRIPT 9
# Diversity, neutrality, haplotype, and differentiation stats
# Native vs Introduced
# ============================================================
#
# This script follows:
#   8_align_seqs_structured_regenerated.R
#
# Core principles:
#   * use script 8's alignment_manifest.tsv as the source of truth;
#   * use the JOINT coverage-filtered alignment for all analyses;
#   * derive Native and Introduced subsets from the SAME alignment;
#   * do NOT realign sequences here;
#   * use complete canonical A/C/G/T sites for haplotype identity,
#     haplotype sharing, Snn, and Phi_ST;
#   * retain pairwise-deletion estimates of pi as descriptive values,
#     but also compute complete-site versions for sensitivity checking;
#   * keep QC as explicit flags rather than silently removing pairs;
#   * add FDR-adjusted p-values across markers.
#
# Required packages:
#   mamba install -c conda-forge -c bioconda \
#       r-readr r-dplyr r-tidyr r-stringr r-purrr r-fs \
#       r-ape r-pegas bioconductor-biostrings
#
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
library(ape)
library(pegas)
library(Biostrings)


# ============================================================
# 2. SETTINGS
# ============================================================

script8_root <- "./outputs/8_align_seqs"

alignment_manifest_file <- file.path(
  script8_root,
  "alignment_manifest.tsv"
)

alignment_qc_file <- file.path(
  script8_root,
  "alignment_qc_summary.tsv"
)

output_root <- "./outputs/9_diversity_stats"

dir_create(output_root, recurse = TRUE)

# ------------------------------------------------------------
# Statistical settings
# ------------------------------------------------------------

# Permutations for Phi_ST and Snn.
n_perm <- 9999L

# Simulations for Ramos-Onsins & Rozas' R2 test.
r2_B <- 9999L

# Reproducibility for permutation/simulation-based tests.
set.seed(1)

# ------------------------------------------------------------
# Optional QC flags
# ------------------------------------------------------------
#
# These are AUDIT FLAGS only. They do not automatically remove rows.
# Script 8 has already done sequence and alignment QC.

qc_pi_abs_threshold <- 0.05
qc_prop_seg_threshold <- 0.20
qc_gap_fraction_threshold <- 0.50
qc_min_fraction_sites_retained <- 0.25

# For per-marker robust outlier detection.
qc_mad_multiplier <- 8

# ------------------------------------------------------------
# Rarefaction
# ------------------------------------------------------------

# Number of random draws when rarefying haplotype richness to the
# smaller of Native / Introduced sample sizes.
haplotype_rarefaction_B <- 1000L


# ============================================================
# 3. HELPERS
# ============================================================

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}


read_dna_safely <- function(path) {
  
  if (
    is.na(path) ||
    path == "" ||
    !file_exists(path)
  ) {
    return(NULL)
  }
  
  tryCatch(
    read.dna(
      path,
      format = "fasta"
    ),
    error = function(e) {
      message(
        "Could not read alignment: ",
        path,
        " -- ",
        conditionMessage(e)
      )
      NULL
    }
  )
}


parse_status_from_names <- function(x) {
  
  case_when(
    str_detect(x, "\\|status=Native(?:\\||$)") ~ "Native",
    str_detect(x, "\\|status=Introduced(?:\\||$)") ~ "Introduced",
    TRUE ~ NA_character_
  )
}


complete_canonical_columns <- function(dnabin) {
  
  if (is.null(dnabin)) {
    return(NULL)
  }
  
  # Convert to ordinary nucleotide characters only to determine
  # which columns contain canonical A/C/G/T in every sequence.
  char_mat <- as.character(dnabin)
  
  if (
    is.null(dim(char_mat)) ||
    nrow(char_mat) == 0L ||
    ncol(char_mat) == 0L
  ) {
    return(NULL)
  }
  
  char_mat <- toupper(char_mat)
  
  valid <- char_mat %in% c(
    "A",
    "C",
    "G",
    "T"
  )
  
  valid <- matrix(
    valid,
    nrow = nrow(char_mat),
    ncol = ncol(char_mat),
    dimnames = dimnames(char_mat)
  )
  
  keep <- colSums(valid) == nrow(char_mat)
  
  if (!any(keep)) {
    return(NULL)
  }
  
  # Preserve DNAbin class by subsetting the original object directly.
  dnabin[
    ,
    keep,
    drop = FALSE
  ]
}


subset_dnabin_rows <- function(
  dnabin,
  idx
) {
  
  if (is.null(dnabin)) {
    return(NULL)
  }
  
  if (length(idx) == 0L) {
    return(NULL)
  }
  
  # Preserve DNAbin class by subsetting the original object directly.
  dnabin[
    idx,
    ,
    drop = FALSE
  ]
}


count_complete_sites <- function(dnabin) {
  
  x <- complete_canonical_columns(dnabin)
  
  if (is.null(x)) return(0L)
  
  ncol(as.matrix(x))
}


# ============================================================
# 4. FU'S Fs
# ============================================================
#
# Point estimate only.
# No rigorous coalescent p-value is produced here.
# Validate against dedicated software before making Fu's Fs a central
# inferential result.
# ============================================================

fu_fs <- function(dnabin) {
  
  mat <- as.matrix(dnabin)
  
  n <- nrow(mat)
  
  if (is.na(n) || n < 2L) {
    return(
      list(
        Fs = NA_real_,
        theta_pi = NA_real_,
        k_observed = NA_integer_
      )
    )
  }
  
  h <- tryCatch(
    haplotype(dnabin),
    error = function(e) NULL
  )
  
  if (is.null(h)) {
    return(
      list(
        Fs = NA_real_,
        theta_pi = NA_real_,
        k_observed = NA_integer_
      )
    )
  }
  
  k_obs <- nrow(h)
  
  d <- tryCatch(
    dist.dna(
      dnabin,
      model = "N",
      pairwise.deletion = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(d)) {
    return(
      list(
        Fs = NA_real_,
        theta_pi = NA_real_,
        k_observed = k_obs
      )
    )
  }
  
  theta_pi <- mean(d)
  
  if (!is.finite(theta_pi) || theta_pi <= 0) {
    return(
      list(
        Fs = NA_real_,
        theta_pi = theta_pi,
        k_observed = k_obs
      )
    )
  }
  
  logS <- matrix(
    -Inf,
    nrow = n + 1L,
    ncol = n + 1L
  )
  
  logS[1L, 1L] <- 0
  
  for (nn in seq_len(n)) {
    
    for (kk in seq_len(nn)) {
      
      term1 <- logS[nn, kk]
      
      term2 <- if (nn - 1L > 0L) {
        log(nn - 1L) + logS[nn, kk + 1L]
      } else {
        -Inf
      }
      
      m <- max(term1, term2)
      
      logS[nn + 1L, kk + 1L] <- if (is.infinite(m)) {
        -Inf
      } else {
        m + log(
          exp(term1 - m) +
            exp(term2 - m)
        )
      }
    }
  }
  
  log_rising_factorial <- sum(
    log(theta_pi + 0:(n - 1L))
  )
  
  log_p <- vapply(
    seq_len(n),
    function(kk) {
      logS[n + 1L, kk + 1L] +
        kk * log(theta_pi) -
        log_rising_factorial
    },
    numeric(1)
  )
  
  vals <- log_p[k_obs:n]
  
  m <- max(vals)
  
  Sprime <- sum(
    exp(vals - m)
  ) * exp(m)
  
  Sprime <- min(
    max(Sprime, 1e-300),
    1 - 1e-300
  )
  
  Fs <- log(
    Sprime /
      (1 - Sprime)
  )
  
  list(
    Fs = Fs,
    theta_pi = theta_pi,
    k_observed = k_obs
  )
}


# ============================================================
# 5. MISMATCH DISTRIBUTION
# ============================================================

compute_mismatch <- function(dnabin) {
  
  d <- tryCatch(
    dist.dna(
      dnabin,
      model = "N",
      pairwise.deletion = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(d)) {
    return(
      list(
        raggedness = NA_real_,
        mean_pairwise_diff = NA_real_
      )
    )
  }
  
  diffs <- as.vector(d)
  
  if (
    length(diffs) == 0L ||
    all(is.na(diffs))
  ) {
    return(
      list(
        raggedness = NA_real_,
        mean_pairwise_diff = NA_real_
      )
    )
  }
  
  diffs <- diffs[is.finite(diffs)]
  
  if (length(diffs) == 0L) {
    return(
      list(
        raggedness = NA_real_,
        mean_pairwise_diff = NA_real_
      )
    )
  }
  
  max_diff <- floor(max(diffs))
  
  freq_table <- table(
    factor(
      floor(diffs),
      levels = 0:max_diff
    )
  )
  
  freq <- as.numeric(freq_table)
  
  if (sum(freq) == 0L) {
    raggedness <- NA_real_
  } else {
    freq <- freq / sum(freq)
    raggedness <- sum(diff(freq)^2)
  }
  
  list(
    raggedness = raggedness,
    mean_pairwise_diff = mean(diffs)
  )
}


# ============================================================
# 6. HAPLOTYPE RICHNESS RAREFACTION
# ============================================================

haplotype_ids <- function(dnabin) {
  
  if (is.null(dnabin)) return(NULL)
  
  h <- tryCatch(
    haplotype(dnabin),
    error = function(e) NULL
  )
  
  if (is.null(h)) return(NULL)
  
  idx <- attr(h, "index")
  
  ids <- rep(
    NA_integer_,
    nrow(as.matrix(dnabin))
  )
  
  for (i in seq_along(idx)) {
    ids[idx[[i]]] <- i
  }
  
  ids
}


rarefied_haplotype_richness <- function(
  hap_ids,
  sample_size,
  B = 1000L
) {
  
  hap_ids <- hap_ids[!is.na(hap_ids)]
  
  n <- length(hap_ids)
  
  if (
    n == 0L ||
    sample_size < 1L ||
    sample_size > n
  ) {
    return(NA_real_)
  }
  
  if (sample_size == n) {
    return(
      length(unique(hap_ids))
    )
  }
  
  draws <- replicate(
    B,
    {
      idx <- sample(
        seq_len(n),
        size = sample_size,
        replace = FALSE
      )
      
      length(
        unique(hap_ids[idx])
      )
    }
  )
  
  mean(draws)
}


# ============================================================
# 7. WITHIN-GROUP DIVERSITY + NEUTRALITY
# ============================================================

compute_diversity_stats <- function(dnabin) {
  
  if (is.null(dnabin)) {
    return(tibble())
  }
  
  mat <- as.matrix(dnabin)
  
  n_seq <- nrow(mat)
  n_sites <- ncol(mat)
  
  if (
    n_seq < 2L ||
    n_sites < 1L
  ) {
    return(
      tibble(
        n_seq = n_seq,
        alignment_sites = n_sites,
        complete_sites = count_complete_sites(dnabin),
        n_seg_sites = NA_integer_,
        n_haplotypes = NA_integer_,
        hap_diversity = NA_real_,
        nuc_diversity_pairwise = NA_real_,
        nuc_diversity_complete = NA_real_,
        theta_w = NA_real_,
        tajima_D = NA_real_,
        tajima_D_pval = NA_real_,
        fu_Fs = NA_real_,
        R2 = NA_real_,
        R2_pval = NA_real_,
        mismatch_raggedness = NA_real_,
        mean_pairwise_diff = NA_real_
      )
    )
  }
  
  complete_dna <- complete_canonical_columns(dnabin)
  
  n_complete <- if (is.null(complete_dna)) {
    0L
  } else {
    ncol(as.matrix(complete_dna))
  }
  
  n_seg_sites <- tryCatch(
    length(seg.sites(dnabin)),
    error = function(e) NA_integer_
  )
  
  hap_dna <- if (!is.null(complete_dna)) {
    complete_dna
  } else {
    dnabin
  }
  
  hap <- tryCatch(
    haplotype(hap_dna),
    error = function(e) NULL
  )
  
  n_haplotypes <- if (!is.null(hap)) {
    nrow(hap)
  } else {
    NA_integer_
  }
  
  hap_diversity <- tryCatch(
    hap.div(
      hap_dna,
      variance = FALSE
    ),
    error = function(e) NA_real_
  )
  
  nuc_div_pairwise <- tryCatch(
    nuc.div(
      dnabin,
      pairwise.deletion = TRUE
    ),
    error = function(e) NA_real_
  )
  
  nuc_div_complete <- if (!is.null(complete_dna)) {
    tryCatch(
      nuc.div(
        complete_dna,
        pairwise.deletion = FALSE
      ),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  
  theta_w <- tryCatch(
    theta.s(dnabin),
    error = function(e) NA_real_
  )
  
  taj <- tryCatch(
    tajima.test(dnabin),
    error = function(e) NULL
  )
  
  tajima_D <- if (!is.null(taj)) {
    as.numeric(taj$D)
  } else {
    NA_real_
  }
  
  tajima_D_pval <- if (!is.null(taj)) {
    as.numeric(taj$Pval.normal)
  } else {
    NA_real_
  }
  
  fs_res <- tryCatch(
    fu_fs(hap_dna),
    error = function(e) {
      list(Fs = NA_real_)
    }
  )
  
  fu_Fs <- fs_res$Fs
  
  mm_res <- tryCatch(
    compute_mismatch(dnabin),
    error = function(e) {
      list(
        raggedness = NA_real_,
        mean_pairwise_diff = NA_real_
      )
    }
  )
  
  r2_res <- tryCatch(
    pegas::R2.test(
      dnabin,
      B = r2_B,
      plot = FALSE
    ),
    error = function(e) NULL
  )
  
  R2 <- NA_real_
  R2_pval <- NA_real_
  
  if (!is.null(r2_res)) {
    
    R2 <- tryCatch(
      as.numeric(r2_res$R2),
      error = function(e) NA_real_
    )
    
    R2_pval <- tryCatch(
      as.numeric(r2_res$P.val),
      error = function(e) NA_real_
    )
  }
  
  tibble(
    n_seq = n_seq,
    alignment_sites = n_sites,
    complete_sites = n_complete,
    n_seg_sites = n_seg_sites,
    n_haplotypes = n_haplotypes,
    hap_diversity = hap_diversity,
    nuc_diversity_pairwise = nuc_div_pairwise,
    nuc_diversity_complete = nuc_div_complete,
    theta_w = theta_w,
    tajima_D = tajima_D,
    tajima_D_pval = tajima_D_pval,
    fu_Fs = fu_Fs,
    R2 = R2,
    R2_pval = R2_pval,
    mismatch_raggedness = mm_res$raggedness,
    mean_pairwise_diff = mm_res$mean_pairwise_diff
  )
}


# ============================================================
# 8. Phi_ST VIA AMOVA
# ============================================================

compute_phi_st <- function(
  dnabin,
  groups,
  n_perm = 9999L
) {
  
  if (
    is.null(dnabin) ||
    length(groups) < 2L ||
    length(unique(groups)) < 2L
  ) {
    return(
      list(
        phi_st = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  d <- tryCatch(
    dist.dna(
      dnabin,
      model = "raw",
      pairwise.deletion = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(d)) {
    return(
      list(
        phi_st = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  df <- data.frame(
    group = factor(groups)
  )
  
  res <- tryCatch(
    pegas::amova(
      d ~ group,
      data = df,
      nperm = n_perm
    ),
    error = function(e) {
      message(
        "AMOVA failed: ",
        conditionMessage(e)
      )
      NULL
    }
  )
  
  if (is.null(res)) {
    return(
      list(
        phi_st = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  phi_st <- NA_real_
  p_value <- NA_real_
  
  try({
    
    vc <- res$varcomp
    
    sigma2 <- vc$sigma2
    names(sigma2) <- rownames(vc)
    
    total_sigma2 <- sum(
      sigma2,
      na.rm = TRUE
    )
    
    if (
      is.finite(total_sigma2) &&
      total_sigma2 != 0
    ) {
      phi_st <- sigma2[["group"]] / total_sigma2
    }
    
    p_value <- vc[
      rownames(vc) == "group",
      "P.value"
    ][1]
  })
  
  list(
    phi_st = safe_numeric(phi_st),
    p_value = safe_numeric(p_value)
  )
}


# ============================================================
# 9. HUDSON'S Snn WITH TIED NEAREST-NEIGHBOUR HANDLING
# ============================================================

compute_snn <- function(
  dnabin,
  groups,
  n_perm = 9999L
) {
  
  if (
    is.null(dnabin) ||
    length(groups) < 2L ||
    length(unique(groups)) < 2L
  ) {
    return(
      list(
        Snn = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  d <- tryCatch(
    as.matrix(
      dist.dna(
        dnabin,
        model = "raw",
        pairwise.deletion = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(d)) {
    return(
      list(
        Snn = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  diag(d) <- NA_real_
  
  snn_stat <- function(g) {
    
    per_sample <- vapply(
      seq_along(g),
      function(i) {
        
        row_d <- d[i, ]
        
        if (all(is.na(row_d))) {
          return(NA_real_)
        }
        
        min_d <- min(
          row_d,
          na.rm = TRUE
        )
        
        nn <- which(
          !is.na(row_d) &
            row_d == min_d
        )
        
        if (length(nn) == 0L) {
          return(NA_real_)
        }
        
        mean(
          g[nn] == g[i]
        )
      },
      numeric(1)
    )
    
    mean(
      per_sample,
      na.rm = TRUE
    )
  }
  
  observed <- snn_stat(groups)
  
  perm_stats <- replicate(
    n_perm,
    snn_stat(
      sample(groups)
    )
  )
  
  p_value <- (
    sum(
      perm_stats >= observed,
      na.rm = TRUE
    ) + 1
  ) / (
    sum(
      is.finite(perm_stats)
    ) + 1
  )
  
  list(
    Snn = observed,
    p_value = p_value
  )
}


# ============================================================
# 10. HAPLOTYPE SHARING
# ============================================================

compute_hap_sharing <- function(
  dnabin,
  groups
) {
  
  if (is.null(dnabin)) {
    return(NULL)
  }
  
  h <- tryCatch(
    haplotype(dnabin),
    error = function(e) NULL
  )
  
  if (is.null(h)) {
    return(NULL)
  }
  
  idx <- attr(h, "index")
  
  hap_groups <- map(
    idx,
    ~ unique(groups[.x])
  )
  
  n_total <- length(idx)
  
  n_native_only <- sum(
    map_lgl(
      hap_groups,
      ~ setequal(.x, "Native")
    )
  )
  
  n_introduced_only <- sum(
    map_lgl(
      hap_groups,
      ~ setequal(.x, "Introduced")
    )
  )
  
  n_shared <- n_total -
    n_native_only -
    n_introduced_only
  
  introduced_haplotypes <- which(
    map_lgl(
      hap_groups,
      ~ "Introduced" %in% .x
    )
  )
  
  n_introduced_haplotypes <- length(
    introduced_haplotypes
  )
  
  n_introduced_shared <- sum(
    map_lgl(
      hap_groups[introduced_haplotypes],
      ~ all(
        c(
          "Native",
          "Introduced"
        ) %in% .x
      )
    )
  )
  
  list(
    n_total = n_total,
    n_native_only = n_native_only,
    n_introduced_only = n_introduced_only,
    n_shared = n_shared,
    
    prop_haplotypes_shared = ifelse(
      n_total > 0L,
      n_shared / n_total,
      NA_real_
    ),
    
    prop_introduced_haplotypes_shared_with_native = ifelse(
      n_introduced_haplotypes > 0L,
      n_introduced_shared / n_introduced_haplotypes,
      NA_real_
    )
  )
}


# ============================================================
# 11. ANALYSE ONE SPECIES × GENOME × MARKER
# ============================================================

compute_pair <- function(row) {
  
  species <- row$taxa_accepted
  genome <- row$genome
  marker <- row$marker
  
  message(
    "Processing: ",
    species,
    " | ",
    genome,
    " | ",
    marker
  )
  
  all_dna <- read_dna_safely(
    row$filtered_all_file
  )
  
  if (is.null(all_dna)) {
    stop(
      "Could not read filtered joint alignment."
    )
  }
  
  mat <- as.matrix(all_dna)
  
  seq_names <- rownames(mat)
  
  groups <- parse_status_from_names(
    seq_names
  )
  
  if (any(is.na(groups))) {
    stop(
      "Could not recover Native/Introduced status from one or more FASTA headers."
    )
  }
  
  native_idx <- which(
    groups == "Native"
  )
  
  introduced_idx <- which(
    groups == "Introduced"
  )
  
  native_dna <- subset_dnabin_rows(
    all_dna,
    native_idx
  )
  
  introduced_dna <- subset_dnabin_rows(
    all_dna,
    introduced_idx
  )
  
  # ----------------------------------------------------------
  # A. Within-group statistics
  # ----------------------------------------------------------
  
  native_stats <- compute_diversity_stats(
    native_dna
  )
  
  introduced_stats <- compute_diversity_stats(
    introduced_dna
  )
  
  names(native_stats) <- paste0(
    "native_",
    names(native_stats)
  )
  
  names(introduced_stats) <- paste0(
    "introduced_",
    names(introduced_stats)
  )
  
  
  # ----------------------------------------------------------
  # B. Complete canonical sites across ALL samples
  # ----------------------------------------------------------
  #
  # Used for:
  #   * haplotype identity
  #   * haplotype sharing
  #   * Snn
  #   * Phi_ST
  #
  # This prevents N / ambiguity / gaps from creating spurious haplotypes.
  # ----------------------------------------------------------
  
  complete_joint <- complete_canonical_columns(
    all_dna
  )
  
  n_complete_joint_sites <- if (is.null(complete_joint)) {
    0L
  } else {
    ncol(as.matrix(complete_joint))
  }
  
  
  # ----------------------------------------------------------
  # C. Between-group differentiation
  # ----------------------------------------------------------
  
  between_stats <- tibble(
    n_complete_joint_sites = n_complete_joint_sites,
    phi_st = NA_real_,
    phi_st_pval = NA_real_,
    snn = NA_real_,
    snn_pval = NA_real_,
    n_haplotypes_combined = NA_integer_,
    n_haplotypes_native_only = NA_integer_,
    n_haplotypes_introduced_only = NA_integer_,
    n_haplotypes_shared = NA_integer_,
    prop_haplotypes_shared = NA_real_,
    prop_introduced_haplotypes_shared_with_native = NA_real_,
    native_rarefied_haplotype_richness = NA_real_,
    introduced_rarefied_haplotype_richness = NA_real_
  )
  
  if (
    !is.null(complete_joint) &&
    n_complete_joint_sites > 0L
  ) {
    
    phi <- compute_phi_st(
      complete_joint,
      groups,
      n_perm = n_perm
    )
    
    between_stats$phi_st <- phi$phi_st
    between_stats$phi_st_pval <- phi$p_value
    
    snn_res <- compute_snn(
      complete_joint,
      groups,
      n_perm = n_perm
    )
    
    between_stats$snn <- snn_res$Snn
    between_stats$snn_pval <- snn_res$p_value
    
    hap_shared <- tryCatch(
      compute_hap_sharing(
        complete_joint,
        groups
      ),
      error = function(e) NULL
    )
    
    if (!is.null(hap_shared)) {
      
      between_stats$n_haplotypes_combined <-
        hap_shared$n_total
      
      between_stats$n_haplotypes_native_only <-
        hap_shared$n_native_only
      
      between_stats$n_haplotypes_introduced_only <-
        hap_shared$n_introduced_only
      
      between_stats$n_haplotypes_shared <-
        hap_shared$n_shared
      
      between_stats$prop_haplotypes_shared <-
        hap_shared$prop_haplotypes_shared
      
      between_stats$prop_introduced_haplotypes_shared_with_native <-
        hap_shared$prop_introduced_haplotypes_shared_with_native
    }
    
    hap_ids <- haplotype_ids(
      complete_joint
    )
    
    if (!is.null(hap_ids)) {
      
      native_hap_ids <- hap_ids[native_idx]
      introduced_hap_ids <- hap_ids[introduced_idx]
      
      common_n <- min(
        length(native_hap_ids),
        length(introduced_hap_ids)
      )
      
      if (common_n >= 1L) {
        
        between_stats$native_rarefied_haplotype_richness <-
          rarefied_haplotype_richness(
            native_hap_ids,
            sample_size = common_n,
            B = haplotype_rarefaction_B
          )
        
        between_stats$introduced_rarefied_haplotype_richness <-
          rarefied_haplotype_richness(
            introduced_hap_ids,
            sample_size = common_n,
            B = haplotype_rarefaction_B
          )
      }
    }
  }
  
  
  # ----------------------------------------------------------
  # D. Pair-level provenance/QC from script 8
  # ----------------------------------------------------------
  
  pair_meta <- tibble(
    taxa_accepted = species,
    genome = genome,
    marker = marker,
    n_native = length(native_idx),
    n_introduced = length(introduced_idx),
    
    filtered_alignment_length =
      safe_numeric(row$filtered_alignment_length),
    
    raw_alignment_length =
      safe_numeric(row$raw_alignment_length),
    
    fraction_sites_retained =
      safe_numeric(row$fraction_sites_retained),
    
    filtered_all_file =
      row$filtered_all_file
  )
  
  
  bind_cols(
    pair_meta,
    native_stats,
    introduced_stats,
    between_stats
  )
}


# ============================================================
# 12. READ SCRIPT-8 MANIFEST
# ============================================================

if (!file_exists(alignment_manifest_file)) {
  stop(
    "Could not find script-8 alignment manifest:\n",
    alignment_manifest_file
  )
}

alignment_manifest <- read_tsv(
  alignment_manifest_file,
  show_col_types = FALSE
)


required_manifest_columns <- c(
  "taxa_accepted",
  "genome",
  "marker",
  "alignment_status",
  "filtered_all_file",
  "raw_alignment_length",
  "filtered_alignment_length",
  "fraction_sites_retained"
)

missing_columns <- setdiff(
  required_manifest_columns,
  names(alignment_manifest)
)

if (length(missing_columns) > 0L) {
  stop(
    "alignment_manifest.tsv is missing required columns: ",
    paste(
      missing_columns,
      collapse = ", "
    )
  )
}


pairs <- alignment_manifest %>%
  
  filter(
    alignment_status == "complete",
    !is.na(filtered_all_file),
    filtered_all_file != ""
  )


message(
  nrow(pairs),
  " species × genome × marker combinations available for diversity analysis."
)


# ============================================================
# 13. RUN ALL PAIRS SAFELY
# ============================================================

safe_compute_pair <- safely(
  compute_pair
)


pair_runs <- map(
  seq_len(nrow(pairs)),
  function(i) {
    
    res <- safe_compute_pair(
      pairs[i, , drop = FALSE]
    )
    
    list(
      row_index = i,
      result = res$result,
      error = res$error
    )
  }
)


pair_results <- map(
  pair_runs,
  "result"
) %>%
  compact()


diversity_results <- bind_rows(
  pair_results
)


analysis_failures <- map_dfr(
  pair_runs,
  function(x) {
    
    if (is.null(x$error)) {
      return(NULL)
    }
    
    i <- x$row_index
    
    tibble(
      taxa_accepted = pairs$taxa_accepted[i],
      genome = pairs$genome[i],
      marker = pairs$marker[i],
      error = conditionMessage(x$error)
    )
  }
)


write_tsv(
  analysis_failures,
  file.path(
    output_root,
    "analysis_failures.tsv"
  )
)


# ============================================================
# 14. DERIVED COMPARATIVE METRICS
# ============================================================

if (nrow(diversity_results) > 0L) {
  
  diversity_results <- diversity_results %>%
    
    mutate(
      delta_nuc_diversity =
        introduced_nuc_diversity_pairwise -
        native_nuc_diversity_pairwise,
      
      delta_nuc_diversity_complete =
        introduced_nuc_diversity_complete -
        native_nuc_diversity_complete,
      
      delta_hap_diversity =
        introduced_hap_diversity -
        native_hap_diversity,
      
      delta_rarefied_haplotype_richness =
        introduced_rarefied_haplotype_richness -
        native_rarefied_haplotype_richness,
      
      pi_ratio = case_when(
        is.na(native_nuc_diversity_pairwise) |
          is.na(introduced_nuc_diversity_pairwise) ~ NA_real_,
        
        native_nuc_diversity_pairwise == 0 &
          introduced_nuc_diversity_pairwise == 0 ~ 1,
        
        native_nuc_diversity_pairwise <= 0 ~ NA_real_,
        
        TRUE ~
          introduced_nuc_diversity_pairwise /
          native_nuc_diversity_pairwise
      ),
      
      log_pi_ratio = case_when(
        is.na(pi_ratio) ~ NA_real_,
        pi_ratio <= 0 ~ NA_real_,
        TRUE ~ log(pi_ratio)
      ),
      
      hap_diversity_ratio = case_when(
        is.na(native_hap_diversity) |
          is.na(introduced_hap_diversity) ~ NA_real_,
        
        native_hap_diversity == 0 &
          introduced_hap_diversity == 0 ~ 1,
        
        native_hap_diversity <= 0 ~ NA_real_,
        
        TRUE ~
          introduced_hap_diversity /
          native_hap_diversity
      ),
      
      diversity_lower_in_introduced =
        delta_nuc_diversity < 0,
      
      tajima_negative_sig_native =
        !is.na(native_tajima_D) &
        !is.na(native_tajima_D_pval) &
        native_tajima_D < 0 &
        native_tajima_D_pval < 0.05,
      
      tajima_negative_sig_introduced =
        !is.na(introduced_tajima_D) &
        !is.na(introduced_tajima_D_pval) &
        introduced_tajima_D < 0 &
        introduced_tajima_D_pval < 0.05,
      
      R2_sig_native =
        !is.na(native_R2_pval) &
        native_R2_pval < 0.05,
      
      R2_sig_introduced =
        !is.na(introduced_R2_pval) &
        introduced_R2_pval < 0.05,
      
      expansion_test_support_native =
        tajima_negative_sig_native |
        R2_sig_native,
      
      expansion_test_support_introduced =
        tajima_negative_sig_introduced |
        R2_sig_introduced
    )
}


# ============================================================
# 15. Phi_ST SINGLE-HAPLOTYPE EDGE CASE
# ============================================================

if (nrow(diversity_results) > 0L) {
  
  diversity_results <- diversity_results %>%
    
    mutate(
      phi_st_adj = if_else(
        is.na(phi_st) &
          n_haplotypes_combined == 1L,
        0,
        phi_st
      )
    )
}


# ============================================================
# 16. MULTIPLE-TESTING CORRECTION
# ============================================================

if (nrow(diversity_results) > 0L) {
  
  diversity_results <- diversity_results %>%
    
    mutate(
      phi_st_padj_BH = p.adjust(
        phi_st_pval,
        method = "BH"
      ),
      
      snn_padj_BH = p.adjust(
        snn_pval,
        method = "BH"
      ),
      
      native_tajima_D_padj_BH = p.adjust(
        native_tajima_D_pval,
        method = "BH"
      ),
      
      introduced_tajima_D_padj_BH = p.adjust(
        introduced_tajima_D_pval,
        method = "BH"
      ),
      
      native_R2_padj_BH = p.adjust(
        native_R2_pval,
        method = "BH"
      ),
      
      introduced_R2_padj_BH = p.adjust(
        introduced_R2_pval,
        method = "BH"
      )
    )
}


# ============================================================
# 17. BRING IN SCRIPT-8 ALIGNMENT QC
# ============================================================

alignment_qc <- tibble()

required_result_cols <- c(
  "taxa_accepted",
  "genome",
  "marker"
)

if (
  nrow(diversity_results) > 0L &&
  all(required_result_cols %in% names(diversity_results)) &&
  file_exists(alignment_qc_file)
) {
  
  alignment_qc <- read_tsv(
    alignment_qc_file,
    show_col_types = FALSE
  )
  
  required_qc_cols <- c(
    "taxa_accepted",
    "genome",
    "marker",
    "alignment_version",
    "gap_fraction",
    "ambiguous_fraction",
    "mean_site_coverage",
    "n_variable_sites",
    "n_parsimony_informative_sites"
  )
  
  if (
    all(required_qc_cols %in% names(alignment_qc))
  ) {
    
    alignment_qc_filtered <- alignment_qc %>%
      
      filter(
        str_detect(
          alignment_version,
          "^coverage_"
        )
      ) %>%
      
      select(
        taxa_accepted,
        genome,
        marker,
        
        alignment_gap_fraction =
          gap_fraction,
        
        alignment_ambiguous_fraction =
          ambiguous_fraction,
        
        alignment_mean_site_coverage =
          mean_site_coverage,
        
        alignment_n_variable_sites =
          n_variable_sites,
        
        alignment_n_parsimony_informative_sites =
          n_parsimony_informative_sites
      )
    
    diversity_results <- diversity_results %>%
      
      left_join(
        alignment_qc_filtered,
        by = c(
          "taxa_accepted",
          "genome",
          "marker"
        )
      )
    
  } else {
    
    warning(
      "alignment_qc_summary.tsv is missing required columns: ",
      paste(
        setdiff(
          required_qc_cols,
          names(alignment_qc)
        ),
        collapse = ", "
      )
    )
  }
  
} else {
  
  if (nrow(diversity_results) == 0L) {
    
    warning(
      "diversity_results is empty. ",
      "Check analysis_failures.tsv for upstream per-marker failures."
    )
    
  } else if (
    !all(required_result_cols %in% names(diversity_results))
  ) {
    
    warning(
      "diversity_results is missing expected columns: ",
      paste(
        setdiff(
          required_result_cols,
          names(diversity_results)
        ),
        collapse = ", "
      )
    )
  }
}


# ============================================================
# 18. QC FLAGS
# ============================================================

if (nrow(diversity_results) > 0L) {
  
  # Ensure optional alignment-QC columns exist even if script 8's
  # alignment_qc_summary.tsv was unavailable.
  if (!"alignment_gap_fraction" %in% names(diversity_results)) {
    diversity_results$alignment_gap_fraction <- NA_real_
  }
  
  if (!"alignment_ambiguous_fraction" %in% names(diversity_results)) {
    diversity_results$alignment_ambiguous_fraction <- NA_real_
  }
  
  if (!"alignment_mean_site_coverage" %in% names(diversity_results)) {
    diversity_results$alignment_mean_site_coverage <- NA_real_
  }
  
  marker_pi_ref <- diversity_results %>%
    
    select(
      marker,
      native_nuc_diversity_pairwise,
      introduced_nuc_diversity_pairwise
    ) %>%
    
    pivot_longer(
      cols = c(
        native_nuc_diversity_pairwise,
        introduced_nuc_diversity_pairwise
      ),
      values_to = "pi"
    ) %>%
    
    filter(
      !is.na(pi)
    ) %>%
    
    group_by(marker) %>%
    
    summarise(
      marker_median_pi = median(
        pi,
        na.rm = TRUE
      ),
      
      marker_mad_pi = mad(
        pi,
        constant = 1.4826,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  diversity_results <- diversity_results %>%
    
    left_join(
      marker_pi_ref,
      by = "marker"
    ) %>%
    
    mutate(
      native_prop_seg_sites = if_else(
        !is.na(native_n_seg_sites) &
          !is.na(native_alignment_sites) &
          native_alignment_sites > 0L,
        native_n_seg_sites /
          native_alignment_sites,
        NA_real_
      ),
      
      introduced_prop_seg_sites = if_else(
        !is.na(introduced_n_seg_sites) &
          !is.na(introduced_alignment_sites) &
          introduced_alignment_sites > 0L,
        introduced_n_seg_sites /
          introduced_alignment_sites,
        NA_real_
      ),
      
      native_modified_z = if_else(
        !is.na(marker_mad_pi) &
          marker_mad_pi > 0,
        0.6745 *
          (
            native_nuc_diversity_pairwise -
              marker_median_pi
          ) /
          marker_mad_pi,
        NA_real_
      ),
      
      introduced_modified_z = if_else(
        !is.na(marker_mad_pi) &
          marker_mad_pi > 0,
        0.6745 *
          (
            introduced_nuc_diversity_pairwise -
              marker_median_pi
          ) /
          marker_mad_pi,
        NA_real_
      ),
      
      qc_flag_pi_native =
        coalesce(
          native_nuc_diversity_pairwise >
            qc_pi_abs_threshold,
          FALSE
        ),
      
      qc_flag_pi_introduced =
        coalesce(
          introduced_nuc_diversity_pairwise >
            qc_pi_abs_threshold,
          FALSE
        ),
      
      qc_flag_seg_native =
        coalesce(
          native_prop_seg_sites >
            qc_prop_seg_threshold,
          FALSE
        ),
      
      qc_flag_seg_introduced =
        coalesce(
          introduced_prop_seg_sites >
            qc_prop_seg_threshold,
          FALSE
        ),
      
      qc_flag_marker_outlier_native =
        coalesce(
          native_modified_z >
            qc_mad_multiplier,
          FALSE
        ),
      
      qc_flag_marker_outlier_introduced =
        coalesce(
          introduced_modified_z >
            qc_mad_multiplier,
          FALSE
        ),
      
      qc_flag_gap_fraction =
        coalesce(
          alignment_gap_fraction >
            qc_gap_fraction_threshold,
          FALSE
        ),
      
      qc_flag_low_site_retention =
        coalesce(
          fraction_sites_retained <
            qc_min_fraction_sites_retained,
          FALSE
        ),
      
      qc_flag_no_complete_joint_sites =
        n_complete_joint_sites == 0L,
      
      qc_flag_pair =
        qc_flag_pi_native |
        qc_flag_pi_introduced |
        qc_flag_seg_native |
        qc_flag_seg_introduced |
        qc_flag_marker_outlier_native |
        qc_flag_marker_outlier_introduced |
        qc_flag_gap_fraction |
        qc_flag_low_site_retention |
        qc_flag_no_complete_joint_sites
    ) %>%
    
    mutate(
      qc_reason = pmap_chr(
        list(
          qc_flag_pi_native,
          qc_flag_pi_introduced,
          qc_flag_seg_native,
          qc_flag_seg_introduced,
          qc_flag_marker_outlier_native,
          qc_flag_marker_outlier_introduced,
          qc_flag_gap_fraction,
          qc_flag_low_site_retention,
          qc_flag_no_complete_joint_sites
        ),
        function(
          pi_n,
          pi_i,
          seg_n,
          seg_i,
          out_n,
          out_i,
          gap,
          retention,
          complete0
        ) {
          
          reasons <- c(
            if (pi_n) "high native pi" else NULL,
            if (pi_i) "high introduced pi" else NULL,
            if (seg_n) "high native segregating-site proportion" else NULL,
            if (seg_i) "high introduced segregating-site proportion" else NULL,
            if (out_n) "native marker-level pi outlier" else NULL,
            if (out_i) "introduced marker-level pi outlier" else NULL,
            if (gap) "high alignment gap fraction" else NULL,
            if (retention) "low fraction of sites retained after filtering" else NULL,
            if (complete0) "no complete canonical joint sites" else NULL
          )
          
          if (length(reasons) == 0L) {
            NA_character_
          } else {
            paste(
              reasons,
              collapse = "; "
            )
          }
        }
      )
    ) %>%
    
    select(
      -marker_median_pi,
      -marker_mad_pi
    )
}


# ============================================================
# 19. ANALYSIS-READY FLAG
# ============================================================
#
# Important:
# "analysis_ready" is conservative, but flagged rows are NOT deleted.
# ============================================================

if (nrow(diversity_results) > 0L) {
  
  diversity_results <- diversity_results %>%
    
    mutate(
      analysis_ready =
        !qc_flag_pair &
        n_native >= 5L &
        n_introduced >= 5L &
        n_complete_joint_sites > 0L
    )
}


# ============================================================
# 20. OUTPUT TABLES
# ============================================================

write_tsv(
  diversity_results,
  file.path(
    output_root,
    "diversity_stats_all.tsv"
  )
)


analysis_ready_results <- diversity_results %>%
  filter(
    analysis_ready
  )

write_tsv(
  analysis_ready_results,
  file.path(
    output_root,
    "diversity_stats_analysis_ready.tsv"
  )
)


qc_flags <- diversity_results %>%
  filter(
    qc_flag_pair
  )

write_tsv(
  qc_flags,
  file.path(
    output_root,
    "qc_flags.tsv"
  )
)


within_population_stats <- diversity_results %>%
  
  select(
    taxa_accepted,
    genome,
    marker,
    
    starts_with("native_"),
    starts_with("introduced_"),
    
    delta_nuc_diversity,
    delta_nuc_diversity_complete,
    delta_hap_diversity,
    delta_rarefied_haplotype_richness,
    pi_ratio,
    log_pi_ratio,
    hap_diversity_ratio
  )

write_tsv(
  within_population_stats,
  file.path(
    output_root,
    "within_population_stats.tsv"
  )
)


between_population_stats <- diversity_results %>%
  
  select(
    taxa_accepted,
    genome,
    marker,
    n_native,
    n_introduced,
    n_complete_joint_sites,
    
    phi_st,
    phi_st_adj,
    phi_st_pval,
    phi_st_padj_BH,
    
    snn,
    snn_pval,
    snn_padj_BH,
    
    n_haplotypes_combined,
    n_haplotypes_native_only,
    n_haplotypes_introduced_only,
    n_haplotypes_shared,
    
    prop_haplotypes_shared,
    prop_introduced_haplotypes_shared_with_native,
    
    native_rarefied_haplotype_richness,
    introduced_rarefied_haplotype_richness,
    delta_rarefied_haplotype_richness
  )

write_tsv(
  between_population_stats,
  file.path(
    output_root,
    "between_population_stats.tsv"
  )
)


neutrality_stats <- diversity_results %>%
  
  select(
    taxa_accepted,
    genome,
    marker,
    
    native_tajima_D,
    native_tajima_D_pval,
    native_tajima_D_padj_BH,
    native_fu_Fs,
    native_R2,
    native_R2_pval,
    native_R2_padj_BH,
    native_mismatch_raggedness,
    native_mean_pairwise_diff,
    expansion_test_support_native,
    
    introduced_tajima_D,
    introduced_tajima_D_pval,
    introduced_tajima_D_padj_BH,
    introduced_fu_Fs,
    introduced_R2,
    introduced_R2_pval,
    introduced_R2_padj_BH,
    introduced_mismatch_raggedness,
    introduced_mean_pairwise_diff,
    expansion_test_support_introduced
  )

write_tsv(
  neutrality_stats,
  file.path(
    output_root,
    "neutrality_stats.tsv"
  )
)


haplotype_stats <- diversity_results %>%
  
  select(
    taxa_accepted,
    genome,
    marker,
    
    native_n_haplotypes,
    native_hap_diversity,
    introduced_n_haplotypes,
    introduced_hap_diversity,
    
    n_haplotypes_combined,
    n_haplotypes_native_only,
    n_haplotypes_introduced_only,
    n_haplotypes_shared,
    
    prop_haplotypes_shared,
    prop_introduced_haplotypes_shared_with_native,
    
    native_rarefied_haplotype_richness,
    introduced_rarefied_haplotype_richness,
    delta_rarefied_haplotype_richness
  )

write_tsv(
  haplotype_stats,
  file.path(
    output_root,
    "haplotype_stats.tsv"
  )
)


multiple_testing_summary <- tibble(
  test = c(
    "Phi_ST",
    "Snn",
    "Native Tajima D",
    "Introduced Tajima D",
    "Native R2",
    "Introduced R2"
  ),
  
  n_raw_p_below_0_05 = c(
    sum(
      diversity_results$phi_st_pval < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$snn_pval < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$native_tajima_D_pval < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$introduced_tajima_D_pval < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$native_R2_pval < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$introduced_R2_pval < 0.05,
      na.rm = TRUE
    )
  ),
  
  n_BH_adjusted_p_below_0_05 = c(
    sum(
      diversity_results$phi_st_padj_BH < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$snn_padj_BH < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$native_tajima_D_padj_BH < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$introduced_tajima_D_padj_BH < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$native_R2_padj_BH < 0.05,
      na.rm = TRUE
    ),
    
    sum(
      diversity_results$introduced_R2_padj_BH < 0.05,
      na.rm = TRUE
    )
  )
)


write_tsv(
  multiple_testing_summary,
  file.path(
    output_root,
    "multiple_testing_summary.tsv"
  )
)


# ============================================================
# 21. QUICK SUMMARY
# ============================================================

message("")
message("Diversity analysis finished.")
message("")

message(
  "Pairs analysed successfully: ",
  nrow(diversity_results)
)

message(
  "Pairs that failed: ",
  nrow(analysis_failures)
)

message(
  "Pairs flagged by QC: ",
  nrow(qc_flags)
)

message(
  "Pairs marked analysis-ready: ",
  nrow(analysis_ready_results)
)


if (nrow(analysis_ready_results) > 0L) {
  
  summary_table <- analysis_ready_results %>%
    
    summarise(
      n_pairs = n(),
      
      mean_native_pi =
        mean(
          native_nuc_diversity_pairwise,
          na.rm = TRUE
        ),
      
      mean_introduced_pi =
        mean(
          introduced_nuc_diversity_pairwise,
          na.rm = TRUE
        ),
      
      n_diversity_lower_in_introduced =
        sum(
          diversity_lower_in_introduced,
          na.rm = TRUE
        ),
      
      n_significant_phi_st_BH =
        sum(
          phi_st_padj_BH < 0.05,
          na.rm = TRUE
        ),
      
      n_significant_snn_BH =
        sum(
          snn_padj_BH < 0.05,
          na.rm = TRUE
        ),
      
      n_native_expansion_test_support =
        sum(
          expansion_test_support_native,
          na.rm = TRUE
        ),
      
      n_introduced_expansion_test_support =
        sum(
          expansion_test_support_introduced,
          na.rm = TRUE
        )
    )
  
  print(
    summary_table
  )
}


message("")
message(
  "Output directory: ",
  path_abs(output_root)
)
