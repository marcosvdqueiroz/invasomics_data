# ============================================================
# GEOLOCATE georeferencing pipeline - chunked HPC version
# with timeout protection
# ------------------------------------------------------------

# ============================================================

library(readr)
library(dplyr)
library(tidyverse)
library(stringr)
library(RJSONIO)
library(RCurl)
library(tidyr)

# ============================================================
# Read command-line arguments
# ------------------------------------------------------------
# Required: START_ROW END_ROW CHUNK_ID
#   - start_row / end_row: the slice of the filtered ("tier 4")
#     table this particular HPC array task is responsible for.
#   - chunk_id: unique ID for this task, used to name its own
#     output file so parallel tasks never write to the same file.
# Optional (new): SLEEP_SECONDS MAX_ATTEMPTS
#   - sleep_seconds: base delay between GEOLocate requests
#     (default 3, matching the original hard-coded value).
#   - max_attempts: how many times to retry a single record's
#     request before giving up and logging it as failed (default 3).
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 4_geolocate_metacentrum.R START_ROW END_ROW CHUNK_ID [SLEEP_SECONDS] [MAX_ATTEMPTS]")
}
start_row <- as.integer(args[1])
end_row   <- as.integer(args[2])
chunk_id  <- as.integer(args[3])

# New, optional tuning knobs - fall back to the original defaults
# if the caller doesn't supply them, so old job scripts still work
# unchanged.
sleep_seconds <- if (length(args) >= 4) as.numeric(args[4]) else 3
max_attempts  <- if (length(args) >= 5) as.integer(args[5]) else 3

message("Running chunk: ", chunk_id)
message("Rows: ", start_row, " to ", end_row)
message("Base delay between requests: ", sleep_seconds,
        "s (+ up to 1s random jitter); max attempts per record: ", max_attempts)

# ------------------------------------------------------------
# IMPORTANT - cross-chunk rate limiting is NOT solved by this
# script alone. Sys.sleep() below only paces THIS process's own
# requests. If the HPC scheduler runs many chunks of this script
# at the same time (which is the whole point of chunking), the
# effective request rate hitting geo-locate.org is
# (sleep_seconds) / (number of chunks running concurrently),
# which can be far more aggressive than intended and risks the
# public service throttling or blocking the cluster's IP - which
# would show up here as a wave of records ending in
# "ERROR GETTING JSON" rather than as an obvious crash.
# Control the *aggregate* rate by limiting how many array tasks
# run concurrently (e.g. a PBS/SLURM array %-limit on
# simultaneous jobs), not just by editing sleep_seconds.
# ------------------------------------------------------------

# ============================================================
# Read input file
# ============================================================
nuc_summary <- read_tsv("nucleotide_entrez_summary_clean_merged.tsv")

# ============================================================
# Create full GEOLOCATE input table (i.e., the tier4 table)
# ------------------------------------------------------------
# Keep only records that:
#   - have no lat_lon yet (nothing to do if we already have coords)
#   - have a "country" field shaped like "Country: locality text",
#     which is the convention used for the free-text locality
#     description that GEOLocate will try to geocode.
# original_row is assigned BEFORE chunking/slicing below, so it's
# a stable ID across the *whole* dataset - this is what lets the
# resume logic (further down) match completed records back to
# input rows even though each chunk only sees a slice of the table.
# ============================================================
glcIn_all <- nuc_summary %>%
  filter(
    is.na(lat_lon),
    str_detect(country, "^[^:]+:")
  ) %>%
  mutate(
    tier = 4,
    country_only = str_extract(country, "^[^:]+"),
    locality = str_trim(str_extract(country, "(?<=:).*")),
    original_row = row_number()
  )

end_row <- min(end_row, nrow(glcIn_all))
if (start_row > nrow(glcIn_all)) {
  message("Chunk ", chunk_id, " starts after the last available row. Nothing to do.")
  quit(save = "no", status = 0)
}

glcIn <- glcIn_all[start_row:end_row, ]
message("Total rows after filtering: ", nrow(glcIn_all))
message("Rows in this chunk (before resume filtering): ", nrow(glcIn))

# ============================================================
# Helper functions
# ============================================================

# Coerce every column to character and blank out NAs, so the
# output TSV always has consistent, quote-safe, non-NA text
# columns regardless of the original column types.
clean_input_row <- function(x) {
  x %>%
    mutate(across(everything(), as.character)) %>%
    mutate(across(everything(), ~ tidyr::replace_na(.x, "")))
}

# Append one row to the chunk's output file. Header + column names
# are written only on the very first write of a genuinely new file
# (recordCounter == 1); every other write - including all writes
# during a resumed run, see below - appends without repeating the
# header.
write_geolocate_row <- function(df, output_file, recordCounter) {
  write.table(
    x = df,
    file = output_file,
    append = recordCounter != 1,
    row.names = FALSE,
    col.names = recordCounter == 1,
    quote = TRUE,
    sep = "\t",
    qmethod = "double"
  )
}

# ------------------------------------------------------------
# Retry wrapper around a single GEOLocate request.
#
# Makes up to `max_attempts` tries. Retries on: connection
# failures, non-200 HTTP status, an empty response body, and a
# response that parses as JSON but is missing the "numResults"
# field we depend on (a sign the service returned something
# other than a normal result payload, e.g. an error page or a
# rate-limit message). Uses exponential backoff (2s, 4s, 8s, ...)
# between attempts so repeated failures back off rather than
# hammering a struggling server harder.
#
# Returns the parsed JSON list on success, or throws (via stop())
# after the final failed attempt - the caller's tryCatch treats
# that as a genuine failure for this record, same as the original
# script did for any single failure.
# ------------------------------------------------------------
fetch_geolocate_response <- function(url, max_attempts = 3, base_backoff = 2) {
  last_error <- "unknown error"
  
  for (attempt in seq_len(max_attempts)) {
    attempt_result <- tryCatch({
      handle <- getCurlHandle()
      response_body <- basicTextGatherer()
      
      curlPerform(
        url = url,
        curl = handle,
        writefunction = response_body$update,
        timeout = 60,
        connecttimeout = 30,
        followlocation = TRUE
      )
      
      status <- getCurlInfo(handle)$response.code
      body_text <- response_body$value()
      
      if (is.null(status) || status != 200) {
        stop(paste0("unexpected HTTP status: ", status))
      }
      if (is.null(body_text) || nchar(body_text) == 0) {
        stop("empty response body")
      }
      
      parsed <- fromJSON(body_text)
      if (is.null(parsed$numResults)) {
        stop("response JSON missing 'numResults' field (unexpected payload shape)")
      }
      
      list(ok = TRUE, data = parsed)
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })
    
    if (isTRUE(attempt_result$ok)) {
      return(attempt_result$data)
    }
    
    last_error <- attempt_result$message
    if (attempt < max_attempts) {
      backoff <- base_backoff * (2 ^ (attempt - 1))
      message("  Attempt ", attempt, "/", max_attempts,
              " failed (", last_error, "). Retrying in ", backoff, "s...")
      Sys.sleep(backoff)
    }
  }
  
  stop(paste0("all ", max_attempts, " attempts failed. Last error: ", last_error))
}

# ============================================================
# Output setup
# ============================================================
output_file <- paste0("output_from_geolocate_", chunk_id, ".tsv")

# ------------------------------------------------------------
# Resume support.
#
# If this chunk's output file already exists - most likely
# because a previous run of this exact chunk got killed partway
# through (HPC walltime limit, node preemption, etc.) - read the
# glcRecNum values already written and skip those input rows this
# time, instead of reprocessing the whole chunk and duplicating
# work (and duplicating rows in the output file).
#
# Safety note: the very last glcRecNum found in the existing file
# is deliberately NOT treated as complete. If the process was
# killed mid-write of a record with multiple GEOLocate candidates,
# some but not all of that record's result rows might have made
# it to disk. Dropping the last one from the "done" set means it
# gets reprocessed - at worst a few duplicate rows for that one
# record, which is a much smaller problem than silently missing
# candidates for it.
# ------------------------------------------------------------
output_exists <- file.exists(output_file)
completed_ids <- integer(0)

if (output_exists) {
  message("Existing output file found for chunk ", chunk_id, " (", output_file, ") - attempting to resume.")
  
  existing <- tryCatch(
    read_tsv(
      output_file,
      col_types = cols(glcRecNum = col_integer(), .default = col_character())
    ),
    error = function(e) {
      message("  Could not parse existing output file (", conditionMessage(e), "). ",
              "Proceeding WITHOUT resume filtering - all rows in this chunk will be ",
              "reprocessed and appended, which may create duplicate lines in ", output_file,
              ". Consider de-duplicating that file afterwards if you see this message.")
      NULL
    }
  )
  
  if (!is.null(existing) && nrow(existing) > 0 && "glcRecNum" %in% names(existing)) {
    all_ids <- unique(existing$glcRecNum)
    all_ids <- all_ids[!is.na(all_ids)]
    
    if (length(all_ids) > 0) {
      last_id <- tail(existing$glcRecNum, 1)
      completed_ids <- setdiff(all_ids, last_id)
    }
    
    message("  Found ", length(completed_ids),
            " already-completed record(s); these will be skipped.")
  }
}

if (length(completed_ids) > 0) {
  glcIn <- glcIn[!(glcIn$original_row %in% completed_ids), ]
  message("Rows remaining in this chunk after resume filtering: ", nrow(glcIn))
}

if (nrow(glcIn) == 0) {
  message("Nothing left to do for chunk ", chunk_id, " (all rows already completed).")
}

# recordCounter drives write_geolocate_row()'s header/append logic.
# - Fresh output file: start at 0, so the first write (recordCounter
#   becomes 1) writes the header and creates the file.
# - Resumed output file: start at 1, so the first write this run
#   (recordCounter becomes 2) always appends without repeating the
#   header - matching write_geolocate_row()'s recordCounter == 1
#   check above.
recordCounter <- if (output_exists) 1 else 0

# ============================================================
# Main GEOLOCATE loop
# ============================================================
for (k in seq_len(nrow(glcIn))) {
  
  message("Processing row ", k, " of ", nrow(glcIn),
          " | chunk ", chunk_id)
  
  # Base delay between requests, plus a little random jitter.
  # The jitter mainly helps desynchronize multiple chunks that
  # were launched at (near) the same instant by the scheduler;
  # it does NOT by itself solve the aggregate cross-chunk rate
  # issue noted above.
  Sys.sleep(sleep_seconds + runif(1, 0, 1))
  
  country  <- glcIn[k, ]$country_only
  locality <- glcIn[k, ]$locality
  
  # Build the query string with proper per-value URL encoding
  # (curlEscape handles spaces, accented characters, "&", "#",
  # "%", apostrophes, etc.) instead of only replacing spaces.
  # This avoids silently truncated/malformed queries for any
  # locality text containing those characters.
  params <- ""
  
  if (!is.na(country) && country != "") {
    params <- paste(
      params,
      paste0("country=", curlEscape(country)),
      sep = "&"
    )
  }
  
  if (!is.na(locality) && locality != "") {
    params <- paste(
      params,
      paste0("locality=", curlEscape(locality)),
      sep = "&"
    )
  }
  
  # Drop the leading "&" left over from the paste() calls above.
  params <- substr(params, 2, nchar(params))
  
  q <- paste0(
    "http://geo-locate.org/webservices/geolocatesvcv2/glcwrap.aspx?",
    params
  )
  # No further gsub() needed here - curlEscape() above already
  # encoded every parameter value correctly.
  
  tryCatch({
    
    glcRecNum <- glcIn[k, ]$original_row
    
    # Request the record with retries/backoff instead of a single
    # unprotected attempt. Throws (caught by the outer tryCatch,
    # same as the original script) only after max_attempts failures.
    glc <- fetch_geolocate_response(q, max_attempts = max_attempts)
    numresults <- glc$numResults
    
    if (numresults > 0) {
      
      for (i in seq_len(numresults)) {
        
        glcRank <- i
        
        glcLongitude <- glc$resultSet$features[[i]]$geometry$coordinates[1]
        glcLatitude  <- glc$resultSet$features[[i]]$geometry$coordinates[2]
        
        glcPrecision <- glc$resultSet$features[[i]]$properties$precision
        glcScore <- glc$resultSet$features[[i]]$properties$score
        glcParsepattern <- glc$resultSet$features[[i]]$properties$parsePattern
        glcUncert <- glc$resultSet$features[[i]]$properties$uncertaintyRadiusMeters
        glcPoly <- glc$resultSet$features[[i]]$properties$uncertaintyPolygon
        
        if ("coordinates" %in% names(glcPoly)) {
          
          # Flatten GEOLocate's nested polygon geometry into a
          # single "lat,lon,lat,lon,..." string so it fits in a
          # flat TSV column.
          sPoly <- ""
          
          for (v in seq_along(glcPoly$coordinates[[1]])) {
            
            vLon <- format(glcPoly$coordinates[[1]][[v]][1])
            vLat <- format(glcPoly$coordinates[[1]][[v]][2])
            
            sPoly <- paste(sPoly, vLat, vLon, sep = ",")
          }
          
          sPoly <- sub("^,+", "", sPoly)
          glcPoly <- sPoly
        }
        
        input_row <- clean_input_row(glcIn[k, ])
        
        # One output row per GEOLocate candidate, combining its
        # geocoding fields with the original input record.
        df <- data.frame(
          glcRecNum = glcRecNum,
          glcRank = glcRank,
          glcLatitude = glcLatitude,
          glcLongitude = glcLongitude,
          glcPrecision = glcPrecision,
          glcScore = glcScore,
          glcParsepattern = glcParsepattern,
          glcUncert = glcUncert,
          glcPoly = glcPoly,
          input_row,
          check.names = FALSE
        )
        
        recordCounter <- recordCounter + 1
        write_geolocate_row(df, output_file, recordCounter)
      }
      
    } else {
      
      # No candidates returned - still write one row so this
      # record is traceable as "attempted, no match" rather than
      # silently missing from the output.
      input_row <- clean_input_row(glcIn[k, ])
      
      df <- data.frame(
        glcRecNum = glcRecNum,
        glcRank = 1,
        glcLatitude = NA,
        glcLongitude = NA,
        glcPrecision = NA,
        glcScore = NA,
        glcParsepattern = NA,
        glcUncert = NA,
        glcPoly = NA,
        input_row,
        check.names = FALSE
      )
      
      recordCounter <- recordCounter + 1
      write_geolocate_row(df, output_file, recordCounter)
    }
    
  }, error = function(err) {
    
    # Reached only after fetch_geolocate_response() has already
    # exhausted all retries (or some other unexpected error
    # occurred). Log a row with glcRank = 0 and the error message
    # in glcPrecision, same convention as the original script, so
    # failed records are visible and filterable in the output.
    glcRecNum <- glcIn[k, ]$original_row
    
    input_row <- clean_input_row(glcIn[k, ])
    
    df <- data.frame(
      glcRecNum = glcRecNum,
      glcRank = 0,
      glcLatitude = NA,
      glcLongitude = NA,
      glcPrecision = paste("ERROR GETTING JSON:", conditionMessage(err)),
      glcScore = 0,
      glcParsepattern = NA,
      glcUncert = NA,
      glcPoly = NA,
      input_row,
      check.names = FALSE
    )
    
    # <<- because this handler is a closure and recordCounter
    # lives in the enclosing loop's scope.
    recordCounter <<- recordCounter + 1
    write_geolocate_row(df, output_file, recordCounter)
  })
}

message("Finished chunk: ", chunk_id)
message("Output written to: ", output_file)