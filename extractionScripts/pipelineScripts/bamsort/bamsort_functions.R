# bamsort_functions.R
# The pipeline-side R library, and the single definition of the functions below.
# It stands alone: the bamsort workers source it via FUNCTIONS_R with no config.R
# and no plotting stack. plotScripts/functions.R, the larger analysis and
# plotting library the .Rmd notebooks source, sources this file rather than
# carrying its own copies.
# R functions for the Ca*Re bamsort pipeline: turning per-region alignment
# tables (stepTwo) into per-UMI category assignments (stepThree) and then into
# per-cell count matrices at a chosen detection threshold (analysis).
#
# Sourced by bamsort_process.R and bamsort_process_host.R via the FUNCTIONS_R
# environment variable, and by downstream analysis scripts that call
# set_minimums() directly.
#
# Requires: data.table, dplyr, tidyr
#
# Viral-target pipeline (SIV mac239 / HIV HXB2):
#   extract_category       category of an alignment CSV, from its filename
#   apply_umi_corrections  reassign categories using confirmed splice junctions
#   resolve_multimap       collapse UMIs assigned to more than one category
#   process_bamsort        stepThree entry point: combine regions, assign a
#                          category per UMI, filter to raw cellIDs, write
#                          <sample>_UMI_read_counts_full.csv
#
# Host-target pipeline (no splice-category resolution):
#   extract_category_host  target of an alignment CSV, from its filename
#   process_bamsort_host   stepThree entry point for host targets
#
# Detection thresholds (called from analysis scripts, not from the pipeline, so
# thresholds can be varied without re-running stepThree):
#   set_minimums               one sample: UMI table -> per-cell counts
#   process_all_set_minimums   set_minimums across a list of samples

# =============================================================================
# Viral targets
# =============================================================================
## CaRe and TenX: extract category from filename
extract_category <- function(filename) {
  # Split the filename into parts
  target_part <- basename(filename)
  target_part <- strsplit(target_part, "_bamsort_alignment_")[[1]][1]
  alignment_part <- strsplit(filename, "_bamsort_alignment_")[[1]][2]
  alignment_part <- sub("\\.csv$", "", alignment_part)
  
  # Probe targets carried in the sample name, e.g. "A_env" -> "env". Single-cell
  # samples carry none ("HIV-A"), leaving this empty: every "%in% targets" test
  # below is then FALSE and the region alone decides the category, which is the
  # correct reading for a TenX library. Only a filename that does not parse at
  # all is uncategorizable.
  targets <- strsplit(target_part, "_")[[1]]
  targets <- targets[-1]
  targets <- targets[!targets %in% c("PL", "CI")]
  if (is.na(target_part) || is.na(alignment_part)) {
    return(NA)
  }

  # universal truths
  if (grepl("KLRB1$", alignment_part)) return("KLRB1")
  if (length(targets) == 1 & "pol" %in% targets) return("US")
  if (grepl("LTR_D1$", alignment_part)) return("any-D1")
  if (grepl("D1_A1$", alignment_part)) return("US")
  
  
  # treat TenX like CaRe env/SSenv/pol
  if (alignment_part == "A1_D4") {
    if ("D1" %in% targets) return("SS-MS4")
    else return("any-D4")
  }
  if (alignment_part == "D4_A7") {
    if ("D1" %in% targets) return("SS")
    else return("US-SS")
  }
  if (alignment_part == "A7_LTR") {
    if ("D1" %in% targets) return("SS-MS7")
    if ("SSenv" %in% targets) return("US-SS")
    if ("tat" %in% targets | "D4" %in% targets) return("MS")
    else return("any-A7")
  }
  
  # Default case
  cat("  Failed to categorize. Assigning as: NA\n")
  return(NA)
}
## CaRe: spliced reads join 
## Load correction CSV and apply category overrides
apply_umi_corrections <- function(umi_read_counts, correction_file, target_category = "MS") {
  # Check if correction file exists
  if (!file.exists(correction_file)) {
    cat("No correction file found at:", correction_file, "\n")
    return(umi_read_counts)
  }
  
  allowed_categories <- c("SS-MS4", "SS-MS7", "any-A7", "any-D4")
  # Load correction data
  corrections <- read.csv(correction_file)
  
  # Validate required columns
  required_cols <- c("CB", "UB")
  if (!all(required_cols %in% colnames(corrections))) {
    cat("Warning: Correction file missing required columns (CB, UB)\n")
    return(umi_read_counts)
  }
  
  # Create composite key for matching
  corrections$key <- paste(corrections$CB, corrections$UB, sep = "_")
  umi_read_counts$key <- paste(umi_read_counts$cellID, umi_read_counts$UMI, sep = "_")
  
  in_correction_file <- umi_read_counts$key %in% corrections$key
  
  # Find matches with allowed categories
  to_correct <- in_correction_file & 
    umi_read_counts$category %in% allowed_categories
  to_toss <- in_correction_file & 
    !umi_read_counts$category %in% allowed_categories
  
  n_corrected <- sum(to_correct)
  n_tossed <- sum(to_toss)
  
  
  if (n_corrected > 0) {
    # Store original categories for reporting
    original_categories <- umi_read_counts$category[to_correct]
    
    cat("Corrected to", target_category, n_corrected, "UMI from:", 
        paste(unique(original_categories), collapse = ", "), "\n")
  } else {
    cat("No matching UMIs found in correction file\n")
  }
  
  if (n_tossed > 0) {
    tossed_categories <- umi_read_counts$category[to_toss]
    cat("Tossed", n_tossed, "UMIs with non-allowed categories:", 
        paste(unique(tossed_categories), collapse = ", "), "\n")
  }
  
  # Apply correction
  umi_read_counts$category[to_correct] <- target_category
  umi_read_counts <- umi_read_counts[!to_toss, ]
  
  # Remove temporary key column
  umi_read_counts$key <- NULL
  
  return(umi_read_counts)
}
## CaRe and TenX: resolve UMI mapping to multiple genes/regions
resolve_multimap <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  valid_pairs <- list(
    # keep the second entry
    c("any-D1", "US"),
    c("any-D4", "US"),
    c("SS-MS4", "US"),
    c("any-D4", "US-SS"), 
    c("any-A7", "US-SS"), 
    c("MS", "US-SS"), 
    c("SS-MS4", "SS"),
    c("SS-MS7", "SS"),
    c("any-D4", "SS-MS4"),
    c("any-A7", "SS-MS7"),
    c("US-SS", "SS"),
    c("SS-MS7", "MS"),
    c("any-A7", "MS")
  )
  
  # Create a composite key
  df$key <- paste(df$cellID, df$UMI, sep = "_")
  
  # Identify multi-category keys
  key_counts <- table(df$key)
  multi_keys <- names(key_counts[key_counts > 1])
  
  # Statistics counters
  total_multimapped <- length(multi_keys)
  resolved_count <- 0
  tossed_count <- 0
  tossed_more_than_2_count <- 0
  consolidated_same_category <- 0 
  
  # Split into single and multi entries
  singles <- df[!df$key %in% multi_keys, ]
  multis <- df[df$key %in% multi_keys, ]
  
  # Initialize result as empty data frame with correct structure
  result <- data.frame(
    cellID = character(0),
    UMI = character(0),
    category = character(0),
    read_count = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Add singles to result
  if (nrow(singles) > 0) {
    singles_clean <- data.frame(
      cellID = as.character(singles$cellID),
      UMI = as.character(singles$UMI),
      category = as.character(singles$category),
      read_count = as.numeric(singles$read_count),
      stringsAsFactors = FALSE
    )
    result <- rbind(result, singles_clean)
  }
  
  # Process each unique key in multis
  if (nrow(multis) > 0) {
    unique_keys <- unique(multis$key)
    
    for (k in unique_keys) {
      # Get all rows for this key
      rows <- multis[multis$key == k, ]
      categories <- unique(rows$category)
      
      # Process based on number of categories
      if (length(categories) == 1) {
        # For keys with 1 category but multiple rows, consolidate by summing read_counts
        total_reads <- sum(as.numeric(rows$read_count))
        keep_row <- data.frame(
          cellID = as.character(rows$cellID[1]),
          UMI = as.character(rows$UMI[1]),
          category = as.character(rows$category[1]),
          read_count = total_reads,
          stringsAsFactors = FALSE
        )
        result <- rbind(result, keep_row)
        consolidated_same_category <- consolidated_same_category + 1
      } else if (length(categories) == 2) {
        # Check against valid pairs
        valid_match <- FALSE
        for (pair in valid_pairs) {
          if (all(sort(categories) == sort(pair))) {
            # Valid pair, keep the second category
            keep_cat <- pair[2]
            keep_rows <- rows[rows$category == keep_cat, ]
            keep_row <- data.frame(
              cellID = as.character(keep_rows$cellID[1]),
              UMI = as.character(keep_rows$UMI[1]),
              category = as.character(keep_cat),
              read_count = sum(as.numeric(rows$read_count)),
              stringsAsFactors = FALSE
            )
            result <- rbind(result, keep_row)
            resolved_count <- resolved_count + 1
            valid_match <- TRUE
            break
          }
        }
        
        # If no valid match found, count as tossed
        if (!valid_match) {
          tossed_count <- tossed_count + 1
        }
      } else if (length(categories) > 2) {
        # Toss UMIs with more than 2 categories
        tossed_more_than_2_count <- tossed_more_than_2_count + 1
      }
    }
  }
  
  # Print statistics
  cat("Total multi-mapped UMIs:", total_multimapped, "\n")
  cat("UMIs successfully resolved:", resolved_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, resolved_count/total_multimapped*100, 0)), "\n")
  cat("UMIs tossed (no valid pair):", tossed_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, tossed_count/total_multimapped*100, 0)), "\n")
  cat("UMIs tossed (>2 categories):", tossed_more_than_2_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, tossed_more_than_2_count/total_multimapped*100, 0)), "\n")
  cat("Total UMIs tossed:", tossed_count + tossed_more_than_2_count, 
      sprintf("(%.1f%%)", ifelse(total_multimapped > 0, 
                                 (tossed_count + tossed_more_than_2_count)/total_multimapped*100, 0)), "\n")
  
  # Clean up row names
  rownames(result) <- NULL
  
  return(result)
}
## CaRe and TenX: read and process BAM
# Read support for a bamsort alignment table. Adjacent regions mean one read can
# clear min_overlap in two neighbouring windows and be written to both files, so
# summing `count` over regions counts it twice. `primary_count` charges each read
# to the single region it overlaps most, which makes the sum a distinct-read
# count. Tables written before that column existed fall back to `count`.
bamsort_read_column <- function(file_data) {
  if ("primary_count" %in% colnames(file_data)) "primary_count" else "count"
}

process_bamsort <- function(samples_list, input.dir, output.dir, raw_cellIDs, 
                              correction_dir = NULL) {
  #
  # optional to read in output from bamsort_splice_sites to update MS
  #
  if (class(samples_list) != "character") {
    stop("samples_list input must be a list of sample names")
  }
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  for (sample in samples_list) {
    cat("\nProcessing sample:", sample, "\n")
    
    # Get all files for current sample
    sample_files <- unlist(lapply(input.dir, function(dir) {
      list.files(dir, 
                 pattern = paste0(".*", sample, ".*_bamsort_alignment_.*\\.csv$"),
                 full.names = TRUE)
    }))
    
    if (length(sample_files) == 0) {
      cat("No files found for sample:", sample, "\n")
      next
    }
    
    # Extract categories for each file
    categories <- sapply(sample_files, extract_category)
    sample_files <- sample_files[!is.na(categories)]
    
    # Create a data frame to store the combined results
    all_data <- data.frame()
    filenames <- basename(sample_files)
    
    # Process each relevant file
    for (j in seq_along(sample_files)) {
      extracted_data <- NULL
      
      tryCatch({
        file_data <- fread(sample_files[j], data.table = FALSE)
        
        if(nrow(file_data) == 0) {
          cat("Warning: File is empty:", filenames[j], "\n")
          next
        }
        
        required_cols <- c("cellID", "UMI", "count")
        missing_cols <- setdiff(required_cols, colnames(file_data))
        
        if (length(missing_cols) > 0) {
          cat("Warning: Missing required columns:", paste(missing_cols, collapse=", "), "\n")
          cat("Available columns:", paste(colnames(file_data), collapse=", "), "\n")
          next
        }
        
        extracted_data <- file_data %>%
          select(cellID, UMI, read = dplyr::all_of(bamsort_read_column(file_data))) %>%
          mutate(category = categories[j])
        
      }, error = function(e) {
        cat("Error reading file:", filenames[j], "\nLikely bad data file\n")
        cat("Error message:", conditionMessage(e), "\n")
      })
      
      if (!is.null(extracted_data) && nrow(extracted_data) > 0) {
        all_data <- rbind(all_data, extracted_data)
      }
    }
    
    # If we have data, process it
    if (nrow(all_data) > 0) {
      all_data <- all_data %>% 
        filter(!is.na(cellID) & !is.na(UMI) & !is.na(read)) %>%
        unique()
      
      tryCatch({
        umi_read_counts <- all_data %>%
          group_by(cellID, UMI, category) %>%
          summarize(
            read_count = sum(as.numeric(read)),
            .groups = 'drop'
          )
        
        cat("Gathered read counts for", nrow(umi_read_counts), "UMIs from", 
            n_distinct(umi_read_counts$cellID), "cells\n")
        
        # Create a wide format with UMI counts per category
        umi_by_category <- umi_read_counts %>%
          group_by(cellID, category) %>%
          summarize(
            category_UMIs = n_distinct(UMI),
            category_reads = sum(read_count),
            .groups = 'drop'
          )
        
        umi_by_category_wide <- tidyr::pivot_wider(
          umi_by_category,
          id_cols = cellID,
          names_from = category,
          names_sep = "_",
          values_from = c(category_UMIs, category_reads),
          values_fill = 0)
        
        # Apply raw_cellIDs filtering
        if (sample %in% names(raw_cellIDs) && !is.null(raw_cellIDs[[sample]])) {
          valid_cells <- raw_cellIDs[[sample]]
          cat("Keeping only cells in", 
              length(valid_cells), "valid cells from 10X data\n")
          
          valid_umi_read_counts <- umi_read_counts %>%
            filter(cellID %in% valid_cells)
          
          # Clean category values
          if (!is.null(names(valid_umi_read_counts$category))) {
            category_values <- as.character(valid_umi_read_counts$category)
            category_values <- gsub('^"(.*)"$', '\\1', category_values)
            valid_umi_read_counts$category <- category_values
          } else {
            valid_umi_read_counts$category <- gsub('^"(.*)"$', '\\1', as.character(valid_umi_read_counts$category))
          }
          
          # Look for splice correction file
          if (!is.null(correction_dir)) {
            correction_file <- file.path(correction_dir, paste0(sample, "_splicesites_filtered.csv"))
            valid_umi_read_counts <- apply_umi_corrections(
              valid_umi_read_counts, 
              correction_file, 
              target_category = "MS"
            )
          }
          
          # Multimap resolution (with corrected categories)
          filtered_umi_read_counts <- resolve_multimap(valid_umi_read_counts)
          
          # Convert to basic types
          filtered_umi_read_counts <- data.frame(
            cellID = as.character(filtered_umi_read_counts$cellID),
            UMI = as.character(filtered_umi_read_counts$UMI),
            category = as.character(filtered_umi_read_counts$category),
            read_count = as.numeric(filtered_umi_read_counts$read_count),
            stringsAsFactors = FALSE
          )
          
          # Consolidate any categories
          filtered_umi_read_counts <- filtered_umi_read_counts %>%
            mutate(category = if_else(grepl("^any", category), "any", category)) %>%
            mutate(category = if_else(grepl("^SS-MS", category), "SS-MS", category))
          
          cat("After filtering: kept", nrow(filtered_umi_read_counts), "UMIs from", 
              n_distinct(filtered_umi_read_counts$cellID), "valid cells\n")
        } else {
          cat("Warning: No 10X data found for sample", sample, "- using unfiltered cell list\n")
          filtered_umi_read_counts <- umi_read_counts
        }
        
        # Convert to basic types
        all_data_clean <- data.frame(
          cellID = as.character(all_data$cellID),
          UMI = as.character(all_data$UMI),
          read = as.numeric(all_data$read),
          category = as.character(all_data$category),
          stringsAsFactors = FALSE
        )
        
        umi_read_counts_clean <- data.frame(
          cellID = as.character(umi_read_counts$cellID),
          UMI = as.character(umi_read_counts$UMI),
          category = as.character(umi_read_counts$category),
          read_count = as.numeric(umi_read_counts$read_count),
          stringsAsFactors = FALSE
        )
        
        # Save outputs
        umi_level_file <- file.path(output.dir, paste0(sample, "_UMI_read_counts_raw.csv"))
        write.csv(umi_read_counts_clean, file = umi_level_file, row.names = FALSE)
        
        filtered_umi_level_file <- file.path(output.dir, paste0(sample, "_UMI_read_counts_full.csv"))
        write.csv(filtered_umi_read_counts, file = filtered_umi_level_file, row.names = FALSE)
        
        cat("Successfully processed\n")
        rm(umi_by_category, umi_by_category_wide, umi_read_counts)
      }, error = function(e) {
        cat("Error processing data for sample", sample, ":", conditionMessage(e), "\n")
        
        tryCatch({
          all_data_clean <- data.frame(
            cellID = as.character(all_data$cellID),
            UMI = as.character(all_data$UMI),
            read = as.numeric(all_data$read),
            category = as.character(all_data$category),
            stringsAsFactors = FALSE
          )
          
          raw_file <- file.path(output.dir, paste0(sample, "_raw_data.csv"))
          write.csv(all_data_clean, file = raw_file, row.names = FALSE)
          cat("Saved raw data to:", raw_file, "\n")
        }, error = function(e2) {
          cat("Could not save raw data:", conditionMessage(e2), "\n")
        })
      })
    } else {
      cat("No data extracted for sample:", sample, "\n")
    }
  }
}
## CaRe and TenX: filtering data for minimums
set_minimums <- function(mode, umi_file, min_reads = 1, min_region_count = 1, 
                           min_umi = 1, min_cells = 0, min_reads_cell = 1,
                           save_valid_umi = TRUE, valid_umi_output_file = NULL,
                           categories = NULL) {
  
  # Define categories based on mode (if not provided explicitly)
  if (is.null(categories)) {
    categories <- switch(mode,
                         "CaRe" = c("US", "US-SS", "SS-MS", "SS", "MS", "any"),
                         # Same six as CaRe, so the two arms are directly
                         # comparable and can be bound without a column
                         # mismatch. SS and SS-MS are currently unreachable for
                         # a TenX library (extract_category needs a target primer
                         # in the filename to resolve them) and MS arrives only
                         # via the splice correction, but an unreachable category
                         # comes back zero-filled rather than absent.
                         "GEX" = c("US", "US-SS", "SS-MS", "SS", "MS", "any"),
                         "KLRB1" = c("KLRB1"),
                         "host" = NULL,
                         stop("Invalid mode. Must be 'CaRe', 'GEX', 'KLRB1', or 'host'")
    )
  }
  
  # Helper function to return empty result. Defined before its first use below,
  # which an all-empty UMI table reaches.
  return_empty <- function(message = "", umi_entries = 0) {
    if (nzchar(message)) cat(message, "\n")
    cats <- if (is.null(categories)) character(0) else categories
    empty_data <- data.frame(cellID = character(0), stringsAsFactors = FALSE)
    empty_data[cats] <- lapply(cats, function(x) numeric(0))
    empty_data$total <- numeric(0)
    cat("Converted", umi_entries, "UMI entries to 0 cells\n")
    return(empty_data)
  }
  
  umi_data <- read.csv(umi_file, stringsAsFactors = FALSE)
  if (nrow(umi_data) == 0) {
    return(return_empty("No UMI data found in file"))
  }
  
  # For host mode (or any mode with NULL categories), infer from data
  if (is.null(categories)) {
    categories <- sort(unique(umi_data$category))
    cat("Inferred categories from data:", paste(categories, collapse = ", "), "\n")
  }
  
  
  original_entries <- nrow(umi_data)
  
  # Apply read count filter and check for empty result
  if (min_reads > 0) {
    umi_data <- umi_data[umi_data$read_count >= min_reads, ]
    if (nrow(umi_data) == 0) {
      return(return_empty("All UMI data filtered out by min_reads threshold", original_entries))
    }
  }
  
  # Calculate cell-level metrics
  cell_stats <- umi_data %>%
    group_by(cellID) %>%
    summarize(
      total_reads = sum(read_count, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Apply reads per cell filter
  if (min_reads_cell > 0) {
    valid_cells <- cell_stats$cellID[cell_stats$total_reads >= min_reads_cell]
    umi_data <- umi_data[umi_data$cellID %in% valid_cells, ]
    if (nrow(umi_data) == 0) {
      return(return_empty("All cells filtered out by min_reads_cell threshold", original_entries))
    }
  }
  
  # Create UMI count matrix
  umi_counts <- umi_data %>%
    count(cellID, category, name = "umi_count") %>%
    pivot_wider(names_from = category, values_from = umi_count, values_fill = 0)
  
  # Initialize result with all categories (including missing ones as 0)
  result <- umi_counts %>%
    select(cellID, any_of(categories)) %>%
    mutate(across(-cellID, ~replace_na(.x, 0)))

  missing_cats <- setdiff(categories, colnames(result))
  if (length(missing_cats) > 0) {
    result[missing_cats] <- 0
  }
  
  result <- result[c("cellID", categories)]
  
  # Apply minimum cells per region filter
  if (min_cells > 0) {
    for (cat in categories) {
      if (sum(result[[cat]] > 0) < min_cells) {
        result[[cat]] <- 0
      }
    }
  }
  
  # Calculate metrics for filtering. `total` counts every category, so for a
  # target with a skipping read-out min_umi is a coverage gate: the cell needs
  # min_umi molecules that report on the region, whichever isoform they carry.
  result <- result %>%
    mutate(
      region_count = rowSums(across(all_of(categories), ~ .x > 0)),
      total = rowSums(across(all_of(categories)))
    )
  
  # Apply filters
  result <- result %>%
    filter(
      region_count >= min_region_count,
      total >= min_umi
    ) %>%
    select(-region_count)
  
  # Check for empty result
  if (nrow(result) == 0) {
    return(return_empty("All cells filtered out by region_count or UMI thresholds", original_entries))
  }
  
  # Save valid cellID-UMI pairs if requested
  if (save_valid_umi && !is.null(valid_umi_output_file) && nrow(result) > 0) {
    valid_cellIDs <- result$cellID
    valid_umi_data <- umi_data %>%
      filter(cellID %in% valid_cellIDs) %>%
      select(cellID, UMI)
    write.csv(valid_umi_data, file = valid_umi_output_file, row.names = FALSE)
    cat("Saved", nrow(valid_umi_data), "valid cellID-UMI pairs\n")
  }
  
  cat("Converted", original_entries, "UMI entries to", nrow(result), "cells\n")
  return(result)
}
## CaRe and TenX: filter multiple samples
process_all_set_minimums <- function(mode, sample_list, input_dir, save_output = TRUE, 
                                       output_dir = "/projects/b1042/GoyalLab/egrody/extractedData/",
                                       min_reads = 1, min_region_count = 1, min_umi = 1, 
                                       min_cells = 0, min_reads_cell = 1,
                                       save_valid_umi = TRUE, valid_umi_dir = NULL,
                                       aggregate_valid_umi = TRUE) {
  
  # Create output directory if needed
  if (save_output == TRUE) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Create separate directory for valid UMI files if needed
  # Use the provided output_dir as the base for valid_umi_dir
  if (save_valid_umi == TRUE) {
    if (is.null(valid_umi_dir)) {
      valid_umi_dir <- file.path(output_dir, "valid_UMI")
    }
    dir.create(valid_umi_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Write filters to file. Only when save_output is TRUE: output_dir is not
  # created otherwise, so writing here would fail on a path that does not exist.
  if (save_output == TRUE) {
    filters_file <- file.path(output_dir, "filters.txt")
    filter_params <- c(
      paste("min_reads =", min_reads),
      paste("min_region_count =", min_region_count),
      paste("min_umi =", min_umi),
      paste("min_cells =", min_cells),
      paste("min_reads_cell =", min_reads_cell),
      paste("mode =", mode),
      paste("processed_date =", Sys.Date())
    )
    writeLines(filter_params, filters_file)
  }
  
  # Define output suffix based on mode
  suffix <- paste0("_", mode, "_filtered.csv")
  umi_suffix <- paste0("_", mode, "_valid_cellID_UMI.csv")
  
  # Initialize list to collect valid UMI data if aggregating
  all_valid_umi <- if (aggregate_valid_umi && save_valid_umi) list() else NULL
  
  # Collect the per-sample filtered tables so they can be returned when the
  # caller is not writing them to disk
  all_results <- list()
  
  # Process samples with better error handling
  results <- sapply(sample_list, function(sample) {
    cat("\nProcessing sample:", sample, "\n")
    
    input_file <- file.path(input_dir, paste0(sample, "_UMI_read_counts_full.csv"))
    output_file <- file.path(output_dir, paste0(sample, suffix))
    
    # Construct valid UMI output filename with sample name and mode included
    valid_umi_file <- if (save_valid_umi) {
      file.path(valid_umi_dir, paste0(sample, umi_suffix))
    } else {
      NULL
    }
    
    # Check if input exists
    if (!file.exists(input_file)) {
      warning("Input file not found: ", input_file)
      return("failed")
    }
    
    # Process with error handling
    tryCatch({
      result <- set_minimums(
        mode = mode, 
        umi_file = input_file, 
        min_reads = min_reads, 
        min_region_count = min_region_count, 
        min_umi = min_umi, 
        min_cells = min_cells,
        min_reads_cell = min_reads_cell,
        save_valid_umi = save_valid_umi,
        valid_umi_output_file = valid_umi_file
      )
      
      all_results[[sample]] <<- result
      
      # Save filtered cell counts
      if (save_output == TRUE) {
        write.csv(result, file = output_file, row.names = FALSE)
      }
      
      # Collect valid UMI data for aggregation
      if (aggregate_valid_umi && save_valid_umi && file.exists(valid_umi_file)) {
        sample_umi_data <- read.csv(valid_umi_file, stringsAsFactors = FALSE)
        sample_umi_data$sample <- sample
        all_valid_umi[[sample]] <<- sample_umi_data
      }
      
      return("success")
      
    }, error = function(e) {
      cat("Error processing sample", sample, ":", conditionMessage(e), "\n")
      return("failed")
    })
  }, USE.NAMES = FALSE)
  
  # Write aggregated valid UMI file if requested
  if (aggregate_valid_umi && save_valid_umi && length(all_valid_umi) > 0) {
    aggregated_file <- file.path(output_dir, paste0("all_samples_", mode, "_valid_cellID_UMI.csv"))
    aggregated_data <- do.call(rbind, all_valid_umi)
    
    # Reorder columns: sample, cellID, UMI
    aggregated_data <- aggregated_data[, c("sample", "cellID", "UMI")]
    
    write.csv(aggregated_data, file = aggregated_file, row.names = FALSE)
    cat("\nSaved aggregated valid UMI data with", nrow(aggregated_data), 
        "entries across", length(all_valid_umi), "samples\n")
  }
  
  # Summary
  successful <- sum(results == "success")
  failed <- sum(results == "failed")
  cat("\nProcessing complete! Successful:", successful, "Failed:", failed, "\n")
  
  if (save_output == TRUE) {
    return(invisible(list(successful = successful, failed = failed)))
  }
  else {
    return(all_results)
  }
}

# =============================================================================
# Host targets
# =============================================================================
## Host: extract target from filename
extract_category_host <- function(filename, targets = NULL) {
  basename_file <- basename(filename)
  parts <- strsplit(basename_file, "_alignment_")[[1]]
  
  if (length(parts) != 2) {
    cat("  Warning: Filename does not match expected format:", basename_file, "\n")
    return(NA)
  }
  alignment_part <- sub("\\.csv$", "", parts[2])
  
  # Extract target_part from the prefix
  # Format: [experiment_name]_[sample_name]_[target_part]
  prefix_parts <- strsplit(parts[1], "_")[[1]]
  
  # test if TenX
  if (!any(c("EGL", "CaRe") %in% prefix_parts)) {
    # process for TenX
    return(alignment_part)
  } else {
    # process for CaRe
    if (length(prefix_parts) < 3) {
      cat("  Warning: Cannot extract target_part from:", parts[1], "\n")
      return(NA)
    }
    
    # target_part is the last element in the prefix
    target_part <- prefix_parts[length(prefix_parts)]
    
    # Filter by allowed targets list (if provided)
    if (!is.null(targets) && !(target_part %in% targets)) {
      return(NA)
    }
    
    # Filter criterion: target_part must match alignment_part
    if (grepl(target_part, alignment_part, fixed = TRUE)) {
      return(target_part)
    }
    
    # remove version tail
    target_base <- sub("v[0-9]+$", "", target_part)
    
    if (target_base != target_part && grepl(target_base, alignment_part, fixed = TRUE)) {
      return(target_part)
    }
    return(NA)
  }
}
## Host: read and process alignment tables
process_bamsort_host <- function(samples_df, input.dir, output.dir, raw_cellIDs, targets = NULL) {
  
  if (!is.data.frame(samples_df) || !all(c("samples", "nicknames") %in% colnames(samples_df))) {
    stop("samples_df must be a data frame with 'samples' and 'nicknames' columns")
  }
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = TRUE)
  }
  
  for (i in seq_len(nrow(samples_df))) {
    sample <- samples_df$samples[i]
    nickname <- samples_df$nicknames[i]
    if (sample == nickname) {
      cat("\nProcessing sample:", sample, "\n")
    } else {
      cat("\nProcessing sample:", sample, "(nickname:", nickname, ")\n")
    }
    
    # Collect files from all input directories using nickname for pattern matching
    sample_files <- unlist(lapply(input.dir, function(dir) {
      list.files(dir, 
                 pattern = paste0(".*", nickname, ".*_alignment_.*\\.csv$"),
                 full.names = TRUE)
    }))
    # Try using full name instead
    if (length(sample_files) == 0) {
      sample_files <- unlist(lapply(input.dir, function(dir) {
        list.files(dir, 
                   pattern = paste0(".*", sample, ".*_alignment_.*\\.csv$"),
                   full.names = TRUE)
      }))
    }
    
    if (length(sample_files) == 0) {
      cat("No files found for sample:", sample, "(nickname:", nickname,")\n")
      next
    }
    
    # Extract categories - NA for files where target_part not in alignment_part
    # or target not in allowed targets list
    categories <- sapply(sample_files, extract_category_host, targets = targets)

    # Filter to only files with valid categories (target in alignment and in allowed list)
    valid_idx <- !is.na(categories)
    n_filtered <- sum(!valid_idx)

    sample_files <- sample_files[valid_idx]
    categories <- categories[valid_idx]

    if (length(sample_files) == 0) {
      cat("No files remaining after target-alignment filtering for sample:", sample,
          "(nickname:", nickname, ")\n")
      next
    }

    cat("Found", length(sample_files), "matching alignment files\n")
    
    # Aggregate data from all valid files
    all_data <- data.frame()
    filenames <- basename(sample_files)
    
    for (j in seq_along(sample_files)) {
      extracted_data <- NULL
      
      tryCatch({
        file_data <- fread(sample_files[j], data.table = FALSE)
        
        if (nrow(file_data) == 0) {
          cat("Warning: File is empty:", filenames[j], "\n")
          next
        }
        
        required_cols <- c("cellID", "UMI", "count")
        missing_cols <- setdiff(required_cols, colnames(file_data))
        
        if (length(missing_cols) > 0) {
          cat("Warning: Missing required columns:", paste(missing_cols, collapse = ", "), 
              "in", filenames[j], "\n")
          cat("Available columns:", paste(colnames(file_data), collapse = ", "), "\n")
          next
        }
        
        extracted_data <- file_data %>%
          select(cellID, UMI, read = dplyr::all_of(bamsort_read_column(file_data))) %>%
          mutate(category = categories[j])
        
      }, error = function(e) {
        cat("Error reading file:", filenames[j], "\n")
        cat("Error message:", conditionMessage(e), "\n")
      })
      
      if (!is.null(extracted_data) && nrow(extracted_data) > 0) {
        all_data <- rbind(all_data, extracted_data)
      }
    }
    
    # Process aggregated data
    if (nrow(all_data) > 0) {
      all_data <- all_data %>%
        filter(!is.na(cellID) & !is.na(UMI) & !is.na(read)) %>%
        unique()
      
      tryCatch({
        # Aggregate UMI read counts (no multimap resolution)
        umi_read_counts <- all_data %>%
          group_by(cellID, UMI, category) %>%
          summarize(
            read_count = sum(as.numeric(read)),
            .groups = 'drop'
          )
        
        cat("Gathered read counts for", nrow(umi_read_counts), "UMIs from", 
            n_distinct(umi_read_counts$cellID), "cells\n")
        
        # Filter to valid 10X cell barcodes
        if (sample %in% names(raw_cellIDs) && !is.null(raw_cellIDs[[sample]])) {
          valid_cells <- raw_cellIDs[[sample]]
          filtered_umi_read_counts <- umi_read_counts %>%
            filter(cellID %in% valid_cells)
          cat("After filtering to", length(valid_cells), "valid 10X cells: kept",
              nrow(filtered_umi_read_counts), "UMIs from",
              n_distinct(filtered_umi_read_counts$cellID), "cells\n")
        } else {
          cat("Warning: No 10X data found for sample", sample, "- using unfiltered cell list\n")
          filtered_umi_read_counts <- umi_read_counts
        }
        
        # Convert to clean data frame for output
        umi_read_counts_clean <- data.frame(
          cellID = as.character(umi_read_counts$cellID),
          UMI = as.character(umi_read_counts$UMI),
          category = as.character(umi_read_counts$category),
          read_count = as.numeric(umi_read_counts$read_count),
          stringsAsFactors = FALSE
        )
        
        filtered_umi_read_counts_clean <- data.frame(
          cellID = as.character(filtered_umi_read_counts$cellID),
          UMI = as.character(filtered_umi_read_counts$UMI),
          category = as.character(filtered_umi_read_counts$category),
          read_count = as.numeric(filtered_umi_read_counts$read_count),
          stringsAsFactors = FALSE
        )
        
        # Save outputs
        umi_level_file <- file.path(output.dir, paste0(sample, "_UMI_read_counts_raw.csv"))
        write.csv(umi_read_counts_clean, file = umi_level_file, row.names = FALSE)
        filtered_umi_level_file <- file.path(output.dir, paste0(sample, "_UMI_read_counts_full.csv"))
        write.csv(filtered_umi_read_counts_clean, file = filtered_umi_level_file, row.names = FALSE)

        rm(umi_read_counts)
        
      }, error = function(e) {
        cat("Error processing data for sample", sample, ":", conditionMessage(e), "\n")
        
        # Attempt to save raw data for debugging
        tryCatch({
          all_data_clean <- data.frame(
            cellID = as.character(all_data$cellID),
            UMI = as.character(all_data$UMI),
            read = as.numeric(all_data$read),
            category = as.character(all_data$category),
            stringsAsFactors = FALSE
          )
          
          raw_file <- file.path(output.dir, paste0(sample, "_raw_data.csv"))
          write.csv(all_data_clean, file = raw_file, row.names = FALSE)
          cat("Saved raw data to:", raw_file, "\n")
        }, error = function(e2) {
          cat("Could not save raw data:", conditionMessage(e2), "\n")
        })
      })
    } else {
      cat("No data extracted for sample:", sample, "(nickname:", nickname, ")\n")
    }
  }
  
  invisible(NULL)
}
