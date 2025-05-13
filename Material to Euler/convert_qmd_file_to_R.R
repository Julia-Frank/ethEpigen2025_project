extract_code_from_qmd <- function(qmd_file, output_file = "extracted_code.R") {
  # Read the .qmd file
  lines <- readLines(qmd_file)
  
  # Identify code chunk start and end
  chunk_starts <- grep("^```\\{.*\\}$", lines)
  chunk_ends <- grep("^```$", lines)
  
  # Ensure matching starts and ends
  if (length(chunk_starts) != length(chunk_ends)) {
    stop("Mismatch between chunk start and end markers.")
  }
  
  # Extract code chunks
  code_chunks <- mapply(
    function(start, end) lines[(start + 1):(end - 1)],
    start = chunk_starts,
    end = chunk_ends
  )
  
  # Flatten and write to an output file
  all_code <- unlist(code_chunks)
  writeLines(all_code, output_file)
  
  message(paste("Code extracted to", output_file))
}

# Example usage
extract_code_from_qmd("data_import_and_clustering.qmd", "opt_data_import_and_clustering.R")
