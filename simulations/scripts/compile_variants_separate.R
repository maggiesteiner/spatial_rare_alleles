# Get a list of all relevant files
file_list <- list.files(pattern = "*.variants.tsv")

# Create a list to store data frames for each 's' value
dfs <- list()

# Iterate through the files
for (file in file_list) {
  # Extract the 's' value from the file name
  s_value <- sub(".*s-([0-9.]+)_.*", "\\1", file)
  
  # Read the data from the current file, skipping the header
  data <- read.table(file, header = TRUE, sep = "\t")
  
  # Add the 's' value as a new column
  data$s_value <- s_value
  
  # Append the data frame to the list
  dfs[[s_value]] <- data
}

# Combine all data frames into a single data frame for each 's' value
for (s_value in names(dfs)) {
  combined_df <- do.call(rbind, dfs[s_value])
  
  # Write the combined data frame to a new file
  output_file <- paste("all_variants_s", s_value, ".tsv", sep = "")
  write.table(combined_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Concatenation completed for 's' value", s_value, ". Output saved as", output_file, "\n")
}

