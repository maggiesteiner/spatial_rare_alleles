library(data.table)
library(dplyr)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
inputdir <- args[1]
# Directory path containing the files
directory <- inputdir
print(inputdir)
# Get list of files in the directory

file_list <- list.files(pattern = "*.sfs.tsv",path = inputdir)
print(paste(length(file_list),"files found."))
# Initialize an empty list to store data frames
dfs <- list()

# Read and process each file
for (file in file_list) {
  # Extract the information from the file name
  file_name <- basename(file)
  info <- strsplit(file_name, "_", fixed = TRUE)[[1]] %>% extract_numeric() 
  info <- info[is.na(info)==F]
  file <- paste0(inputdir,"/",file)
  # Read the file, skipping the first row (commented out line)
  data <- fread(file, skip = 1)
  
  # Extract the metadata from the first line
  metadata <- readLines(file, n = 1)
  num_inds_gt_1 <- as.numeric(str_extract(metadata, "\\d+\\.\\d+"))
  
  # Add columns for the file-specific information
  data$rep <- info[1]
  data$width <- info[2]
  data$s <- info[3] %>% as.numeric()
  data$mu <- info[4] %>% as.numeric()
  data$K <- info[5] %>% as.numeric()
  data$sigma <- info[6] %>% as.numeric()
  data$f_geq_1 <- num_inds_gt_1
  
  # Append the processed data frame to the list
  dfs[[file]] <- data
}



# Initialize empty data frame to store the data
df <- data.frame()

# Loop over the file names
for (file_name in file_list) {
  # Extract the parameters from the file name
  params <- strsplit(file_name, "_")[[1]] %>% extract_numeric() 
  params <- params[is.na(params)==F]
  
  replicate <- as.integer(params[1])
  width <- as.numeric(params[2])
  selection_coefficient <- as.numeric(params[3])
  mutation_rate <- as.numeric(params[4])
  pop_density <- as.numeric(params[5])
  sigma <- as.numeric(params[6])

  num_samples <- as.integer(sub("\\.sfs\\.tsv", "", params[7]))
  
  # Read the file into a data frame
  file_path <- paste(inputdir, file_name,sep = "/")  # Replace with the actual file path
  data <- read.table(file_path, skip = 1 , header= F,na.strings = "NA") 
  # filter on row sums
  header <- str_split(read_lines(file_path, skip = 0)[1],"\t")[[1]]
  
  colnames(data) <- header
  
  # Add the parameters as columns
  data$rep <- replicate
  data$w <- width
  data$sigma <- sigma
  data$mu <- mutation_rate
  data$s <- selection_coefficient
  data$K <- pop_density
  data$sampnum <- num_samples
  
  # Append the data to the main data frame
  df <- rbind(df, data)
}

# Assuming your data frame is named 'df'
df_summary <- df %>%
  pivot_longer(cols=all_of(header[-length(header)]),names_to = "width",values_to = "f") %>%
  group_by(s,width,`allele counts`) %>% # mean across replicates
  summarise(se_f=sd(f,na.rm=T),f=mean(f,na.rm = T),n=n()) %>%
  mutate(f_per_kb = (f/1e8) * 1e3) %>%
  ungroup()

df_summary <- df_summary %>% 
  filter(`allele counts`>0,width>0) 

df_summary %>% write_delim(delim="\t",file.path(inputdir,"summary_df.tsv"))

