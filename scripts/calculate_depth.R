suppressPackageStartupMessages(library(tidyverse))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  input_file <- input_args[1]
  output_file <- input_args[2]
}

# need to get list of all unique gene_exon values with their exon sizes
all_ids <- read_tsv("files/UCSC_hg38_exon_sizes.tsv",
                    col_types = cols()) %>%
  mutate(id = str_c(gene, exon, sep = "_")) %>%
  dplyr::select(id, size)

# function to take depth data and normalize it by exon length and reads per sample
normalize_depth_data <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt", "")
  # sample_name <- str_replace(sample_name, "output/fragstats/depth_files/", "")
  
  depth_data <- read_tsv(input_file,
                        col_names = c("gene", "exon", "count"),
                        col_types = cols())
  
  num_reads <- sum(depth_data$count)
  
  output <- depth_data %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    left_join(all_ids, by = "id") %>% 
    mutate(norm_depth = (count / size) / (num_reads / 1e6)) %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, norm_depth)
  
  # output <- all_ids %>% 
  #   left_join(depth_data_norm, by = "id") %>% 
  #   dplyr::select(id, norm_depth) %>% 
  #   pivot_wider(names_from = id, values_from = norm_depth) %>% 
  #   mutate(sample = sample_name) %>% 
  #   relocate(sample)
  
  return(output)
}

# normalize sample
sample_normalized <- normalize_depth_data(input_file)

# set output file and write to file
# output_file <- str_replace(input_file, ".txt", ".normlong.tsv")
write_tsv(sample_normalized, output_file)
