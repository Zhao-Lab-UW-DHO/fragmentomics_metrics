suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  input_file <- input_args[1]
  output_file <- input_args[2]
}

# need to get list of all unique gene_exon values with their exon sizes
# all_ids <- read_tsv("/mnt/Scratch02/helzer_temp/fragmentomics_mixing/grail/unique_grail_gene_exons.txt",
#                     col_types = cols()) %>% 
#   mutate(id = str_c(gene, exon, sep = "_")) %>% 
#   dplyr::select(id)

# function to read in file and calculate Shannon entropy
calculate_SE <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt.gz", "")
  # sample_name <- str_replace(sample_name, "output/fragstats/SE_files/", "")
  
  read_data <- read_tsv(input_file,
                         col_names = c("gene", "exon", "fragsize", "count"),
                         col_types = cols())
  
  output <- read_data %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    group_by(id) %>% 
    summarise(se = entropy(count)) %>% 
    ungroup() %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, se)
    
  # output <- all_ids %>% 
  #   left_join(se_data, by = "id") %>% 
  #   dplyr::select(id, se) %>% 
  #   pivot_wider(names_from = id, values_from = se) %>% 
  #   mutate(sample = sample_name) %>% 
  #   relocate(sample)
  
  return(output)
}

# function to take depth data and normalize it by exon length and reads per sample
# normalize_depth_data <- function(input_file) {
#   sample_name <- str_replace(input_file, ".txt", "")
#   
#   depth_data <- read_tsv(input_file,
#                         col_names = c("gene", "exon", "count"),
#                         col_types = cols())
#   
#   num_reads <- sum(depth_data$count)
#   
#   depth_data_norm <- depth_data %>% 
#     mutate(id = str_c(gene, exon, sep = "_")) %>% 
#     left_join(all_ids, by = "id") %>% 
#     mutate(norm_depth = (count / size) / (num_reads / 1e6))
#   
#   output <- all_ids %>% 
#     left_join(depth_data_norm, by = "id") %>% 
#     dplyr::select(id, norm_depth) %>% 
#     pivot_wider(names_from = id, values_from = norm_depth) %>% 
#     mutate(sample = sample_name) %>% 
#     relocate(sample)
#   
#   return(output)
# }

# normalize sample
sample_SE <- calculate_SE(input_file)

# set output file and write to file
# output_file <- str_replace(input_file, ".txt.gz", ".SElong.tsv")
write_tsv(sample_SE, output_file)
