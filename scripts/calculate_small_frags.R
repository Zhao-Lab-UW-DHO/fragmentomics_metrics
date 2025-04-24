suppressPackageStartupMessages(library(tidyverse))

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
calculate_small_frag_prop <- function(input_file) {
  sample_name <- str_replace(basename(input_file), ".txt.gz", "")
  # sample_name <- str_replace(sample_name, "output/fragstats/SE_files/", "")
  
  read_data <- read_tsv(input_file,
                         col_names = c("gene", "exon", "fragsize", "count"),
                         col_types = cols())
  
  output <- read_data %>% 
    mutate(bin = cut(fragsize, breaks = c(0, 150, 500), labels = c("small", "large"))) %>% 
    group_by(gene, exon, bin) %>% 
    summarise(bin_count = sum(count)) %>% 
    ungroup() %>% 
    group_by(gene, exon) %>% 
    mutate(bin_prop = bin_count / sum(bin_count)) %>% 
    ungroup() %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    filter(bin == "small") %>% 
    mutate(sample = sample_name) %>% 
    dplyr::select(sample, id, bin_prop)
    
  # output <- all_ids %>% 
  #   left_join(small_frag_data, by = "id") %>% 
  #   dplyr::select(id, bin_prop) %>% 
  #   pivot_wider(names_from = id, values_from = bin_prop) %>% 
  #   mutate(sample = sample_name) %>% 
  #   relocate(sample)
  
  return(output)
}

# normalize sample
sample_output <- calculate_small_frag_prop(input_file)

# set output file and write to file
# output_file <- str_replace(input_file, ".txt.gz", ".smallfrag.tsv")
# output_file <- str_replace(output_file, "mixed_SE_files/", "mixed_ratio_metrics/small_frags/")
write_tsv(sample_output, output_file)
