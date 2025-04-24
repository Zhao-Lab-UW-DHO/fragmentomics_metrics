suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  sample_file <- input_args[1]
  output_file <- input_args[2]
}

# all_atac <- read_tsv("unique_ATAC_cancers.txt", col_types = cols())

calculate_atac_entropy <- function(input_file) {
  
  sample_name <- str_replace(basename(input_file), "_ATAC_frag_count.txt.gz", "")
  # sample_name <- str_replace(sample_name, "output/metrics/ATAC_entropy/", "")
  
  output <- read_table(input_file,
                             col_names = c("count", "cancer", "fragsize"), 
                             col_types = cols()) %>% 
    group_by(cancer) %>% 
    summarise(ATAC_entropy = entropy(count)) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, cancer, ATAC_entropy)
  
  # output <- all_atac %>% 
  #   left_join(atac_entropy, by = "cancer") %>% 
  #   pivot_wider(names_from = cancer, values_from = ATAC_entropy) %>% 
  #   mutate(sample = sample_name) %>% 
  #   relocate(sample)
  
  return(output)
}

# calculate entropy by TF
sample_ATAC_entropy <- calculate_atac_entropy(sample_file)

# write to file
# output_file <- str_replace(sample_file, "_ATAC_frag_count.txt.gz", "_ATAC_entropy.txt")
write_tsv(sample_ATAC_entropy, output_file)
