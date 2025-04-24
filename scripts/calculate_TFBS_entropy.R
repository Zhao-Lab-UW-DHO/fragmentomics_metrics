suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  sample_file <- input_args[1]
  output_file <- input_args[2]
}

# all_tfbs <- read_tsv("unique_TFs_top5k.txt", col_types = cols())

calculate_tfbs_entropy <- function(input_file) {
  
  sample_name <- str_replace(basename(input_file), "_TFBS_frag_count.txt.gz", "")
  # sample_name <- str_replace(sample_name, "output/metrics/TFBS_entropy/", "")
  
  output <- read_table(input_file,
                             col_names = c("count", "frag_size", "tfbs"),
                             col_types = cols()) %>% 
    group_by(tfbs) %>% 
    summarise(TF_entropy = entropy(count)) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, tfbs, TF_entropy)
  
  # output <- all_tfbs %>% 
  #   left_join(tfbs_entropy, by = "tfbs") %>% 
  #   pivot_wider(names_from = tfbs, values_from = TF_entropy) %>% 
  #   mutate(sample = sample_name) %>% 
  #   relocate(sample)
  
  return(output)
}

# calculate entropy by TF
sample_TFBS_entropy <- calculate_tfbs_entropy(sample_file)

# write to file
# output_file <- str_replace(sample_file, "_TFBS_frag_count.txt.gz", "_TFBS_entropy.txt")
write_tsv(sample_TFBS_entropy, output_file)