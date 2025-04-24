suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(entropy))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  left_file <- input_args[1]
  right_file <- input_args[2]
  output_file <- input_args[3]
}

# need to get list of all unique gene_exon values with their exon sizes
# all_ids <- read_tsv("/mnt/Scratch02/helzer_temp/fragmentomics_mixing/grail/unique_grail_gene_exons.txt",
#                     col_types = cols()) %>% 
#   mutate(id = str_c(gene, exon, sep = "_")) %>% 
#   dplyr::select(id)

calculate_mds <- function(left_data, right_data) {
  
  sample_name <- str_replace(basename(left_data), "_left_4mers.txt.gz", "")
  # sample_name <- str_replace(sample_name, "output/metrics/mds/", "")
  
  left_motifs <- read_table(left_data,
                          col_names = c("count", "gene", "exon", "motif"),
                          col_types = cols())
  right_motifs <- read_table(right_data,
                           col_names = c("count", "gene", "exon", "motif"),
                           col_types = cols())
  
  combined_data <- bind_rows(left_motifs, right_motifs) %>% 
    group_by(gene, exon, motif) %>% 
    summarise(
      count = sum(count)
    ) %>% 
    ungroup() %>% 
    mutate(id = str_c(gene, exon, sep = "_")) %>% 
    dplyr::select(id, motif, count)
  
  output <- combined_data %>% 
    group_by(id) %>% 
    summarise(
      mds = entropy(count)
    ) %>% 
    mutate(sample = sample_name) %>% 
    relocate(sample, id, mds)
  
  # output <- all_ids %>% 
  #   left_join(mds, by = "id") %>% 
  #   pivot_wider(names_from = id, values_from = mds) %>% 
  #   mutate(sample = sample_name) %>% 
  #   relocate(sample)
  
  return(output)
}

sample_mds <- calculate_mds(left_data = left_file, right_data = right_file)

# output_file <- str_replace(left_file, "_left_4mers.txt.gz", "_mds.txt")
write_tsv(sample_mds, output_file)
