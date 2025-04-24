suppressPackageStartupMessages(library(tidyverse))

downsample_levels <- c("100M", "50M", "25M", "10M", "5M", "1M")

feature_tables <- c("ATAC_entropy", "depth", "frag_bins", "full_gene_depth", "mds", "se", "small_frags", "TFBS_entropy")

for (ft in feature_tables) {
  
  file_prefix <- case_when(
    ft == "ATAC_entropy" ~ "*_ATAC_entropy.txt",
    ft == "depth" ~ "*.depth.tsv",
    ft == "frag_bins" ~ "*.fragbins.tsv",
    ft == "full_gene_depth" ~ "*.fullgenedepth.tsv",
    ft == "mds" ~ "*_mds.txt",
    ft == "se" ~ "*.SE.tsv", 
    ft == "small_frags" ~ "*.smallfrag.tsv", 
    ft == "TFBS_entropy" ~ "*_TFBS_entropy.txt"
  )
  
  for (dl in downsample_levels) {
    print(str_c("Building feature table for", ft, dl, sep = " "))
    
    file_path <- str_c("output/downsample/metrics/", ft, "/", dl, "/")
    
    out_path <- str_c("output/downsample/feature_tables/grail_", ft, "_", dl, ".rds")
    
    file_list <- list.files(path = file_path, pattern = file_prefix)
    
    all_data <- purrr::map_dfr(str_c(file_path, file_list), ~ read_tsv(.x, show_col_types = FALSE))
    
    
    all_data <- all_data %>% 
      dplyr::rename("sample" = 1, "feature" = 2, "value" = 3) %>% 
      filter(str_detect(feature, "CRLF2", negate = TRUE)) # this gene ambiguously maps to both X and Y chr so is removed
    
    output_data <- all_data %>% 
      pivot_wider(names_from = feature, values_from = value)
    
    # remove or replace NAs with 0?
    output_data <- output_data %>% 
      mutate(across(everything(), ~replace_na(., 0)))
    
    saveRDS(output_data, out_path)
  }
}