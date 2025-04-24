suppressPackageStartupMessages(library(tidyverse))

get_top_15k_by_variance <- function(df) {
  output <- sapply(df %>% dplyr::select(-c(sample)), var) %>% 
    as.data.frame() %>% 
    as_tibble(rownames = "feature") %>% 
    dplyr::rename(variance = 2) %>% 
    arrange(desc(variance)) %>% 
    head(n = 15000) %>% 
    pull(feature)
  return(output)
}

append_suffix_to_colnames <- function(df, suffix) {
  output <- df %>% 
    pivot_longer(cols = -c(sample), names_to = "feature", values_to = "value") %>% 
    mutate(feature = str_c(feature, suffix, sep = "__"))
  return(output)
}

gene_panels <- c("default", "tempus", "guardant", "foundationOne")

for (gp in gene_panels) {
  print(str_c("Generating Combined Feature Table for", gp, sep = " "))
  
  outpath <- str_c("output/feature_tables/", gp, "/all_combined.rds")
  
  df1 <- read_rds(str_c("output/feature_tables/", gp, "/ATAC_entropy.rds")) %>% append_suffix_to_colnames(suffix = "ATAC_entropy")
  df2 <- read_rds(str_c("output/feature_tables/", gp, "/depth.rds")) %>% append_suffix_to_colnames(suffix = "depth")
  df3 <- read_rds(str_c("output/feature_tables/", gp, "/frag_bins.rds")) %>% append_suffix_to_colnames(suffix = "frag_bins")
  df4 <- read_rds(str_c("output/feature_tables/", gp, "/full_gene_depth.rds")) %>% append_suffix_to_colnames(suffix = "full_gene_depth")
  df5 <- read_rds(str_c("output/feature_tables/", gp, "/mds.rds")) %>% append_suffix_to_colnames(suffix = "mds")
  df6 <- read_rds(str_c("output/feature_tables/", gp, "/se.rds")) %>% append_suffix_to_colnames(suffix = "se")
  df7 <- read_rds(str_c("output/feature_tables/", gp, "/small_frags.rds")) %>% append_suffix_to_colnames(suffix = "small_frags")
  df8 <- read_rds(str_c("output/feature_tables/", gp, "/TFBS_entropy.rds")) %>% append_suffix_to_colnames(suffix = "TFBS_entropy")
  
  combined_df <- bind_rows(df1, df2, df3, df4, df5, df6, df7, df8) %>% 
    pivot_wider(names_from = feature, values_from = value)
  
  top_features <- get_top_15k_by_variance(combined_df)
  
  combined_df_filtered <- combined_df %>% dplyr::select(sample, all_of(top_features))
  
  saveRDS(combined_df_filtered, outpath)
}