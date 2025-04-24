suppressPackageStartupMessages(library(tidyverse))

# Get Strand Info
gene_strands <- read_tsv("files/UCSC_hg38_canonical_exons.bed", show_col_types = F) %>% 
  dplyr::select(gene = gene_symbol, strand) %>% 
  unique()


feature_tables <- c("depth","mds", "se", "small_frags")
gene_panels <- c("default", "tempus", "guardant", "foundationOne")

# Extract First Exon
for (ft in feature_tables) {
  for (gp in gene_panels) {
    print(str_c("Filtering E1 for", gp, ft, sep = " "))
    
    rds_inpath <- str_c("output/feature_tables/", gp, "/", ft, ".rds")
    
    outpath <- str_c("output/feature_tables/", gp, "/", ft, "_E1.rds")
    
    full_exon_data <- read_rds(rds_inpath)
    
    first_exon_data <- full_exon_data %>% 
      pivot_longer(cols = -sample, names_to = "feature", values_to = "value") %>% 
      separate(feature, into = c("gene", "exon"), sep = "_", remove = FALSE) %>% 
      left_join(gene_strands, by = "gene") %>% 
      group_by(gene) %>% 
      mutate(
        first_exon = ifelse(strand == "+" & exon == min(exon), TRUE, FALSE),
        first_exon = ifelse(strand == "-" & exon == max(exon), TRUE, first_exon)
      ) %>% 
      filter(first_exon) %>% 
      ungroup() %>% 
      dplyr::select(sample, feature, value) %>% 
      pivot_wider(names_from = feature, values_from = value)
    
    saveRDS(first_exon_data, outpath)
  }
}