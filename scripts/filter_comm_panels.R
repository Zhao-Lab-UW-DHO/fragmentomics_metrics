suppressPackageStartupMessages(library(tidyverse))

input_args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  input_file <- input_args[1]
}



output_file <- str_replace(basename(input_file), ".gz", "")

# Load gene lists
tempus_genes <- read_tsv("files/manuscript_data/exon_panels/tempus_xF.tsv", show_col_types = F) %>% pull(gene)
guardant_genes <- read_tsv("files/manuscript_data/exon_panels/guardant360CDx.tsv", show_col_types = F) %>% pull(gene)
f1_genes <- read_tsv("files/manuscript_data/exon_panels/foundationOneCDx.tsv", show_col_types = F) %>% pull(gene)

default_reads <- read_tsv(input_file, show_col_types = F,
         col_names = c("chr", "start", "stop", "ens_id", "refseq", "gene", "exon", "strand"))


# filter reads for genes in each panel

tempus_reads <- default_reads %>% filter(gene %in% tempus_genes)
guardant_reads <- default_reads %>% filter(gene %in% guardant_genes)
f1_reads <- default_reads %>% filter(gene %in% f1_genes)


write_tsv(default_reads, str_c("data/default/", output_file), col_names = FALSE)
write_tsv(tempus_reads, str_c("data/tempus/", output_file), col_names = FALSE)
write_tsv(guardant_reads, str_c("data/guardant/", output_file), col_names = FALSE)
write_tsv(f1_reads, str_c("data/foundationOne/", output_file), col_names = FALSE)

# system2("gzip", args = str_c("data/default/", output_file))
# system2("gzip", args = str_c("data/tempus/", output_file))
# system2("gzip", args = str_c("data/guardant/", output_file))
# system2("gzip", args = str_c("data/foundationOne/", output_file))