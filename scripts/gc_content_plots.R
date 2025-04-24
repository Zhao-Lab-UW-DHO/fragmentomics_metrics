library(tidyverse)
library(janitor)
library(patchwork)

# metadata
# uw_metadata <- read_rds("/mnt/Data01/projects/fragmentomics/data/clean_rds/zhao_metadata_full_n449_updated_ctfrac_breast_subtype.rds")
uw_metadata <- read_rds("files/uw_ml_metadata.rds")

samples_in_manuscript <- uw_metadata %>% 
  filter(phenotype %in% c("Cancer prostate", "Cancer breast", "Cancer lung", "Cancer rcc", "Non cancer", "Cancer bladder")) %>% 
  pull(sample)

# read in feature tables
uw_atac <- read_rds("output/feature_tables/default/ATAC_entropy.rds")
uw_depth <- read_rds("output/feature_tables/default/depth.rds")
uw_frag_bins <- read_rds("output/feature_tables/default/frag_bins.rds")
uw_small_frags <- read_rds("output/feature_tables/default/small_frags.rds")
uw_mds <- read_rds("output/feature_tables/default/mds.rds")
uw_se <- read_rds("output/feature_tables/default/se.rds")
uw_tfbs <- read_rds("output/feature_tables/default/TFBS_entropy.rds")
uw_full_gene_depth <- read_rds("output/feature_tables/default/full_gene_depth.rds")

grail_atac <- read_rds("output/feature_tables/grail/default/ATAC_entropy.rds")
grail_depth <- read_rds("output/feature_tables/grail/default/depth.rds")
grail_frag_bins <- read_rds("output/feature_tables/grail/default/frag_bins.rds")
grail_small_frags <- read_rds("output/feature_tables/grail/default/small_frags.rds")
grail_mds <- read_rds("output/feature_tables/grail/default/mds.rds")
grail_se <- read_rds("output/feature_tables/grail/default/se.rds")
grail_tfbs <- read_rds("output/feature_tables/grail/default/TFBS_entropy.rds")
grail_full_gene_depth <- read_rds("output/feature_tables/grail/default/full_gene_depth.rds")

# read in GC content info
exons_gc_content <- read_tsv("files/gc_content/UCSC_hg38_canonical_exons_gc_content.bed") %>% 
  clean_names() %>% 
  mutate(exon_id = str_c(x6_usercol, x7_usercol, sep = "_")) %>% 
  dplyr::select(exon_id, gc_content = x10_pct_gc)

# first_exon_gc_content <- read_tsv("files/gc_content/UCSC_hg38_canonical_exons_gc_content.bed") %>% 
#   clean_names()

full_genes_gc_content <- read_tsv("files/gc_content/UCSC_hg38_canonical_exons_gc_content.bed") %>% 
  clean_names() %>% 
  group_by(x6_usercol) %>% 
  summarise(
    a_total = sum(x11_num_a),
    c_total = sum(x12_num_c),
    g_total = sum(x13_num_g),
    t_total = sum(x14_num_t),
    bases_total = sum(x17_seq_len)
  ) %>% 
  mutate(gc_content = (c_total + g_total) / bases_total) %>% 
  dplyr::select(gene = x6_usercol, gc_content)

atac_gc_uw <- read_tsv("files/gc_content/ATAC_UW_v1_panel_overlap_gc_content.bed") %>% 
  clean_names() %>% 
  separate(x4_usercol, into = c("cancer_type", "region_id"), remove = FALSE, sep = "_") %>% 
  group_by(cancer_type) %>% 
  summarise(
    a_total = sum(x7_num_a),
    c_total = sum(x8_num_c),
    g_total = sum(x9_num_g),
    t_total = sum(x10_num_t),
    bases_total = sum(x13_seq_len)
  ) %>% 
  mutate(gc_content = (c_total + g_total) / bases_total) %>% 
  dplyr::select(cancer_type, gc_content)

atac_gc_grail <- read_tsv("files/gc_content/ATAC_grail_panel_overlap_gc_content.bed") %>% 
  clean_names() %>% 
  separate(x4_usercol, into = c("cancer_type", "region_id"), remove = FALSE, sep = "_") %>% 
  group_by(cancer_type) %>% 
  summarise(
    a_total = sum(x7_num_a),
    c_total = sum(x8_num_c),
    g_total = sum(x9_num_g),
    t_total = sum(x10_num_t),
    bases_total = sum(x13_seq_len)
  ) %>% 
  mutate(gc_content = (c_total + g_total) / bases_total) %>% 
  dplyr::select(cancer_type, gc_content)

tfbs_gc_uw <- read_tsv("files/gc_content/TFBS_overlaps_UW_v1_panel_hg38_gc_content.bed") %>% 
  clean_names() %>% 
  group_by(x4_usercol) %>% 
  summarise(
    a_total = sum(x7_num_a),
    c_total = sum(x8_num_c),
    g_total = sum(x9_num_g),
    t_total = sum(x10_num_t),
    bases_total = sum(x13_seq_len)
  ) %>% 
  mutate(gc_content = (c_total + g_total) / bases_total) %>% 
  dplyr::select(TF = x4_usercol, gc_content)

tfbs_gc_grail <- read_tsv("files/gc_content/TFBS_overlaps_grail_panel_hg19_gc_content.bed") %>% 
  clean_names() %>% 
  group_by(x4_usercol) %>% 
  summarise(
    a_total = sum(x7_num_a),
    c_total = sum(x8_num_c),
    g_total = sum(x9_num_g),
    t_total = sum(x10_num_t),
    bases_total = sum(x13_seq_len)
  ) %>% 
  mutate(gc_content = (c_total + g_total) / bases_total) %>% 
  dplyr::select(TF = x4_usercol, gc_content)

##############################
### Make Correlation Plots ###
##############################

make_cor_text <- function(cor_output) {
  rsquared <- signif((cor_output$estimate)^2, 3)
  output <- str_c("R2 = ", rsquared)
  return(output)
}

color_scale <- scale_fill_gradient(low = "#FBFBFB", high = "blue")

# ATAC Entropy
uw_atac_long <- uw_atac %>% 
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:trial), names_to = "cancer_type", values_to = "metric") %>% 
  left_join(atac_gc_uw, by = "cancer_type")
uw_atac_cor <- cor.test(uw_atac_long$metric, uw_atac_long$gc_content)
uw_atac_cor_text <- make_cor_text(cor_output = uw_atac_cor)

uw_atac_gc_plot <- uw_atac_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("ATAC Entropy") +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: ATAC Entropy") +
  annotate("text", label = uw_atac_cor_text, x = 0.63, y = 6.0) +
  theme(
    aspect.ratio = 1
  )

grail_atac_long <- grail_atac %>% 
  pivot_longer(cols = -c(sample:phenotype, ctdna_fraction, gdna_fraction), names_to = "cancer_type", values_to = "metric") %>% 
  left_join(atac_gc_grail, by = "cancer_type")
grail_atac_cor <- cor.test(grail_atac_long$metric, grail_atac_long$gc_content)
grail_atac_cor_text <- make_cor_text(cor_output = grail_atac_cor)

grail_atac_gc_plot <- grail_atac_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("ATAC Entropy") +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  annotate("text", label = grail_atac_cor_text, x = 0.6, y = 6.0) +
  ggtitle("Cohort: GRAIL\nFeature: ATAC Entropy") +
  theme(
    aspect.ratio = 1
  )


# Depth
uw_depth_long <- uw_depth %>% 
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:trial), names_to = "exon_id", values_to = "metric") %>% 
  left_join(exons_gc_content, by = "exon_id")
uw_depth_cor <- cor.test(uw_depth_long$metric, uw_depth_long$gc_content)
uw_depth_cor_text <- make_cor_text(cor_output = uw_depth_cor)

uw_depth_gc_plot <- uw_depth_long %>% 
  ggplot(aes(gc_content, log10(metric))) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("log10(Normalized Exon Depth)") +
  annotate("text", label = uw_depth_cor_text, x = 0.75, y = 2.0) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: Depth") +
  theme(
    aspect.ratio = 1
  )

grail_depth_long <- grail_depth %>% 
  pivot_longer(cols = -c(sample:phenotype), names_to = "exon_id", values_to = "metric") %>% 
  left_join(exons_gc_content, by = "exon_id")
grail_depth_cor <- cor.test(grail_depth_long$metric, grail_depth_long$gc_content)
grail_depth_cor_text <- make_cor_text(cor_output = grail_depth_cor)

grail_depth_gc_plot <- grail_depth_long %>% 
  ggplot(aes(gc_content, log10(metric))) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("log10(Normalized Exon Depth)") +
  annotate("text", label = grail_depth_cor_text, x = 0.75, y = 2.0) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: Depth") +
  theme(
    aspect.ratio = 1
  )

# Frag Bins
starttime <- Sys.time() # ~ 10 min
uw_frag_bins_long <- uw_frag_bins %>%
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample), names_to = "exon_bin_id", values_to = "metric") %>%
  separate(exon_bin_id, into = c("exon_id", "bin"), sep = "_bin")
uw_frag_bins_long <- uw_frag_bins_long %>% 
  left_join(exons_gc_content, by = "exon_id")
elapsedtime <- Sys.time() - starttime
print(elapsedtime)

uw_frag_bins_cor <- cor.test(uw_frag_bins_long$metric, uw_frag_bins_long$gc_content)
uw_frag_bins_cor_text <- make_cor_text(cor_output = uw_frag_bins_cor)
uw_frag_bins_gc_plot <- uw_frag_bins_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("Read Proportion in Bin") +
  annotate("text", label = uw_frag_bins_cor_text, x = 0.75, y = 0.7) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: Fragment Bins") +
  theme(
    aspect.ratio = 1
  )
  
starttime <- Sys.time()
grail_frag_bins_long <- grail_frag_bins %>%
  pivot_longer(cols = -c(sample), names_to = "exon_bin_id", values_to = "metric") %>%
  separate(exon_bin_id, into = c("exon_id", "bin"), sep = "_bin") %>% 
  left_join(exons_gc_content, by = "exon_id")
elapsedtime <- Sys.time() - starttime
print(elapsedtime)

grail_frag_bins_cor <- cor.test(grail_frag_bins_long$metric, grail_frag_bins_long$gc_content)
grail_frag_bins_cor_text <- make_cor_text(cor_output = grail_frag_bins_cor)
grail_frag_bins_gc_plot <- grail_frag_bins_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("Read Proportion in Bin") +
  annotate("text", label = grail_frag_bins_cor_text, x = 0.75, y = 0.7) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: Fragment Bins") +
  theme(
    aspect.ratio = 1
  )


# small frags
uw_small_frags_long <- uw_small_frags %>%
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:trial), names_to = "exon_id", values_to = "metric") %>%
  left_join(exons_gc_content, by = "exon_id")

uw_small_frags_cor <- cor.test(uw_small_frags_long$metric, uw_small_frags_long$gc_content)
uw_small_frags_cor_text <- make_cor_text(cor_output = uw_small_frags_cor)
  
uw_small_frags_gc_plot <- uw_small_frags_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("Proportion Small Reads") +
  annotate("text", label = uw_small_frags_cor_text, x = 0.75, y = 0.95) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: Small Fragments") +
  theme(
    aspect.ratio = 1
  )


grail_small_frags_long <- grail_small_frags %>%
  pivot_longer(cols = -c(sample:phenotype), names_to = "exon_id", values_to = "metric") %>%
  left_join(exons_gc_content, by = "exon_id")

grail_small_frags_cor <- cor.test(grail_small_frags_long$metric, grail_small_frags_long$gc_content)
grail_small_frags_cor_text <- make_cor_text(cor_output = grail_small_frags_cor)

grail_small_frags_gc_plot <- grail_small_frags_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("Proportion Small Reads") +
  annotate("text", label = grail_small_frags_cor_text, x = 0.75, y = 0.95) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: Small Fragments") +
  theme(
    aspect.ratio = 1
  )

# MDS
uw_mds_long <- uw_mds %>%
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:trial), names_to = "exon_id", values_to = "metric") %>%
  left_join(exons_gc_content, by = "exon_id")

uw_mds_cor <- cor.test(uw_mds_long$metric, uw_mds_long$gc_content)
uw_mds_cor_text <- make_cor_text(cor_output = uw_mds_cor)

uw_mds_gc_plot <- uw_mds_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("MDS") +
  annotate("text", label = uw_mds_cor_text, x = 0.75, y = 1.0) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: MDS") +
  theme(
    aspect.ratio = 1
  )

grail_mds_long <- grail_mds %>%
  pivot_longer(cols = -c(sample:phenotype), names_to = "exon_id", values_to = "metric") %>%
  left_join(exons_gc_content, by = "exon_id")

grail_mds_cor <- cor.test(grail_mds_long$metric, grail_mds_long$gc_content)
grail_mds_cor_text <- make_cor_text(cor_output = grail_mds_cor)

grail_mds_gc_plot <- grail_mds_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("MDS") +
  annotate("text", label = grail_mds_cor_text, x = 0.75, y = 0.75) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: MDS") +
  theme(
    aspect.ratio = 1
  )


# SE
uw_se_long <- uw_se %>%
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:subtype), names_to = "exon_id", values_to = "metric") %>%
  left_join(exons_gc_content, by = "exon_id")

uw_se_cor <- cor.test(uw_se_long$metric, uw_se_long$gc_content)
uw_se_cor_text <- make_cor_text(cor_output = uw_se_cor)

uw_se_gc_plot <- uw_se_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("Shannon Entropy") +
  annotate("text", label = uw_se_cor_text, x = 0.75, y = 0.75) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: Shannon Entropy") +
  theme(
    aspect.ratio = 1
  )

grail_se_long <- grail_se %>%
  pivot_longer(cols = -c(sample:phenotype), names_to = "exon_id", values_to = "metric") %>%
  left_join(exons_gc_content, by = "exon_id")

grail_se_cor <- cor.test(grail_se_long$metric, grail_se_long$gc_content)
grail_se_cor_text <- make_cor_text(cor_output = grail_se_cor)

grail_se_gc_plot <- grail_se_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("Shannon Entropy") +
  annotate("text", label = grail_se_cor_text, x = 0.75, y = 0.75) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: Shannon Entropy") +
  theme(
    aspect.ratio = 1
  )

# TFBS
uw_tfbs_long <- uw_tfbs %>% 
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:trial), names_to = "TF", values_to = "metric") %>% 
  left_join(tfbs_gc_uw, by = c("TF"))

uw_tfbs_cor <- cor.test(uw_tfbs_long$metric, uw_tfbs_long$gc_content)
uw_tfbs_cor_text <- make_cor_text(cor_output = uw_tfbs_cor)


uw_tfbs_gc_plot <- uw_tfbs_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("TFBS Entropy") +
  annotate("text", label = uw_tfbs_cor_text, x = 0.75, y = 3.8) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  # theme_classic() +
  ggtitle("Cohort: UW\nFeature: TFBS Entropy") +
  theme(
    aspect.ratio = 1
  )

grail_tfbs_long <- grail_tfbs %>% 
  pivot_longer(cols = -c(sample:phenotype), names_to = "TF", values_to = "metric") %>% 
  left_join(tfbs_gc_uw, by = c("TF"))

grail_tfbs_cor <- cor.test(grail_tfbs_long$metric, grail_tfbs_long$gc_content)
grail_tfbs_cor_text <- make_cor_text(cor_output = grail_tfbs_cor)

grail_tfbs_gc_plot <- grail_tfbs_long %>% 
  ggplot(aes(gc_content, metric)) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("TFBS Entropy") +
  annotate("text", label = grail_tfbs_cor_text, x = 0.75, y = 4) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: TFBS Entropy") +
  theme(
    aspect.ratio = 1
  )

# Full Genes Depth
uw_full_gene_depth_long <- uw_full_gene_depth %>% 
  filter(sample %in% samples_in_manuscript) %>% 
  pivot_longer(cols = -c(sample:subtype), names_to = "gene", values_to = "metric") %>% 
  left_join(full_genes_gc_content, by = "gene")

uw_full_gene_depth_cor <- cor.test(uw_full_gene_depth_long$metric, uw_full_gene_depth_long$gc_content)
uw_full_gene_depth_cor_text <- make_cor_text(cor_output = uw_full_gene_depth_cor)

uw_full_gene_depth_gc_plot <- uw_full_gene_depth_long %>% 
  ggplot(aes(gc_content, log10(metric))) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("log10(Full Gene Depth)") +
  annotate("text", label = uw_full_gene_depth_cor_text, x = 0.65, y = -1.3) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: UW\nFeature: Full Gene Depth") +
  theme(
    aspect.ratio = 1
  )


grail_full_gene_depth_long <- grail_full_gene_depth %>% 
  pivot_longer(cols = -c(sample:phenotype), names_to = "gene", values_to = "metric") %>% 
  left_join(full_genes_gc_content, by = "gene")

grail_full_gene_depth_cor <- cor.test(grail_full_gene_depth_long$metric, grail_full_gene_depth_long$gc_content)
grail_full_gene_depth_cor_text <- make_cor_text(cor_output = grail_full_gene_depth_cor)

grail_full_gene_depth_gc_plot <- grail_full_gene_depth_long %>% 
  ggplot(aes(gc_content, log10(metric))) + 
  geom_bin_2d(bins = 100) + 
  geom_smooth(method = "lm", se = FALSE) +
  xlab("GC Content") +
  ylab("log10(Full Gene Depth)") +
  annotate("text", label = grail_full_gene_depth_cor_text, x = 0.65, y = -1.7) +
  # scale_fill_continuous(type = "viridis") + 
  color_scale +
  ggtitle("Cohort: GRAIL\nFeature: Full Gene Depth") +
  theme(
    aspect.ratio = 1
  )


grail_gc_plots <- (grail_atac_gc_plot + grail_depth_gc_plot) / (grail_frag_bins_gc_plot + grail_full_gene_depth_gc_plot) / (grail_mds_gc_plot + grail_se_gc_plot) / (grail_small_frags_gc_plot + grail_tfbs_gc_plot)
ggsave(filename = "output/figures/Figure_S3_GC_Content_GRAIL.pdf", plot = grail_gc_plots, height = 16, width = 10)

uw_gc_plots <- (uw_atac_gc_plot + uw_depth_gc_plot) / (uw_frag_bins_gc_plot + uw_full_gene_depth_gc_plot) / (uw_mds_gc_plot + uw_se_gc_plot) / (uw_small_frags_gc_plot + uw_tfbs_gc_plot)
ggsave(filename = "output/figures/Figure_S3_GC_Content_UW.pdf", plot = uw_gc_plots, height = 16, width = 10)
