library(tidyverse)
library(RColorBrewer)
library(pals)
library(patchwork)

setwd("/mnt/Data01/projects/fragmentomics/scripts/paper_2_code_upload/fragmentomics_metrics/")

# Load in Data

# Machine Learning Outputs
grail_comparison_preds <- read_rds("files/manuscript_data/grail_predictions.rds")
grail_comparison_ROCs <- read_rds("files/manuscript_data/grail_AUROCs.rds")

uw_comparison_preds <- read_rds("files/manuscript_data/uw_predictions.rds")
uw_comparison_ROCs <- read_rds("files/manuscript_data/uw_AUROCs.rds")

mixed_sample_ROCs <- read_rds("files/manuscript_data/grail_mixing_AUROCs.rds")

grail_downsampling_ROCs <- read_rds("files/manuscript_data/grail_downsampling_AUROCs.rds")

# metadata
uw_sample_phenotypes <- read_rds("files/manuscript_data/uw_sample_phenotypes.rds")
grail_sample_phenotypes <- read_rds("files/manuscript_data/grail_sample_phenotypes.rds")

# Fragment Counts and Stats
uw_fragstats <- read_rds("files/manuscript_data/uw_global_fragstats.rds")
grail_fragstats <- read_rds("files/manuscript_data/grail_global_fragstats.rds")


# Data Cleanup
remove_fts <- TRUE

remove_and_rename_fts <- function(df, fts_remove) {
  output <- df %>% 
    filter(!(feature_table %in% fts_remove)) %>% 
    mutate(feature_table = str_replace_all(feature_table, "_", " ")) %>% 
    mutate(
      feature_table = case_when(
        feature_table == "E1SE" ~ "E1 SE",
        feature_table == "E1depth" ~ "E1 Depth",
        feature_table == "MDS" ~ "E1 MDS",
        feature_table == "MDS all exons" ~ "All Exons MDS",
        feature_table == "small fragments" ~ "E1 Small Fragments",
        feature_table == "small fragments all exons" ~ "All Exons Small Fragments",
        feature_table == "all exons SE" ~ "All Exons SE",
        feature_table == "all exons depth" ~ "All Exons Depth",
        feature_table == "fragment bins" ~ "Fragment Bins",
        feature_table == "full gene depth" ~ "Full Gene Depth",
        feature_table == "ATAC entropy" ~ "ATAC Entropy",
        feature_table == "TFBS entropy" ~ "TFBS Entropy",
        feature_table == "all combined" ~ "All Combined",
        TRUE ~ feature_table
      )
    )
  return(output)
}

rename_mixed_data_fts <- function(df) {
  output <- df %>% 
    mutate(feature_table = str_replace_all(feature_table, "_", " ")) %>% 
    mutate(
      feature_table = case_when(
        feature_table == "depth" ~ "Depth",
        feature_table == "frag bins" ~ "Fragment Bins",
        feature_table == "full gene depth" ~ "Full Gene Depth",
        feature_table == "mds" ~ "MDS",
        feature_table == "small frags" ~ "Small Fragments",
        feature_table == "ATAC entropy" ~ "ATAC Entropy",
        feature_table == "TFBS entropy" ~ "TFBS Entropy",
        feature_table == "se" ~ "SE",
        TRUE ~ feature_table
      )
    )
  return(output)
}

if (remove_fts) {
  fts_to_remove <- c("E1SE_and_E1depth", "all_exons_SE_and_depth", "last_exon_depth")
  
  grail_comparison_preds <- remove_and_rename_fts(df = grail_comparison_preds, fts_remove = fts_to_remove)
  grail_comparison_ROCs <- remove_and_rename_fts(df = grail_comparison_ROCs, fts_remove = fts_to_remove)
  
  uw_comparison_preds <- remove_and_rename_fts(df = uw_comparison_preds, fts_remove = fts_to_remove)
  uw_comparison_ROCs <- remove_and_rename_fts(df = uw_comparison_ROCs, fts_remove = fts_to_remove)
  
  mixed_sample_ROCs <- rename_mixed_data_fts(df = mixed_sample_ROCs)
  
  grail_downsampling_ROCs <- rename_mixed_data_fts(df = grail_downsampling_ROCs)
}


### Setting up global color vars

colorpal <- c('black',
              '#1B9E77',
              '#D95F02',
              '#7570B3',
              '#E7298A',
              '#666666',
              '#66A61E',
              '#A6761D',
              "deepskyblue2",
              "red2",
              "gold",
              "darkorange"
              )
names(colorpal) <- c('All',
                     'Bladder',
                     'Breast',
                     'Lung',
                     'NEPC',
                     'Normal',
                     'Prostate',
                     'RCC',
                     "ERpos",
                     "ERneg",
                     "NSCLC",
                     "SCLC")


feature_colorpal <- c("#1F601A", "#33A02C", "#B2DF8A", 
                      "#1F78B4", "#A6CEE3", 
                      "#E31A1C", "#FB9A99", 
                      "#6A3D9A", "#CAB2D6", 
                      "gold",
                      # 'grey40',
                      'darkorange',
                      "#765341",
                      "black",
                      "grey30")
names(feature_colorpal) <- c("All Exons Depth", "Full Gene Depth", "E1 Depth", 
                             "All Exons SE", "E1 SE", 
                             "All Exons MDS", "E1 MDS", 
                             "All Exons Small Fragments", "E1 Small Fragments", 
                             "Fragment Bins",
                             "TFBS Entropy",
                             "ATAC Entropy",
                             "Griffin",
                             "All Combined")

### Plotting

panel2plot <- "default"
# panel2plot <- "tempus"
# panel2plot <- "guardant"
# panel2plot <- "foundationOne"

grail_median_ROCs <- grail_comparison_ROCs %>% 
  group_by(feature_table, gene_panel, phenotype) %>% 
  summarise(med_roc = median(roc)) %>% 
  ungroup()

### TABLE S1 ###

grail_large_table_summary_medians <- grail_median_ROCs %>% 
  group_by(feature_table, gene_panel) %>% 
  summarise(mean_roc = mean(med_roc)) %>% 
  ungroup()

grail_large_table_summary <- grail_median_ROCs %>% 
  pivot_wider(names_from = phenotype, values_from = med_roc) %>% 
  left_join(grail_large_table_summary_medians, by = c("feature_table", "gene_panel")) %>% 
  mutate(gene_panel = case_when(
    gene_panel == "default" ~ "GRAIL",
    gene_panel == "foundationOne" ~ "FoundationOne Liquid CDx",
    gene_panel == "guardant" ~ "Guardant360 CDx",
    gene_panel == "tempus" ~ "Tempus xF"
  )) %>% 
  mutate(gene_panel = factor(gene_panel, level = c("GRAIL", "FoundationOne Liquid CDx", "Guardant360 CDx", "Tempus xF"))) %>% 
  arrange(gene_panel, desc(mean_roc))

write_tsv(grail_large_table_summary, "output/figures/Table_S1_GRAIL_large_table_summary.tsv")

### TABLE S2 ###

uw_median_ROCs <- uw_comparison_ROCs %>% 
  group_by(feature_table, gene_panel, phenotype) %>% 
  summarise(med_roc = median(roc)) %>% 
  ungroup()

uw_large_table_summary_medians <- uw_median_ROCs %>% 
  group_by(feature_table, gene_panel) %>% 
  summarise(mean_roc = mean(med_roc)) %>% 
  ungroup()

uw_large_table_summary <- uw_median_ROCs %>% 
  pivot_wider(names_from = phenotype, values_from = med_roc) %>% 
  left_join(uw_large_table_summary_medians, by = c("feature_table", "gene_panel")) %>% 
  mutate(gene_panel = case_when(
    gene_panel == "default" ~ "UW",
    gene_panel == "foundationOne" ~ "FoundationOne Liquid CDx",
    gene_panel == "guardant" ~ "Guardant360 CDx",
    gene_panel == "tempus" ~ "Tempus xF"
  )) %>% 
  mutate(gene_panel = factor(gene_panel, level = c("UW", "FoundationOne Liquid CDx", "Guardant360 CDx", "Tempus xF"))) %>% 
  arrange(gene_panel, desc(mean_roc))

write_tsv(uw_large_table_summary, "output/figures/Table_S2_UW_large_table_summary.tsv")


### FIGURE 2 ###

uw_metric_levels <- names(feature_colorpal)
uw_phenotype_levels <- c("Normal", "Prostate", "NEPC", "ERpos", "ERneg", "NSCLC", "SCLC", "Bladder", "RCC")

fig2 <- ggplot() + 
  geom_boxplot(data = uw_comparison_ROCs %>% 
                 filter(feature_table != "Griffin") %>% 
                 filter(gene_panel == panel2plot) %>% 
                 mutate(feature_table = factor(feature_table, levels = uw_metric_levels)) %>% 
                 mutate(phenotype = factor(phenotype, levels = uw_phenotype_levels)),
               aes(feature_table, roc, color = feature_table), outlier.shape = NA) +
  geom_point(data = uw_comparison_ROCs %>% 
               filter(gene_panel == panel2plot) %>% 
               filter(feature_table != "Griffin") %>% 
               mutate(feature_table = factor(feature_table, levels = uw_metric_levels)) %>% 
               mutate(phenotype = factor(phenotype, levels = uw_phenotype_levels)),
             aes(feature_table, roc, color = feature_table), size = 2, alpha = 0.2) + 
  facet_wrap(~ phenotype) + 
  ggtitle("UW") +
  xlab("Feature Used") + 
  ylab("AUROC") +
  labs(color = "Phenotype") +
  scale_color_manual(values = feature_colorpal) +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    aspect.ratio = 1,
    legend.position = "none"
  )

ggsave(plot = fig2, filename = "output/figures/Figure_2_UW_AUROCs.pdf", height = 10, width = 10)

### FIGURE 3 ###

grail_metric_levels <- names(feature_colorpal)
grail_phenotype_levels <- c("Normal", "Prostate", "Breast", "Lung")

fig3 <- ggplot() + 
  geom_boxplot(data = grail_comparison_ROCs %>% 
                 filter(gene_panel == panel2plot) %>% 
                 mutate(feature_table = factor(feature_table, levels = grail_metric_levels)) %>% 
                 mutate(phenotype = factor(phenotype, levels = grail_phenotype_levels)),
               aes(feature_table, roc, color = feature_table), outlier.shape = NA) +
  geom_point(data = grail_comparison_ROCs %>% 
               filter(gene_panel == panel2plot) %>% 
               mutate(feature_table = factor(feature_table, levels = grail_metric_levels)) %>% 
               mutate(phenotype = factor(phenotype, levels = grail_phenotype_levels)),
             aes(feature_table, roc, color = feature_table), size = 2, alpha = 0.2) + 
  facet_wrap(~ phenotype) + 
  ggtitle("GRAIL") +
  xlab("Feature Used") + 
  ylab("AUROC") +
  labs(color = "Phenotype") +
  scale_color_manual(values = feature_colorpal, breaks = names(feature_colorpal)) +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    aspect.ratio = 1,
    legend.position = "none"
  )

ggsave(plot = fig3, filename = "output/figures/Figure_3_GRAIL_AUROCs.pdf", height = 8, width = 8)

### Comparison of panels ###

gene_panel_colors <- brewer.pal(n = 5, name = "Set1")
names(gene_panel_colors) <- c("UW", "GRAIL", "FoundationOne Liquid CDx", "Guardant360 CDx", "Tempus xF")

### FIGURE 4A ###

fig4A <- uw_median_ROCs %>%
  filter(feature_table != "Griffin") %>% 
  mutate(feature_table = factor(feature_table, levels = uw_metric_levels)) %>% 
  mutate(phenotype = factor(phenotype, levels = uw_phenotype_levels)) %>% 
  mutate(gene_panel = case_when(
    gene_panel == "default" ~ "UW",
    gene_panel == "foundationOne" ~ "FoundationOne Liquid CDx",
    gene_panel == "guardant" ~ "Guardant360 CDx",
    gene_panel == "tempus" ~ "Tempus xF"
  )) %>% 
  mutate(gene_panel = factor(gene_panel, levels = c("UW", "FoundationOne Liquid CDx", "Guardant360 CDx", "Tempus xF"))) %>% 
  ggplot(aes(feature_table, med_roc, color = feature_table, shape = gene_panel)) +
  geom_point(size = 2.5) + 
  facet_wrap(~ phenotype) + 
  xlab("Feature Used") + 
  ylab("Median AUROC") +
  labs(color = "Feature Table", shape = "Gene Panel") +
  scale_color_manual(values = feature_colorpal, breaks = names(feature_colorpal), guide = "none") +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  )

ggsave(plot = fig4A, filename = "output/figures/Figure_4A_UW_Commercial_Panel_AUROCs.pdf", height = 12, width = 14)

### FIGURE 4B ###

fig4B <- uw_median_ROCs %>%
  filter(feature_table != "Griffin") %>% 
  pivot_wider(names_from = gene_panel, values_from = med_roc) %>%
  mutate(
    FoundationOne = foundationOne - default,
    Guardant = guardant - default,
    Tempus = tempus - default,
  ) %>% 
  dplyr::select(feature_table, phenotype, FoundationOne:Tempus) %>% 
  pivot_longer(cols = c(FoundationOne, Guardant, Tempus), names_to = "gene_panel", values_to = "roc_diff") %>% 
  mutate(feature_table = factor(feature_table, levels = names(feature_colorpal))) %>% 
  ggplot(aes(gene_panel, roc_diff, fill = feature_table)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.8) + 
  geom_point(position = position_dodge(0.8), size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = feature_colorpal, breaks = names(feature_colorpal)) +
  scale_color_manual(values = feature_colorpal, breaks = names(feature_colorpal)) +
  xlab("Gene Panel") + 
  ylab("Change in median AUROC") + 
  labs(fill = "Feature") +
  ggtitle("UW Cohort") +
  theme_bw() + 
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(plot = fig4B, filename = "output/figures/Figure_4B_UW_Commercial_Panel_Diffs.pdf", height = 8, width = 12)


### FIGURE 5A ###

fig5A <- grail_median_ROCs %>%
  mutate(feature_table = factor(feature_table, levels = grail_metric_levels)) %>% 
  mutate(phenotype = factor(phenotype, levels = grail_phenotype_levels)) %>% 
  mutate(gene_panel = case_when(
    gene_panel == "default" ~ "GRAIL",
    gene_panel == "foundationOne" ~ "FoundationOne Liquid CDx",
    gene_panel == "guardant" ~ "Guardant360 CDx",
    gene_panel == "tempus" ~ "Tempus xF"
  )) %>% 
  mutate(gene_panel = factor(gene_panel, levels = c("GRAIL", "FoundationOne Liquid CDx", "Guardant360 CDx", "Tempus xF"))) %>% 
  ggplot(aes(feature_table, med_roc, color = feature_table, shape = gene_panel)) + 
  geom_point(size = 3) + 
  facet_wrap(~ phenotype) + 
  xlab("Feature Used") + 
  ylab("Median AUROC") +
  labs(color = "Feature Table", shape = "Gene Panel") +
  scale_color_manual(values = feature_colorpal, breaks = names(feature_colorpal), guide = "none") +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    strip.text = element_text(size = 12, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  )

ggsave(plot = fig5A, filename = "output/figures/Figure_5A_GRAIL_Commercial_Panel_AUROCs.pdf", height = 8, width = 10)


### FIGURE 5B ###

fig5B <- grail_median_ROCs %>% 
  pivot_wider(names_from = gene_panel, values_from = med_roc) %>%
  mutate(
    FoundationOne = foundationOne - default,
    Guardant = guardant - default,
    Tempus = tempus - default,
  ) %>% 
  dplyr::select(feature_table, phenotype, FoundationOne:Tempus) %>% 
  pivot_longer(cols = c(FoundationOne, Guardant, Tempus), names_to = "gene_panel", values_to = "roc_diff") %>% 
  mutate(feature_table = factor(feature_table, levels = names(feature_colorpal))) %>% 
  ggplot(aes(gene_panel, roc_diff, fill = feature_table)) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.8) +
  geom_point(position = position_dodge(0.8), size = 1.5, shape = 21, color = "black") +
  scale_fill_manual(values = feature_colorpal, breaks = names(feature_colorpal)) +
  scale_color_manual(values = feature_colorpal, breaks = names(feature_colorpal)) +
  xlab("Gene Panel") + 
  ylab("Change in median AUROC") + 
  labs(fill = "Feature") +
  ggtitle("GRAIL Cohort") +
  theme_bw() + 
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )

ggsave(plot = fig5B, filename = "output/figures/Figure_5B_Commercial_Panel_Diffs.pdf", height = 8, width = 12)





### FIGURE 6B ###

ctdna_levels <- c("0.1 - 1",
                  "0.01 - 0.1",
                  "0.001 - 0.01",
                  "0.0001 - 0.001",
                  "0.00001 - 0.0001")

mixed_sample_ROCs_ctdna_leves <- mixed_sample_ROCs %>% 
  filter(ctdna_group != "all") %>% 
  mutate(ctdna_group = factor(ctdna_group, levels = ctdna_levels)) %>% 
  mutate(phenotype = factor(phenotype, levels = grail_phenotype_levels))

fig6B <- mixed_sample_ROCs_ctdna_leves %>% 
  ggplot(aes(ctdna_group, roc, color = phenotype)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(), size = 0.25) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.25, size = 0.75) + 
  scale_color_manual(values = colorpal) +
  facet_wrap(~ feature_table) + 
  xlab("ctDNA Bin") + 
  ylab("AUROC") +
  labs(color = "Phenotype") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    aspect.ratio = 1,
    legend.position = c(5/6, 0.12)
  )

ggsave(plot = fig6B, filename = "output/figures/Figure_6B_GRAIL_Mixing_ctDNA_AUROCs.pdf", height = 8, width = 8)


### FIGURE S1 ###

uw_genes <- read_tsv("files/manuscript_data/exon_panels/uw_genes.tsv", col_names = "gene") %>% pull(gene)
grail_genes <- read_tsv("files/manuscript_data/exon_panels/grail_genes.tsv") %>% pull(gene)
tempus_genes <- read_tsv("files/manuscript_data/exon_panels/tempus_xF.tsv") %>% pull(gene)
f1_genes <- read_tsv("files/manuscript_data/exon_panels/foundationOneCDx.tsv") %>% pull(gene)
guardant_genes <- read_tsv("files/manuscript_data/exon_panels/guardant360CDx.tsv") %>% pull(gene)

base_gene_sets <- list(uw_genes, grail_genes)
subset_gene_sets <- list(tempus_genes, f1_genes, guardant_genes)

base_vs_subset_data <- matrix(0, nrow = length(base_gene_sets), ncol = length(subset_gene_sets))
for (i in 1:length(base_gene_sets)) {
  for (j in 1:length(subset_gene_sets)) {
    base_vs_subset_data[i, j] <- length(intersect(base_gene_sets[[i]], subset_gene_sets[[j]]))
  }
}

base_vs_subset_data_df <- as.data.frame(base_vs_subset_data)
rownames(base_vs_subset_data_df) <- c("uw_genes", "grail_genes")
colnames(base_vs_subset_data_df) <- c("tempus_genes", "f1_genes", "guardant_genes")

base_vs_subset_data_tibble <- as_tibble(base_vs_subset_data_df, rownames = "base_panel") %>% 
  pivot_longer(cols = -base_panel, names_to = "subset_panel", values_to = "count") %>% 
  mutate(
    total_in_subset = case_when(
      subset_panel == "tempus_genes" ~ 105,
      subset_panel == "f1_genes" ~ 309,
      subset_panel == "guardant_genes" ~ 55
    )
  ) %>% 
  mutate(fraction_in_subset = count / total_in_subset) %>% 
  mutate(tile_text = str_c(count, "/", total_in_subset))

fig_s1 <- base_vs_subset_data_tibble %>% 
  mutate(
    base_panel = case_when(
      base_panel == "uw_genes" ~ "UW Panel\n(Gene N = 821)",
      base_panel == "grail_genes" ~ "GRAIL Panel\n(Gene N = 508)"
    )
  ) %>% 
  mutate(
    subset_panel = case_when(
      subset_panel == "tempus_genes" ~ "Tempus xF",
      subset_panel == "guardant_genes" ~ "Guardant360 CDx", 
      subset_panel == "f1_genes" ~ "FoundationOne Liquid CDx"
    )
  ) %>% 
  ggplot(aes(base_panel, subset_panel, fill = fraction_in_subset)) + 
  geom_tile() + 
  geom_text(aes(label = tile_text), size = 6) + 
  ylab("Number of genes from commercial panel [A]...") +
  xlab("...in analysis panel [B]") +
  labs(fill = "Fraction in Analysis Panel") +
  scale_fill_viridis_c(limits = c(0,1)) +
  scale_x_discrete(expand = c(0,0)) + 
  scale_y_discrete(expand = c(0,0)) + 
  theme_bw() + 
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.title.x = element_text(margin = margin(t = 20)),
    axis.title.y = element_text(margin = margin(r = 20)),
    axis.text = element_text(size = 12, color = "black"),
    aspect.ratio = 1
  )

ggsave(plot = fig_s1, filename = "output/figures/Figure_S1_Panel_Overlaps.pdf", height = 8, width = 8)


### FIGURE S2 ###

# UW
uw_read_count_data <- uw_fragstats %>% 
  group_by(sample) %>% 
  summarise(total_reads = sum(count)) %>% 
  ungroup() %>% 
  left_join(uw_sample_phenotypes, by = "sample")

pheno_order <- c("Prostate", "NEPC", "ERpos", "ERneg", "SCLC", "NSCLC", "Bladder", "RCC", "Normal")

pheno_colors <- c(rev(RColorBrewer::brewer.pal(n = 6, name = "Paired")), RColorBrewer::brewer.pal(n = 5, name = "Set1")[4:5], "black")
names(pheno_colors) <- pheno_order

uw_read_count_plot <- uw_read_count_data %>% 
  mutate(pheno_simple = factor(pheno_simple, levels = pheno_order)) %>% 
  ggplot(aes(pheno_simple, total_reads / 1e6, color = pheno_simple)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  geom_jitter(width = 0.1) + 
  scale_y_continuous(limits = c(0,62), expand = c(0,0)) +
  scale_color_manual(values = pheno_colors) +
  xlab("Phenotype") + 
  ylab("Total Reads Post-Deduplication (Millions)") + 
  ggtitle("UW Cohort") +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 20, color = "black", hjust = 0.5),
    legend.position = "none"
  )

# GRAIL
grail_read_count_data <- grail_fragstats %>% 
  group_by(sample) %>% 
  summarise(total_reads = sum(count)) %>% 
  ungroup() %>% 
  left_join(grail_sample_phenotypes, by = "sample")

grail_pheno_order <- c("Prostate", "Breast", "Lung", "Normal")

grail_pheno_colors <- c("#E31A1C", "#33A02C", "#1F78B4", "#000000")
names(grail_pheno_colors) <- grail_pheno_order

grail_read_count_plot <- grail_read_count_data %>% 
  mutate(phenotype = factor(phenotype, levels = grail_pheno_order)) %>% 
  ggplot(aes(phenotype, total_reads / 1e6, color = phenotype)) + 
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  geom_jitter(width = 0.1) + 
  scale_color_manual(values = grail_pheno_colors) +
  xlab("Phenotype") + 
  ylab("Total Reads Post-Deduplication (Millions)") + 
  ggtitle("GRAIL Cohort") +
  theme_classic() + 
  theme(
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 20, color = "black", hjust = 0.5),
    legend.position = "none"
  )

# Combine into one plot
fig_s2 <- uw_read_count_plot + grail_read_count_plot

ggsave(plot = fig_s2, filename = "output/figures/Figure_S2_Reads_Post_Dedup.pdf", height = 10, width = 16)


### FIGURE S3 ###

# Check if feature table files are present
gc_files_exist_uw <- all(file.exists("output/feature_tables/default/ATAC_entropy.rds",
                                 "output/feature_tables/default/depth.rds",
                                 "output/feature_tables/default/frag_bins.rds",
                                 "output/feature_tables/default/small_frags.rds",
                                 "output/feature_tables/default/mds.rds",
                                 "output/feature_tables/default/se.rds",
                                 "output/feature_tables/default/TFBS_entropy.rds",
                                 "output/feature_tables/default/full_gene_depth.rds"))
gc_files_exist_grail <- all(file.exists("output/feature_tables/grail/default/ATAC_entropy.rds",
                                    "output/feature_tables/grail/default/depth.rds",
                                    "output/feature_tables/grail/default/frag_bins.rds",
                                    "output/feature_tables/grail/default/small_frags.rds",
                                    "output/feature_tables/grail/default/mds.rds",
                                    "output/feature_tables/grail/default/se.rds",
                                    "output/feature_tables/grail/default/TFBS_entropy.rds",
                                    "output/feature_tables/grail/default/full_gene_depth.rds"))


if (gc_files_exist_uw & gc_files_exist_grail) {
  source("scripts/gc_content_plots.R")
}

### FIGURE S4 ###

fig_s4 <- ggplot() + 
  geom_boxplot(data = uw_comparison_ROCs %>% 
                 filter(feature_table != "All Combined") %>% 
                 filter(gene_panel == panel2plot) %>% 
                 mutate(feature_table = factor(feature_table, levels = uw_metric_levels)) %>% 
                 mutate(phenotype = factor(phenotype, levels = uw_phenotype_levels)),
               aes(feature_table, roc, color = feature_table), outlier.shape = NA) +
  geom_point(data = uw_comparison_ROCs %>% 
               filter(gene_panel == panel2plot) %>% 
               filter(feature_table != "All Combined") %>% 
               mutate(feature_table = factor(feature_table, levels = uw_metric_levels)) %>% 
               mutate(phenotype = factor(phenotype, levels = uw_phenotype_levels)),
             aes(feature_table, roc, color = feature_table), size = 2, alpha = 0.2) + 
  facet_wrap(~ phenotype) + 
  ggtitle("UW") +
  xlab("Feature Used") + 
  ylab("AUROC") +
  labs(color = "Phenotype") +
  scale_color_manual(values = feature_colorpal) +
  scale_y_continuous(limits = c(0.5, 1.0)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    aspect.ratio = 1,
    legend.position = "none"
  )

ggsave(plot = fig_s4, filename = "output/figures/Figure_S4_UW_Griffin_AUROCs.pdf", height = 10, width = 10)

### FIGURES S5 + S6 ###

plot_conf_matrix_facets_from_preds <- function(preds, gp, plot_title = gp) {
  
  preds <- preds %>% filter(gene_panel == gp)
  preds <- preds %>% filter(!(feature_table %in% c("All Combined", "Griffin")))
  
  feature_accs <- preds %>% 
    dplyr::count(feature_table, actual, pred_class) %>% 
    mutate(correct_pred = actual == pred_class) %>% 
    group_by(feature_table, correct_pred) %>% 
    summarise(count = sum(n)) %>% 
    ungroup() %>% 
    group_by(feature_table) %>% 
    mutate(total_n = sum(count)) %>% 
    filter(correct_pred) %>% 
    mutate(accuracy = count / total_n) %>% 
    mutate(acc_label = str_c(feature_table, "\nAccuracy: ", round(accuracy, 3))) %>% 
    dplyr::select(feature_table, acc_label) %>% 
    ungroup() %>% 
    deframe()
  
  out_plot <- preds %>% 
    dplyr::count(feature_table, actual, pred_class) %>% 
    group_by(feature_table, actual) %>% 
    mutate(total_actual = sum(n)) %>% 
    ungroup() %>% 
    mutate(prop = n / total_actual) %>% 
    mutate(tile_text = round(prop, digits = 3)) %>% 
    complete(feature_table, actual, pred_class, fill = list(prop = 0)) %>% 
    ggplot(aes(actual, pred_class, fill = prop)) + 
    geom_tile(color = "black") + 
    scale_fill_gradient(low = "white", high = "green4", limits = c(0,1)) + 
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) + 
    xlab("Actual") + 
    ylab("Predicted") +
    ggtitle(str_c("Panel: ", plot_title)) +
    theme_classic() + 
    facet_wrap(~ feature_table) +
    theme(
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 90)
    )
  
  return(out_plot)
}


uw_default_conf_plot <- plot_conf_matrix_facets_from_preds(preds = uw_comparison_preds, gp = "default", plot_title = "UW")
uw_tempus_conf_plot <- plot_conf_matrix_facets_from_preds(preds = uw_comparison_preds, gp = "tempus", plot_title = "Tempus")
uw_guardant_conf_plot <- plot_conf_matrix_facets_from_preds(preds = uw_comparison_preds, gp = "guardant", plot_title = "Guardant")
uw_foundationone_conf_plot <- plot_conf_matrix_facets_from_preds(preds = uw_comparison_preds, gp = "foundationOne", plot_title = "Foundation One")

fig_s5 <- (uw_default_conf_plot + uw_tempus_conf_plot) / (uw_guardant_conf_plot + uw_foundationone_conf_plot)
ggsave(filename = "output/figures/Figure_S5_UW_Conf_Plot.pdf", plot = fig_s5, height = 10, width = 16)

grail_default_conf_plot <- plot_conf_matrix_facets_from_preds(preds = grail_comparison_preds, gp = "default", plot_title = "GRAIL")
grail_tempus_conf_plot <- plot_conf_matrix_facets_from_preds(preds = grail_comparison_preds, gp = "tempus", plot_title = "Tempus")
grail_guardant_conf_plot <- plot_conf_matrix_facets_from_preds(preds = grail_comparison_preds, gp = "guardant", plot_title = "Guardant")
grail_foundationone_conf_plot <- plot_conf_matrix_facets_from_preds(preds = grail_comparison_preds, gp = "foundationOne", plot_title = "Foundation One")

fig_s6 <- (grail_default_conf_plot + grail_tempus_conf_plot) / (grail_guardant_conf_plot + grail_foundationone_conf_plot)
ggsave(filename = "output/figures/Figure_S6_GRAIL_Conf_Plot.pdf", plot = fig_s6, height = 10, width = 16)

### FIGURE S7 ###

fig_s7 <- grail_downsampling_ROCs %>% 
  mutate(downsample = factor(downsample, levels = c("100M", "50M", "25M", "10M", "5M", "1M"))) %>% 
  mutate(phenotype = factor(phenotype, levels = c("Normal", "Breast", "Lung", "Prostate"))) %>% 
  ggplot(aes(downsample, roc, color = phenotype)) + 
  geom_jitter(alpha = 0.25, position = position_dodge(width = 0.75)) +
  scale_color_manual(values = colorpal) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  facet_wrap(~ feature_table) + 
  xlab("Downsample Level") + 
  ylab("AUROC") +
  labs(color = "Phenotype") +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    aspect.ratio = 1,
    legend.position = c(5/6, 0.12)
  )

ggsave(plot = fig_s7, filename = "output/figures/Figure_S7_GRAIL_Downsampling.pdf", height = 10, width = 10)
