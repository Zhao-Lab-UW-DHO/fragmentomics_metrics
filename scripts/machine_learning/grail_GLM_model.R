library(tidyverse)
library(rsample)
library(tidymodels)
library(parsnip)
library(glmnet)
library(fmsb)
library(pROC)
library(doSNOW)
library(ggplot2)
library(dplyr)
library(themis)
library(patchwork)

cohort <- "GRAIL"

feature_tables <- c(
  "E1SE",
  "E1depth",
  "all_exons_SE",
  "all_exons_depth",
  "full_gene_depth",
  "fragment_bins",
  "small_fragments",
  "small_fragments_all_exons",
  "MDS",
  "MDS_all_exons",
  "TFBS_entropy",
  "ATAC_entropy",
  "all_combined"
)

gene_panels <- c("default", "tempus", "guardant", "foundationOne")
# gene_panels <- "default"

starttime <- Sys.time()
GRAIL_comparison_model_results <- tibble()
GRAIL_comparison_ROCs <- tibble()
for (ft in feature_tables) {
  for (gp in gene_panels) {
    # Set variable for file names
    out_text <- str_c(ft, "_", gp)
    
    ft_path <- case_when(
      ft == "E1SE" ~ "se_E1.rds",
      ft == "E1depth" ~ "depth_E1.rds",
      ft == "all_exons_SE" ~ "se.rds",
      ft == "all_exons_depth" ~ "depth.rds",
      ft == "full_gene_depth" ~ "full_gene_depth.rds",
      ft == "fragment_bins" ~ "frag_bins.rds",
      ft == "small_fragments" ~ "small_frags_E1.rds",
      ft == "small_fragments_all_exons" ~ "small_frags.rds",
      ft == "MDS" ~ "mds_E1.rds",
      ft == "MDS_all_exons" ~ "mds.rds",
      ft == "TFBS_entropy" ~ "TFBS_entropy.rds",
      ft == "ATAC_entropy" ~ "ATAC_entropy.rds",
      ft == "all_combined" ~ "all_combined.rds"
    )
    
    feature_table_to_use <- str_c("output/feature_tables/grail/", gp, "/", ft_path)
    
    # Load in correct helper function
    print("Loading Helper Functions and Data")
    
    source("scripts/machine_learning/helper_load_grail.R")
    
    # set variables
    usegrail <- T
    
    if(usegrail) {
      df2use <- grail_df
      colors2use <- colorsGRAIL
    } else {
      df2use <- internal_df
      colors2use <- colorsUW
    }
    
    print("All samples table")
    print(table(df2use$phenotype))
  
    
    # TRAINING ----------------------------------------------------------------
    
    for (seed in seq(1,25)) {
    # for (seed in seq(25,25)) {
      set.seed(seed)
      folds <- vfold_cv(df2use, v=nfold, strata=phenotype, repeats=1)
      
      print('looping through folds to train models')
      print(out_text)
      print(str_c("Seed: ", seed))
      grail_training_preds <- NULL
      model_coefs <- tibble()
      for(i in 1:nfold) {
        print(i)
        split <- folds$splits[[i]]
        traindatafold <- split %>% analysis() %>% ungroup()
        traindatafold <- traindatafold %>% dplyr::select(-c(ctdna_fraction, sample))
        fit <- trainmodel(traindatafold)
        
        validdatafold <- split %>% assessment()
        pred <- predictmodel(fit,validdatafold)
        grail_training_preds <- bind_rows(grail_training_preds,pred)
        
        split_coefs <- extract_fit_parsnip(fit) %>% broom::tidy() %>%
          mutate(fold = i)
        model_coefs <- bind_rows(model_coefs, split_coefs)
      }
      
      # Data Collecting ---------------------------------------------------------
      
      # source("feature_comparison_scripts/roc_extraction_functions.R")
      
      temp_model_output <- grail_training_preds %>% mutate(feature_table = ft, gene_panel = gp, seed = seed)
      GRAIL_comparison_model_results <- bind_rows(GRAIL_comparison_model_results, temp_model_output)
      
      temp_ROC_results <- extract_ROC(predictions = grail_training_preds, feature_table = ft, gene_panel = gp, seedused = seed)
      GRAIL_comparison_ROCs <- bind_rows(GRAIL_comparison_ROCs, temp_ROC_results)
      
      # gets timestamp and appends to file output
      file_timestamp <- str_replace_all(string = Sys.time(), pattern = "[-:. ]", replacement = "_") %>% str_sub(1, 19)
      model_output_filename <- str_c("GRAIL_comparison_model_results_", file_timestamp, ".rds")
      roc_output_filename <- str_c("GRAIL_comparison_ROCs_", file_timestamp, ".rds")
      
      # saves after every iteration, in case something goes wrong down the line and needs to restart
      if (seed == 25) {
        saveRDS(GRAIL_comparison_model_results, str_c("output/machine_learning/", model_output_filename))
        saveRDS(GRAIL_comparison_ROCs, str_c("output/machine_learning/", roc_output_filename))
      }
      
      if (seed == 25 & ft == "all_combined") {
        saveRDS(GRAIL_comparison_model_results, "output/machine_learning/grail_predictions.rds")
        saveRDS(GRAIL_comparison_ROCs, "output/machine_learning/grail_AUROCs.rds")
      }
    }
  }
}

elapsedtime <- Sys.time() - starttime
print(elapsedtime)
